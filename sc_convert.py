#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# python sc_convert.py /path/to/jiaData.rds --assay-name RNA
# python sc_convert.py /path/to/jiaData.h5ad --assay-name RNA
"""
Python-first .rds <-> .h5ad converter
- Python is the main runtime
- R is embedded via rpy2
- strict mode: fail instead of silently dropping unsupported data
"""

import argparse
import gzip
import json
import re
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread, mmwrite
from rpy2 import robjects as ro


# =========================================================
# Logging
# =========================================================
def log(msg: str) -> None:
    print(f"[sc_convert] {msg}", flush=True)


# =========================================================
# Embedded R helpers
# =========================================================
R_HELPERS = r'''
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(jsonlite)
})

rlog <- function(msg) {
  message(sprintf("[R] %s", msg))
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]", "_", x)
}

make_key <- function(name) {
  clean <- gsub("[^A-Za-z0-9]", "", toupper(name))
  if (clean == "") clean <- "RED"
  paste0(clean, "_")
}

pick_assay <- function(obj, assay_name = NULL) {
  assays <- tryCatch(as.character(SeuratObject::Assays(obj)), error = function(e) names(obj@assays))
  if (!is.null(assay_name) && nzchar(assay_name) && assay_name %in% assays) {
    return(assay_name)
  }
  if ("RNA" %in% assays) return("RNA")
  da <- tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) NULL)
  if (!is.null(da) && nzchar(da)) return(da)
  assays[[1]]
}

write_lines_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(as.character(x), con = con, sep = "\n")
}

read_lines_gz <- function(path) {
  readLines(gzfile(path), warn = FALSE)
}

write_table_gz <- function(df, path) {
  if (is.null(df)) {
    df <- data.frame()
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  idx <- rownames(df)
  if (is.null(idx)) idx <- as.character(seq_len(nrow(df)))

  out <- data.frame(`__index__` = idx, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) > 0) {
    out <- cbind(out, df)
  }

  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.csv(out, con, row.names = FALSE, quote = TRUE, na = "")
}

read_table_gz <- function(path) {
  df <- read.csv(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  if (!"__index__" %in% colnames(df)) {
    stop("Missing __index__ column in: ", path)
  }
  rn <- as.character(df$`__index__`)
  df$`__index__` <- NULL
  rownames(df) <- rn
  df
}

write_matrix_mtx <- function(mat, path) {
  if (inherits(mat, "Matrix")) {
    Matrix::writeMM(mat, path)
  } else {
    Matrix::writeMM(Matrix::Matrix(mat, sparse = FALSE), path)
  }
}

safe_feature_meta <- function(obj, assay_name) {
  out <- tryCatch(obj[[assay_name]][[]], error = function(e) NULL)
  if (is.null(out)) {
    out <- data.frame(row.names = rownames(obj[[assay_name]]))
  }
  out
}

safe_layers <- function(obj, assay_name) {
  out <- tryCatch(SeuratObject::Layers(obj[[assay_name]]), error = function(e) character(0))
  unique(as.character(out))
}

safe_reductions <- function(obj) {
  out <- tryCatch(SeuratObject::Reductions(obj), error = function(e) names(obj@reductions))
  unique(as.character(out))
}

safe_key <- function(obj, red_name) {
  tryCatch(as.character(Key(obj[[red_name]])), error = function(e) make_key(red_name))
}

safe_loadings <- function(obj, red_name) {
  tryCatch(Loadings(obj[[red_name]]), error = function(e) NULL)
}

dump_rds_bundle <- function(input_file, out_dir, assay_name = "", allow_loss = FALSE) {
  rlog(paste0("reading: ", input_file))
  obj <- readRDS(input_file)

  info <- list(
    input = normalizePath(input_file, winslash = "/", mustWork = TRUE),
    source_class = class(obj),
    selected_assay = NULL,
    all_assays = NULL,
    project = NULL,
    x_source = NULL,
    idents_levels = character(),
    layers = list(),
    reductions = list(),
    notes = character()
  )

  if (inherits(obj, "Seurat")) {
    use_assay <- pick_assay(obj, assay_name)
    info$selected_assay <- use_assay
    info$all_assays <- tryCatch(as.character(SeuratObject::Assays(obj)), error = function(e) names(obj@assays))
    info$project <- tryCatch(obj@project.name, error = function(e) NULL)
    info$idents_levels <- tryCatch(as.character(levels(SeuratObject::Idents(obj))), error = function(e) character())

    other_assays <- setdiff(info$all_assays, use_assay)
    if (length(other_assays) > 0) {
      msg <- paste0("other assays not exported: ", paste(other_assays, collapse = ", "))
      if (!allow_loss) stop(msg)
      info$notes <- c(info$notes, msg)
    }

    SeuratObject::DefaultAssay(obj) <- use_assay

    obj <- tryCatch(
      JoinLayers(object = obj, assay = use_assay),
      error = function(e) obj
    )

    cells <- colnames(obj)
    features <- rownames(obj[[use_assay]])

    meta <- obj[[]]
    meta <- meta[cells, , drop = FALSE]
    if (!"seurat_ident" %in% colnames(meta)) {
      meta$seurat_ident <- as.character(SeuratObject::Idents(obj))
    }

    feat <- safe_feature_meta(obj, use_assay)
    feat <- feat[features, , drop = FALSE]

    write_lines_gz(cells, file.path(out_dir, "cells.tsv.gz"))
    write_lines_gz(features, file.path(out_dir, "features.tsv.gz"))
    write_table_gz(meta, file.path(out_dir, "obs.csv.gz"))
    write_table_gz(feat, file.path(out_dir, "var.csv.gz"))

    layer_names <- safe_layers(obj, use_assay)
    preferred <- unique(c("counts", "data", "scale.data", layer_names))

    rlog("exporting layers")
    for (layer_name in preferred) {
      mat <- tryCatch(
        SeuratObject::LayerData(object = obj, assay = use_assay, layer = layer_name),
        error = function(e) NULL
      )
      if (is.null(mat)) next
      if (!is.matrix(mat) && !inherits(mat, "Matrix")) next
      if (nrow(mat) == 0 || ncol(mat) == 0) next

      # AnnData.layers must have the same shape as X
      if (!identical(rownames(mat), features) || !identical(colnames(mat), cells)) {
        msg <- paste0(
          "layer '", layer_name,
          "' is not full feature x cell shape and cannot be stored in AnnData.layers losslessly"
        )
        if (!allow_loss) stop(msg)
        info$notes <- c(info$notes, msg)
        next
      }

      layer_file <- paste0("layer__", sanitize_name(layer_name), ".mtx")
      write_matrix_mtx(mat, file.path(out_dir, layer_file))
      info$layers[[layer_name]] <- layer_file
    }

    if (length(info$layers) == 0) {
      stop("no exportable full-size layers found")
    }

    if ("data" %in% names(info$layers)) {
      info$x_source <- "data"
    } else if ("counts" %in% names(info$layers)) {
      info$x_source <- "counts"
    } else {
      info$x_source <- names(info$layers)[[1]]
    }

    red_names <- safe_reductions(obj)
    if (length(red_names) > 0) rlog("exporting reductions")

    for (red in red_names) {
      emb <- tryCatch(Embeddings(object = obj, reduction = red), error = function(e) NULL)
      if (is.null(emb) || nrow(emb) == 0 || ncol(emb) == 0) next

      emb <- emb[cells, , drop = FALSE]
      emb_file <- paste0("obsm__", sanitize_name(red), ".csv.gz")
      write_table_gz(as.data.frame(emb), file.path(out_dir, emb_file))

      red_info <- list(
        embeddings = emb_file,
        key = safe_key(obj, red),
        loadings = NULL
      )

      load <- safe_loadings(obj, red)
      if (!is.null(load) && nrow(load) > 0 && ncol(load) > 0) {
        load_file <- paste0("varm__", sanitize_name(red), "_loadings.csv.gz")
        write_table_gz(as.data.frame(load), file.path(out_dir, load_file))
        red_info$loadings <- load_file
      }

      info$reductions[[red]] <- red_info
    }

  } else if (
    inherits(obj, "dgCMatrix") || inherits(obj, "dgTMatrix") ||
    inherits(obj, "matrix") || inherits(obj, "data.frame")
  ) {
    mat <- obj
    if (inherits(mat, "data.frame")) {
      mat <- as.matrix(mat)
    }

    if (is.null(rownames(mat))) rownames(mat) <- paste0("feature_", seq_len(nrow(mat)))
    if (is.null(colnames(mat))) colnames(mat) <- paste0("cell_", seq_len(ncol(mat)))

    cells <- colnames(mat)
    features <- rownames(mat)

    write_lines_gz(cells, file.path(out_dir, "cells.tsv.gz"))
    write_lines_gz(features, file.path(out_dir, "features.tsv.gz"))
    write_table_gz(data.frame(row.names = cells), file.path(out_dir, "obs.csv.gz"))
    write_table_gz(data.frame(row.names = features), file.path(out_dir, "var.csv.gz"))

    layer_file <- "layer__counts.mtx"
    write_matrix_mtx(mat, file.path(out_dir, layer_file))

    info$selected_assay <- ifelse(nzchar(assay_name), assay_name, "RNA")
    info$all_assays <- info$selected_assay
    info$x_source <- "counts"
    info$layers[["counts"]] <- layer_file
    info$notes <- c(info$notes, "input was not Seurat; exported as counts only")

  } else {
    stop(
      "unsupported RDS class: ", paste(class(obj), collapse = "/"),
      ". Supported: Seurat, matrix, data.frame, dgCMatrix, dgTMatrix"
    )
  }

  write(
    jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE, null = "null"),
    file = file.path(out_dir, "manifest.json")
  )
  invisible(out_dir)
}

inspect_rds <- function(input_file, assay_name = "") {
  obj <- readRDS(input_file)
  info <- list(class = class(obj))

  if (inherits(obj, "Seurat")) {
    use_assay <- pick_assay(obj, assay_name)
    info$dims <- unname(dim(obj))
    info$assays <- tryCatch(as.character(SeuratObject::Assays(obj)), error = function(e) names(obj@assays))
    info$default_assay <- tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) NULL)
    info$layers <- safe_layers(obj, use_assay)
    info$reductions <- safe_reductions(obj)
    info$obs_cols <- colnames(obj[[]])
  }

  jsonlite::toJSON(info, auto_unbox = TRUE, pretty = FALSE, null = "null")
}

build_rds_from_bundle <- function(bundle_dir, output_file, assay_name = "RNA", project = "converted") {
  manifest <- jsonlite::fromJSON(file.path(bundle_dir, "manifest.json"), simplifyVector = FALSE)

  cells <- read_lines_gz(file.path(bundle_dir, "cells.tsv.gz"))
  features <- read_lines_gz(file.path(bundle_dir, "features.tsv.gz"))

  counts_path <- NULL
  if (!is.null(manifest$layers$counts)) {
    counts_path <- file.path(bundle_dir, manifest$layers$counts)
  } else if (!is.null(manifest$x_source) && manifest$x_source %in% names(manifest$layers)) {
    counts_path <- file.path(bundle_dir, manifest$layers[[manifest$x_source]])
  } else if (length(manifest$layers) > 0) {
    counts_path <- file.path(bundle_dir, manifest$layers[[1]])
  }

  if (is.null(counts_path) || !file.exists(counts_path)) {
    stop("no usable matrix found in bundle for Seurat creation")
  }

  counts <- Matrix::readMM(counts_path)
  counts <- methods::as(counts, "dgCMatrix")
  rownames(counts) <- features
  colnames(counts) <- cells

  obs <- read_table_gz(file.path(bundle_dir, "obs.csv.gz"))
  obs <- obs[cells, , drop = FALSE]

  project_name <- project
  if (!is.null(manifest$project) && is.character(manifest$project) && nzchar(manifest$project)) {
    project_name <- manifest$project
  }

  obj <- CreateSeuratObject(
    counts = counts,
    assay = assay_name,
    meta.data = obs,
    project = project_name
  )

  var_path <- file.path(bundle_dir, "var.csv.gz")
  if (file.exists(var_path)) {
    var <- read_table_gz(var_path)
    var <- var[features, , drop = FALSE]
    if (ncol(var) > 0) {
      obj[[assay_name]] <- AddMetaData(object = obj[[assay_name]], metadata = var)
    }
  }

  if (length(manifest$layers) > 0) {
    for (layer_name in names(manifest$layers)) {
      if (layer_name == "counts") next

      layer_file <- file.path(bundle_dir, manifest$layers[[layer_name]])
      if (!file.exists(layer_file)) next

      mat <- Matrix::readMM(layer_file)
      mat <- methods::as(mat, "dgCMatrix")
      rownames(mat) <- features
      colnames(mat) <- cells

      SeuratObject::LayerData(object = obj, assay = assay_name, layer = layer_name) <- mat
    }
  }

  if (length(manifest$reductions) > 0) {
    for (red_name in names(manifest$reductions)) {
      red_info <- manifest$reductions[[red_name]]

      if (is.character(red_info)) {
        emb_file <- red_info
        load_file <- NULL
        red_key <- make_key(red_name)
      } else {
        emb_file <- red_info$embeddings
        load_file <- red_info$loadings
        red_key <- if (!is.null(red_info$key) && nzchar(red_info$key)) red_info$key else make_key(red_name)
      }

      emb_path <- file.path(bundle_dir, emb_file)
      if (!file.exists(emb_path)) next

      emb <- read_table_gz(emb_path)
      emb <- as.matrix(emb[cells, , drop = FALSE])

      loadings <- NULL
      if (!is.null(load_file) && is.character(load_file) && nzchar(load_file)) {
        load_path <- file.path(bundle_dir, load_file)
        if (file.exists(load_path)) {
          ldf <- read_table_gz(load_path)
          if (ncol(ldf) > 0) {
            full <- matrix(
              0,
              nrow = length(features),
              ncol = ncol(ldf),
              dimnames = list(features, colnames(ldf))
            )
            common <- intersect(features, rownames(ldf))
            if (length(common) > 0) {
              full[common, ] <- as.matrix(ldf[common, , drop = FALSE])
            }
            loadings <- full
          }
        }
      }

      if (is.null(loadings)) {
        dr <- CreateDimReducObject(
          embeddings = emb,
          assay = assay_name,
          key = red_key
        )
      } else {
        dr <- CreateDimReducObject(
          embeddings = emb,
          loadings = loadings,
          assay = assay_name,
          key = red_key
        )
      }

      obj[[red_name]] <- dr
    }
  }

  if ("seurat_ident" %in% colnames(obj[[]])) {
    lv <- tryCatch(unlist(manifest$idents_levels, use.names = FALSE), error = function(e) character())
    if (length(lv) > 0) {
      SeuratObject::Idents(obj) <- factor(obj$seurat_ident, levels = lv)
    } else {
      SeuratObject::Idents(obj) <- obj$seurat_ident
    }
  }

  obj@misc$rds_bridge <- list(
    source_class = manifest$source_class,
    selected_assay = manifest$selected_assay,
    x_source = manifest$x_source,
    notes = manifest$notes
  )

  uns_path <- file.path(bundle_dir, "uns.json.gz")
  if (file.exists(uns_path)) {
    uns_json <- paste(readLines(gzfile(uns_path), warn = FALSE), collapse = "\n")
    obj@misc$anndata_uns_json <- uns_json
  }

  saveRDS(obj, file = output_file)
  invisible(output_file)
}
'''

ro.r(R_HELPERS)
R_DUMP = ro.globalenv["dump_rds_bundle"]
R_BUILD = ro.globalenv["build_rds_from_bundle"]
R_INSPECT = ro.globalenv["inspect_rds"]


# =========================================================
# Python helpers
# =========================================================
def infer_output_path(input_path: Path, output_path: str | None = None) -> Path:
    input_path = input_path.resolve()
    ext = input_path.suffix.lower()
    if ext not in {".rds", ".h5ad"}:
        raise ValueError("Only .rds and .h5ad are supported")

    if output_path is None:
        return input_path.with_suffix(".h5ad" if ext == ".rds" else ".rds")

    out = Path(output_path)
    if out.suffix == "":
        out = out.with_suffix(".h5ad" if ext == ".rds" else ".rds")
    return out.resolve()


def sanitize_key(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def reduction_to_obsm_key(name: str) -> str:
    return name if name.startswith("X_") else f"X_{name}"


def obsm_to_reduction_name(name: str) -> str:
    return name[2:] if name.startswith("X_") else name


def write_lines_gz_py(items, path: Path) -> None:
    with gzip.open(path, "wt") as fh:
        for x in items:
            fh.write(f"{x}\n")


def read_lines_gz_py(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n") for line in fh]


def write_table_gz_py(df: pd.DataFrame, path: Path) -> None:
    out = df.copy()
    out.index = out.index.astype(str)
    out.insert(0, "__index__", out.index)
    out.to_csv(path, index=False, compression="gzip")


def read_table_gz_py(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, compression="gzip")
    if "__index__" not in df.columns:
        raise ValueError(f"Missing __index__ column: {path}")
    df = df.set_index("__index__")
    df.index = df.index.astype(str)
    return df


def read_mtx_transposed(mtx_path: Path):
    x = mmread(str(mtx_path))
    if sp.issparse(x):
        return x.T.tocsr()
    return np.asarray(x).T


def write_mtx_from_anndata(x, mtx_path: Path) -> None:
    if sp.issparse(x):
        mmwrite(str(mtx_path), x.T.tocoo())
    else:
        mmwrite(str(mtx_path), np.asarray(x).T)


def looks_like_counts(x, sample_n: int = 10000) -> bool:
    if sp.issparse(x):
        vals = x.data
    else:
        vals = np.asarray(x).ravel()

    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return True

    if vals.size > sample_n:
        rng = np.random.default_rng(1)
        vals = rng.choice(vals, size=sample_n, replace=False)

    return bool(np.all(vals >= 0) and np.mean(np.abs(vals - np.round(vals)) < 1e-8) > 0.98)


def to_jsonable(x):
    if x is None or isinstance(x, (str, int, float, bool)):
        return x
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, pd.Index):
        return [str(v) for v in x.tolist()]
    if isinstance(x, pd.Series):
        return to_jsonable(x.tolist())
    if isinstance(x, pd.DataFrame):
        return {
            "__dataframe__": True,
            "shape": list(x.shape),
            "columns": [str(c) for c in x.columns],
        }
    if isinstance(x, np.ndarray):
        return x.tolist()
    if sp.issparse(x):
        return {"__sparse__": True, "shape": list(x.shape), "nnz": int(x.nnz)}
    if isinstance(x, dict):
        return {str(k): to_jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [to_jsonable(v) for v in x]
    return repr(x)


def choose_x_name(manifest: dict) -> str:
    layers = list((manifest.get("layers") or {}).keys())
    if not layers:
        raise RuntimeError("No layers found in manifest")

    x_source = manifest.get("x_source")
    if isinstance(x_source, str) and x_source in layers:
        return x_source
    if "data" in layers:
        return "data"
    if "counts" in layers:
        return "counts"
    return layers[0]


def inspect_rds_json(rds_file: Path, assay_name: str) -> dict:
    txt = str(R_INSPECT(str(rds_file), assay_name or "")[0])
    return json.loads(txt)


# =========================================================
# Validation
# =========================================================
def validate_rds_to_h5ad(output_file: Path, manifest: dict, cells: list[str], features: list[str]) -> None:
    adata = ad.read_h5ad(str(output_file))

    errors = []
    if adata.n_obs != len(cells) or adata.n_vars != len(features):
        errors.append("dimension mismatch")
    if list(map(str, adata.obs_names.tolist())) != list(map(str, cells)):
        errors.append("obs_names mismatch")
    if list(map(str, adata.var_names.tolist())) != list(map(str, features)):
        errors.append("var_names mismatch")

    expected_layers = set((manifest.get("layers") or {}).keys())
    missing_layers = sorted(expected_layers - set(adata.layers.keys()))
    if missing_layers:
        errors.append(f"missing layers: {missing_layers}")

    expected_obsm = {reduction_to_obsm_key(k) for k in (manifest.get("reductions") or {}).keys()}
    missing_obsm = sorted(expected_obsm - set(adata.obsm.keys()))
    if missing_obsm:
        errors.append(f"missing obsm: {missing_obsm}")

    if errors:
        raise RuntimeError("validation failed: " + "; ".join(errors))

    log(
        "validate ok | "
        f"cells={adata.n_obs} genes={adata.n_vars} | "
        f"layers={list(adata.layers.keys())} | "
        f"obsm={list(adata.obsm.keys())}"
    )


def validate_h5ad_to_rds(output_file: Path, manifest: dict, cells: list[str], features: list[str], assay_name: str) -> None:
    info = inspect_rds_json(output_file, assay_name)

    errors = []
    dims = info.get("dims") or []
    if list(map(int, dims)) != [len(features), len(cells)]:
        errors.append("dimension mismatch")

    expected_layers = set((manifest.get("layers") or {}).keys())
    have_layers = set(info.get("layers") or [])
    missing_layers = sorted(expected_layers - have_layers)
    if missing_layers:
        errors.append(f"missing layers: {missing_layers}")

    expected_reductions = set((manifest.get("reductions") or {}).keys())
    have_reductions = set(info.get("reductions") or [])
    missing_reductions = sorted(expected_reductions - have_reductions)
    if missing_reductions:
        errors.append(f"missing reductions: {missing_reductions}")

    if errors:
        raise RuntimeError("validation failed: " + "; ".join(errors))

    log(
        "validate ok | "
        f"features={dims[0]} cells={dims[1]} | "
        f"layers={info.get('layers', [])} | "
        f"reductions={info.get('reductions', [])}"
    )


# =========================================================
# Conversion directions
# =========================================================
def rds_to_h5ad(
    input_file: Path,
    output_file: Path,
    assay_name: str | None,
    compression: str | None = "gzip",
    allow_loss: bool = False,
    validate: bool = True,
) -> None:
    with tempfile.TemporaryDirectory(prefix="rds_to_h5ad_") as tmpdir:
        tmpdir_p = Path(tmpdir)
        assay_arg = "" if assay_name is None else assay_name

        log("1/5 export bundle from embedded R")
        R_DUMP(str(input_file), str(tmpdir_p), assay_arg, allow_loss)

        log("2/5 read exported bundle")
        manifest = json.loads((tmpdir_p / "manifest.json").read_text())
        cells = read_lines_gz_py(tmpdir_p / "cells.tsv.gz")
        features = read_lines_gz_py(tmpdir_p / "features.tsv.gz")

        obs = read_table_gz_py(tmpdir_p / "obs.csv.gz").reindex(cells)
        var = read_table_gz_py(tmpdir_p / "var.csv.gz").reindex(features)

        layer_map = {}
        for layer_name, rel_path in (manifest.get("layers") or {}).items():
            layer_map[layer_name] = read_mtx_transposed(tmpdir_p / rel_path)

        x_name = choose_x_name(manifest)

        log(f"3/5 build AnnData (X={x_name})")
        adata = ad.AnnData(X=layer_map[x_name], obs=obs, var=var)

        # keep all original layers explicitly, including x_source
        for layer_name, mat in layer_map.items():
            adata.layers[layer_name] = mat

        for red_name, red_info in (manifest.get("reductions") or {}).items():
            if isinstance(red_info, str):
                emb_file = red_info
                load_file = None
            else:
                emb_file = red_info.get("embeddings")
                load_file = red_info.get("loadings")

            if emb_file:
                emb = read_table_gz_py(tmpdir_p / emb_file).reindex(cells)
                adata.obsm[reduction_to_obsm_key(red_name)] = emb.to_numpy()

            if load_file:
                load = read_table_gz_py(tmpdir_p / load_file).reindex(features)
                if load.shape[1] > 0:
                    adata.varm[f"{red_name}_loadings"] = load.fillna(0.0).to_numpy()

        adata.uns["rds_bridge"] = {
            "source_class": manifest.get("source_class", []),
            "selected_assay": manifest.get("selected_assay"),
            "all_assays": manifest.get("all_assays", []),
            "project": manifest.get("project"),
            "x_source": x_name,
            "idents_levels": manifest.get("idents_levels", []),
            "reductions": manifest.get("reductions", {}),
            "notes": manifest.get("notes", []),
        }

        log("4/5 write h5ad")
        adata.write_h5ad(str(output_file), compression=compression)

        if validate:
            log("5/5 validate output")
            validate_rds_to_h5ad(output_file, manifest, cells, features)


def h5ad_to_rds(
    input_file: Path,
    output_file: Path,
    assay_name: str,
    project: str,
    allow_loss: bool = False,
    validate: bool = True,
) -> None:
    log("1/5 read h5ad")
    adata = ad.read_h5ad(str(input_file))

    bridge = adata.uns.get("rds_bridge", {})
    if not isinstance(bridge, dict):
        bridge = {}

    x_source = bridge.get("x_source")
    if not isinstance(x_source, str) or not x_source:
        if "data" in adata.layers:
            x_source = "data"
        elif "counts" in adata.layers:
            x_source = "counts"
        else:
            x_source = "counts" if looks_like_counts(adata.X) else "data"

    with tempfile.TemporaryDirectory(prefix="h5ad_to_rds_") as tmpdir:
        tmpdir_p = Path(tmpdir)

        log("2/5 write temporary bundle")
        cells = [str(x) for x in adata.obs_names]
        features = [str(x) for x in adata.var_names]

        write_lines_gz_py(cells, tmpdir_p / "cells.tsv.gz")
        write_lines_gz_py(features, tmpdir_p / "features.tsv.gz")

        obs = adata.obs.copy()
        var = adata.var.copy()
        obs.index = obs.index.astype(str)
        var.index = var.index.astype(str)
        write_table_gz_py(obs, tmpdir_p / "obs.csv.gz")
        write_table_gz_py(var, tmpdir_p / "var.csv.gz")

        manifest = {
            "input": str(input_file),
            "source_class": ["AnnData"],
            "selected_assay": assay_name,
            "all_assays": [assay_name],
            "project": bridge.get("project", project),
            "x_source": x_source,
            "idents_levels": bridge.get("idents_levels", []),
            "layers": {},
            "reductions": {},
            "notes": [],
        }

        # write all named layers
        for layer_name in adata.layers.keys():
            file_name = f"layer__{sanitize_key(layer_name)}.mtx"
            write_mtx_from_anndata(adata.layers[layer_name], tmpdir_p / file_name)
            manifest["layers"][layer_name] = file_name

        # ensure X also has a named layer for round-trip
        if x_source not in manifest["layers"]:
            file_name = f"layer__{sanitize_key(x_source)}.mtx"
            write_mtx_from_anndata(adata.X, tmpdir_p / file_name)
            manifest["layers"][x_source] = file_name

        # if no counts layer, try to recover one from raw or X
        if "counts" not in manifest["layers"]:
            if adata.raw is not None:
                file_name = "layer__counts.mtx"
                write_mtx_from_anndata(adata.raw.X, tmpdir_p / file_name)
                manifest["layers"]["counts"] = file_name
                manifest["notes"].append("counts recovered from AnnData.raw.X")
            elif looks_like_counts(adata.X):
                file_name = "layer__counts.mtx"
                write_mtx_from_anndata(adata.X, tmpdir_p / file_name)
                manifest["layers"]["counts"] = file_name
                manifest["notes"].append("counts recovered from AnnData.X")
            else:
                msg = "no counts-like matrix found; Seurat counts will fall back to X"
                if not allow_loss:
                    raise RuntimeError(msg)
                manifest["notes"].append(msg)

        bridge_reductions = bridge.get("reductions", {}) if isinstance(bridge, dict) else {}

        for key in adata.obsm.keys():
            red_name = obsm_to_reduction_name(key)
            emb = pd.DataFrame(np.asarray(adata.obsm[key]), index=adata.obs_names)
            emb_file = f"obsm__{sanitize_key(red_name)}.csv.gz"
            write_table_gz_py(emb, tmpdir_p / emb_file)

            red_meta = bridge_reductions.get(red_name, {}) if isinstance(bridge_reductions, dict) else {}
            red_info = {
                "embeddings": emb_file,
                "key": red_meta.get("key", None),
                "loadings": None,
            }

            varm_key = f"{red_name}_loadings"
            if varm_key in adata.varm.keys():
                load = pd.DataFrame(np.asarray(adata.varm[varm_key]), index=adata.var_names)
                load_file = f"varm__{sanitize_key(red_name)}_loadings.csv.gz"
                write_table_gz_py(load, tmpdir_p / load_file)
                red_info["loadings"] = load_file

            manifest["reductions"][red_name] = red_info

        with gzip.open(tmpdir_p / "uns.json.gz", "wt") as fh:
            fh.write(json.dumps(to_jsonable(adata.uns), ensure_ascii=False, indent=2))

        (tmpdir_p / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False))

        log("3/5 build Seurat in embedded R")
        R_BUILD(str(tmpdir_p), str(output_file), assay_name, str(manifest["project"]))

        if validate:
            log("4/5 validate output")
            validate_h5ad_to_rds(output_file, manifest, cells, features, assay_name)

        log("5/5 done")


# =========================================================
# CLI
# =========================================================
def main():
    parser = argparse.ArgumentParser(
        description="Convert .rds <-> .h5ad with Python as primary runtime and embedded R."
    )
    parser.add_argument("input_file", help="Input .rds or .h5ad")
    parser.add_argument("-o", "--output-file", default=None)
    parser.add_argument("--assay-name", default="RNA")
    parser.add_argument("--project", default="converted")
    parser.add_argument("--compression", default="gzip", choices=["gzip", "lzf", "none"])
    parser.add_argument("--allow-loss", action="store_true", help="Allow unsupported structures to be skipped with notes.")
    parser.add_argument("--no-validate", action="store_true", help="Skip post-write validation.")
    args = parser.parse_args()

    input_file = Path(args.input_file).resolve()
    output_file = infer_output_path(input_file, args.output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    validate = not args.no_validate

    if input_file.suffix.lower() == ".rds":
        compression = None if args.compression == "none" else args.compression
        rds_to_h5ad(
            input_file=input_file,
            output_file=output_file,
            assay_name=args.assay_name,
            compression=compression,
            allow_loss=args.allow_loss,
            validate=validate,
        )
    elif input_file.suffix.lower() == ".h5ad":
        h5ad_to_rds(
            input_file=input_file,
            output_file=output_file,
            assay_name=args.assay_name,
            project=args.project,
            allow_loss=args.allow_loss,
            validate=validate,
        )
    else:
        raise ValueError("Only .rds and .h5ad are supported")

    print(json.dumps({
        "input": str(input_file),
        "output": str(output_file),
        "status": "ok",
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
