#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rds_to_h5ad.py

Tolerant Seurat/RDS -> AnnData/H5AD converter.

Design:
- Core data must succeed: cells, features, obs, var, and at least one usable matrix layer.
- Optional data are attempted one by one: layers, reductions, loadings.
- Optional parts that cannot be converted are skipped instead of killing the whole job.
- A terminal summary and a sidecar JSON report are written, listing converted and skipped parts.

Examples:
  python rds_to_h5ad.py input.rds --assay-name RNA
  python rds_to_h5ad.py input.rds -o output.h5ad --assay-name RNA
  python rds_to_h5ad.py input.rds --assay-name RNA --strict
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
from rpy2 import robjects as ro


SCRIPT_VERSION = "2026-05-04-tolerant-v3"


def log(msg: str) -> None:
    print(f"[rds_to_h5ad] {msg}", flush=True)


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
  if (!is.null(assay_name) && nzchar(assay_name) && assay_name %in% assays) return(assay_name)
  if ("RNA" %in% assays) return("RNA")
  da <- tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) NULL)
  if (!is.null(da) && nzchar(da)) return(da)
  assays[[1]]
}

append_item <- function(x, value) {
  c(x, as.character(value))
}

write_lines_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(as.character(x), con, sep = "\n")
}

write_table_gz <- function(df, path) {
  if (is.null(df)) df <- data.frame()
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  idx <- rownames(df)
  if (is.null(idx)) idx <- as.character(seq_len(nrow(df)))

  out <- data.frame(`__index__` = idx, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) > 0) out <- cbind(out, df)

  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.csv(out, con, row.names = FALSE, quote = TRUE, na = "")
}

write_matrix_mtx <- function(mat, path) {
  if (inherits(mat, "Matrix")) {
    mat <- methods::as(mat, "dgCMatrix")
  } else {
    mat <- Matrix::Matrix(mat, sparse = TRUE)
  }
  Matrix::writeMM(mat, path)
}

safe_feature_meta <- function(obj, assay_name) {
  features <- tryCatch(rownames(obj[[assay_name]]), error = function(e) NULL)
  out <- tryCatch(obj[[assay_name]][[]], error = function(e) NULL)
  if (is.null(out)) out <- data.frame(row.names = features)
  out
}

safe_layers <- function(obj, assay_name) {
  out <- tryCatch(SeuratObject::Layers(obj[[assay_name]]), error = function(e) character(0))
  out <- unique(as.character(out))
  if (length(out) == 0) {
    candidates <- c("counts", "data", "scale.data")
    have <- character(0)
    for (ly in candidates) {
      m <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay_name, slot = ly), error = function(e) NULL)
      if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) have <- c(have, ly)
    }
    out <- have
  }
  unique(out)
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

get_layer_matrix <- function(obj, assay_name, layer_name) {
  m <- tryCatch(
    SeuratObject::LayerData(object = obj, assay = assay_name, layer = layer_name, fast = FALSE),
    error = function(e) NULL
  )
  if (!is.null(m)) return(m)

  m <- tryCatch(
    SeuratObject::GetAssayData(object = obj, assay = assay_name, layer = layer_name),
    error = function(e) NULL
  )
  if (!is.null(m)) return(m)

  tryCatch(
    SeuratObject::GetAssayData(object = obj, assay = assay_name, slot = layer_name),
    error = function(e) NULL
  )
}

standardize_layer <- function(mat, features, cells, layer_name, strict = FALSE) {
  if (is.null(mat)) return(list(mat = NULL, notes = character(), reason = "not found"))
  if (inherits(mat, "data.frame")) mat <- as.matrix(mat)
  if (!is.matrix(mat) && !inherits(mat, "Matrix")) {
    return(list(mat = NULL, notes = character(), reason = "not matrix-like"))
  }
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(list(mat = NULL, notes = character(), reason = "empty matrix"))
  }

  notes <- character()
  rn <- rownames(mat)
  cn <- colnames(mat)
  has_rn <- !is.null(rn) && length(rn) == nrow(mat) && !any(is.na(rn))
  has_cn <- !is.null(cn) && length(cn) == ncol(mat) && !any(is.na(cn))
  orientation <- NA_character_

  if (has_rn && has_cn) {
    if (all(rn %in% features) && all(cn %in% cells)) {
      orientation <- "feature_by_cell"
    } else if (all(rn %in% cells) && all(cn %in% features)) {
      mat <- Matrix::t(mat)
      rn <- rownames(mat); cn <- colnames(mat)
      orientation <- "feature_by_cell"
      notes <- c(notes, paste0("layer '", layer_name, "' was transposed to feature x cell"))
    }
  }

  if (is.na(orientation)) {
    if (nrow(mat) == length(features) && ncol(mat) == length(cells)) {
      rownames(mat) <- features; colnames(mat) <- cells
      rn <- features; cn <- cells
      orientation <- "feature_by_cell"
      notes <- c(notes, paste0("layer '", layer_name, "' had missing/unusable dimnames; assigned assay feature/cell names"))
    } else if (nrow(mat) == length(cells) && ncol(mat) == length(features)) {
      mat <- Matrix::t(mat)
      rownames(mat) <- features; colnames(mat) <- cells
      rn <- features; cn <- cells
      orientation <- "feature_by_cell"
      notes <- c(notes, paste0("layer '", layer_name, "' looked cell x feature; transposed and assigned dimnames"))
    }
  }

  if (is.na(orientation)) {
    msg <- paste0(
      "cannot align to selected assay features/cells; matrix dim=", paste(dim(mat), collapse = "x"),
      ", target dim=", length(features), "x", length(cells)
    )
    if (strict) stop(paste0("layer '", layer_name, "' ", msg))
    return(list(mat = NULL, notes = character(), reason = msg))
  }

  rn <- rownames(mat); cn <- colnames(mat)
  if (anyDuplicated(rn) || anyDuplicated(cn)) {
    msg <- "duplicated rownames or colnames"
    if (strict) stop(paste0("layer '", layer_name, "' has ", msg))
    return(list(mat = NULL, notes = character(), reason = msg))
  }
  if (!all(rn %in% features) || !all(cn %in% cells)) {
    msg <- "contains names outside selected assay features/cells"
    if (strict) stop(paste0("layer '", layer_name, "' ", msg))
    return(list(mat = NULL, notes = character(), reason = msg))
  }

  if (identical(rn, features) && identical(cn, cells)) {
    return(list(mat = mat, notes = notes, reason = NULL))
  }

  if (setequal(rn, features) && setequal(cn, cells)) {
    mat <- mat[features, cells, drop = FALSE]
    notes <- c(notes, paste0("layer '", layer_name, "' reordered to assay feature/cell order"))
    return(list(mat = mat, notes = notes, reason = NULL))
  }

  msg <- paste0("partial layer: features ", length(rn), "/", length(features), ", cells ", length(cn), "/", length(cells))
  if (strict) stop(paste0("layer '", layer_name, "' is ", msg, "; AnnData.layers requires full shape"))

  out <- Matrix::Matrix(0, nrow = length(features), ncol = length(cells), sparse = TRUE)
  rownames(out) <- features
  colnames(out) <- cells
  out[match(rn, features), match(cn, cells)] <- mat
  notes <- c(notes, paste0("layer '", layer_name, "' is ", msg, "; zero-padded for AnnData.layers compatibility"))
  list(mat = out, notes = notes, reason = NULL)
}

dump_rds_bundle <- function(input_file, out_dir, assay_name = "", strict = FALSE) {
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
    notes = character(),
    strict = strict,
    converted = list(layers = character(), reductions = character(), loadings = character()),
    skipped = list(layers = character(), reductions = character(), loadings = character(), assays = character())
  )

  if (inherits(obj, "Seurat")) {
    use_assay <- pick_assay(obj, assay_name)
    info$selected_assay <- use_assay
    info$all_assays <- tryCatch(as.character(SeuratObject::Assays(obj)), error = function(e) names(obj@assays))
    info$project <- tryCatch(obj@project.name, error = function(e) NULL)
    info$idents_levels <- tryCatch(as.character(levels(SeuratObject::Idents(obj))), error = function(e) character())

    other_assays <- setdiff(info$all_assays, use_assay)
    if (length(other_assays) > 0) {
      msg <- paste0("other assays not exported because this script exports one assay per run: ", paste(other_assays, collapse = ", "))
      info$notes <- c(info$notes, msg)
      info$skipped$assays <- c(info$skipped$assays, other_assays)
    }

    SeuratObject::DefaultAssay(obj) <- use_assay
    obj <- tryCatch(JoinLayers(object = obj, assay = use_assay), error = function(e) obj)

    cells <- colnames(obj)
    features <- rownames(obj[[use_assay]])
    if (is.null(cells) || is.null(features)) stop("could not determine cells/features from selected assay")

    meta <- obj[[]]
    meta <- meta[cells, , drop = FALSE]
    if (!"seurat_ident" %in% colnames(meta)) meta$seurat_ident <- as.character(SeuratObject::Idents(obj))

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
      mat0 <- get_layer_matrix(obj, use_assay, layer_name)
      std <- standardize_layer(mat0, features = features, cells = cells, layer_name = layer_name, strict = strict)
      if (length(std$notes) > 0) info$notes <- c(info$notes, std$notes)
      if (is.null(std$mat)) {
        reason <- ifelse(is.null(std$reason), "unknown reason", std$reason)
        info$skipped$layers <- c(info$skipped$layers, paste0(layer_name, ": ", reason))
        rlog(paste0("  skip layer ", layer_name, ": ", reason))
        next
      }

      layer_file <- paste0("layer__", sanitize_name(layer_name), ".mtx")
      ok <- tryCatch({
        write_matrix_mtx(std$mat, file.path(out_dir, layer_file))
        TRUE
      }, error = function(e) {
        info$skipped$layers <<- c(info$skipped$layers, paste0(layer_name, ": writeMM failed: ", conditionMessage(e)))
        FALSE
      })
      if (!ok) next

      info$layers[[layer_name]] <- layer_file
      info$converted$layers <- c(info$converted$layers, layer_name)
      rlog(paste0("  layer ", layer_name, " -> ", layer_file, " dim=", paste(dim(std$mat), collapse = "x")))
    }

    if (length(info$layers) == 0) stop("no exportable matrix layers found; cannot build AnnData")
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
      if (is.null(emb) || nrow(emb) == 0 || ncol(emb) == 0) {
        info$skipped$reductions <- c(info$skipped$reductions, paste0(red, ": empty or unreadable embeddings"))
        next
      }
      emb <- as.data.frame(emb)
      if (is.null(rownames(emb))) {
        info$skipped$reductions <- c(info$skipped$reductions, paste0(red, ": embeddings have no cell names"))
        next
      }
      common <- intersect(cells, rownames(emb))
      if (length(common) == 0) {
        info$skipped$reductions <- c(info$skipped$reductions, paste0(red, ": embeddings do not match any selected cells"))
        next
      }

      full_emb <- matrix(NA_real_, nrow = length(cells), ncol = ncol(emb), dimnames = list(cells, colnames(emb)))
      full_emb[common, ] <- as.matrix(emb[common, , drop = FALSE])
      if (length(common) < length(cells)) {
        info$notes <- c(info$notes, paste0("reduction '", red, "' missing ", length(cells) - length(common), " cells; filled missing embeddings with NA"))
      }

      emb_file <- paste0("obsm__", sanitize_name(red), ".csv.gz")
      ok <- tryCatch({
        write_table_gz(as.data.frame(full_emb), file.path(out_dir, emb_file))
        TRUE
      }, error = function(e) {
        info$skipped$reductions <<- c(info$skipped$reductions, paste0(red, ": writing embeddings failed: ", conditionMessage(e)))
        FALSE
      })
      if (!ok) next

      red_info <- list(embeddings = emb_file, key = safe_key(obj, red))
      info$converted$reductions <- c(info$converted$reductions, red)

      load <- safe_loadings(obj, red)
      if (!is.null(load) && nrow(load) > 0 && ncol(load) > 0) {
        load_file <- paste0("varm__", sanitize_name(red), "_loadings.csv.gz")
        ok_load <- tryCatch({
          write_table_gz(as.data.frame(load), file.path(out_dir, load_file))
          TRUE
        }, error = function(e) {
          info$skipped$loadings <<- c(info$skipped$loadings, paste0(red, ": writing loadings failed: ", conditionMessage(e)))
          FALSE
        })
        if (ok_load) {
          red_info$loadings <- load_file
          info$converted$loadings <- c(info$converted$loadings, red)
        }
      } else {
        info$skipped$loadings <- c(info$skipped$loadings, paste0(red, ": no loadings available"))
      }

      info$reductions[[red]] <- red_info
    }

  } else if (
    inherits(obj, "dgCMatrix") || inherits(obj, "dgTMatrix") ||
    inherits(obj, "matrix") || inherits(obj, "data.frame")
  ) {
    mat <- obj
    if (inherits(mat, "data.frame")) mat <- as.matrix(mat)
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
    info$converted$layers <- c(info$converted$layers, "counts")
    info$notes <- c(info$notes, "input was not Seurat; exported as counts only")
  } else {
    stop("unsupported RDS class: ", paste(class(obj), collapse = "/"), ". Supported: Seurat, matrix, data.frame, dgCMatrix, dgTMatrix")
  }

  write(jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE, null = "null"), file = file.path(out_dir, "manifest.json"))
  invisible(out_dir)
}
'''

ro.r(R_HELPERS)
R_DUMP = ro.globalenv["dump_rds_bundle"]


def infer_output_path(input_path: Path, output_path: str | None = None) -> Path:
    input_path = input_path.resolve()
    if input_path.suffix.lower() != ".rds":
        raise ValueError("Input must be .rds for rds_to_h5ad.py")
    if output_path is None:
        return input_path.with_suffix(".h5ad")
    out = Path(output_path)
    if out.suffix == "":
        out = out.with_suffix(".h5ad")
    return out.resolve()


def report_path_for(output_file: Path) -> Path:
    return Path(str(output_file) + ".conversion_report.json")


def reduction_to_obsm_key(name: str) -> str:
    return name if name.startswith("X_") else f"X_{name}"


def read_lines_gz_py(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n") for line in fh]


def read_table_gz_py(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, compression="gzip", low_memory=False)
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


def choose_x_name(manifest: dict, available_layers: dict[str, Any]) -> str:
    layers = list(available_layers.keys())
    if not layers:
        raise RuntimeError("No readable layers found after export; cannot build AnnData")
    x_source = manifest.get("x_source")
    if isinstance(x_source, str) and x_source in available_layers:
        return x_source
    if "data" in available_layers:
        return "data"
    if "counts" in available_layers:
        return "counts"
    return layers[0]


def drop_none_for_h5ad_uns(x):
    if x is None:
        return None
    if isinstance(x, dict):
        out = {}
        for k, v in x.items():
            vv = drop_none_for_h5ad_uns(v)
            if vv is not None:
                out[str(k)] = vv
        return out
    if isinstance(x, (list, tuple)):
        out = []
        for v in x:
            vv = drop_none_for_h5ad_uns(v)
            if vv is not None:
                out.append(vv)
        return out
    if isinstance(x, np.generic):
        return x.item()
    return x


def init_report(input_file: Path, output_file: Path, direction: str) -> dict[str, Any]:
    return {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "direction": direction,
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "input": str(input_file),
        "output": str(output_file),
        "status": "running",
        "converted": {"layers": [], "obsm": [], "varm": [], "obs": False, "var": False},
        "skipped": {"layers": [], "reductions": [], "loadings": [], "assays": [], "other": []},
        "notes": [],
    }


def write_report(report: dict[str, Any], output_file: Path) -> Path:
    report["finished_at"] = datetime.now().isoformat(timespec="seconds")
    path = report_path_for(output_file)
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


def validate_rds_to_h5ad(output_file: Path, cells: list[str], features: list[str]) -> None:
    adata = ad.read_h5ad(str(output_file))
    errors: list[str] = []
    if adata.n_obs != len(cells) or adata.n_vars != len(features):
        errors.append("dimension mismatch")
    if list(map(str, adata.obs_names.tolist())) != list(map(str, cells)):
        errors.append("obs_names mismatch")
    if list(map(str, adata.var_names.tolist())) != list(map(str, features)):
        errors.append("var_names mismatch")
    if errors:
        raise RuntimeError("validation failed: " + "; ".join(errors))
    log(f"validate ok | cells={adata.n_obs} genes={adata.n_vars} | layers={list(adata.layers.keys())} | obsm={list(adata.obsm.keys())}")


def rds_to_h5ad(
    input_file: Path,
    output_file: Path,
    assay_name: str | None,
    compression: str | None = "gzip",
    strict: bool = False,
    validate: bool = True,
) -> dict[str, Any]:
    report = init_report(input_file, output_file, "rds_to_h5ad")

    with tempfile.TemporaryDirectory(prefix="rds_to_h5ad_") as tmpdir:
        tmpdir_p = Path(tmpdir)
        assay_arg = "" if assay_name is None else assay_name

        log("1/5 export bundle from embedded R")
        R_DUMP(str(input_file), str(tmpdir_p), assay_arg, strict)

        log("2/5 read exported bundle")
        manifest = json.loads((tmpdir_p / "manifest.json").read_text())
        cells = read_lines_gz_py(tmpdir_p / "cells.tsv.gz")
        features = read_lines_gz_py(tmpdir_p / "features.tsv.gz")

        report["selected_assay"] = manifest.get("selected_assay")
        report["all_assays"] = manifest.get("all_assays")
        report["notes"].extend(manifest.get("notes") or [])
        skipped_manifest = manifest.get("skipped") or {}
        for k, vals in skipped_manifest.items():
            if k in report["skipped"] and vals:
                if isinstance(vals, list):
                    report["skipped"][k].extend(map(str, vals))
                else:
                    report["skipped"][k].append(str(vals))

        obs = read_table_gz_py(tmpdir_p / "obs.csv.gz").reindex(cells)
        var = read_table_gz_py(tmpdir_p / "var.csv.gz").reindex(features)
        report["converted"]["obs"] = True
        report["converted"]["var"] = True

        layer_map: dict[str, Any] = {}
        for layer_name, rel_path in (manifest.get("layers") or {}).items():
            try:
                mat = read_mtx_transposed(tmpdir_p / rel_path)
                if mat.shape != (len(cells), len(features)):
                    report["skipped"]["layers"].append(f"{layer_name}: read shape {mat.shape} != expected {(len(cells), len(features))}")
                    continue
                layer_map[layer_name] = mat
                report["converted"]["layers"].append(layer_name)
            except Exception as e:
                report["skipped"]["layers"].append(f"{layer_name}: Python read failed: {e}")

        x_name = choose_x_name(manifest, layer_map)
        report["x_source"] = x_name

        log(f"3/5 build AnnData (X={x_name})")
        adata = ad.AnnData(X=layer_map[x_name], obs=obs, var=var)
        for layer_name, mat in layer_map.items():
            try:
                adata.layers[layer_name] = mat
            except Exception as e:
                report["skipped"]["layers"].append(f"{layer_name}: could not assign to adata.layers: {e}")
                if layer_name in report["converted"]["layers"]:
                    report["converted"]["layers"].remove(layer_name)

        successful_reductions: dict[str, Any] = {}
        for red_name, red_info in (manifest.get("reductions") or {}).items():
            try:
                emb_file = red_info if isinstance(red_info, str) else red_info.get("embeddings")
                if not emb_file:
                    report["skipped"]["reductions"].append(f"{red_name}: no embeddings file")
                    continue
                emb = read_table_gz_py(tmpdir_p / emb_file).reindex(cells)
                emb_arr = emb.to_numpy(dtype=float)
                if emb_arr.shape[0] != len(cells) or emb_arr.ndim != 2 or emb_arr.shape[1] == 0:
                    report["skipped"]["reductions"].append(f"{red_name}: invalid embedding shape {emb_arr.shape}")
                    continue
                if np.all(~np.isfinite(emb_arr)):
                    report["skipped"]["reductions"].append(f"{red_name}: all embedding values are NA/NaN")
                    continue
                obsm_key = reduction_to_obsm_key(red_name)
                adata.obsm[obsm_key] = emb_arr
                report["converted"]["obsm"].append(obsm_key)
                keep_info = {"embeddings": emb_file}
                if isinstance(red_info, dict):
                    if red_info.get("key"):
                        keep_info["key"] = red_info.get("key")
                successful_reductions[red_name] = keep_info
            except Exception as e:
                report["skipped"]["reductions"].append(f"{red_name}: {e}")
                continue

            try:
                load_file = None if isinstance(red_info, str) else red_info.get("loadings")
                if not load_file:
                    continue
                load = read_table_gz_py(tmpdir_p / load_file).reindex(features)
                if load.shape[1] == 0:
                    report["skipped"]["loadings"].append(f"{red_name}: empty loadings table")
                    continue
                arr = load.fillna(0.0).to_numpy(dtype=float)
                adata.varm[f"{red_name}_loadings"] = arr
                report["converted"]["varm"].append(f"{red_name}_loadings")
                successful_reductions[red_name]["loadings"] = load_file
            except Exception as e:
                report["skipped"]["loadings"].append(f"{red_name}: {e}")

        bridge = {
            "source_class": manifest.get("source_class", []),
            "selected_assay": manifest.get("selected_assay"),
            "all_assays": manifest.get("all_assays", []),
            "project": manifest.get("project"),
            "x_source": x_name,
            "idents_levels": manifest.get("idents_levels", []),
            "reductions": successful_reductions,
            "notes": report["notes"],
            "skipped": report["skipped"],
            "script_version": SCRIPT_VERSION,
        }
        adata.uns["rds_bridge"] = drop_none_for_h5ad_uns(bridge)

        log("4/5 write h5ad")
        adata.write_h5ad(str(output_file), compression=compression)

        if validate:
            log("5/5 validate output")
            validate_rds_to_h5ad(output_file, cells, features)

    report["status"] = "ok"
    report["n_obs"] = len(cells)
    report["n_vars"] = len(features)
    report_path = write_report(report, output_file)
    print_summary(report, report_path)
    return report


def print_summary(report: dict[str, Any], report_path: Path) -> None:
    log("SUCCESS: h5ad written")
    log(f"  output: {report['output']}")
    log(f"  report: {report_path}")
    log(f"  cells={report.get('n_obs')} genes={report.get('n_vars')} X={report.get('x_source')}")
    log("converted:")
    log(f"  layers: {report['converted'].get('layers') or []}")
    log(f"  obsm:   {report['converted'].get('obsm') or []}")
    log(f"  varm:   {report['converted'].get('varm') or []}")
    skipped_any = any(report["skipped"].get(k) for k in report["skipped"])
    if skipped_any:
        log("skipped / unavailable parts:")
        for k, vals in report["skipped"].items():
            if vals:
                log(f"  {k}:")
                for v in vals:
                    log(f"    - {v}")
    else:
        log("skipped / unavailable parts: none")


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert Seurat/RDS to h5ad with embedded R via rpy2. Optional parts that fail are skipped and reported.")
    parser.add_argument("input_file", help="Input .rds file")
    parser.add_argument("-o", "--output-file", default=None, help="Output .h5ad path")
    parser.add_argument("--assay-name", default="RNA", help="Seurat assay to export, default: RNA")
    parser.add_argument("--compression", default="gzip", choices=["gzip", "lzf", "none"], help="h5ad compression")
    parser.add_argument("--strict", action="store_true", help="Fail on partial/misaligned layers instead of reordering/padding")
    parser.add_argument("--no-validate", action="store_true", help="Skip post-write validation")
    args = parser.parse_args()

    input_file = Path(args.input_file).resolve()
    if not input_file.exists():
        raise FileNotFoundError(input_file)
    output_file = infer_output_path(input_file, args.output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    compression = None if args.compression == "none" else args.compression
    report = rds_to_h5ad(
        input_file=input_file,
        output_file=output_file,
        assay_name=args.assay_name,
        compression=compression,
        strict=args.strict,
        validate=not args.no_validate,
    )

    print(json.dumps({"input": str(input_file), "output": str(output_file), "status": report["status"], "report": str(report_path_for(output_file))}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
