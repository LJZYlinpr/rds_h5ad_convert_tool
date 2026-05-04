#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
h5ad_to_rds.py

Tolerant AnnData/H5AD -> Seurat/RDS converter.

Design:
- Core object must succeed: cells, features, metadata, and one matrix usable for Seurat creation.
- Optional data are attempted one by one: layers, reductions, loadings, var metadata, uns metadata.
- Optional parts that cannot be converted are skipped instead of killing the whole job.
- A terminal summary and a sidecar JSON report are written, listing converted and skipped parts.

Examples:
  python h5ad_to_rds.py input.h5ad --assay-name RNA
  python h5ad_to_rds.py input.h5ad -o output.rds --assay-name RNA --project converted
  python h5ad_to_rds.py input.h5ad --assay-name RNA --strict-counts
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
from scipy.io import mmwrite
from rpy2 import robjects as ro


SCRIPT_VERSION = "2026-05-04-tolerant-v3"


def log(msg: str) -> None:
    print(f"[h5ad_to_rds] {msg}", flush=True)


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

make_key <- function(name) {
  clean <- gsub("[^A-Za-z0-9]", "", toupper(name))
  if (clean == "") clean <- "RED"
  paste0(clean, "_")
}

append_item <- function(x, value) {
  c(x, as.character(value))
}

read_lines_gz <- function(path) {
  readLines(gzfile(path), warn = FALSE)
}

read_table_gz <- function(path) {
  df <- read.csv(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  if (!"__index__" %in% colnames(df)) stop("Missing __index__ column in: ", path)
  rn <- as.character(df$`__index__`)
  df$`__index__` <- NULL
  rownames(df) <- rn
  df
}

safe_layers <- function(obj, assay_name) {
  out <- tryCatch(SeuratObject::Layers(obj[[assay_name]]), error = function(e) character(0))
  unique(as.character(out))
}

safe_reductions <- function(obj) {
  out <- tryCatch(SeuratObject::Reductions(obj), error = function(e) names(obj@reductions))
  unique(as.character(out))
}

inspect_rds <- function(input_file, assay_name = "") {
  obj <- readRDS(input_file)
  info <- list(class = class(obj))
  if (inherits(obj, "Seurat")) {
    info$dims <- unname(dim(obj))
    info$assays <- tryCatch(as.character(SeuratObject::Assays(obj)), error = function(e) names(obj@assays))
    info$default_assay <- tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) NULL)
    if (!nzchar(assay_name)) assay_name <- info$default_assay
    info$layers <- safe_layers(obj, assay_name)
    info$reductions <- safe_reductions(obj)
    info$obs_cols <- colnames(obj[[]])
  }
  jsonlite::toJSON(info, auto_unbox = TRUE, pretty = FALSE, null = "null")
}

build_rds_from_bundle <- function(bundle_dir, output_file, assay_name = "RNA", project = "converted") {
  manifest <- jsonlite::fromJSON(file.path(bundle_dir, "manifest.json"), simplifyVector = FALSE)

  build_report <- list(
    converted = list(layers = character(), reductions = character(), loadings = character(), var_metadata = FALSE, uns = FALSE, idents = FALSE),
    skipped = list(layers = character(), reductions = character(), loadings = character(), var_metadata = character(), uns = character(), idents = character(), other = character()),
    notes = character()
  )

  cells <- read_lines_gz(file.path(bundle_dir, "cells.tsv.gz"))
  features <- read_lines_gz(file.path(bundle_dir, "features.tsv.gz"))

  counts_path <- NULL
  counts_name <- NULL
  if (!is.null(manifest$layers$counts)) {
    counts_path <- file.path(bundle_dir, manifest$layers$counts)
    counts_name <- "counts"
  } else if (!is.null(manifest$x_source) && manifest$x_source %in% names(manifest$layers)) {
    counts_path <- file.path(bundle_dir, manifest$layers[[manifest$x_source]])
    counts_name <- manifest$x_source
  } else if (length(manifest$layers) > 0) {
    counts_name <- names(manifest$layers)[[1]]
    counts_path <- file.path(bundle_dir, manifest$layers[[1]])
  }
  if (is.null(counts_path) || !file.exists(counts_path)) stop("no usable matrix found in bundle for Seurat creation")

  counts <- Matrix::readMM(counts_path)
  counts <- methods::as(counts, "dgCMatrix")
  if (nrow(counts) != length(features) || ncol(counts) != length(cells)) {
    stop("core matrix dimension mismatch: ", paste(dim(counts), collapse = "x"), " vs ", length(features), "x", length(cells))
  }
  rownames(counts) <- features
  colnames(counts) <- cells

  obs <- read_table_gz(file.path(bundle_dir, "obs.csv.gz"))
  obs <- obs[cells, , drop = FALSE]

  project_name <- project
  if (!is.null(manifest$project) && is.character(manifest$project) && nzchar(manifest$project)) project_name <- manifest$project

  obj <- CreateSeuratObject(counts = counts, assay = assay_name, meta.data = obs, project = project_name)
  build_report$converted$layers <- c(build_report$converted$layers, "counts")
  if (!identical(counts_name, "counts")) {
    build_report$notes <- c(build_report$notes, paste0("Seurat counts layer was created from ", counts_name))
  }

  var_path <- file.path(bundle_dir, "var.csv.gz")
  if (file.exists(var_path)) {
    ok_var <- tryCatch({
      var <- read_table_gz(var_path)
      var <- var[features, , drop = FALSE]
      if (ncol(var) > 0) {
        obj[[assay_name]] <- AddMetaData(object = obj[[assay_name]], metadata = var)
        build_report$converted$var_metadata <- TRUE
      }
      TRUE
    }, error = function(e) {
      build_report$skipped$var_metadata <<- c(build_report$skipped$var_metadata, conditionMessage(e))
      FALSE
    })
  }

  if (length(manifest$layers) > 0) {
    for (layer_name in names(manifest$layers)) {
      if (layer_name == "counts") next
      layer_file <- file.path(bundle_dir, manifest$layers[[layer_name]])
      if (!file.exists(layer_file)) {
        build_report$skipped$layers <- c(build_report$skipped$layers, paste0(layer_name, ": file not found"))
        next
      }
      ok_layer <- tryCatch({
        mat <- Matrix::readMM(layer_file)
        mat <- methods::as(mat, "dgCMatrix")
        if (nrow(mat) != length(features) || ncol(mat) != length(cells)) stop("shape mismatch")
        rownames(mat) <- features
        colnames(mat) <- cells
        SeuratObject::LayerData(object = obj, assay = assay_name, layer = layer_name) <- mat
        TRUE
      }, error = function(e) {
        build_report$skipped$layers <<- c(build_report$skipped$layers, paste0(layer_name, ": ", conditionMessage(e)))
        FALSE
      })
      if (ok_layer) build_report$converted$layers <- c(build_report$converted$layers, layer_name)
    }
  }

  if (length(manifest$reductions) > 0) {
    for (red_name in names(manifest$reductions)) {
      red_info <- manifest$reductions[[red_name]]
      emb_file <- if (is.character(red_info)) red_info else red_info$embeddings
      load_file <- if (!is.character(red_info) && !is.null(red_info$loadings) && nzchar(red_info$loadings)) red_info$loadings else NULL
      red_key <- if (!is.character(red_info) && !is.null(red_info$key) && nzchar(red_info$key)) red_info$key else make_key(red_name)

      emb_path <- file.path(bundle_dir, emb_file)
      if (!file.exists(emb_path)) {
        build_report$skipped$reductions <- c(build_report$skipped$reductions, paste0(red_name, ": embeddings file not found"))
        next
      }

      ok_red <- tryCatch({
        emb <- read_table_gz(emb_path)
        emb <- as.matrix(emb[cells, , drop = FALSE])
        storage.mode(emb) <- "numeric"
        if (nrow(emb) != length(cells) || ncol(emb) == 0) stop("invalid embedding shape")
        emb[!is.finite(emb)] <- 0
        colnames(emb) <- paste0(red_key, seq_len(ncol(emb)))

        loadings <- NULL
        if (!is.null(load_file)) {
          load_path <- file.path(bundle_dir, load_file)
          if (file.exists(load_path)) {
            loadings <- tryCatch({
              ldf <- read_table_gz(load_path)
              if (ncol(ldf) == 0) stop("empty loadings")
              full <- matrix(0, nrow = length(features), ncol = ncol(ldf), dimnames = list(features, colnames(ldf)))
              common <- intersect(features, rownames(ldf))
              if (length(common) > 0) full[common, ] <- as.matrix(ldf[common, , drop = FALSE])
              storage.mode(full) <- "numeric"
              if (ncol(full) == ncol(emb)) colnames(full) <- colnames(emb) else colnames(full) <- paste0(red_key, seq_len(ncol(full)))
              full
            }, error = function(e) {
              build_report$skipped$loadings <<- c(build_report$skipped$loadings, paste0(red_name, ": ", conditionMessage(e)))
              NULL
            })
          } else {
            build_report$skipped$loadings <- c(build_report$skipped$loadings, paste0(red_name, ": loadings file not found"))
          }
        }

        if (is.null(loadings)) {
          dr <- CreateDimReducObject(embeddings = emb, assay = assay_name, key = red_key)
        } else {
          dr <- CreateDimReducObject(embeddings = emb, loadings = loadings, assay = assay_name, key = red_key)
          build_report$converted$loadings <<- c(build_report$converted$loadings, red_name)
        }
        obj[[red_name]] <- dr
        TRUE
      }, error = function(e) {
        build_report$skipped$reductions <<- c(build_report$skipped$reductions, paste0(red_name, ": ", conditionMessage(e)))
        FALSE
      })

      if (ok_red) build_report$converted$reductions <- c(build_report$converted$reductions, red_name)
    }
  }

  if ("seurat_ident" %in% colnames(obj[[]])) {
    ok_id <- tryCatch({
      lv <- tryCatch(unlist(manifest$idents_levels, use.names = FALSE), error = function(e) character())
      if (length(lv) > 0) SeuratObject::Idents(obj) <- factor(obj$seurat_ident, levels = lv) else SeuratObject::Idents(obj) <- obj$seurat_ident
      TRUE
    }, error = function(e) {
      build_report$skipped$idents <<- c(build_report$skipped$idents, conditionMessage(e))
      FALSE
    })
    if (ok_id) build_report$converted$idents <- TRUE
  }

  ok_misc <- tryCatch({
    obj@misc$rds_bridge <- list(
      source_class = manifest$source_class,
      selected_assay = manifest$selected_assay,
      x_source = manifest$x_source,
      notes = manifest$notes,
      skipped = manifest$skipped,
      converted = manifest$converted
    )
    uns_path <- file.path(bundle_dir, "uns.json.gz")
    if (file.exists(uns_path)) {
      uns_json <- paste(readLines(gzfile(uns_path), warn = FALSE), collapse = "\n")
      obj@misc$anndata_uns_json <- uns_json
    }
    TRUE
  }, error = function(e) {
    build_report$skipped$uns <<- c(build_report$skipped$uns, conditionMessage(e))
    FALSE
  })
  if (ok_misc) build_report$converted$uns <- TRUE

  saveRDS(obj, file = output_file)
  write(jsonlite::toJSON(build_report, auto_unbox = TRUE, pretty = TRUE, null = "null"), file = file.path(bundle_dir, "build_report.json"))
  invisible(output_file)
}
'''

ro.r(R_HELPERS)
R_BUILD = ro.globalenv["build_rds_from_bundle"]
R_INSPECT = ro.globalenv["inspect_rds"]


def infer_output_path(input_path: Path, output_path: str | None = None) -> Path:
    input_path = input_path.resolve()
    if input_path.suffix.lower() != ".h5ad":
        raise ValueError("Input must be .h5ad for h5ad_to_rds.py")
    if output_path is None:
        return input_path.with_suffix(".rds")
    out = Path(output_path)
    if out.suffix == "":
        out = out.with_suffix(".rds")
    return out.resolve()


def report_path_for(output_file: Path) -> Path:
    return Path(str(output_file) + ".conversion_report.json")


def sanitize_key(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))


def obsm_to_reduction_name(name: str) -> str:
    name = str(name)
    return name[2:] if name.startswith("X_") else name


def reduction_key(name: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9]", "", str(name).upper())
    return f"{clean or 'RED'}_"


def make_unique_str(names: list[Any]) -> tuple[list[str], list[str]]:
    seen: dict[str, int] = {}
    out: list[str] = []
    notes: list[str] = []
    changed = 0
    for x in names:
        base = str(x)
        if base not in seen:
            seen[base] = 0
            out.append(base)
        else:
            seen[base] += 1
            new = f"{base}-{seen[base]}"
            while new in seen:
                seen[base] += 1
                new = f"{base}-{seen[base]}"
            seen[new] = 0
            out.append(new)
            changed += 1
    if changed:
        notes.append(f"made {changed} duplicated names unique using '-N' suffix")
    return out, notes


def write_lines_gz_py(items, path: Path) -> None:
    with gzip.open(path, "wt") as fh:
        for x in items:
            fh.write(f"{x}\n")


def write_table_gz_py(df: pd.DataFrame, path: Path) -> None:
    out = df.copy()
    out.index = out.index.astype(str)
    if "__index__" in out.columns:
        out = out.rename(columns={"__index__": "__index___column"})
    out.insert(0, "__index__", out.index)
    out.to_csv(path, index=False, compression="gzip")


def write_mtx_from_anndata(x, mtx_path: Path) -> None:
    """Write AnnData-style cell x gene matrix as sparse MatrixMarket gene x cell."""
    if sp.issparse(x):
        mat = x.T.tocoo()
    else:
        arr = np.asarray(x)
        if arr.ndim != 2:
            raise ValueError(f"Expected 2D matrix, got shape={arr.shape}")
        mat = sp.coo_matrix(arr.T)
    mmwrite(str(mtx_path), mat)


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
        return {"__dataframe__": True, "shape": list(x.shape), "columns": [str(c) for c in x.columns]}
    if isinstance(x, np.ndarray):
        if x.size > 1000:
            return {"__ndarray__": True, "shape": list(x.shape), "dtype": str(x.dtype)}
        return x.tolist()
    if sp.issparse(x):
        return {"__sparse__": True, "shape": list(x.shape), "nnz": int(x.nnz)}
    if isinstance(x, dict):
        return {str(k): to_jsonable(v) for k, v in x.items() if v is not None}
    if isinstance(x, (list, tuple)):
        if len(x) > 1000:
            return {"__sequence__": True, "length": len(x)}
        return [to_jsonable(v) for v in x if v is not None]
    return repr(x)


def init_report(input_file: Path, output_file: Path, direction: str) -> dict[str, Any]:
    return {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "direction": direction,
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "input": str(input_file),
        "output": str(output_file),
        "status": "running",
        "converted": {"layers": [], "reductions": [], "loadings": [], "obs": False, "var": False, "var_metadata": False, "uns": False, "idents": False},
        "skipped": {"layers": [], "reductions": [], "loadings": [], "var_metadata": [], "uns": [], "idents": [], "other": []},
        "notes": [],
    }


def write_report(report: dict[str, Any], output_file: Path) -> Path:
    report["finished_at"] = datetime.now().isoformat(timespec="seconds")
    path = report_path_for(output_file)
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


def inspect_rds_json(rds_file: Path, assay_name: str) -> dict:
    txt = str(R_INSPECT(str(rds_file), assay_name or "")[0])
    return json.loads(txt)


def validate_h5ad_to_rds(output_file: Path, cells: list[str], features: list[str], assay_name: str) -> None:
    info = inspect_rds_json(output_file, assay_name)
    dims = info.get("dims") or []
    if list(map(int, dims)) != [len(features), len(cells)]:
        raise RuntimeError(f"validation failed: dimension mismatch; RDS dims={dims}, expected={[len(features), len(cells)]}")
    log(f"validate ok | features={dims[0]} cells={dims[1]} | layers={info.get('layers', [])} | reductions={info.get('reductions', [])}")


def write_matrix_layer_checked(name: str, x: Any, expected_shape: tuple[int, int], out_path: Path) -> None:
    shape = x.shape
    if tuple(shape) != tuple(expected_shape):
        raise ValueError(f"shape {shape} != expected {expected_shape}")
    write_mtx_from_anndata(x, out_path)


def h5ad_to_rds(
    input_file: Path,
    output_file: Path,
    assay_name: str,
    project: str,
    strict_counts: bool = False,
    validate: bool = True,
) -> dict[str, Any]:
    report = init_report(input_file, output_file, "h5ad_to_rds")

    log("1/5 read h5ad")
    adata = ad.read_h5ad(str(input_file))
    n_obs, n_vars = adata.n_obs, adata.n_vars
    report["n_obs"] = n_obs
    report["n_vars"] = n_vars

    cells, cell_notes = make_unique_str(list(adata.obs_names))
    features, feature_notes = make_unique_str(list(adata.var_names))
    report["notes"].extend(["obs_names: " + n for n in cell_notes])
    report["notes"].extend(["var_names: " + n for n in feature_notes])

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
    report["x_source"] = x_source

    with tempfile.TemporaryDirectory(prefix="h5ad_to_rds_") as tmpdir:
        tmpdir_p = Path(tmpdir)
        log("2/5 write temporary bundle")

        write_lines_gz_py(cells, tmpdir_p / "cells.tsv.gz")
        write_lines_gz_py(features, tmpdir_p / "features.tsv.gz")

        obs = adata.obs.copy()
        var = adata.var.copy()
        obs.index = cells
        var.index = features
        write_table_gz_py(obs, tmpdir_p / "obs.csv.gz")
        write_table_gz_py(var, tmpdir_p / "var.csv.gz")
        report["converted"]["obs"] = True
        report["converted"]["var"] = True

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
            "notes": report["notes"],
            "converted": report["converted"],
            "skipped": report["skipped"],
            "script_version": SCRIPT_VERSION,
        }

        expected_shape = (n_obs, n_vars)

        # Write named layers first. Any bad optional layer is skipped.
        for layer_name in list(adata.layers.keys()):
            try:
                file_name = f"layer__{sanitize_key(layer_name)}.mtx"
                write_matrix_layer_checked(layer_name, adata.layers[layer_name], expected_shape, tmpdir_p / file_name)
                manifest["layers"][str(layer_name)] = file_name
                report["converted"]["layers"].append(str(layer_name))
            except Exception as e:
                report["skipped"]["layers"].append(f"{layer_name}: {e}")

        # Ensure X has a named layer for round-trip. If x_source already exists, do not duplicate.
        if x_source not in manifest["layers"]:
            try:
                file_name = f"layer__{sanitize_key(x_source)}.mtx"
                write_matrix_layer_checked(x_source, adata.X, expected_shape, tmpdir_p / file_name)
                manifest["layers"][x_source] = file_name
                report["converted"]["layers"].append(x_source)
            except Exception as e:
                raise RuntimeError(f"core X matrix could not be written as layer '{x_source}': {e}") from e

        # Seurat creation requires a counts matrix. Prefer real counts; otherwise fall back and say so.
        if "counts" not in manifest["layers"]:
            recovered = False
            if adata.raw is not None:
                try:
                    raw_var_names = [str(x) for x in adata.raw.var_names]
                    if adata.raw.X.shape == expected_shape and raw_var_names == features:
                        file_name = "layer__counts.mtx"
                        write_matrix_layer_checked("counts", adata.raw.X, expected_shape, tmpdir_p / file_name)
                        manifest["layers"]["counts"] = file_name
                        report["converted"]["layers"].append("counts")
                        report["notes"].append("counts recovered from AnnData.raw.X")
                        recovered = True
                    else:
                        report["skipped"]["layers"].append("raw.X: shape or raw.var_names do not match adata.X/var_names; not used as counts")
                except Exception as e:
                    report["skipped"]["layers"].append(f"raw.X: could not use as counts: {e}")
            if not recovered and looks_like_counts(adata.X):
                file_name = "layer__counts.mtx"
                write_matrix_layer_checked("counts", adata.X, expected_shape, tmpdir_p / file_name)
                manifest["layers"]["counts"] = file_name
                report["converted"]["layers"].append("counts")
                report["notes"].append("counts recovered from AnnData.X because X looks integer-like")
                recovered = True
            if not recovered:
                msg = "no counts-like matrix found; Seurat counts will be created from x_source/X, so counts are not guaranteed raw counts"
                if strict_counts:
                    raise RuntimeError(msg + "; rerun without --strict-counts to allow fallback")
                report["notes"].append(msg)

        bridge_reductions = bridge.get("reductions", {}) if isinstance(bridge, dict) else {}

        for key in list(adata.obsm.keys()):
            red_name = obsm_to_reduction_name(key)
            try:
                arr = np.asarray(adata.obsm[key])
                if arr.ndim != 2 or arr.shape[0] != n_obs or arr.shape[1] == 0:
                    raise ValueError(f"invalid shape {arr.shape}, expected ({n_obs}, n_components)")
                arr = arr.astype(float, copy=False)
                if np.all(~np.isfinite(arr)):
                    raise ValueError("all values are NA/NaN/Inf")
                arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
                emb = pd.DataFrame(arr, index=cells)
                emb_file = f"obsm__{sanitize_key(red_name)}.csv.gz"
                write_table_gz_py(emb, tmpdir_p / emb_file)

                red_meta = bridge_reductions.get(red_name, {}) if isinstance(bridge_reductions, dict) else {}
                red_info = {"embeddings": emb_file, "key": red_meta.get("key") or reduction_key(red_name)}

                varm_key = f"{red_name}_loadings"
                if varm_key in adata.varm.keys():
                    try:
                        load_arr = np.asarray(adata.varm[varm_key])
                        if load_arr.ndim != 2 or load_arr.shape[0] != n_vars or load_arr.shape[1] == 0:
                            raise ValueError(f"invalid shape {load_arr.shape}, expected ({n_vars}, n_components)")
                        load_arr = np.nan_to_num(load_arr.astype(float, copy=False), nan=0.0, posinf=0.0, neginf=0.0)
                        load = pd.DataFrame(load_arr, index=features)
                        load_file = f"varm__{sanitize_key(red_name)}_loadings.csv.gz"
                        write_table_gz_py(load, tmpdir_p / load_file)
                        red_info["loadings"] = load_file
                        report["converted"]["loadings"].append(red_name)
                    except Exception as e:
                        report["skipped"]["loadings"].append(f"{red_name}: {e}")

                manifest["reductions"][red_name] = red_info
                report["converted"]["reductions"].append(red_name)
            except Exception as e:
                report["skipped"]["reductions"].append(f"{key}: {e}")

        # .uns is optional. Store a compact JSON representation in Seurat misc if possible.
        try:
            with gzip.open(tmpdir_p / "uns.json.gz", "wt") as fh:
                fh.write(json.dumps(to_jsonable(adata.uns), ensure_ascii=False, indent=2))
        except Exception as e:
            report["skipped"]["uns"].append(f"uns: could not serialize: {e}")

        # Update manifest after optional writes.
        manifest["notes"] = report["notes"]
        manifest["converted"] = report["converted"]
        manifest["skipped"] = report["skipped"]
        (tmpdir_p / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

        log("3/5 build Seurat in embedded R")
        R_BUILD(str(tmpdir_p), str(output_file), assay_name, str(manifest["project"]))

        # Merge R-side build report, because R may skip optional components while creating the Seurat object.
        build_report_path = tmpdir_p / "build_report.json"
        if build_report_path.exists():
            try:
                build_report = json.loads(build_report_path.read_text())
                for k, vals in (build_report.get("skipped") or {}).items():
                    if k in report["skipped"]:
                        if isinstance(vals, list):
                            report["skipped"][k].extend([str(v) for v in vals if str(v)])
                        elif vals:
                            report["skipped"][k].append(str(vals))
                # R-side converted is the source of truth for actual RDS content.
                rconv = build_report.get("converted") or {}
                for k in ["layers", "reductions", "loadings"]:
                    if k in rconv:
                        report["converted"][k] = list(dict.fromkeys(map(str, rconv.get(k) or [])))
                for k in ["var_metadata", "uns", "idents"]:
                    if k in rconv:
                        report["converted"][k] = bool(rconv.get(k))
                for note in build_report.get("notes") or []:
                    if str(note) and str(note) not in report["notes"]:
                        report["notes"].append(str(note))
            except Exception as e:
                report["skipped"]["other"].append(f"could not read R build report: {e}")

        if validate:
            log("4/5 validate output")
            validate_h5ad_to_rds(output_file, cells, features, assay_name)

    log("5/5 done")
    report["status"] = "ok"
    report_path = write_report(report, output_file)
    print_summary(report, report_path)
    return report


def print_summary(report: dict[str, Any], report_path: Path) -> None:
    log("SUCCESS: rds written")
    log(f"  output: {report['output']}")
    log(f"  report: {report_path}")
    log(f"  cells={report.get('n_obs')} genes={report.get('n_vars')} x_source={report.get('x_source')}")
    log("converted:")
    log(f"  layers:     {report['converted'].get('layers') or []}")
    log(f"  reductions: {report['converted'].get('reductions') or []}")
    log(f"  loadings:   {report['converted'].get('loadings') or []}")
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
    if report.get("notes"):
        log("notes:")
        for note in report["notes"]:
            log(f"  - {note}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert h5ad to Seurat/RDS with embedded R via rpy2. Optional parts that fail are skipped and reported.")
    parser.add_argument("input_file", help="Input .h5ad file")
    parser.add_argument("-o", "--output-file", default=None, help="Output .rds path")
    parser.add_argument("--assay-name", default="RNA", help="Seurat assay name, default: RNA")
    parser.add_argument("--project", default="converted", help="Seurat project name")
    parser.add_argument("--strict-counts", action="store_true", help="Fail if no counts-like matrix is available instead of falling back to X/x_source")
    parser.add_argument("--no-validate", action="store_true", help="Skip post-write validation")
    args = parser.parse_args()

    input_file = Path(args.input_file).resolve()
    if not input_file.exists():
        raise FileNotFoundError(input_file)
    output_file = infer_output_path(input_file, args.output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    report = h5ad_to_rds(
        input_file=input_file,
        output_file=output_file,
        assay_name=args.assay_name,
        project=args.project,
        strict_counts=args.strict_counts,
        validate=not args.no_validate,
    )

    print(json.dumps({"input": str(input_file), "output": str(output_file), "status": report["status"], "report": str(report_path_for(output_file))}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
