# Species-level parallelism that is robust on Windows: PSOCK workers (future::multisession, no forking), each
# running whole units of work (e.g. one species: all draws sequentially) and returning one small result row.

#' Number of parallel workers
#'
#' Physical cores minus `reserve`, capped by available RAM (`mem_per_worker_gb` each) and `max_workers`.
#' @param max_workers Upper limit (NULL = none).
#' @param reserve Cores left free (default 2).
#' @param mem_per_worker_gb RAM per worker in GB (default 0.5, the kernel model at 5 km).
#' @param avail_gb Available RAM in GB; default from `ps::ps_system_memory()` when installed, else no RAM cap.
#' @return Integer >= 1.
#' @export
sg_workers <- function(max_workers = NULL, reserve = 2, mem_per_worker_gb = 0.5, avail_gb = NULL) {
  cores <- parallel::detectCores(logical = FALSE)
  if (is.na(cores)) cores <- 1L
  n <- cores - reserve
  if (is.null(avail_gb) && requireNamespace("ps", quietly = TRUE)) {
    avail_gb <- tryCatch(ps::ps_system_memory()$avail / 2^30, error = function(e) NULL)
  }
  if (!is.null(avail_gb)) n <- min(n, floor(avail_gb / mem_per_worker_gb))
  if (!is.null(max_workers)) n <- min(n, max_workers)
  max(1L, as.integer(n))
}

#' Prepare an R process for one worker
#'
#' One thread for terra, data.table and OpenMP (parallelism is across workers), and a short local temp directory for
#' terra (avoids OneDrive sync and Windows MAX_PATH problems).
#' @param dir Local directory for terra's temporary files (created if missing).
#' @return Invisibly, the temp directory.
#' @export
sg_worker_setup <- function(dir = file.path(Sys.getenv("TEMP", tempdir()), "spreadgen")) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  terra::terraOptions(tempdir = dir, progress = 0)
  tryCatch(terra::terraOptions(threads = 1), error = function(e) NULL)
  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
  invisible(dir)
}

#' Run units of work in batches on PSOCK workers
#'
#' Batches of `workers` items run in parallel with `future.apply::future_lapply()` on `future::multisession`
#' (Windows-safe; no forking). Every item is wrapped in `tryCatch`: an error becomes a row with
#' `status = "FAILED: <message>"`. If a whole batch fails (e.g. a worker process dies), it is retried once
#' sequentially in the main process. After each batch, `on_batch(rows)` runs in the main process (write the master
#' table, copy outputs), so an interrupted run keeps every finished batch.
#'
#' Workers load the INSTALLED package (not `devtools::load_all()`); the run stops if a worker's spreadgen version
#' differs from `expected_version`.
#'
#' @param items List (or vector) of work items.
#' @param fun Function `fun(item, ...)` returning a one-row data frame (or named list). Runs in a worker: it must read
#'   its own inputs (no terra objects in `items` or `...`).
#' @param ... Further arguments for `fun` (plain R objects only).
#' @param workers Number of workers; 1 = sequential in the main process.
#' @param on_batch Function called with the data frame of the batch's rows (main process); NULL = none.
#' @param key_fun Function giving an item's key for FAILED rows (default `as.character`).
#' @param packages Packages attached in the workers (default `"spreadgen"`).
#' @param expected_version spreadgen version the workers must have (default: the version loaded in the main process).
#' @param setup Arguments for [sg_worker_setup()] in each worker.
#' @param log_fun Function for progress messages (default `message`).
#' @return Data frame of all rows (`data.table::rbindlist(fill = TRUE)` semantics).
#' @export
run_batches <- function(items, fun, ..., workers = 1, on_batch = NULL, key_fun = function(x) as.character(x)[1],
                        packages = "spreadgen", expected_version = as.character(utils::packageVersion("spreadgen")),
                        setup = list(), log_fun = message) {
  items <- as.list(items)
  args <- list(...)
  safe <- function(item) {
    t0 <- Sys.time()
    out <- tryCatch({
      do.call(sg_worker_setup, setup)
      r <- do.call(fun, c(list(item), args))
      as.data.frame(r, stringsAsFactors = FALSE)
    }, error = function(e) data.frame(key = key_fun(item), status = paste("FAILED:", conditionMessage(e)),
                                      stringsAsFactors = FALSE))
    out$worker_pid <- Sys.getpid()
    out$worker_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
    out
  }
  bind <- function(rows) {
    if (requireNamespace("data.table", quietly = TRUE)) as.data.frame(data.table::rbindlist(rows, fill = TRUE))
    else do.call(rbind, rows)
  }
  if (workers > 1) {
    old <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(old), add = TRUE)
    v <- future::value(future::future(as.character(utils::packageVersion("spreadgen")), packages = packages,
                                      seed = NULL))
    if (!identical(v, expected_version)) {
      stop("spreadgen version in the workers (", v, ") differs from the main process (", expected_version,
           "): reinstall the package (R CMD INSTALL / remotes::install_local) before the run")
    }
  }
  batches <- split(seq_along(items), ceiling(seq_along(items) / max(1, workers)))
  all_rows <- list()
  for (b in seq_along(batches)) {
    ix <- batches[[b]]
    rows <- if (workers > 1) {
      tryCatch(future.apply::future_lapply(items[ix], safe, future.packages = packages, future.seed = NULL),
               error = function(e) {
                 log_fun(sprintf("batch %d failed in the workers (%s): retrying sequentially", b, conditionMessage(e)))
                 lapply(items[ix], safe)
               })
    } else lapply(items[ix], safe)
    rows <- bind(rows)
    if (!is.null(on_batch)) on_batch(rows)
    all_rows[[b]] <- rows
    log_fun(sprintf("batch %d/%d done (%d item(s))", b, length(batches), length(ix)))
  }
  bind(all_rows)
}
