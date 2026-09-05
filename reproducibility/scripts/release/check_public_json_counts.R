#!/usr/bin/env Rscript
expected <- c(
  "reproducibility/primary/Census_immune/final" = 16L,
  "reproducibility/primary/Sikkema_lung/final" = 11L,
  "reproducibility/primary/TS_kidney/final" = 4L,
  "reproducibility/primary/TS_pancreas/final" = 10L,
  "reproducibility/primary/Zheng_blood/final" = 9L,
  "reproducibility/external/Anderson_DLPFC/final" = 18L,
  "reproducibility/external/Zha_AD_mouse/final" = 9L,
  "reproducibility/external/S6K1_organoid/final" = 15L
)
root <- Sys.getenv("TRIAGE_HOME", unset = getwd())
failed <- FALSE
for (rel in names(expected)) {
  d <- file.path(root, rel)
  n <- if (dir.exists(d)) length(list.files(d, pattern = "^cluster_.*\\.json$")) else 0L
  ok <- identical(as.integer(n), as.integer(expected[[rel]]))
  cat(sprintf("%-50s %3d / %3d  %s\n", rel, n, expected[[rel]], if (ok) "OK" else "MISMATCH"))
  if (!ok) failed <- TRUE
}
if (failed) quit(status = 1L)
