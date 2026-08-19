#!/usr/bin/env Rscript
expected <- c(
  "data/primary/Census_immune/final" = 16L,
  "data/primary/Sikkema_lung/final" = 11L,
  "data/primary/TS_kidney/final" = 4L,
  "data/primary/TS_pancreas/final" = 10L,
  "data/primary/Zheng_blood/final" = 9L,
  "data/validation/Anderson_DLPFC/final" = 18L,
  "data/validation/Zha_AD_mouse/final" = 9L,
  "data/validation/S6K1_organoid/final" = 15L
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
