test_that("manifest run directories preserve absolute paths", {
  triage_home <- "/repo/Triage"
  dataset <- "Census_immune"

  relative <- Triage:::resolve_manifest_run_dir(triage_home, dataset, "20260906_033548")
  expect_equal(relative, file.path(triage_home, "outputs", dataset, "20260906_033548"))

  unix_absolute <- "/tmp/triage_census_smoke3/Census_immune/20260906_033548"
  expect_identical(
    Triage:::resolve_manifest_run_dir(triage_home, dataset, unix_absolute),
    unix_absolute
  )

  windows_absolute <- "C:\\triage\\runs\\Census_immune\\20260906_033548"
  expect_identical(
    Triage:::resolve_manifest_run_dir(triage_home, dataset, windows_absolute),
    windows_absolute
  )
  expect_false(startsWith(
    Triage:::resolve_manifest_run_dir(triage_home, dataset, windows_absolute),
    file.path(triage_home, "outputs", dataset)
  ))
})
