test_that("installed package ships the mapper/verifier prompt resources", {
  prompts_path <- system.file("prompts", "cl_mapper_prompts.R", package = "Triage")
  expect_true(nzchar(prompts_path), info = "inst/prompts/cl_mapper_prompts.R is missing from the installed package")
  expect_true(file.exists(prompts_path))

  env <- new.env(parent = baseenv())
  source(prompts_path, local = env)

  expect_type(env$MAPPER_SYSTEM_PROMPT, "character")
  expect_true(nzchar(env$MAPPER_SYSTEM_PROMPT))
  expect_type(env$MAPPER_USER_TEMPLATE, "closure")
  expect_type(env$MAPPER_USER_TEMPLATE_EVIDENCE, "closure")
  expect_type(env$VERIFIER_SYSTEM_PROMPT, "character")
  expect_true(nzchar(env$VERIFIER_SYSTEM_PROMPT))
  expect_type(env$VERIFIER_USER_TEMPLATE, "closure")

  expect_match(env$MAPPER_SYSTEM_PROMPT, "Cell Ontology", fixed = TRUE)
  expect_match(env$VERIFIER_SYSTEM_PROMPT, "ontology mapping verifier", fixed = TRUE)
})
