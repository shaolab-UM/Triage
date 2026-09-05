#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(jsonlite); library(optparse) })
opts <- parse_args(OptionParser(option_list = list(
  make_option("--final_dir", type="character", help="Dataset final/ directory containing public cluster JSON files"),
  make_option("--output", type="character", help="Output .jsonl file")
)))
if (is.null(opts$final_dir) || !dir.exists(opts$final_dir)) stop("--final_dir must exist")
if (is.null(opts$output) || !nzchar(opts$output)) stop("--output is required")
files <- list.files(opts$final_dir, pattern="^cluster_.*\\.json$", full.names=TRUE)
if (!length(files)) stop("No public cluster JSON files found")
con <- file(opts$output, open="wt", encoding="UTF-8")
on.exit(close(con), add=TRUE)
for (fp in files) {
  x <- fromJSON(fp, simplifyVector=FALSE)
  writeLines(toJSON(x, auto_unbox=TRUE, null="null", na="null"), con)
}
message("Wrote ", length(files), " records to ", opts$output)
