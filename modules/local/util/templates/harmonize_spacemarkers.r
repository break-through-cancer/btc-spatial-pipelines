#!/usr/bin/env Rscript

scores_path  <- "${scores}"
sample <- "${prefix}"
source <- "${source}"

#process imscores - can be 2 flavors: csv (old) or rds (new)
if (length(grep(".csv", scores_path)) > 0) {
    imscores <- read.csv(scores_path, row.names = 1)
    } else {
    imscores <- readRDS(scores_path)
}

dir.create(file.path(sample, source), showWarnings = FALSE, recursive = TRUE)

write.csv(imscores,
          file = gzfile(paste0(file.path(sample, source, scores_path), ".gz")))
