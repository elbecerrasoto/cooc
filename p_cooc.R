#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(DataExplorer)

PAIRS <- "results/pairs.tsv"
META <- "results/genomes_metadata.tsv"
HITS <- "results/hits_query.tsv"

pairs <- read_tsv(PAIRS)
meta <- read_tsv(META)
hits <- read_tsv(HITS)


ordered_DE <- function(s1, s2, p1, p2, deam = "YwqJ-deaminase", endo = "Endonuclease_5") {
  out <- s1 == s2

  for (i in seq_along(out)) {
    if (out[i]) {
      pair <- c(p1[i], p2[i])

      if (s1[i] == "+") {
        valid <- all(pair == c(deam, endo))
        out[i] <- valid
      } else if (s1[i] == "-") {
        valid <- all(pair == c(endo, deam))
        out[i] <- valid
      } else {
        print("error")
      }
    }
  }
  out
}


meta_bad <- meta |>
  filter(!str_starts(org, "[A-Z][a-z]+") |
    is.na(tax_id) |
    status != "current")
print(glue("Excluded: {nrow(meta_bad)} from meta"))


meta <- meta |>
  filter(
    str_starts(org, "[A-Z][a-z]+"),
    !is.na(tax_id),
    status == "current"
  )


pairs_bad <- pairs |>
  filter(!ordered_DE(strand_1, strand_2, query_1, query_2))
print(glue("Excluded: {nrow(pairs_bad)} from meta"))



pairs <- pairs |>
  filter(ordered_DE(strand_1, strand_2, query_1, query_2))
