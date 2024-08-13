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



# EDA ---------------------------------------------------------------------

plot_missing(meta)
plot_missing(pairs)

L_LEVEL <- c(
  "Complete Genome",
  "Chromosome",
  "Scaffold",
  "Contig"
)

L_STRAND <- c(
  "+",
  "-"
)

L_QUERY <- c(
  "YwqJ-deaminase",
  "Endonuclease_5"
)


meta <- meta |>
  select(-status) |>
  mutate(
    tax_id = as_factor(tax_id),
    level = factor(level,
      levels = L_LEVEL,
      ordered = TRUE
    ),
    cds = as.integer(cds)
  )


pairs <- pairs |>
  mutate(
    distance = as.integer(distance),
    genes_inbet = as.integer(genes_inbet),
    order_1 = as.integer(order_1),
    order_2 = as.integer(order_2),
    start_1 = as.integer(start_1),
    start_2 = as.integer(start_2),
    end_1 = as.integer(end_1),
    end_2 = as.integer(end_2),
    strand_1 = factor(strand_1,
      levels = L_STRAND
    ),
    strand_2 = factor(strand_2,
      levels = L_STRAND
    ),
    query_1 = factor(query_1,
      levels = L_QUERY
    ),
    query_2 = factor(query_2,
      levels = L_QUERY
    )
  )

create_report(meta,
  output_file = "meta.html",
  report_title = "Genome Metadata Bacillota Report"
)

create_report(pairs,
  output_file = "pairs.html",
  report_title = "Deaminase Endonuclease Bacillota Pairs Report"
)
