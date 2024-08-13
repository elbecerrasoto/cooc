#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(stringr)
library(skimr)
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



# DataExplorer ----

# plot_missing(meta)
# plot_missing(pairs)

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

# create_report(meta,
#   output_file = "meta.html",
#   report_title = "Genome Metadata Bacillota Report"
# )
#
# create_report(pairs,
#   output_file = "pairs.html",
#   report_title = "Deaminase Endonuclease Bacillota Pairs Report"
# )

# Genomes with more than 1 one pair ----

g_count <- pairs |>
  group_by(genome) |>
  summarise(n = n()) |>
  arrange(desc(n))

g_count |>
  pull(n) |>
  boxplot()

g_count |>
  pull(n) |>
  qplot()

table(g_count$n)

# Outlier on distance pairs ----


qplot(pairs$distance)

boxplot(pairs$distance)

x <- c(
  seq(0, 1000, 100),
  seq(1000, 10000, 1000)
) |>
  unique()

Pdistance <- cut(pairs$distance,
  breaks = x
)

boxplot(Pdistance)

qplot(Pdistance) +
  scale_x_discrete(labels = str_c("+", x))

skim(pairs)
skim(pairs$distance)

distanceBP <- boxplot(pairs$distance)
length(distanceBP$out)

q95 <- quantile(pairs$distance, seq(0, 1, 0.01))["95%"]

pairs_distance_outs <- pairs |>
  filter(distance > q95)

pairs <- pairs |>
  filter(distance <= q95)


# Pair set ----

nrow(pairs)

qpair_set <- unique(pairs$genome)
length(qpair_set)

all_set <- unique(meta$genome)
length(all_set)

map_chr(hits, class)
hits <- hits |>
  mutate(
    query_description = as_factor(query_description),
    order = as.integer(order),
    start = as.integer(start),
    end = as.integer(end),
    strand = as_factor(strand),
    query = as_factor(query)
  )


# create_report(hits,
#               output_file = "hits.html",
#               report_title = "PFAM hits Bacillota Report"
# )


hits <- hits |>
  filter(!is.na(query))

deam_set <- hits |>
  filter(query_description == "YwqJ-deaminase") |>
  pull(genome) |>
  unique()
length(deam_set)


endo_set <- hits |>
  filter(query_description == "Endonuclease_5") |>
  pull(genome) |>
  unique()
length(endo_set)

l <- length

intersect(endo_set, deam_set) |>
  l()

endo_set |> l()
deam_set |> l()


# pair <- intersect(endo_set, deam_set)

# assert equal
meta |> nrow() == unique(meta$genome) |> l()

all_set <- meta$genome




# set ari ----

# genome taxid set deam endo pair qpair qset r1 ... rn

# 1. all_set
# 2. endo_set
# 3. deam_set
# 4. qpair_set

granks <- tibble(
  genome = all_set,
  deam = all_set %in% deam_set,
  endo = all_set %in% endo_set,
  qpair = all_set %in% qpair_set
)

granks <- granks |>
  mutate(
    pair = deam & endo
  )

Egset <- rep(".", length(all_set))

attach(granks)
Egset[deam] <- "Deam"
Egset[endo] <- "Endo"
Egset[pair] <- "Pair"
detach(granks)

granks <- granks |>
  mutate(
    gset = Egset
  )

Eqgset <- rep(".", length(all_set))

attach(granks)
Eqgset[deam] <- "Deam"
Eqgset[endo] <- "Endo"
Eqgset[pair] <- "Pair"
Eqgset[qpair] <- "Qpair"
detach(granks)

granks <- granks |>
  mutate(
    qgset = Eqgset
  )

COLS <- c(
  "genome", "org", "tax_id",
  "deam", "endo", "pair",
  "qpair", "gset", "qgset"
)

granks <- granks |>
  left_join(meta, join_by(genome)) |>
  select(COLS)


# write_tsv(granks, "granks_bacillota.tsv")
