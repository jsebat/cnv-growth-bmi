#!/usr/bin/env Rscript

# Reverse pathway analysis:
# 1. Discover exact pathways enriched among Kosmicki height genes.
# 2. Test negative-effect bias in those same pathways.
# 3. Test those same pathways in asymmetric and mirror CNV genes.
# 4. Prune only highly overlapping Kosmicki-enriched pathways.
#
# Required input files:
#   kosmicki_207_genes_annotated.csv
#   cnv_genes_annotated.csv
#   GO_Biological_Process_2023.gmt
#   KEGG_2021_Human.gmt
#
# The GMT-like files distributed with this analysis have:
#   pathway_name <tab> blank_description <tab> gene1 <tab> gene2 ...

options(stringsAsFactors = FALSE)

required_packages <- c("data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install required R package(s) before running: ",
    paste(missing_packages, collapse = ", ")
  )
}
library(data.table)

# ----------------------------
# Configuration
# ----------------------------

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else getwd()
output_dir <- if (length(args) >= 2) args[2] else file.path(input_dir, "kosmicki_reverse_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

kosmicki_file <- file.path(input_dir, "kosmicki_207_genes_annotated.csv")
cnv_file <- file.path(input_dir, "cnv_genes_annotated.csv")
go_file <- file.path(input_dir, "GO_Biological_Process_2023.gmt")
kegg_file <- file.path(input_dir, "KEGG_2021_Human.gmt")

min_set_size <- 5L
max_set_size <- 500L
min_kosmicki_overlap <- 3L
jaccard_threshold <- 0.50
n_locus_permutations <- 5000L
random_seed <- 20260719L

for (f in c(kosmicki_file, cnv_file, go_file, kegg_file)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

clean_gene <- function(x) toupper(trimws(as.character(x)))

safe_or <- function(a, b, c, d) {
  vals <- c(a, b, c, d)
  if (any(vals == 0)) {
    return(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
  }
  (a * d) / (b * c)
}

one_sided_hypergeom <- function(overlap, universe_n, target_n, set_n) {
  if (set_n <= 0 || target_n <= 0 || universe_n <= 0) return(1)
  phyper(
    q = overlap - 1,
    m = target_n,
    n = universe_n - target_n,
    k = set_n,
    lower.tail = FALSE
  )
}

read_gmt_like <- function(path, library_name) {
  lines <- readLines(path, warn = FALSE)
  parsed <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3) return(NULL)
    term <- trimws(fields[1])
    genes <- unique(clean_gene(fields[3:length(fields)]))
    genes <- genes[nzchar(genes)]
    if (length(genes) < min_set_size || length(genes) > max_set_size) return(NULL)
    list(
      pathway_key = paste(library_name, term, sep = "::"),
      library = library_name,
      pathway = term,
      genes = genes
    )
  })
  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  names(parsed) <- vapply(parsed, `[[`, character(1), "pathway_key")
  parsed
}

ora_one <- function(target, universe, geneset) {
  target <- intersect(target, universe)
  geneset <- intersect(geneset, universe)

  a <- length(intersect(target, geneset))
  b <- length(target) - a
  c <- length(geneset) - a
  d <- length(universe) - a - b - c

  list(
    overlap = a,
    target_n = length(target),
    set_n = length(geneset),
    universe_n = length(universe),
    odds_ratio = safe_or(a, b, c, d),
    fold_enrichment = if (length(target) > 0 && length(geneset) > 0) {
      (a / length(target)) / (length(geneset) / length(universe))
    } else {
      NA_real_
    },
    p = one_sided_hypergeom(a, length(universe), length(target), length(geneset))
  )
}

# Connected components among enriched pathways using exact full-set Jaccard overlap.
prune_pathways <- function(enriched, membership_list, threshold = 0.50) {
  keys <- enriched$pathway_key
  n <- length(keys)

  adjacency <- vector("list", n)
  for (i in seq_len(n)) adjacency[[i]] <- integer()

  if (n >= 2) {
    for (i in seq_len(n - 1)) {
      si <- membership_list[[keys[i]]]
      for (j in (i + 1):n) {
        sj <- membership_list[[keys[j]]]
        inter_n <- length(intersect(si, sj))
        if (inter_n == 0) next
        jaccard <- inter_n / length(union(si, sj))
        if (jaccard >= threshold) {
          adjacency[[i]] <- c(adjacency[[i]], j)
          adjacency[[j]] <- c(adjacency[[j]], i)
        }
      }
    }
  }

  cluster <- rep(NA_integer_, n)
  cluster_id <- 0L
  for (i in seq_len(n)) {
    if (!is.na(cluster[i])) next
    cluster_id <- cluster_id + 1L
    queue <- i
    cluster[i] <- cluster_id
    while (length(queue) > 0) {
      current <- queue[1]
      queue <- queue[-1]
      neighbors <- adjacency[[current]]
      new_nodes <- neighbors[is.na(cluster[neighbors])]
      if (length(new_nodes) > 0) {
        cluster[new_nodes] <- cluster_id
        queue <- c(queue, new_nodes)
      }
    }
  }

  enriched[, redundancy_cluster := cluster]
  enriched[, redundancy_cluster_size := .N, by = redundancy_cluster]

  representatives <- enriched[
    order(
      kosmicki_q_within_library,
      kosmicki_p,
      -kosmicki_fold_enrichment,
      pathway_size
    ),
    .SD[1],
    by = redundancy_cluster
  ][, .(redundancy_cluster, representative_pathway_key = pathway_key)]

  enriched <- merge(
    enriched,
    representatives,
    by = "redundancy_cluster",
    all.x = TRUE,
    sort = FALSE
  )
  enriched[, retained_after_pruning := pathway_key == representative_pathway_key]
  enriched
}

# ----------------------------
# Load genes and exact gene sets
# ----------------------------

kos <- fread(kosmicki_file)
cnv <- fread(cnv_file)

kos[, Gene_std := clean_gene(Gene)]
kos[, Negative := tolower(Effect_direction) == "negative"]
kos[, Inside_CNV := tolower(as.character(Inside_any_analyzed_CNV_interval)) == "yes"]
kos[, Effect_SD := as.numeric(Effect_SD)]

cnv[, Gene_std := clean_gene(Gene)]
cnv[, Dose_response_class := tolower(Dose_response_class)]

go_sets <- read_gmt_like(go_file, "GO_Biological_Process_2023")
kegg_sets <- read_gmt_like(kegg_file, "KEGG_2021_Human")
all_sets <- c(go_sets, kegg_sets)

library_names <- c("GO_Biological_Process_2023", "KEGG_2021_Human")
library_universes <- setNames(vector("list", length(library_names)), library_names)
for (lib in library_names) {
  lib_sets <- all_sets[vapply(all_sets, function(x) x$library == lib, logical(1))]
  library_universes[[lib]] <- unique(unlist(lapply(lib_sets, `[[`, "genes")))
}

membership_list <- lapply(all_sets, `[[`, "genes")

# ----------------------------
# 1. Systematic Kosmicki pathway discovery
# ----------------------------

discovery_rows <- lapply(all_sets, function(x) {
  universe <- library_universes[[x$library]]
  target <- intersect(kos$Gene_std, universe)
  st <- ora_one(target, universe, x$genes)

  data.table(
    pathway_key = x$pathway_key,
    library = x$library,
    pathway = x$pathway,
    pathway_size = st$set_n,
    library_universe_n = st$universe_n,
    kosmicki_mapped_n = st$target_n,
    kosmicki_overlap_n = st$overlap,
    kosmicki_odds_ratio = st$odds_ratio,
    kosmicki_fold_enrichment = st$fold_enrichment,
    kosmicki_p = st$p,
    overlap_genes = paste(sort(intersect(target, x$genes)), collapse = ";")
  )
})
results <- rbindlist(discovery_rows, fill = TRUE)

results[
  ,
  kosmicki_q_within_library := p.adjust(kosmicki_p, method = "BH"),
  by = library
]
results[, kosmicki_q_global := p.adjust(kosmicki_p, method = "BH")]
results[
  ,
  kosmicki_enriched_primary :=
    kosmicki_q_within_library < 0.05 &
    kosmicki_overlap_n >= min_kosmicki_overlap &
    kosmicki_odds_ratio > 1
]
results[
  ,
  kosmicki_enriched_global_sensitivity :=
    kosmicki_q_global < 0.05 &
    kosmicki_overlap_n >= min_kosmicki_overlap &
    kosmicki_odds_ratio > 1
]

# ----------------------------
# 2. Directional tests and exact CNV tests
#    These are calculated for every tested pathway in the master table.
# ----------------------------

all_kosmicki <- unique(kos$Gene_std)
all_negative <- unique(kos[Negative == TRUE, Gene_std])
independent_kosmicki <- unique(kos[Inside_CNV == FALSE, Gene_std])
independent_negative <- intersect(all_negative, independent_kosmicki)

all_cnv_raw <- unique(cnv$Gene_std)
asymmetric_raw <- unique(cnv[Dose_response_class == "asymmetric", Gene_std])
mirror_raw <- unique(cnv[Dose_response_class == "mirror", Gene_std])

followup_rows <- lapply(all_sets, function(x) {
  geneset <- x$genes
  universe <- library_universes[[x$library]]

  # Primary directional background: Kosmicki genes mapped to this library.
  mapped_kosmicki <- intersect(all_kosmicki, universe)
  mapped_negative <- intersect(all_negative, mapped_kosmicki)

  n <- length(intersect(geneset, mapped_kosmicki))
  a <- length(intersect(geneset, mapped_negative))
  b <- n - a
  c <- length(mapped_negative) - a
  d <- length(mapped_kosmicki) - length(mapped_negative) - b

  # All-207 sensitivity.
  n207 <- length(intersect(geneset, all_kosmicki))
  a207 <- length(intersect(geneset, all_negative))
  b207 <- n207 - a207
  c207 <- length(all_negative) - a207
  d207 <- length(all_kosmicki) - length(all_negative) - b207

  # Exclude four Kosmicki genes lying inside analyzed CNV intervals.
  mapped_independent <- intersect(independent_kosmicki, universe)
  mapped_independent_negative <- intersect(independent_negative, mapped_independent)
  ni <- length(intersect(geneset, mapped_independent))
  ai <- length(intersect(geneset, mapped_independent_negative))
  bi <- ni - ai
  ci <- length(mapped_independent_negative) - ai
  di <- length(mapped_independent) - length(mapped_independent_negative) - bi

  # CNV tests use the same pathway library's mapped-gene universe.
  cnv_universe <- intersect(all_cnv_raw, universe)
  asymmetric <- intersect(asymmetric_raw, cnv_universe)
  mirror <- intersect(mirror_raw, cnv_universe)
  asymmetric_only <- setdiff(asymmetric, mirror)
  mirror_only <- setdiff(mirror, asymmetric)

  set_cnv_n <- length(intersect(geneset, cnv_universe))

  asym_overlap <- length(intersect(geneset, asymmetric))
  asym_b <- set_cnv_n - asym_overlap
  asym_c <- length(asymmetric) - asym_overlap
  asym_d <- length(cnv_universe) - length(asymmetric) - asym_b

  mirror_overlap <- length(intersect(geneset, mirror))
  mirror_b <- set_cnv_n - mirror_overlap
  mirror_c <- length(mirror) - mirror_overlap
  mirror_d <- length(cnv_universe) - length(mirror) - mirror_b

  ao <- length(intersect(geneset, asymmetric_only))
  mo <- length(intersect(geneset, mirror_only))

  # Circularity sensitivity: remove all Kosmicki genes from CNV tests.
  cnv_u2 <- setdiff(cnv_universe, all_kosmicki)
  asym2 <- setdiff(asymmetric, all_kosmicki)
  mirror2 <- setdiff(mirror, all_kosmicki)
  set2_n <- length(intersect(geneset, cnv_u2))
  aa2 <- length(intersect(geneset, asym2))
  ma2 <- length(intersect(geneset, mirror2))

  data.table(
    pathway_key = x$pathway_key,

    direction_overlap_n = n,
    direction_negative_n = a,
    direction_positive_n = b,
    direction_negative_fraction = if (n > 0) a / n else NA_real_,
    direction_background_negative_fraction =
      if ((c + d) > 0) c / (c + d) else NA_real_,
    direction_odds_ratio = safe_or(a, b, c, d),
    direction_fisher_p_one_sided =
      one_sided_hypergeom(a, length(mapped_kosmicki), length(mapped_negative), n),

    all207_direction_overlap_n = n207,
    all207_direction_negative_n = a207,
    all207_direction_positive_n = b207,
    all207_direction_negative_fraction = if (n207 > 0) a207 / n207 else NA_real_,
    all207_direction_odds_ratio = safe_or(a207, b207, c207, d207),
    all207_direction_fisher_p_one_sided =
      one_sided_hypergeom(a207, length(all_kosmicki), length(all_negative), n207),

    independent_direction_overlap_n = ni,
    independent_direction_negative_n = ai,
    independent_direction_positive_n = bi,
    independent_direction_negative_fraction = if (ni > 0) ai / ni else NA_real_,
    independent_direction_odds_ratio = safe_or(ai, bi, ci, di),
    independent_direction_fisher_p_one_sided =
      one_sided_hypergeom(
        ai,
        length(mapped_independent),
        length(mapped_independent_negative),
        ni
      ),

    cnv_library_mapped_universe_n = length(cnv_universe),
    asymmetric_cnv_target_n = length(asymmetric),
    asymmetric_cnv_pathway_n_in_universe = set_cnv_n,
    asymmetric_cnv_overlap_n = asym_overlap,
    asymmetric_cnv_odds_ratio = safe_or(asym_overlap, asym_b, asym_c, asym_d),
    asymmetric_cnv_fold_enrichment =
      if (length(asymmetric) > 0 && set_cnv_n > 0) {
        (asym_overlap / length(asymmetric)) /
          (set_cnv_n / length(cnv_universe))
      } else {
        NA_real_
      },
    asymmetric_cnv_p =
      one_sided_hypergeom(
        asym_overlap,
        length(cnv_universe),
        length(asymmetric),
        set_cnv_n
      ),

    mirror_cnv_target_n = length(mirror),
    mirror_cnv_pathway_n_in_universe = set_cnv_n,
    mirror_cnv_overlap_n = mirror_overlap,
    mirror_cnv_odds_ratio =
      safe_or(mirror_overlap, mirror_b, mirror_c, mirror_d),
    mirror_cnv_fold_enrichment =
      if (length(mirror) > 0 && set_cnv_n > 0) {
        (mirror_overlap / length(mirror)) /
          (set_cnv_n / length(cnv_universe))
      } else {
        NA_real_
      },
    mirror_cnv_p =
      one_sided_hypergeom(
        mirror_overlap,
        length(cnv_universe),
        length(mirror),
        set_cnv_n
      ),

    asymmetric_only_overlap_n = ao,
    asymmetric_only_n = length(asymmetric_only),
    mirror_only_overlap_n = mo,
    mirror_only_n = length(mirror_only),
    asymmetric_vs_mirror_odds_ratio =
      safe_or(
        ao,
        length(asymmetric_only) - ao,
        mo,
        length(mirror_only) - mo
      ),
    asymmetric_vs_mirror_p_one_sided =
      one_sided_hypergeom(
        ao,
        length(asymmetric_only) + length(mirror_only),
        ao + mo,
        length(asymmetric_only)
      ),

    asymmetric_cnv_or_excluding_kosmicki_genes =
      safe_or(
        aa2,
        set2_n - aa2,
        length(asym2) - aa2,
        length(cnv_u2) - length(asym2) - (set2_n - aa2)
      ),
    asymmetric_cnv_p_excluding_kosmicki_genes =
      one_sided_hypergeom(aa2, length(cnv_u2), length(asym2), set2_n),
    mirror_cnv_or_excluding_kosmicki_genes =
      safe_or(
        ma2,
        set2_n - ma2,
        length(mirror2) - ma2,
        length(cnv_u2) - length(mirror2) - (set2_n - ma2)
      ),
    mirror_cnv_p_excluding_kosmicki_genes =
      one_sided_hypergeom(ma2, length(cnv_u2), length(mirror2), set2_n)
  )
})

followup <- rbindlist(followup_rows, fill = TRUE)
results <- merge(results, followup, by = "pathway_key", all.x = TRUE, sort = FALSE)

results[, direction_q_all_tested :=
          p.adjust(direction_fisher_p_one_sided, method = "BH")]
results[, asymmetric_cnv_q_all_tested :=
          p.adjust(asymmetric_cnv_p, method = "BH")]
results[, mirror_cnv_q_all_tested :=
          p.adjust(mirror_cnv_p, method = "BH")]
results[, asymmetric_vs_mirror_q_all_tested :=
          p.adjust(asymmetric_vs_mirror_p_one_sided, method = "BH")]

results[
  kosmicki_enriched_primary == TRUE,
  direction_q_among_kosmicki_enriched :=
    p.adjust(direction_fisher_p_one_sided, method = "BH")
]
results[
  kosmicki_enriched_primary == TRUE,
  independent_direction_q_among_kosmicki_enriched :=
    p.adjust(independent_direction_fisher_p_one_sided, method = "BH")
]
results[
  kosmicki_enriched_primary == TRUE,
  asymmetric_cnv_q_among_kosmicki_enriched :=
    p.adjust(asymmetric_cnv_p, method = "BH")
]
results[
  kosmicki_enriched_primary == TRUE,
  mirror_cnv_q_among_kosmicki_enriched :=
    p.adjust(mirror_cnv_p, method = "BH")
]
results[
  kosmicki_enriched_primary == TRUE,
  asymmetric_vs_mirror_q_among_kosmicki_enriched :=
    p.adjust(asymmetric_vs_mirror_p_one_sided, method = "BH")
]

# ----------------------------
# 3. Prune only exact pathways already enriched in Kosmicki
# ----------------------------

enriched <- results[kosmicki_enriched_primary == TRUE]
enriched_pruning <- prune_pathways(
  copy(enriched),
  membership_list,
  threshold = jaccard_threshold
)

results <- merge(
  results,
  enriched_pruning[
    ,
    .(
      pathway_key,
      redundancy_cluster,
      redundancy_cluster_size,
      representative_pathway_key,
      retained_after_pruning
    )
  ],
  by = "pathway_key",
  all.x = TRUE,
  sort = FALSE
)
results[is.na(retained_after_pruning), retained_after_pruning := FALSE]

results[
  retained_after_pruning == TRUE,
  direction_q_pruned := p.adjust(direction_fisher_p_one_sided, method = "BH")
]
results[
  retained_after_pruning == TRUE,
  independent_direction_q_pruned :=
    p.adjust(independent_direction_fisher_p_one_sided, method = "BH")
]
results[
  retained_after_pruning == TRUE,
  asymmetric_cnv_q_pruned := p.adjust(asymmetric_cnv_p, method = "BH")
]
results[
  retained_after_pruning == TRUE,
  mirror_cnv_q_pruned := p.adjust(mirror_cnv_p, method = "BH")
]
results[
  retained_after_pruning == TRUE,
  asymmetric_vs_mirror_q_pruned :=
    p.adjust(asymmetric_vs_mirror_p_one_sided, method = "BH")
]

# ----------------------------
# 4. Locus-block permutation for retained exact pathways
# ----------------------------

loci <- sort(unique(cnv$CNV_ID))
locus_class <- setNames(
  vapply(
    loci,
    function(id) cnv[CNV_ID == id, Dose_response_class][1],
    character(1)
  ),
  loci
)
locus_genes <- setNames(
  lapply(loci, function(id) unique(cnv[CNV_ID == id, Gene_std])),
  loci
)
locus_sizes <- vapply(locus_genes, length, integer(1))

# Gene-count tertiles; preserve the observed number of asymmetric loci per stratum.
breaks <- unique(quantile(locus_sizes, probs = c(0, 1/3, 2/3, 1), type = 1))
if (length(breaks) < 3) {
  strata <- rep(1L, length(loci))
} else {
  strata <- cut(
    locus_sizes,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
}

set.seed(random_seed)
retained_keys <- results[retained_after_pruning == TRUE, pathway_key]
perm_results <- vector("list", length(retained_keys))

for (path_i in seq_along(retained_keys)) {
  key <- retained_keys[path_i]
  row <- results[pathway_key == key][1]
  geneset <- membership_list[[key]]
  lib_universe <- library_universes[[row$library]]

  observed_asym_genes <- unique(unlist(locus_genes[locus_class == "asymmetric"]))
  observed_mirror_genes <- unique(unlist(locus_genes[locus_class == "mirror"]))
  observed_asym_genes <- intersect(observed_asym_genes, lib_universe)
  observed_mirror_genes <- intersect(observed_mirror_genes, lib_universe)

  observed_difference <-
    length(intersect(geneset, observed_asym_genes)) /
      max(length(observed_asym_genes), 1) -
    length(intersect(geneset, observed_mirror_genes)) /
      max(length(observed_mirror_genes), 1)

  exceed <- 0L
  null_values <- numeric(n_locus_permutations)

  for (perm_i in seq_len(n_locus_permutations)) {
    perm_asym <- character()

    for (stratum in sort(unique(strata))) {
      ids <- loci[strata == stratum]
      observed_n <- sum(locus_class[ids] == "asymmetric")
      if (observed_n > 0) {
        perm_asym <- c(perm_asym, sample(ids, observed_n, replace = FALSE))
      }
    }

    perm_mirror <- setdiff(loci, perm_asym)
    perm_asym_genes <- intersect(
      unique(unlist(locus_genes[perm_asym])),
      lib_universe
    )
    perm_mirror_genes <- intersect(
      unique(unlist(locus_genes[perm_mirror])),
      lib_universe
    )

    value <-
      length(intersect(geneset, perm_asym_genes)) /
        max(length(perm_asym_genes), 1) -
      length(intersect(geneset, perm_mirror_genes)) /
        max(length(perm_mirror_genes), 1)

    null_values[perm_i] <- value
    if (value >= observed_difference - 1e-15) exceed <- exceed + 1L
  }

  perm_results[[path_i]] <- data.table(
    pathway_key = key,
    cnv_locus_perm_observed_fraction_difference = observed_difference,
    cnv_locus_perm_p_asymmetric_greater =
      (exceed + 1) / (n_locus_permutations + 1),
    cnv_locus_perm_null_mean = mean(null_values),
    cnv_locus_perm_null_sd = sd(null_values)
  )
}

perm_results <- rbindlist(perm_results, fill = TRUE)
perm_results[
  ,
  cnv_locus_perm_q_pruned :=
    p.adjust(cnv_locus_perm_p_asymmetric_greater, method = "BH")
]
results <- merge(results, perm_results, by = "pathway_key", all.x = TRUE, sort = FALSE)

# ----------------------------
# 5. Evidence definitions
# ----------------------------

results[
  ,
  negative_bias_nominal :=
    kosmicki_enriched_primary &
    direction_odds_ratio > 1 &
    direction_fisher_p_one_sided < 0.05
]
results[
  ,
  asymmetric_enrichment_nominal :=
    kosmicki_enriched_primary &
    asymmetric_cnv_odds_ratio > 1 &
    asymmetric_cnv_p < 0.05
]
results[
  ,
  asymmetric_specific_locus_nominal :=
    kosmicki_enriched_primary &
    cnv_locus_perm_observed_fraction_difference > 0 &
    cnv_locus_perm_p_asymmetric_greater < 0.05
]
results[
  ,
  triple_evidence_nominal_gene_level :=
    kosmicki_enriched_primary &
    negative_bias_nominal &
    asymmetric_enrichment_nominal
]
results[
  ,
  triple_evidence_nominal_locus_specific :=
    kosmicki_enriched_primary &
    negative_bias_nominal &
    asymmetric_specific_locus_nominal
]
results[
  ,
  triple_evidence_strict_pruned :=
    retained_after_pruning &
    direction_q_pruned < 0.05 &
    asymmetric_cnv_q_pruned < 0.05
]
results[
  ,
  triple_evidence_strict_locus_pruned :=
    retained_after_pruning &
    direction_q_pruned < 0.05 &
    cnv_locus_perm_q_pruned < 0.05
]

results[
  ,
  evidence_count_nominal :=
    as.integer(kosmicki_enriched_primary) +
    as.integer(negative_bias_nominal) +
    as.integer(asymmetric_enrichment_nominal)
]
results[
  ,
  priority_score :=
    -log10(pmax(kosmicki_q_within_library, 1e-300)) -
    log10(pmax(direction_fisher_p_one_sided, 1e-300)) -
    log10(pmax(asymmetric_cnv_p, 1e-300))
]

setorder(
  results,
  -triple_evidence_strict_pruned,
  -triple_evidence_nominal_gene_level,
  -evidence_count_nominal,
  kosmicki_q_within_library,
  kosmicki_p
)

# Add gene lists for the Kosmicki-enriched pathways only.
results[, `:=`(
  negative_genes = NA_character_,
  positive_genes = NA_character_,
  asymmetric_cnv_overlap_genes = NA_character_,
  mirror_cnv_overlap_genes = NA_character_
)]

for (key in results[kosmicki_enriched_primary == TRUE, pathway_key]) {
  geneset <- membership_list[[key]]
  results[
    pathway_key == key,
    negative_genes := paste(sort(intersect(geneset, all_negative)), collapse = ";")
  ]
  results[
    pathway_key == key,
    positive_genes :=
      paste(sort(intersect(geneset, setdiff(all_kosmicki, all_negative))), collapse = ";")
  ]
  results[
    pathway_key == key,
    asymmetric_cnv_overlap_genes :=
      paste(sort(intersect(geneset, asymmetric_raw)), collapse = ";")
  ]
  results[
    pathway_key == key,
    mirror_cnv_overlap_genes :=
      paste(sort(intersect(geneset, mirror_raw)), collapse = ";")
  ]
}

# ----------------------------
# Outputs
# ----------------------------

fwrite(
  results,
  file.path(output_dir, "kosmicki_pathway_cnv_master_results.csv")
)
fwrite(
  results[kosmicki_enriched_primary == TRUE],
  file.path(output_dir, "kosmicki_enriched_pathways_full_results.csv")
)
fwrite(
  results[retained_after_pruning == TRUE],
  file.path(output_dir, "kosmicki_enriched_pathways_pruned_results.csv")
)
fwrite(
  results[triple_evidence_nominal_gene_level == TRUE],
  file.path(output_dir, "kosmicki_three_way_candidates.csv")
)

membership_rows <- rbindlist(
  lapply(
    results[kosmicki_enriched_primary == TRUE, pathway_key],
    function(key) {
      row <- results[pathway_key == key][1]
      genes <- membership_list[[key]]
      data.table(
        pathway_key = key,
        library = row$library,
        pathway = row$pathway,
        gene = genes,
        is_kosmicki = genes %in% all_kosmicki,
        kosmicki_effect_direction =
          kos$Effect_direction[match(genes, kos$Gene_std)],
        in_asymmetric_cnv = genes %in% asymmetric_raw,
        in_mirror_cnv = genes %in% mirror_raw,
        retained_after_pruning = row$retained_after_pruning
      )
    }
  ),
  fill = TRUE
)
fwrite(
  membership_rows,
  file.path(output_dir, "kosmicki_enriched_pathway_gene_memberships.csv")
)

summary <- data.table(
  metric = c(
    "Kosmicki genes",
    "Kosmicki negative genes",
    "Kosmicki positive genes",
    "GO BP pathways tested",
    "KEGG pathways tested",
    "Kosmicki-enriched pathways",
    "Retained pathways after Jaccard pruning",
    "Nominal three-way candidates",
    "Strict FDR three-way candidates",
    "Asymmetric CNV loci",
    "Mirror CNV loci"
  ),
  value = c(
    nrow(kos),
    sum(kos$Negative),
    sum(!kos$Negative),
    sum(results$library == "GO_Biological_Process_2023"),
    sum(results$library == "KEGG_2021_Human"),
    sum(results$kosmicki_enriched_primary),
    sum(results$retained_after_pruning),
    sum(results$triple_evidence_nominal_gene_level),
    sum(results$triple_evidence_strict_pruned),
    uniqueN(cnv[Dose_response_class == "asymmetric", CNV_ID]),
    uniqueN(cnv[Dose_response_class == "mirror", CNV_ID])
  )
)
fwrite(summary, file.path(output_dir, "analysis_summary.csv"))

message("Analysis complete. Results written to: ", normalizePath(output_dir))
