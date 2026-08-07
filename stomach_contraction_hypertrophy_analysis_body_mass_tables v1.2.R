# Stomach contraction, relaxation, AUC, and GAMM analysis
# Author: Jack Auty
# Version: 4.0
#
# Naming convention used throughout:
#   bird_body_mass                  = whole-bird body mass
#   proventriculus_tissue_mass      = proventriculus/stomach tissue mass
#   proventriculus_plastic_mass     = mass of plastic in the proventriculus
#   total_plastic_mass              = total plastic mass across sampled gut regions
#
# This naming is deliberate. The original source column "Stomach_mass" in
# masses.csv was used as PLASTIC mass, so that ambiguous name is removed here.
#
# WC, HB, and Cul are preserved exactly as written in Bird_Body_Morphometrics.csv
# because their full names/units are not specified in the source file.

options(scipen = 999)

# ============================================================
# PACKAGES
# ============================================================

packages <- c(
  "ggplot2",
  "cowplot",
  "viridisLite",
  "sjPlot",
  "grid",
  "usethis",
  "mgcv"
)

invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))
  pdf(
    file = file.path(
      DATA_DIR,
      "proventriculus_plastic_mass_distribution_relaxation.pdf"
    ),
    width = 7,
    height = 5
  )

  h <- hist(
    relaxation_bird_check$proventriculus_plastic_mass,
    breaks = c(0, 0.5, 1, 5, 10, 60),
    plot = FALSE
  )
  plot(
    h,
    log = "x",
    main = "Relaxation birds: proventriculus plastic mass",
    xlab = "Proventriculus plastic mass"
  )

  dev.off()
  txt <- gsub("[_*\"“”]", "", txt)
  txt <- gsub(", \\.$|\\s*\\.$", "", txt)
  trimws(txt)
}

for (i in pkg) {
  raw <- paste(
    capture.output(print(citation(i), style = "text", bibtex = FALSE)),
    collapse = "\n"
  )
  citation_list <- c(citation_list, clean_citation(raw))
}

sjPlot::tab_df(
  data.frame(References = citation_list),
  col.header = "References",
  sep = " "
)


# ============================================================
# LOAD DATA
# ============================================================

DATA_DIR <- usethis::proj_get()

ts_path           <- file.path(DATA_DIR, "final_trim.csv")
map_path          <- file.path(DATA_DIR, "video_map.csv")
mass_path         <- file.path(DATA_DIR, "masses.csv")
morphometrics_path <- file.path(DATA_DIR, "Bird_Body_Morphometrics.csv")

required_files <- c(
  ts_path,
  map_path,
  mass_path,
  morphometrics_path
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files)) {
  stop(
    "The following required data files were not found: ",
    paste(basename(missing_files), collapse = ", ")
  )
}

ts_df <- read.csv(
  ts_path,
  check.names = TRUE
)

map_df <- read.csv(
  map_path,
  check.names = TRUE,
  strip.white = TRUE
)

mass_df <- read.csv(
  mass_path,
  check.names = TRUE,
  strip.white = TRUE
)

morphometrics_df <- read.csv(
  morphometrics_path,
  check.names = TRUE,
  strip.white = TRUE
)


# ============================================================
# CHECK AND RENAME MORPHOMETRIC COLUMNS
# ============================================================
#
# Bird_Body_Morphometrics.csv is expected to contain:
#   ID
#   Bird_Body_Weight
#   WC
#   HB
#   Cul
#   Stom_tiss_mass
#   Prov_plastic_mass
#
# Only the ambiguous mass columns are renamed.
# WC, HB, and Cul are retained exactly as supplied.

required_morphometric_columns <- c(
  "ID",
  "Bird_Body_Weight",
  "WC",
  "HB",
  "Cul",
  "Stom_tiss_mass",
  "Prov_plastic_mass"
)

missing_morphometric_columns <- setdiff(
  required_morphometric_columns,
  names(morphometrics_df)
)

if (length(missing_morphometric_columns)) {
  stop(
    "Bird_Body_Morphometrics.csv is missing required column(s): ",
    paste(missing_morphometric_columns, collapse = ", ")
  )
}

names(morphometrics_df)[
  names(morphometrics_df) == "Bird_Body_Weight"
] <- "bird_body_mass"

names(morphometrics_df)[
  names(morphometrics_df) == "Stom_tiss_mass"
] <- "proventriculus_tissue_mass"

names(morphometrics_df)[
  names(morphometrics_df) == "Prov_plastic_mass"
] <- "proventriculus_plastic_mass"

# Make the three mass variables explicitly numeric.
morphometrics_df$bird_body_mass <- suppressWarnings(
  as.numeric(morphometrics_df$bird_body_mass)
)

morphometrics_df$proventriculus_tissue_mass <- suppressWarnings(
  as.numeric(morphometrics_df$proventriculus_tissue_mass)
)

morphometrics_df$proventriculus_plastic_mass <- suppressWarnings(
  as.numeric(morphometrics_df$proventriculus_plastic_mass)
)


# ============================================================
# NORMALISE BIRD IDS
# ============================================================
#
# This is one of the few places where a helper function is useful because
# exactly the same ID rule must be applied to three separate source files.

names(map_df) <- sub(
  "BIRD[._-]?ID",
  "BIRD_ID",
  names(map_df),
  ignore.case = TRUE
)

normalise_bird_key <- function(x) {

  x <- trimws(as.character(x))
  out <- x

  # Example: NL-2025-10 -> NL-10
  m <- grepl("^[A-Z]+-20[0-9]{2}-[0-9]+$", x)
  out[m] <- sub(
    "^([A-Z]+)-20[0-9]{2}-([0-9]+)$",
    "\\1-\\2",
    x[m]
  )

  # Example: 380-02015 -> 02015
  m <- grepl("^[0-9]{3}-[0-9]+$", out)
  out[m] <- sub(
    "^[0-9]{3}-([0-9]+)$",
    "\\1",
    out[m]
  )

  # Example: BW-1 -> BW-01
  out <- sub(
    "^BW-([0-9])$",
    "BW-0\\1",
    out
  )

  out
}

map_df$BIRD_KEY <- normalise_bird_key(map_df$BIRD_ID)
mass_df$BIRD_KEY <- normalise_bird_key(mass_df$ID)
morphometrics_df$BIRD_KEY <- normalise_bird_key(morphometrics_df$ID)


# ============================================================
# CLEAN masses.csv NAMES
# ============================================================
#
# IMPORTANT:
# In the original script, "Stomach_mass" was used as the amount of PLASTIC
# in the proventriculus, whereas "Wt" was the stomach/proventriculus tissue
# weight. Those names are too easy to confuse, so they are made explicit here.

# Collapse repeated mass rows first, as in the original analysis.
is_numeric_mass_column <- sapply(mass_df, is.numeric)

sum_keep_na <- function(x) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    sum(x, na.rm = TRUE)
  }
}

mass_df <- aggregate(
  mass_df[, is_numeric_mass_column, drop = FALSE],
  by = list(
    ID = mass_df$ID,
    BIRD_KEY = mass_df$BIRD_KEY
  ),
  FUN = sum_keep_na
)

# If two different source IDs collapse to the same normalised BIRD_KEY,
# flag that explicitly before any cross-file comparison.
duplicate_mass_keys <- unique(
  mass_df$BIRD_KEY[
    duplicated(mass_df$BIRD_KEY) |
      duplicated(mass_df$BIRD_KEY, fromLast = TRUE)
  ]
)

if (length(duplicate_mass_keys)) {
  warning(
    "More than one masses.csv row remains for these normalised BIRD_KEY values: ",
    paste(duplicate_mass_keys, collapse = ", "),
    ". Inspect masses.csv before relying on the source-comparison table."
  )
}

if ("Stomach_mass" %in% names(mass_df)) {
  names(mass_df)[names(mass_df) == "Stomach_mass"] <-
    "proventriculus_plastic_mass"
}

if ("Wt" %in% names(mass_df)) {
  names(mass_df)[names(mass_df) == "Wt"] <-
    "proventriculus_tissue_mass"
}

if ("Stomach_number" %in% names(mass_df)) {
  names(mass_df)[names(mass_df) == "Stomach_number"] <-
    "proventriculus_plastic_count"
}

if ("Total_mass" %in% names(mass_df)) {
  names(mass_df)[names(mass_df) == "Total_mass"] <-
    "total_plastic_mass"
}

if ("Total_number" %in% names(mass_df)) {
  names(mass_df)[names(mass_df) == "Total_number"] <-
    "total_plastic_count"
}


# ============================================================
# ONE-ROW-PER-BIRD MORPHOMETRIC TABLE
# ============================================================
#
# Unlike masses.csv, morphometrics should already contain one biological row
# per bird. We DO NOT silently average duplicate morphometric rows.
# If normalisation creates a duplicate BIRD_KEY, the script stops so it can
# be inspected rather than hidden.

duplicate_morphometric_keys <- unique(
  morphometrics_df$BIRD_KEY[
    duplicated(morphometrics_df$BIRD_KEY) |
      duplicated(morphometrics_df$BIRD_KEY, fromLast = TRUE)
  ]
)

if (length(duplicate_morphometric_keys)) {
  stop(
    "More than one row in Bird_Body_Morphometrics.csv maps to these BIRD_KEY values: ",
    paste(duplicate_morphometric_keys, collapse = ", "),
    ". Inspect these rows before modelling."
  )
}

bird_morphometrics <- morphometrics_df[, c(
  "ID",
  "BIRD_KEY",
  "bird_body_mass",
  "WC",
  "HB",
  "Cul",
  "proventriculus_tissue_mass",
  "proventriculus_plastic_mass"
)]

names(bird_morphometrics)[
  names(bird_morphometrics) == "ID"
] <- "morphometrics_ID"

bird_morphometrics <- bird_morphometrics[
  order(bird_morphometrics$BIRD_KEY),
]


# ============================================================
# CHECK PLASTIC/TISSUE MASS ACROSS THE TWO SOURCE FILES
# ============================================================
#
# masses.csv and Bird_Body_Morphometrics.csv both contain measures relating
# to proventriculus plastic/tissue mass. The analysis below uses the clearly
# named columns from Bird_Body_Morphometrics.csv as the bird-level metadata.
#
# This table is exported so the two sources can be checked side by side.

mass_source_check <- mass_df[, intersect(
  c(
    "BIRD_KEY",
    "proventriculus_plastic_mass",
    "proventriculus_tissue_mass"
  ),
  names(mass_df)
)]

if ("proventriculus_plastic_mass" %in% names(mass_source_check)) {
  names(mass_source_check)[
    names(mass_source_check) == "proventriculus_plastic_mass"
  ] <- "proventriculus_plastic_mass_masses_csv"
}

if ("proventriculus_tissue_mass" %in% names(mass_source_check)) {
  names(mass_source_check)[
    names(mass_source_check) == "proventriculus_tissue_mass"
  ] <- "proventriculus_tissue_mass_masses_csv"
}

mass_alignment_check <- merge(
  bird_morphometrics,
  mass_source_check,
  by = "BIRD_KEY",
  all = TRUE
)

if (
  all(c(
    "proventriculus_plastic_mass",
    "proventriculus_plastic_mass_masses_csv"
  ) %in% names(mass_alignment_check))
) {
  mass_alignment_check$plastic_mass_difference <-
    mass_alignment_check$proventriculus_plastic_mass -
    mass_alignment_check$proventriculus_plastic_mass_masses_csv
}

if (
  all(c(
    "proventriculus_tissue_mass",
    "proventriculus_tissue_mass_masses_csv"
  ) %in% names(mass_alignment_check))
) {
  mass_alignment_check$tissue_mass_difference <-
    mass_alignment_check$proventriculus_tissue_mass -
    mass_alignment_check$proventriculus_tissue_mass_masses_csv
}

cat("\n============================================================\n")
cat("PROVENTRICULUS MASS ALIGNMENT CHECK\n")
cat("============================================================\n")
print(mass_alignment_check)

write.csv(
  mass_alignment_check,
  "proventriculus_mass_source_alignment_check.csv",
  row.names = FALSE
)


# ============================================================
# HELPERS FOR THE TIME-SERIES DATA
# ============================================================

# Convert the wide force traces to long format.
make_long <- function(df, suffix) {

  keep <- names(df)[
    grepl(
      paste0("(", suffix, ")$"),
      names(df),
      ignore.case = TRUE
    )
  ]

  long_df <- stack(df[keep])
  long_df$time <- rep.int(
    seq_len(nrow(df)),
    times = length(keep)
  )

  names(long_df) <- c(
    "values",
    "ind",
    "time"
  )

  long_df$VIDEO <- sub(
    paste0("_(?i:", suffix, ")$"),
    "",
    long_df$ind,
    perl = TRUE
  )

  long_df
}


# Baseline correction is kept as a helper because it must be applied
# identically to every video.
baseline_correct <- function(df) {

  df <- df[
    order(df$VIDEO, df$time),
  ]

  df$baseline0 <- NA_real_

  for (video in unique(df$VIDEO)) {

    idx <- which(df$VIDEO == video)
    video_data <- df[idx, ]

    baseline_index <- which(
      video_data$time == 0 &
        !is.na(video_data$values)
    )

    if (!length(baseline_index)) {
      baseline_index <- which(
        !is.na(video_data$values)
      )[1]
    }

    if (length(baseline_index)) {
      df$baseline0[idx] <-
        video_data$values[baseline_index[1]]
    }
  }

  df$value_delta <- df$values - df$baseline0

  df
}


# Used only for the raw mean ± SEM plots.
mean_sem <- function(x) {

  mean_value <- mean(x, na.rm = TRUE)
  sd_value <- stats::sd(x, na.rm = TRUE)
  n_value <- sum(!is.na(x))

  c(
    mean = mean_value,
    sd = sd_value,
    n = n_value
  )
}


# Trapezoidal integration is a distinct operation and is therefore retained
# as a small helper. Contraction and relaxation AUC calculations themselves
# are written out separately below.
trapz_auc <- function(time, force) {

  ok <- is.finite(time) & is.finite(force)

  time <- time[ok]
  force <- force[ok]

  if (length(time) < 2) {
    return(NA_real_)
  }

  order_index <- order(time)
  time <- time[order_index]
  force <- force[order_index]

  sum(
    (force[-length(force)] + force[-1]) *
      0.5 *
      diff(time),
    na.rm = TRUE
  )
}


# ============================================================
# PLOT SETTINGS
# ============================================================

pal_group <- c(
  "Low plastic (≤0.5)" = "#0072B2",
  "High plastic (>0.5)" = "#D55E00"
)

na_col <- "#999999"

theme_small <- theme_minimal(base_size = 7) +
  theme(
    panel.grid.minor = element_blank()
  )

# 120 frames = 60 minutes because the trace is sampled twice per minute.
t_max <- 120


# ============================================================
# BUILD LONG TIME-SERIES DATA AND ATTACH BIRD METADATA
# ============================================================

raw_long <- make_long(
  ts_df,
  "raw"
)

sm_long <- make_long(
  ts_df,
  "smoothed"
)

map_keep <- c(
  "VIDEO",
  "BIRD_ID",
  "BIRD_KEY",
  "PHASE"
)

raw_m <- merge(
  raw_long,
  map_df[map_keep],
  by = "VIDEO",
  all.x = TRUE
)

sm_m <- merge(
  sm_long,
  map_df[map_keep],
  by = "VIDEO",
  all.x = TRUE
)

# Attach the same clearly named bird-level morphometric data to every trace.
merged_raw <- merge(
  raw_m,
  bird_morphometrics,
  by = "BIRD_KEY",
  all.x = TRUE
)

merged_sm <- merge(
  sm_m,
  bird_morphometrics,
  by = "BIRD_KEY",
  all.x = TRUE
)

cat("\n============================================================\n")
cat("TIME-SERIES / MORPHOMETRIC ID CHECK\n")
cat("============================================================\n")

time_series_bird_check <- unique(
  merged_sm[, c(
    "BIRD_KEY",
    "BIRD_ID",
    "morphometrics_ID",
    "bird_body_mass",
    "WC",
    "HB",
    "Cul",
    "proventriculus_tissue_mass",
    "proventriculus_plastic_mass"
  )]
)

time_series_bird_check <- time_series_bird_check[
  order(time_series_bird_check$BIRD_KEY),
]

print(time_series_bird_check)

write.csv(
  time_series_bird_check,
  "time_series_bird_morphometric_alignment_check.csv",
  row.names = FALSE
)

missing_morphometric_ids <- unique(
  merged_sm$BIRD_KEY[
    is.na(merged_sm$morphometrics_ID)
  ]
)

if (length(missing_morphometric_ids)) {
  warning(
    "No Bird_Body_Morphometrics.csv row was matched for: ",
    paste(missing_morphometric_ids, collapse = ", ")
  )
}


# ============================================================
# PREPARE CONTRACTION DATA
# ============================================================

sm_contract <- subset(
  merged_sm,
  grepl(
    "^Contract",
    PHASE,
    ignore.case = TRUE
  )
)

sm_contract$VIDEO <- factor(
  sm_contract$VIDEO,
  levels = unique(sm_contract$VIDEO)
)

sm_contract <- baseline_correct(
  sm_contract
)

sm_contract <- subset(
  sm_contract,
  time <= t_max
)

sm_contract$plastic_group <- factor(
  ifelse(
    is.na(sm_contract$proventriculus_plastic_mass),
    NA,
    ifelse(
      sm_contract$proventriculus_plastic_mass > 0.5,
      "High plastic (>0.5)",
      "Low plastic (≤0.5)"
    )
  ),
  levels = names(pal_group)
)


# ============================================================
# PREPARE RELAXATION DATA
# ============================================================

relax_sm <- subset(
  merged_sm,
  grepl(
    "^(Relax|Relaxation)",
    PHASE,
    ignore.case = TRUE
  )
)

relax_sm$VIDEO <- factor(
  relax_sm$VIDEO,
  levels = unique(relax_sm$VIDEO)
)

relax_sm <- baseline_correct(
  relax_sm
)

relax_sm <- subset(
  relax_sm,
  time <= t_max
)

relax_sm$plastic_group <- factor(
  ifelse(
    is.na(relax_sm$proventriculus_plastic_mass),
    NA,
    ifelse(
      relax_sm$proventriculus_plastic_mass > 0.5,
      "High plastic (>0.5)",
      "Low plastic (≤0.5)"
    )
  ),
  levels = names(pal_group)
)



# Figure 1: proventriculus tissue mass vs proventriculus plastic mass, plus contraction mean ± SEM

# Scatter (A)
df_scatter <- subset(bird_morphometrics, is.finite(proventriculus_plastic_mass) & is.finite(proventriculus_tissue_mass))
pear <- suppressWarnings(cor(df_scatter$proventriculus_plastic_mass, df_scatter$proventriculus_tissue_mass, use = "complete.obs", method = "pearson"))
pA <- ggplot(df_scatter, aes(x = proventriculus_plastic_mass, y = proventriculus_tissue_mass)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, colour = "grey55", fill = "grey85", na.rm = TRUE) +
  geom_point(shape = 21, size = 2.4, stroke = 0.3, colour = "grey25",
             fill = scales::alpha("#0072B2", 0.35), na.rm = TRUE) +
  labs(x = "Proventriculus plastic mass", y = "Proventriculus tissue mass",
       subtitle = sprintf("Linear fit • Pearson r = %.2f", pear)) +
  theme_classic(base_size = 11)

# Contract (smoothed) mean ± SEM (B)
# sm_contract was prepared explicitly above.

agg_s <- aggregate(value_delta ~ time + plastic_group,
                   data = subset(sm_contract, !is.na(plastic_group)),
                   FUN = mean_sem)
agg_s <- data.frame(time = agg_s$time,
                    plastic_group = agg_s$plastic_group,
                    mean = agg_s$value_delta[, "mean"],
                    sd   = agg_s$value_delta[, "sd"],
                    n    = agg_s$value_delta[, "n"])
agg_s$sem   <- with(agg_s, sd / sqrt(pmax(n, 1)))
agg_s$lower <- agg_s$mean - agg_s$sem
agg_s$upper <- agg_s$mean + agg_s$sem
agg_s$time_min <- agg_s$time / 2

pB_core <- ggplot(agg_s, aes(x = time_min, y = mean, colour = plastic_group, fill = plastic_group)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = pal_group, guide = guide_legend(title = "Plastic load")) +
  scale_fill_manual(values = pal_group, guide = "none") +
  labs(x = "Time (minutes)", y = "Mean contraction force (Δ g)") +
  coord_cartesian(xlim = c(0, t_max/2)) +
  theme_classic(base_size = 11)

leg_B <- cowplot::get_legend(pB_core + theme(legend.position = "bottom"))
pB    <- cowplot::plot_grid(pB_core + theme(legend.position = "none"),
                            leg_B, ncol = 1, rel_heights = c(1, 0.14))

# Combine A and B
fig1 <- cowplot::plot_grid(pA + labs(subtitle = NULL), pB,
                           labels = c("A","B"),
                           label_fontface = "bold", label_size = 14,
                           nrow = 1, rel_widths = c(1, 1.25))


print(fig1)

ggsave("Figure 1.pdf", fig1,
       width = 270, height = 95, units = "mm")

# Figure 2: ABCD square – A scatter, B contract, C relax, D legend (270 × 180 mm)

# A
dfA <- transform(df_scatter,
                 plastic_group = factor(ifelse(proventriculus_plastic_mass > 0.5, "High plastic (>0.5)", "Low plastic (≤0.5)"),
                                        levels = names(pal_group)))
pA2 <- ggplot(dfA, aes(proventriculus_plastic_mass, proventriculus_tissue_mass, colour = plastic_group, fill = plastic_group)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40", fill = scales::alpha("grey40", 0.125), linewidth = 0.5) +
  geom_point(shape = 21, size = 2.4, stroke = 0.35, alpha = 0.4) +
  scale_colour_manual(values = pal_group, name = "Plastic load") +
  scale_fill_manual(values = pal_group, guide = "none") +
  guides(colour = "none") +
  labs(x = "Proventriculus plastic mass", y = "Proventriculus tissue mass") +
  theme_classic(base_size = 11)

# B
agg_s2 <- transform(agg_s, time_min = time/2)
pB2 <- ggplot(agg_s2, aes(x = time_min, y = mean, colour = plastic_group, fill = plastic_group)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = pal_group, name = "Plastic load") +
  scale_fill_manual(values = pal_group, guide = "none") +
  labs(x = "Time (minutes)", y = "Contraction force (Δ g)") +
  coord_cartesian(xlim = c(0, t_max/2)) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")

# C
# relax_sm was prepared explicitly above.

agg_relax <- aggregate(value_delta ~ time + plastic_group,
                       data = subset(relax_sm, !is.na(plastic_group)),
                       FUN = mean_sem)
agg_relax <- data.frame(time = agg_relax$time,
                        plastic_group = agg_relax$plastic_group,
                        mean = agg_relax$value_delta[, "mean"],
                        sd   = agg_relax$value_delta[, "sd"],
                        n    = agg_relax$value_delta[, "n"])
agg_relax$sem   <- with(agg_relax, sd / sqrt(pmax(n, 1)))
agg_relax$lower <- agg_relax$mean - agg_relax$sem
agg_relax$upper <- agg_relax$mean + agg_relax$sem
agg_relax2 <- transform(agg_relax, time_min = time / 2)

pC2 <- ggplot(agg_relax2, aes(x = time_min, y = mean, colour = plastic_group, fill = plastic_group)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = pal_group, name = "Plastic load") +
  scale_fill_manual(values = pal_group, guide = "none") +
  labs(x = "Time (minutes)", y = "Relaxation (Δ g)") +
  coord_cartesian(xlim = c(0, t_max/2)) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")

# D: centred legend
legend_df <- data.frame(plastic_group = factor(names(pal_group), levels = names(pal_group)), x = 1, y = 1)
pD <- ggplot(legend_df, aes(x, y, colour = plastic_group)) +
  geom_line(size = 1.1, alpha = 0) +
  scale_colour_manual(values = pal_group, name = "Plastic load") +
  guides(colour = guide_legend(title.position = "top", nrow = 1, override.aes = list(alpha = 1, size = 1.1))) +
  theme_void(base_size = 11) +
  theme(legend.position = c(0.5, 0.5),
        legend.justification = c(0.5, 0.5),
        legend.direction = "horizontal",
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10))

top_row <- cowplot::plot_grid(pA2, pB2, ncol = 2, labels = c("A","B"), label_fontface = "bold", label_size = 14)
bot_row <- cowplot::plot_grid(pC2, pD,  ncol = 2, labels = c("C","D"), label_fontface = "bold", label_size = 14)
fig2 <- cowplot::plot_grid(top_row, bot_row, nrow = 2, rel_heights = c(1, 1))
print(fig2)

ggsave("Figure 2.pdf", fig2,
       width = 270, height = 180, units = "mm")




# ============================================================
# EXPORT CLEAN BIRD-LEVEL MASS DATA FROM masses.csv
# ============================================================

mass_export_columns <- intersect(
  c(
    "BIRD_KEY",
    "proventriculus_plastic_mass",
    "proventriculus_plastic_count",
    "total_plastic_mass",
    "total_plastic_count",
    "proventriculus_tissue_mass"
  ),
  names(mass_df)
)

mass_df_export <- mass_df[
  ,
  mass_export_columns,
  drop = FALSE
]

write.csv(
  mass_df_export,
  "stomach_hypertrophy_and_plastic_mass.csv",
  row.names = FALSE
)


# ============================================================
# EXPORT CONTRACTION TRACE AT 20 MINUTES
# ============================================================

contraction_data_20minutes <- sm_contract[
  sm_contract$time == 120,
]

write.csv(
  contraction_data_20minutes,
  "stomach_contraction_20mins.csv",
  row.names = FALSE
)


# ============================================================
# BIRD-LEVEL AUC ANALYSIS
# ============================================================
#
# The AUC calculation is deliberately written out separately for
# contraction and relaxation.
#
# bird_auc becomes the central ONE-ROW-PER-BIRD analysis/audit dataframe.
# It contains:
#   - original bird IDs
#   - all requested morphometric fields
#   - bird body mass
#   - proventriculus tissue mass
#   - proventriculus plastic mass
#   - plastic group
#   - contraction AUC
#   - relaxation AUC


# ------------------------------------------------------------
# CONTRACTION AUC
# ------------------------------------------------------------

contraction_auc_time <- sm_contract[, c(
  "BIRD_KEY",
  "time",
  "value_delta"
)]

contraction_auc_time$time_min <-
  contraction_auc_time$time / 2

# If more than one contraction trace exists for a bird, average the force
# at each time point before integrating so each bird contributes one curve.
contraction_auc_time <- aggregate(
  value_delta ~ BIRD_KEY + time_min,
  data = contraction_auc_time,
  FUN = mean,
  na.rm = TRUE
)

contraction_bird_keys <- sort(
  unique(contraction_auc_time$BIRD_KEY)
)

contraction_auc <- data.frame(
  BIRD_KEY = contraction_bird_keys,
  auc_contraction = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(contraction_bird_keys)) {

  this_bird <- contraction_bird_keys[i]

  this_curve <- contraction_auc_time[
    contraction_auc_time$BIRD_KEY == this_bird,
  ]

  contraction_auc$auc_contraction[i] <- trapz_auc(
    this_curve$time_min,
    this_curve$value_delta
  )
}


# ------------------------------------------------------------
# RELAXATION AUC
# ------------------------------------------------------------

relaxation_auc_time <- relax_sm[, c(
  "BIRD_KEY",
  "time",
  "value_delta"
)]

relaxation_auc_time$time_min <-
  relaxation_auc_time$time / 2

# Again, each bird contributes one mean curve before integration.
relaxation_auc_time <- aggregate(
  value_delta ~ BIRD_KEY + time_min,
  data = relaxation_auc_time,
  FUN = mean,
  na.rm = TRUE
)

relaxation_bird_keys <- sort(
  unique(relaxation_auc_time$BIRD_KEY)
)

relaxation_auc <- data.frame(
  BIRD_KEY = relaxation_bird_keys,
  auc_relaxation = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(relaxation_bird_keys)) {

  this_bird <- relaxation_bird_keys[i]

  this_curve <- relaxation_auc_time[
    relaxation_auc_time$BIRD_KEY == this_bird,
  ]

  relaxation_auc$auc_relaxation[i] <- trapz_auc(
    this_curve$time_min,
    this_curve$value_delta
  )
}


# ------------------------------------------------------------
# BIRD ID LOOKUP
# ------------------------------------------------------------

bird_id_lookup <- unique(
  map_df[, c(
    "BIRD_KEY",
    "BIRD_ID"
  )]
)

# A normalised BIRD_KEY should correspond to only one original BIRD_ID.
duplicate_map_keys <- unique(
  bird_id_lookup$BIRD_KEY[
    duplicated(bird_id_lookup$BIRD_KEY) |
      duplicated(bird_id_lookup$BIRD_KEY, fromLast = TRUE)
  ]
)

if (length(duplicate_map_keys)) {
  stop(
    "These BIRD_KEY values map to more than one BIRD_ID in video_map.csv: ",
    paste(duplicate_map_keys, collapse = ", "),
    ". Inspect the ID mapping before continuing."
  )
}


# ------------------------------------------------------------
# BUILD bird_auc
# ------------------------------------------------------------
#
# Start with experimental birds from video_map.csv.
# Then attach the morphometrics and both independently calculated AUCs.

bird_auc <- merge(
  bird_id_lookup,
  bird_morphometrics,
  by = "BIRD_KEY",
  all.x = TRUE
)

bird_auc <- merge(
  bird_auc,
  contraction_auc,
  by = "BIRD_KEY",
  all.x = TRUE
)

bird_auc <- merge(
  bird_auc,
  relaxation_auc,
  by = "BIRD_KEY",
  all.x = TRUE
)

bird_auc$plastic_group <- factor(
  ifelse(
    is.na(bird_auc$proventriculus_plastic_mass),
    NA,
    ifelse(
      bird_auc$proventriculus_plastic_mass > 0.5,
      "High plastic (>0.5)",
      "Low plastic (≤0.5)"
    )
  ),
  levels = names(pal_group)
)

# Put the columns into an intentionally readable order.
bird_auc <- bird_auc[, c(
  "BIRD_KEY",
  "BIRD_ID",
  "morphometrics_ID",
  "bird_body_mass",
  "WC",
  "HB",
  "Cul",
  "proventriculus_tissue_mass",
  "proventriculus_plastic_mass",
  "plastic_group",
  "auc_contraction",
  "auc_relaxation"
)]

bird_auc <- bird_auc[
  order(bird_auc$BIRD_KEY),
]


# ------------------------------------------------------------
# AUDIT bird_auc BEFORE STATISTICS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL ONE-ROW-PER-BIRD AUC DATAFRAME\n")
cat("============================================================\n")

print(
  bird_auc,
  row.names = FALSE
)

cat("\nBirds by plastic group:\n")
print(
  table(
    bird_auc$plastic_group,
    useNA = "ifany"
  )
)

cat("\nMissing bird body mass:\n")
print(
  bird_auc[
    is.na(bird_auc$bird_body_mass),
    c(
      "BIRD_KEY",
      "BIRD_ID",
      "morphometrics_ID"
    )
  ]
)

cat("\nMissing proventriculus plastic mass:\n")
print(
  bird_auc[
    is.na(bird_auc$proventriculus_plastic_mass),
    c(
      "BIRD_KEY",
      "BIRD_ID",
      "morphometrics_ID"
    )
  ]
)

write.csv(
  bird_auc,
  "bird_level_auc_with_morphometrics.csv",
  row.names = FALSE
)


# ============================================================
# ORIGINAL AUC GROUP COMPARISONS
# ============================================================
#
# These retain the original Welch t-test comparison between low- and
# high-plastic birds.

auc_contraction_ttest_data <- subset(
  bird_auc,
  !is.na(plastic_group) &
    is.finite(auc_contraction)
)

tt_contraction_auc <- t.test(
  auc_contraction ~ plastic_group,
  data = auc_contraction_ttest_data
)

cat("\n============================================================\n")
cat("CONTRACTION AUC: UNADJUSTED PLASTIC-GROUP T-TEST\n")
cat("============================================================\n")
print(tt_contraction_auc)

stripchart(
  auc_contraction ~ plastic_group,
  data = auc_contraction_ttest_data,
  vertical = TRUE,
  method = "jitter",
  ylab = "Contraction AUC",
  xlab = "Plastic group"
)


auc_relaxation_ttest_data <- subset(
  bird_auc,
  !is.na(plastic_group) &
    is.finite(auc_relaxation)
)

tt_relaxation_auc <- t.test(
  auc_relaxation ~ plastic_group,
  data = auc_relaxation_ttest_data
)

cat("\n============================================================\n")
cat("RELAXATION AUC: UNADJUSTED PLASTIC-GROUP T-TEST\n")
cat("============================================================\n")
print(tt_relaxation_auc)


# ============================================================
# BODY-MASS-ADJUSTED AUC MODELS
# ============================================================
#
# Two transparent versions are fitted for each phase:
#
# 1. Group model:
#      AUC ~ plastic_group + bird body mass
#
# 2. Continuous plastic-mass model:
#      AUC ~ proventriculus plastic mass + bird body mass
#
# Body mass is centred and expressed per 100 source units only to make
# the coefficient readable. If body mass is recorded in grams, the
# coefficient is the AUC change per 100 g.


# ------------------------------------------------------------
# CONTRACTION AUC MODELS
# ------------------------------------------------------------

auc_contraction_model_data <- subset(
  bird_auc,
  is.finite(auc_contraction) &
    is.finite(bird_body_mass) &
    is.finite(proventriculus_plastic_mass) &
    !is.na(plastic_group)
)

auc_contraction_body_mass_mean <- mean(
  auc_contraction_model_data$bird_body_mass,
  na.rm = TRUE
)

auc_contraction_model_data$bird_body_mass_100_c <-
  (
    auc_contraction_model_data$bird_body_mass -
      auc_contraction_body_mass_mean
  ) / 100

cat("\n============================================================\n")
cat("CONTRACTION AUC: BIRDS ENTERING ADJUSTED MODELS\n")
cat("============================================================\n")

print(
  auc_contraction_model_data[, c(
    "BIRD_KEY",
    "BIRD_ID",
    "bird_body_mass",
    "proventriculus_tissue_mass",
    "proventriculus_plastic_mass",
    "plastic_group",
    "auc_contraction"
  )],
  row.names = FALSE
)

auc_contraction_group_body_model <- lm(
  auc_contraction ~
    plastic_group +
    bird_body_mass_100_c,
  data = auc_contraction_model_data
)

auc_contraction_plastic_mass_body_model <- lm(
  auc_contraction ~
    proventriculus_plastic_mass +
    bird_body_mass_100_c,
  data = auc_contraction_model_data
)

cat("\nContraction AUC model: plastic GROUP + bird body mass\n")
print(
  summary(
    auc_contraction_group_body_model
  )
)

cat("\nContraction AUC model: continuous plastic MASS + bird body mass\n")
print(
  summary(
    auc_contraction_plastic_mass_body_model
  )
)


# ------------------------------------------------------------
# RELAXATION AUC MODELS
# ------------------------------------------------------------

auc_relaxation_model_data <- subset(
  bird_auc,
  is.finite(auc_relaxation) &
    is.finite(bird_body_mass) &
    is.finite(proventriculus_plastic_mass) &
    !is.na(plastic_group)
)

auc_relaxation_body_mass_mean <- mean(
  auc_relaxation_model_data$bird_body_mass,
  na.rm = TRUE
)

auc_relaxation_model_data$bird_body_mass_100_c <-
  (
    auc_relaxation_model_data$bird_body_mass -
      auc_relaxation_body_mass_mean
  ) / 100

cat("\n============================================================\n")
cat("RELAXATION AUC: BIRDS ENTERING ADJUSTED MODELS\n")
cat("============================================================\n")

print(
  auc_relaxation_model_data[, c(
    "BIRD_KEY",
    "BIRD_ID",
    "bird_body_mass",
    "proventriculus_tissue_mass",
    "proventriculus_plastic_mass",
    "plastic_group",
    "auc_relaxation"
  )],
  row.names = FALSE
)

auc_relaxation_group_body_model <- lm(
  auc_relaxation ~
    plastic_group +
    bird_body_mass_100_c,
  data = auc_relaxation_model_data
)

auc_relaxation_plastic_mass_body_model <- lm(
  auc_relaxation ~
    proventriculus_plastic_mass +
    bird_body_mass_100_c,
  data = auc_relaxation_model_data
)

cat("\nRelaxation AUC model: plastic GROUP + bird body mass\n")
print(
  summary(
    auc_relaxation_group_body_model
  )
)

cat("\nRelaxation AUC model: continuous plastic MASS + bird body mass\n")
print(
  summary(
    auc_relaxation_plastic_mass_body_model
  )
)



# ============================================================
# CONTRACTION GAMM ANALYSIS
# ============================================================

# Start from the contraction dataset already baseline-corrected above.
contraction_gamm <- subset(
  sm_contract,
  !is.na(plastic_group) &
    !is.na(value_delta) &
    !is.na(BIRD_KEY) &
    is.finite(bird_body_mass)
)

contraction_gamm$time_min <- contraction_gamm$time / 2
contraction_gamm$BIRD_KEY <- factor(contraction_gamm$BIRD_KEY)
contraction_gamm$plastic_group <- droplevels(contraction_gamm$plastic_group)

# One-row-per-bird audit table showing EXACTLY which birds enter the model.
contraction_bird_check <- unique(
  contraction_gamm[, c(
    "BIRD_KEY",
    "BIRD_ID",
    "morphometrics_ID",
    "bird_body_mass",
    "WC",
    "HB",
    "Cul",
    "proventriculus_tissue_mass",
    "proventriculus_plastic_mass",
    "plastic_group"
  )]
)
contraction_bird_check <- contraction_bird_check[order(contraction_bird_check$BIRD_KEY), ]

cat("\n=========================\n")
cat("CONTRACTION: BIRDS ENTERING MODEL\n")
cat("=========================\n")
print(contraction_bird_check)
print(table(contraction_bird_check$plastic_group, useNA = "ifany"))

write.csv(
  contraction_bird_check,
  "contraction_gamm_bird_alignment_check.csv",
  row.names = FALSE
)

# Centre body mass using exactly one weight per bird, not one value per time point.
contraction_body_mass_by_bird <- unique(
  contraction_gamm[, c("BIRD_KEY", "bird_body_mass")]
)
contraction_body_mass_mean <- mean(
  contraction_body_mass_by_bird$bird_body_mass,
  na.rm = TRUE
)

# Express the coefficient per 100 units of bird_body_mass.
# If bird_body_mass is stored in grams, this is the effect per 100 g.
# This changes only the coefficient scale, not the fitted model or p-value.
contraction_gamm$bird_body_mass_100_c <-
  (contraction_gamm$bird_body_mass - contraction_body_mass_mean) / 100

cat("Mean contraction bird body mass:", contraction_body_mass_mean, "\n")
cat(
  "Contraction body-mass range:",
  paste(range(contraction_body_mass_by_bird$bird_body_mass, na.rm = TRUE), collapse = " to "),
  "\n"
)

# ------------------------------------------------------------
# Contraction model 0: original plastic model WITHOUT body mass
# This is useful for checking whether adding body mass changes anything.
# ------------------------------------------------------------
contraction_no_body <- mgcv::gamm(
  value_delta ~ plastic_group +
    s(time_min, by = plastic_group, k = 10),
  random = list(BIRD_KEY = ~1),
  data = contraction_gamm,
  method = "ML"
)

# ------------------------------------------------------------
# Contraction model 1: reduced model WITH body mass
# Common time curve, adjusted for body mass.
# ------------------------------------------------------------
contraction_reduced_body <- mgcv::gamm(
  value_delta ~ bird_body_mass_100_c +
    s(time_min, k = 10),
  random = list(BIRD_KEY = ~1),
  data = contraction_gamm,
  method = "ML"
)

# ------------------------------------------------------------
# Contraction model 2: full model WITH body mass
# Bird body mass + plastic group + different smooth time curve by plastic group.
# ------------------------------------------------------------
contraction_full_body <- mgcv::gamm(
  value_delta ~ bird_body_mass_100_c +
    plastic_group +
    s(time_min, by = plastic_group, k = 10),
  random = list(BIRD_KEY = ~1),
  data = contraction_gamm,
  method = "ML"
)

cat("\n=========================\n")
cat("CONTRACTION MODEL SUMMARIES\n")
cat("=========================\n")

cat("\nOriginal plastic model, no body mass\n")
print(summary(contraction_no_body$gam))

cat("\nReduced model, adjusted for body mass\n")
print(summary(contraction_reduced_body$gam))

cat("\nFull plastic model, adjusted for body mass\n")
print(summary(contraction_full_body$gam))

# Print the body-mass coefficient with more precision than sjPlot's default.
cat("\nContraction body-mass coefficient, effect per 100 body-mass units\n")
print(
  summary(contraction_full_body$gam)$p.table[
    "bird_body_mass_100_c",
    ,
    drop = FALSE
  ],
  digits = 8
)

# Does adding body mass improve the full plastic model?
cat("\nContraction: full plastic model WITHOUT vs WITH body mass\n")
print(AIC(contraction_no_body$lme, contraction_full_body$lme))
print(anova(contraction_no_body$lme, contraction_full_body$lme))

# Does plastic group improve fit AFTER adjustment for body mass?
cat("\nContraction: reduced vs full model, both adjusted for body mass\n")
print(AIC(contraction_reduced_body$lme, contraction_full_body$lme))
print(anova(contraction_reduced_body$lme, contraction_full_body$lme))

# ============================================================
# RELAXATION GAMM ANALYSIS
# ============================================================

relaxation_gamm <- subset(
  relax_sm,
  !is.na(plastic_group) &
    !is.na(value_delta) &
    !is.na(BIRD_KEY) &
    is.finite(bird_body_mass)
)

relaxation_gamm$time_min <- relaxation_gamm$time / 2
relaxation_gamm$BIRD_KEY <- factor(relaxation_gamm$BIRD_KEY)
relaxation_gamm$plastic_group <- droplevels(relaxation_gamm$plastic_group)

# One-row-per-bird audit table showing EXACTLY which birds enter the model.
relaxation_bird_check <- unique(
  relaxation_gamm[, c(
    "BIRD_KEY",
    "BIRD_ID",
    "morphometrics_ID",
    "bird_body_mass",
    "WC",
    "HB",
    "Cul",
    "proventriculus_tissue_mass",
    "proventriculus_plastic_mass",
    "plastic_group"
  )]
)
relaxation_bird_check <- relaxation_bird_check[order(relaxation_bird_check$BIRD_KEY), ]

cat("\n=========================\n")
cat("RELAXATION: BIRDS ENTERING MODEL\n")
cat("=========================\n")
print(relaxation_bird_check)
print(table(relaxation_bird_check$plastic_group, useNA = "ifany"))

write.csv(
  relaxation_bird_check,
  "relaxation_gamm_bird_alignment_check.csv",
  row.names = FALSE
)

# Centre body mass using exactly one weight per bird.
relaxation_body_mass_by_bird <- unique(
  relaxation_gamm[, c("BIRD_KEY", "bird_body_mass")]
)
relaxation_body_mass_mean <- mean(
  relaxation_body_mass_by_bird$bird_body_mass,
  na.rm = TRUE
)

# Express the coefficient per 100 body-mass units for readability.
relaxation_gamm$bird_body_mass_100_c <-
  (relaxation_gamm$bird_body_mass - relaxation_body_mass_mean) / 100

cat("Mean relaxation bird body mass:", relaxation_body_mass_mean, "\n")
cat(
  "Relaxation body-mass range:",
  paste(range(relaxation_body_mass_by_bird$bird_body_mass, na.rm = TRUE), collapse = " to "),
  "\n"
)

# ------------------------------------------------------------
# Relaxation model 0: original plastic model WITHOUT body mass
# ------------------------------------------------------------
relaxation_no_body <- mgcv::gamm(
  value_delta ~ plastic_group +
    s(time_min, by = plastic_group, k = 10),
  random = list(BIRD_KEY = ~1),
  data = relaxation_gamm,
  method = "ML"
)

# ------------------------------------------------------------
# Relaxation model 1: reduced model WITH body mass
# ------------------------------------------------------------
relaxation_reduced_body <- mgcv::gamm(
  value_delta ~ bird_body_mass_100_c +
    s(time_min, k = 10),
  random = list(BIRD_KEY = ~1),
  data = relaxation_gamm,
  method = "ML"
)

# ------------------------------------------------------------
# Relaxation model 2: full model WITH body mass
# ------------------------------------------------------------
relaxation_full_body <- mgcv::gamm(
  value_delta ~ bird_body_mass_100_c +
    plastic_group +
    s(time_min, by = plastic_group, k = 10),
  random = list(BIRD_KEY = ~1),
  data = relaxation_gamm,
  method = "ML"
)

cat("\n=========================\n")
cat("RELAXATION MODEL SUMMARIES\n")
cat("=========================\n")

cat("\nOriginal plastic model, no body mass\n")
print(summary(relaxation_no_body$gam))

cat("\nReduced model, adjusted for body mass\n")
print(summary(relaxation_reduced_body$gam))

cat("\nFull plastic model, adjusted for body mass\n")
print(summary(relaxation_full_body$gam))

cat("\nRelaxation body-mass coefficient, effect per 100 body-mass units\n")
print(
  summary(relaxation_full_body$gam)$p.table[
    "bird_body_mass_100_c",
    ,
    drop = FALSE
  ],
  digits = 8
)

# Does adding body mass improve the full plastic model?
cat("\nRelaxation: full plastic model WITHOUT vs WITH body mass\n")
print(AIC(relaxation_no_body$lme, relaxation_full_body$lme))
print(anova(relaxation_no_body$lme, relaxation_full_body$lme))

# Does plastic group improve fit AFTER adjustment for body mass?
cat("\nRelaxation: reduced vs full model, both adjusted for body mass\n")
print(AIC(relaxation_reduced_body$lme, relaxation_full_body$lme))
print(anova(relaxation_reduced_body$lme, relaxation_full_body$lme))

# ============================================================
# BODY-MASS-ADJUSTED PREDICTED CURVES
# ============================================================
# Predictions are made at centred body mass = 0, i.e. mean body mass.

prediction_time <- seq(0, t_max / 2, length.out = 300)

# Contraction predictions
contraction_prediction_data <- expand.grid(
  time_min = prediction_time,
  plastic_group = levels(contraction_gamm$plastic_group),
  KEEP.OUT.ATTRS = FALSE
)
contraction_prediction_data$plastic_group <- factor(
  contraction_prediction_data$plastic_group,
  levels = levels(contraction_gamm$plastic_group)
)
contraction_prediction_data$bird_body_mass_100_c <- 0

contraction_prediction <- predict(
  contraction_full_body$gam,
  newdata = contraction_prediction_data,
  type = "response",
  se.fit = TRUE
)

contraction_prediction_data$predicted <- as.numeric(contraction_prediction$fit)
contraction_prediction_data$se <- as.numeric(contraction_prediction$se.fit)
contraction_prediction_data$lower_95 <-
  contraction_prediction_data$predicted - 1.96 * contraction_prediction_data$se
contraction_prediction_data$upper_95 <-
  contraction_prediction_data$predicted + 1.96 * contraction_prediction_data$se

p_contraction_adjusted <- ggplot(
  contraction_prediction_data,
  aes(
    x = time_min,
    y = predicted,
    colour = plastic_group,
    fill = plastic_group
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_ribbon(
    aes(ymin = lower_95, ymax = upper_95),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = pal_group, name = "Plastic load") +
  scale_fill_manual(values = pal_group, guide = "none") +
  labs(
    x = "Time (minutes)",
    y = "Body-mass-adjusted contraction force (Δ g)",
    subtitle = paste0(
      "Predicted at mean bird body mass: ",
      round(contraction_body_mass_mean, 1)
    )
  ) +
  coord_cartesian(xlim = c(0, t_max / 2)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

print(p_contraction_adjusted)

# Relaxation predictions
relaxation_prediction_data <- expand.grid(
  time_min = prediction_time,
  plastic_group = levels(relaxation_gamm$plastic_group),
  KEEP.OUT.ATTRS = FALSE
)
relaxation_prediction_data$plastic_group <- factor(
  relaxation_prediction_data$plastic_group,
  levels = levels(relaxation_gamm$plastic_group)
)
relaxation_prediction_data$bird_body_mass_100_c <- 0

relaxation_prediction <- predict(
  relaxation_full_body$gam,
  newdata = relaxation_prediction_data,
  type = "response",
  se.fit = TRUE
)

relaxation_prediction_data$predicted <- as.numeric(relaxation_prediction$fit)
relaxation_prediction_data$se <- as.numeric(relaxation_prediction$se.fit)
relaxation_prediction_data$lower_95 <-
  relaxation_prediction_data$predicted - 1.96 * relaxation_prediction_data$se
relaxation_prediction_data$upper_95 <-
  relaxation_prediction_data$predicted + 1.96 * relaxation_prediction_data$se

p_relaxation_adjusted <- ggplot(
  relaxation_prediction_data,
  aes(
    x = time_min,
    y = predicted,
    colour = plastic_group,
    fill = plastic_group
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_ribbon(
    aes(ymin = lower_95, ymax = upper_95),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = pal_group, name = "Plastic load") +
  scale_fill_manual(values = pal_group, guide = "none") +
  labs(
    x = "Time (minutes)",
    y = "Body-mass-adjusted relaxation (Δ g)",
    subtitle = paste0(
      "Predicted at mean bird body mass: ",
      round(relaxation_body_mass_mean, 1)
    )
  ) +
  coord_cartesian(xlim = c(0, t_max / 2)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

print(p_relaxation_adjusted)

# ============================================================
# HTML MODEL TABLES
# ============================================================

model_table_dir <- file.path(DATA_DIR, "model_tables")
dir.create(model_table_dir, showWarnings = FALSE, recursive = TRUE)

# Contraction HTML table
sjPlot::tab_model(
  contraction_no_body$gam,
  contraction_reduced_body$gam,
  contraction_full_body$gam,
  title = "Contraction GAMM models",
  dv.labels = c(
    "Plastic model without body mass",
    "Reduced model with body mass",
    "Full plastic model with body mass"
  ),
  auto.label = FALSE,
  show.intercept = TRUE,
  show.est = TRUE,
  show.se = TRUE,
  show.ci = 0.95,
  show.stat = TRUE,
  show.p = TRUE,
  show.r2 = TRUE,
  show.aic = TRUE,
  show.dev = TRUE,
  show.loglik = TRUE,
  show.icc = FALSE,
  show.re.var = FALSE,
  show.ngroups = FALSE,
  digits = 5,
  digits.p = 4,
  emph.p = TRUE,
  file = file.path(model_table_dir, "GAMM_contraction_models.html")
)

# Relaxation HTML table
sjPlot::tab_model(
  relaxation_no_body$gam,
  relaxation_reduced_body$gam,
  relaxation_full_body$gam,
  title = "Relaxation GAMM models",
  dv.labels = c(
    "Plastic model without body mass",
    "Reduced model with body mass",
    "Full plastic model with body mass"
  ),
  auto.label = FALSE,
  show.intercept = TRUE,
  show.est = TRUE,
  show.se = TRUE,
  show.ci = 0.95,
  show.stat = TRUE,
  show.p = TRUE,
  show.r2 = TRUE,
  show.aic = TRUE,
  show.dev = TRUE,
  show.loglik = TRUE,
  show.icc = FALSE,
  show.re.var = FALSE,
  show.ngroups = FALSE,
  digits = 5,
  digits.p = 4,
  emph.p = TRUE,
  file = file.path(model_table_dir, "GAMM_relaxation_models.html")
)


# AUC model tables
sjPlot::tab_model(
  auc_contraction_group_body_model,
  auc_contraction_plastic_mass_body_model,
  title = "Contraction AUC models adjusted for bird body mass",
  dv.labels = c(
    "Plastic group + bird body mass",
    "Continuous proventriculus plastic mass + bird body mass"
  ),
  auto.label = FALSE,
  show.est = TRUE,
  show.se = TRUE,
  show.ci = 0.95,
  show.stat = TRUE,
  show.p = TRUE,
  show.r2 = TRUE,
  show.aic = TRUE,
  digits = 5,
  digits.p = 4,
  emph.p = TRUE,
  file = file.path(
    model_table_dir,
    "AUC_contraction_models.html"
  )
)

sjPlot::tab_model(
  auc_relaxation_group_body_model,
  auc_relaxation_plastic_mass_body_model,
  title = "Relaxation AUC models adjusted for bird body mass",
  dv.labels = c(
    "Plastic group + bird body mass",
    "Continuous proventriculus plastic mass + bird body mass"
  ),
  auto.label = FALSE,
  show.est = TRUE,
  show.se = TRUE,
  show.ci = 0.95,
  show.stat = TRUE,
  show.p = TRUE,
  show.r2 = TRUE,
  show.aic = TRUE,
  digits = 5,
  digits.p = 4,
  emph.p = TRUE,
  file = file.path(
    model_table_dir,
    "AUC_relaxation_models.html"
  )
)

cat(
  "\nHTML model tables written to:\n",
  normalizePath(model_table_dir),
  "\n"
)

graph_dat<-relaxation_bird_check[-4,]

ggplot(
  graph_dat,
  aes(x = proventriculus_plastic_mass)
) +
  geom_histogram(
    bins = 15,
    colour = "black",
    fill = "grey70"
  ) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(
      sigma = 0.05,
      base = 10
    ),
    breaks = c(0, 0.5, 1, 3, 10, 30, 60),
    labels = c("0", "0.5", "1", "3", "10", "30", "60")
  ) +
  labs(
    x = "Proventriculus plastic mass (g)",
    y = "Number of birds",
    title = "Relaxation birds: proventriculus plastic mass"
  ) +
  theme_classic()

Distribution <- ggplot(
  graph_dat,
  aes(x = proventriculus_plastic_mass)
) +
  geom_histogram(
    bins = 15,
    colour = "black",
    fill = "grey70"
  ) +
  scale_x_continuous(
    trans = scales::sqrt_trans(),

    breaks = c(0, 0.5, 1, 3, 10, 30, 60),
    labels = c("0", "0.5", "1", "3", "10", "30", "60")
  ) +
  labs(
    x = "Proventriculus plastic mass (g)",
    y = "Number of birds",
    title = "Relaxation birds: proventriculus plastic mass"
  ) +
  theme_classic()

