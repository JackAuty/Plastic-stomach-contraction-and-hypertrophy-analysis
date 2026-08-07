# Stomach contraction and hypertrophy analysis
# Author: Jack Auty
# Version: 3.2
# Converted from R Markdown to a standard R script.

options(scipen = 999)

# Load packages once

packages <- c("ggplot2", "cowplot", "viridisLite", "sjPlot", "grid", "usethis", "mgcv")
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# Citations

pkg <- packages
citation_list <- character()
cln <- function(txt){
  txt <- gsub("\\\\texttt\\{(.*?)\\}", "\\1", txt)
  txt <- gsub("[_*\"“”]", "", txt)
  txt <- gsub(", \\.$|\\s*\\.$", "", txt)
  trimws(txt)
}
for (i in pkg) {
  raw <- paste(capture.output(print(citation(i), style = "text", bibtex = FALSE)), collapse = "\n")
  citation_list <- c(citation_list, cln(raw))
}
sjPlot::tab_df(data.frame(References = citation_list), col.header = "References", sep = " ")

# Load data

DATA_DIR <- usethis::proj_get()
ts_path    <- file.path(DATA_DIR, "final_trim.csv")
map_path   <- file.path(DATA_DIR, "video_map.csv")
mass_path  <- file.path(DATA_DIR, "masses.csv")
body_path  <- file.path(DATA_DIR, "Bird_Body_Morphometrics.csv")

required_files <- c(ts_path, map_path, mass_path, body_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("The following required data files were not found: ",
       paste(basename(missing_files), collapse = ", "))
}

ts_df   <- read.csv(ts_path,   check.names = TRUE)
map_df  <- read.csv(map_path,  check.names = TRUE, strip.white = TRUE)
mass_df <- read.csv(mass_path, check.names = TRUE, strip.white = TRUE)
body_df <- read.csv(body_path, check.names = TRUE, strip.white = TRUE)

required_body_columns <- c("ID", "Bird_Body_Weight")
missing_body_columns <- setdiff(required_body_columns, names(body_df))
if (length(missing_body_columns)) {
  stop("Bird_Body_Morphometrics.csv is missing required column(s): ",
       paste(missing_body_columns, collapse = ", "))
}

# Prepare keys and collapse duplicates

# Common BIRD_KEY
names(map_df) <- sub("BIRD[._-]?ID", "BIRD_ID", names(map_df), ignore.case = TRUE)

norm_id_mass <- function(x){
  x <- trimws(as.character(x))
  out <- x
  # NL-2025-10 -> NL-10
  m <- grepl("^[A-Z]+-20[0-9]{2}-[0-9]+$", x)
  out[m] <- sub("^([A-Z]+)-20[0-9]{2}-([0-9]+)$", "\\1-\\2", x[m])
  # 380-02015 -> 02015
  m <- grepl("^[0-9]{3}-[0-9]+$", x)
  out[m] <- sub("^[0-9]{3}-([0-9]+)$", "\\1", x[m])
  out
}

# Normalise IDs across all files
normalise_bird_key <- function(x){
  out <- norm_id_mass(x)
  sub("^BW-([0-9])$", "BW-0\\1", out)
}

map_df$BIRD_KEY  <- normalise_bird_key(map_df$BIRD_ID)
mass_df$BIRD_KEY <- normalise_bird_key(mass_df$ID)
body_df$BIRD_KEY <- normalise_bird_key(body_df$ID)

# Collapse duplicate rows in mass_df by ID + BIRD_KEY (sum numeric cols)
is_num <- sapply(mass_df, is.numeric)
sum_keep_na <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)

mass_df <- aggregate(
  mass_df[, is_num, drop = FALSE],
  by = list(ID = mass_df$ID, BIRD_KEY = mass_df$BIRD_KEY),
  FUN = sum_keep_na
)

# Keep one body-weight value per bird. If duplicate rows exist, use their mean.
mean_keep_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
body_df$Bird_Body_Weight <- suppressWarnings(as.numeric(body_df$Bird_Body_Weight))
body_mass <- aggregate(
  Bird_Body_Weight ~ BIRD_KEY,
  data = body_df,
  FUN = mean_keep_na,
  na.action = na.pass
)

# Helpers

# Long format for columns ending with the given suffix
make_long <- function(df, suffix){
  keep <- names(df)[grepl(paste0("(", suffix, ")$"), names(df), ignore.case = TRUE)]
  L <- stack(df[keep])
  L$time <- rep.int(seq_len(nrow(df)), times = length(keep))
  names(L) <- c("values", "ind", "time")
  L$VIDEO <- sub(paste0("_(?i:", suffix, ")$"), "", L$ind, perl = TRUE)
  L
}

# Baseline: subtract value at time == 0 (or first non-NA if no 0 present)
baseline_correct <- function(df){
  df <- df[order(df$VIDEO, df$time), ]
  df$baseline0 <- NA_real_
  for (v in unique(df$VIDEO)){
    idx <- which(df$VIDEO == v)
    d   <- df[idx, ]
    i0  <- which(d$time == 0 & !is.na(d$values))
    if (!length(i0)) i0 <- which(!is.na(d$values))[1]
    if (length(i0)) df$baseline0[idx] <- d$values[i0[1]]
  }
  df$value_delta <- df$values - df$baseline0
  df
}

# Simple summary
mean_sem <- function(x){
  m  <- mean(x, na.rm = TRUE)
  sd <- stats::sd(x, na.rm = TRUE)
  n  <- sum(!is.na(x))
  c(mean = m, sd = sd, n = n)
}

# Palette and settings
pal_group <- c("Low plastic (≤0.5)" = "#0072B2", "High plastic (>0.5)" = "#D55E00")
na_col    <- "#999999"
theme_small <- theme_minimal(base_size = 7) +
  theme(panel.grid.minor = element_blank())
t_max <- 120

# Two-level grouping by Stomach_mass
grp2 <- function(x){
  out <- ifelse(is.na(x), NA,
         ifelse(x > 0.5, "High plastic (>0.5)", "Low plastic (≤0.5)"))
  factor(out, levels = names(pal_group))
}

# Build long data and merge metadata

raw_long <- make_long(ts_df, "raw")
sm_long  <- make_long(ts_df, "smoothed")

keep_cols <- c("VIDEO","BIRD_ID","BIRD_KEY","PHASE")
raw_m <- merge(raw_long, map_df[keep_cols], by = "VIDEO", all.x = TRUE)
sm_m  <- merge(sm_long,  map_df[keep_cols], by = "VIDEO", all.x = TRUE)

# Add Stomach_mass etc. so we can form plastic groups
mass_keep  <- intersect(c("BIRD_KEY","Total_mass","Stomach_mass","Wt","Stomach_mass","Total_number"), names(mass_df))
merged_raw <- merge(raw_m, mass_df[mass_keep], by = "BIRD_KEY", all.x = TRUE)
merged_sm  <- merge(sm_m, mass_df[mass_keep], by = "BIRD_KEY", all.x = TRUE)

# Add bird body mass to both time-series datasets
merged_raw <- merge(merged_raw, body_mass, by = "BIRD_KEY", all.x = TRUE)
merged_sm  <- merge(merged_sm,  body_mass, by = "BIRD_KEY", all.x = TRUE)

# Quick sanity checks
unique(merged_sm$BIRD_KEY)
cat("Birds with body-mass data:",
    length(unique(merged_sm$BIRD_KEY[!is.na(merged_sm$Bird_Body_Weight)])),
    "of", length(unique(merged_sm$BIRD_KEY)), "\n")

missing_body_ids <- setdiff(unique(merged_sm$BIRD_KEY), body_mass$BIRD_KEY)
if (length(missing_body_ids)) {
  warning("No Bird_Body_Weight was matched for: ",
          paste(missing_body_ids, collapse = ", "))
}

# Figure 1: AB – Wt vs plastic mass and Contract SM mean ± SEM (legend below)

# Scatter (A)
df_scatter <- subset(mass_df, is.finite(Stomach_mass) & is.finite(Wt))
pear <- suppressWarnings(cor(df_scatter$Stomach_mass, df_scatter$Wt, use = "complete.obs", method = "pearson"))
pA <- ggplot(df_scatter, aes(x = Stomach_mass, y = Wt)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, colour = "grey55", fill = "grey85", na.rm = TRUE) +
  geom_point(shape = 21, size = 2.4, stroke = 0.3, colour = "grey25",
             fill = scales::alpha("#0072B2", 0.35), na.rm = TRUE) +
  labs(x = "Plastic mass (g)", y = "Stomach weight (g)",
       subtitle = sprintf("Linear fit • Pearson r = %.2f", pear)) +
  theme_classic(base_size = 11)

# Contract (smoothed) mean ± SEM (B)
sm_contract <- subset(merged_sm, grepl("^Contract", PHASE, ignore.case = TRUE))
sm_contract$VIDEO <- factor(sm_contract$VIDEO, levels = unique(sm_contract$VIDEO))
sm_contract <- baseline_correct(sm_contract)
sm_contract <- subset(sm_contract, time <= t_max)
sm_contract$plastic_group <- grp2(sm_contract$Stomach_mass)

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
                 plastic_group = factor(ifelse(Stomach_mass > 0.5, "High plastic (>0.5)", "Low plastic (≤0.5)"),
                                        levels = names(pal_group)))
pA2 <- ggplot(dfA, aes(Stomach_mass, Wt, colour = plastic_group, fill = plastic_group)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40", fill = scales::alpha("grey40", 0.125), linewidth = 0.5) +
  geom_point(shape = 21, size = 2.4, stroke = 0.35, alpha = 0.4) +
  scale_colour_manual(values = pal_group, name = "Plastic load") +
  scale_fill_manual(values = pal_group, guide = "none") +
  guides(colour = "none") +
  labs(x = "Plastic mass (g)", y = "Stomach weight (g)") +
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
relax_sm <- subset(merged_sm, grepl("^(Relax|Relaxation)", PHASE, ignore.case = TRUE))
relax_sm$VIDEO <- factor(relax_sm$VIDEO, levels = unique(relax_sm$VIDEO))
relax_sm <- baseline_correct(relax_sm)
relax_sm <- subset(relax_sm, time <= t_max)
relax_sm$plastic_group <- grp2(relax_sm$Stomach_mass)

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


# Export hypertrophy CSV

mass_df_new <- mass_df[, c("BIRD_KEY", "Stomach_mass", "Stomach_number", "Total_mass", "Total_number", "Wt")]
colnames(mass_df_new) <- c("ID", "Prov_plastic_mass", "Prov_plasic_number", "Total_plastic_mass", "Total_plastic_number", "Stomach_weight")
mass_df_new <- mass_df_new[!is.na(mass_df_new$Stomach_weight),]
write.csv(mass_df_new, "stomach_hypertrophy_prov_mass.csv", row.names = FALSE)

# Export contraction


contraction_data_20minutes <- sm_contract[sm_contract$time==120,]
write.csv(contraction_data_20minutes, "stomach_contraction_20mins.csv", row.names = FALSE)


# Statistical inference

# ==========================================
# Bird-level AUC for contraction + relaxation
# (23 birds, ~half/half; 1 contract + 1 relax video per bird)
# ==========================================

trapz_auc <- function(t, y){
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]; y <- y[ok]
  if (length(t) < 2) return(NA_real_)
  o <- order(t)
  t <- t[o]; y <- y[o]
  sum((y[-length(y)] + y[-1]) * 0.5 * diff(t), na.rm = TRUE)
}

calc_phase_auc <- function(df, phase_pattern, t_max = 120){
  d <- subset(df, grepl(phase_pattern, PHASE, ignore.case = TRUE))
  if (!nrow(d)) return(data.frame())

  # baseline per VIDEO (your function loops VIDEO internally)
  d <- baseline_correct(d)

  # restrict time window
  d <- subset(d, time <= t_max)

  # convert to minutes (your plots use time/2)
  d$time_min <- d$time / 2

  # ensure bird-level time series (safe even if only 1 video per bird/phase)
  bird_time <- aggregate(value_delta ~ BIRD_KEY + time_min, data = d, FUN = mean, na.rm = TRUE)

  # attach metadata (one row per bird)
  meta <- unique(d[, c("BIRD_KEY","BIRD_ID","Stomach_mass","Wt")])
  meta$plastic_group <- grp2(meta$Stomach_mass)

  bird_time <- merge(bird_time, meta, by = "BIRD_KEY", all.x = TRUE)

  # AUC per bird
  sp <- split(bird_time, bird_time$BIRD_KEY)
  out <- lapply(sp, function(x){
    data.frame(
      BIRD_KEY = x$BIRD_KEY[1],
      BIRD_ID  = x$BIRD_ID[1],
      Stomach_mass = x$Stomach_mass[1],
      Wt = x$Wt[1],
      plastic_group = x$plastic_group[1],
      AUC = trapz_auc(x$time_min, x$value_delta),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

# ---- Calculate AUCs ----
auc_contract <- calc_phase_auc(merged_sm, "^Contract", t_max = t_max)
names(auc_contract)[names(auc_contract) == "AUC"] <- "AUC_contraction"

auc_relax <- calc_phase_auc(merged_sm, "^(Relax|Relaxation)", t_max = t_max)
names(auc_relax)[names(auc_relax) == "AUC"] <- "AUC_relaxation"

# ---- Merge to one bird-level dataframe ----
bird_auc <- merge(
  auc_contract,
  auc_relax[, c("BIRD_KEY","AUC_relaxation")],
  by = "BIRD_KEY",
  all = TRUE
)

bird_auc

# ---- Sanity checks ----
# birds per group
table(bird_auc$plastic_group, useNA = "ifany")

# ---- t-tests (Welch by default) ----
tt_con <- t.test(AUC_contraction ~ plastic_group,
                 data = subset(bird_auc, !is.na(plastic_group) & !is.na(AUC_contraction)))

stripchart(AUC_contraction ~ plastic_group,
                 data = subset(bird_auc, !is.na(plastic_group) & !is.na(AUC_contraction)))
tt_rel <- t.test(AUC_relaxation ~ plastic_group,
                 data = subset(bird_auc, !is.na(plastic_group) & !is.na(AUC_relaxation)))

tt_con
tt_rel

# Optional export
write.csv(bird_auc, "bird_level_auc_contraction_relaxation.csv", row.names = FALSE)


# =========================
# GAMM model comparison, adjusted for bird body mass
# =========================
# Both models include bird body mass. The comparison therefore tests whether
# plastic group and group-specific time curves improve the fit after adjustment
# for body mass.
fit_gamm_compare_body_mass <- function(dat, phase_label = "phase"){
  dat <- subset(
    dat,
    !is.na(plastic_group) &
      !is.na(value_delta) &
      !is.na(BIRD_KEY) &
      is.finite(Bird_Body_Weight)
  )

  if (!nrow(dat)) {
    stop("No complete observations remain for ", phase_label,
         " after requiring Bird_Body_Weight.")
  }

  dat$time_min <- dat$time / 2
  dat$BIRD_KEY <- factor(dat$BIRD_KEY)
  dat$plastic_group <- droplevels(dat$plastic_group)

  # Centre body mass so the intercept represents a bird of average body mass.
  body_mass_mean <- mean(dat$Bird_Body_Weight, na.rm = TRUE)
  dat$Bird_Body_Weight_c <- dat$Bird_Body_Weight - body_mass_mean

  if (nlevels(dat$plastic_group) < 2) {
    stop("Fewer than two plastic groups remain for ", phase_label,
         " after matching body-mass data.")
  }

  # Reduced model: common time curve plus a linear body-mass effect.
  m0 <- mgcv::gamm(
    value_delta ~ Bird_Body_Weight_c + s(time_min, k = 10),
    random = list(BIRD_KEY = ~1),
    data = dat,
    method = "ML"
  )

  # Full model: body mass plus plastic-group main effect and separate time curves.
  m1 <- mgcv::gamm(
    value_delta ~ Bird_Body_Weight_c + plastic_group +
      s(time_min, by = plastic_group, k = 10),
    random = list(BIRD_KEY = ~1),
    data = dat,
    method = "ML"
  )

  cat("\n=========================\n", phase_label,
      " adjusted for bird body mass\n=========================\n")
  cat("Birds included:", nlevels(dat$BIRD_KEY), "\n")
  cat("Mean bird body mass:", round(body_mass_mean, 3), "\n")
  cat("Body-mass range:",
      paste(round(range(dat$Bird_Body_Weight), 3), collapse = " to "), "\n")

  cat("\nAIC comparison using maximum likelihood\n")
  print(AIC(m0$lme, m1$lme))

  cat("\nLikelihood-ratio comparison using maximum likelihood\n")
  print(anova(m0$lme, m1$lme))

  cat("\nGAM summary: reduced model\n")
  print(summary(m0$gam))
  cat("\nGAM summary: full model\n")
  print(summary(m1$gam))

  list(
    m0 = m0,
    m1 = m1,
    dat = dat,
    body_mass_mean = body_mass_mean
  )
}

# =========================
# Contraction
# =========================
d_con <- sm_contract
res_con_body_mass <- fit_gamm_compare_body_mass(
  d_con,
  "Contraction"
)

# =========================
# Relaxation
# =========================
d_rel <- relax_sm
res_rel_body_mass <- fit_gamm_compare_body_mass(
  d_rel,
  "Relaxation"
)

# ============================================================
# Plot GAMM curves adjusted for bird body weight
# ============================================================

make_body_mass_adjusted_predictions <- function(
    model_result,
    n_time_points = 300
) {

  # Full GAMM containing body mass, plastic group and
  # group-specific smooths over time
  model <- model_result$m1$gam
  model_data <- model_result$dat

  # Prediction range
  time_sequence <- seq(
    min(model_data$time_min, na.rm = TRUE),
    max(model_data$time_min, na.rm = TRUE),
    length.out = n_time_points
  )

  # Preserve the factor levels used when fitting the model
  plastic_levels <- levels(model_data$plastic_group)

  prediction_data <- expand.grid(
    time_min = time_sequence,
    plastic_group = plastic_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  prediction_data$plastic_group <- factor(
    prediction_data$plastic_group,
    levels = plastic_levels
  )

  # Body mass was centred before fitting the GAMM.
  # A centred value of zero represents the mean bird body mass.
  prediction_data$Bird_Body_Weight_c <- 0

  predictions <- predict(
    model,
    newdata = prediction_data,
    type = "response",
    se.fit = TRUE
  )

  prediction_data$predicted <- as.numeric(predictions$fit)
  prediction_data$se <- as.numeric(predictions$se.fit)

  prediction_data$lower <- (
    prediction_data$predicted -
      prediction_data$se
  )

  prediction_data$upper <- (
    prediction_data$predicted +
      prediction_data$se
  )

  prediction_data$Bird_Body_Weight <- model_result$body_mass_mean

  prediction_data
}


# Generate adjusted predictions
pred_contraction <- make_body_mass_adjusted_predictions(
  res_con_body_mass
)

pred_relaxation <- make_body_mass_adjusted_predictions(
  res_rel_body_mass
)


# ============================================================
# Contraction graph
# ============================================================

p_contraction_adjusted <- ggplot(
  pred_contraction,
  aes(
    x = time_min,
    y = predicted,
    colour = plastic_group,
    fill = plastic_group
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    colour = "grey70"
  ) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_colour_manual(
    values = pal_group,
    name = "Plastic load"
  ) +
  scale_fill_manual(
    values = pal_group,
    guide = "none"
  ) +
  labs(
    x = "Time (minutes)",
    y = "Body-weight-adjusted contraction force (Δ g)",
    )
  ) +
  coord_cartesian(
    xlim = c(0, t_max / 2)
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "bottom"
  )


# ============================================================
# Relaxation graph
# ============================================================

p_relaxation_adjusted <- ggplot(
  pred_relaxation,
  aes(
    x = time_min,
    y = predicted,
    colour = plastic_group,
    fill = plastic_group
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    colour = "grey70"
  ) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_colour_manual(
    values = pal_group,
    name = "Plastic load"
  ) +
  scale_fill_manual(
    values = pal_group,
    guide = "none"
  ) +
  labs(
    x = "Time (minutes)",
    y = "Body-weight-adjusted relaxation (Δ g)",

    )
  ) +
  coord_cartesian(
    xlim = c(0, t_max / 2)
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "bottom"
  )


# Display separately
print(p_contraction_adjusted)
print(p_relaxation_adjusted)


# ============================================================
# Combine into one figure
# ============================================================

body_mass_adjusted_figure <- cowplot::plot_grid(
  p_contraction_adjusted +
    theme(legend.position = "none"),
  p_relaxation_adjusted +
    theme(legend.position = "none"),
  labels = c("A", "B"),
  label_fontface = "bold",
  label_size = 14,
  nrow = 1
)

adjusted_legend <- cowplot::get_legend(
  p_contraction_adjusted +
    theme(legend.position = "bottom")
)

body_mass_adjusted_figure <- cowplot::plot_grid(
  body_mass_adjusted_figure,
  adjusted_legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

print(body_mass_adjusted_figure)


# Save figure
ggsave(
  filename = "Figure_GAMM_body_weight_adjusted.pdf",
  plot = body_mass_adjusted_figure,
  width = 270,
  height = 120,
  units = "mm"
)


# Optional: export the plotted predictions
write.csv(
  pred_contraction,
  "GAMM_body_weight_adjusted_contraction_predictions.csv",
  row.names = FALSE
)

write.csv(
  pred_relaxation,
  "GAMM_body_weight_adjusted_relaxation_predictions.csv",
  row.names = FALSE
)