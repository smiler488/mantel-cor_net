#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Correlation heatmap with Mantel-test links
#
# Upper triangle: Pearson r between explanatory variables, with significance
# stars. Lines from each target (e.g. Yield, PUE) to every explanatory
# variable show a Mantel test between the two distance matrices:
#   colour = Mantel p, width = |Mantel r|, solid / dashed = sign of r.
#
# Usage:
#   Rscript cor_mantel_overlay.R [input.csv] [output_prefix]
# Defaults: demo.csv and cor_upper_heatmap_mantel (writes .pdf and .svg).
# Edit the settings block below for other column names.
#
# Liangchao Deng (邓良超), MIT License
# ---------------------------------------------------------------------------

# ---- settings -------------------------------------------------------------
vars <- c("Year", "cultivar", "Treatment",
          "0-10cm AP", "10-20cm AP", "20-30cm AP", "30-40cm AP",
          "P-aMT", "P-aMTR", "P-aMCR", "A-aMCR", "PA", "PT",
          "P-aPCR", "A-aPCR", "PHI", "PEP", "BR", "SP")
targets      <- c("Yield", "PUE", "SB")
permutations <- 999
seed         <- 123
width_in     <- 10
height_in    <- 8

# ---- packages -------------------------------------------------------------
for (pkg in c("ggplot2", "vegan", "ggnewscale", "svglite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(vegan)
  library(ggnewscale)
})

# ---- input ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[1] else "demo.csv"
out_prefix <- if (length(args) >= 2) args[2] else "cor_upper_heatmap_mantel"

df <- read.csv(input_path, check.names = FALSE)
missing_cols <- setdiff(c(vars, targets), names(df))
if (length(missing_cols)) {
  stop("Columns not found in ", input_path, ": ", paste(missing_cols, collapse = ", "))
}

df_sel <- df[, vars, drop = FALSE]
is_cat <- !vapply(df_sel, is.numeric, logical(1))
if (any(is_cat)) {
  message("Categorical variables (no Pearson r; Mantel uses a 0/1 mismatch distance): ",
          paste(names(df_sel)[is_cat], collapse = ", "))
}
n <- length(vars)

stars <- function(p) {
  if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
}

# ---- Pearson correlations, upper triangle --------------------------------
# Pearson r on arbitrary factor codes has no meaning, so cells that involve a
# categorical variable are left empty (grey).
grid <- subset(expand.grid(i = seq_len(n), j = seq_len(n)), j > i)
cor_pair <- function(i, j) {
  if (is_cat[i] || is_cat[j]) return(c(NA_real_, NA_real_))
  x <- df_sel[[i]]; y <- df_sel[[j]]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(c(NA_real_, NA_real_))
  ct <- cor.test(x[ok], y[ok], method = "pearson")
  c(unname(ct$estimate), ct$p.value)
}
rp <- t(mapply(cor_pair, grid$i, grid$j))

df_plot <- data.frame(
  xi = grid$j,
  yi = n - grid$i + 1,
  r  = rp[, 1],
  p  = rp[, 2]
)
df_plot$label   <- vapply(df_plot$p, stars, character(1))
df_plot$sig_cat <- factor(ifelse(df_plot$label == "", "ns", df_plot$label),
                          levels = c("***", "**", "*", "ns"))
diag_pos <- data.frame(var = vars, x = seq_len(n), y = n - seq_len(n) + 1)

p <- ggplot(df_plot, aes(x = xi, y = yi, fill = r)) +
  geom_tile(width = 1, height = 1, color = "grey80", linewidth = 0.3) +
  geom_text(aes(label = label, colour = sig_cat), size = 5, fontface = "bold",
            na.rm = TRUE, show.legend = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1), na.value = "grey95", name = "Pearson r",
                       guide = guide_colorbar(order = 1)) +
  scale_colour_manual(
    values = c("***" = "black", "**" = "black", "*" = "black", "ns" = "black"),
    breaks = c("***", "**", "*"),
    labels = c("p < 0.001", "p < 0.01", "p < 0.05"),
    guide  = guide_legend(title = "Significance",
                          override.aes = list(label = c("***", "**", "*")), order = 2)
  ) +
  geom_point(data = diag_pos, aes(x = x, y = y), shape = 16, size = 5,
             inherit.aes = FALSE, show.legend = FALSE) +
  scale_x_continuous(position = "top", breaks = seq_len(n), labels = vars,
                     expand = c(0, 0), limits = c(-0.2, n + 0.6)) +
  scale_y_continuous(position = "right", breaks = seq_len(n), labels = rev(vars),
                     expand = c(0, 0), limits = c(-0.4, n + 0.6)) +
  coord_fixed(clip = "off") +
  theme_minimal() +
  theme(
    axis.title        = element_blank(),
    panel.grid        = element_blank(),
    axis.text.x.top   = element_text(angle = 45, hjust = 0, vjust = 0, margin = margin(b = 8)),
    axis.text.y.right = element_text(hjust = 0, margin = margin(l = 6)),
    plot.margin       = margin(t = 20, r = 20, b = 40, l = 40)
  )

# ---- Mantel tests ----------------------------------------------------------
# Distance for one variable: |a - b| for numbers, 0/1 mismatch for categories.
var_dist <- function(x) {
  if (is.numeric(x)) {
    dist(x, method = "euclidean")
  } else {
    f <- as.factor(x)
    as.dist(outer(f, f, function(a, b) as.numeric(a != b)))
  }
}

set.seed(seed)
mantel_rows <- list()
for (tg in targets) {
  xt <- df[[tg]]
  for (v in vars) {
    yv <- df_sel[[v]]
    ok <- !is.na(xt) & !is.na(yv)
    if (sum(ok) < 3 || length(unique(xt[ok])) < 2 || length(unique(yv[ok])) < 2) next
    mt <- tryCatch(
      vegan::mantel(var_dist(xt[ok]), var_dist(yv[ok]), method = "pearson",
                    permutations = permutations),
      error = function(e) NULL)
    if (!is.null(mt)) {
      mantel_rows[[length(mantel_rows) + 1]] <-
        data.frame(target = tg, var = v, r = unname(mt$statistic), p = unname(mt$signif))
    }
  }
}
mantel_df <- if (length(mantel_rows)) do.call(rbind, mantel_rows) else NULL

p_all <- p
if (!is.null(mantel_df) && nrow(mantel_df) > 0) {
  write.csv(mantel_df, paste0(out_prefix, "_mantel.csv"), row.names = FALSE)

  mantel_df$p_cat <- cut(mantel_df$p, c(-Inf, 0.01, 0.05, Inf),
                         labels = c("p_le_01", "p01_05", "p_gt_05"))
  mantel_df$r_bin <- cut(abs(mantel_df$r), c(-Inf, 0.25, 0.5, Inf), right = FALSE,
                         labels = c("r_lt_0.25", "r_0.25_0.5", "r_ge_0.5"))
  mantel_df$dir   <- ifelse(mantel_df$r >= 0, "pos", "neg")

  # Target nodes sit in the empty lower-left triangle, spread along its diagonal
  k <- length(targets)
  s <- if (k == 1) 0.5 else seq(0, 1, length.out = k)
  t_pos <- data.frame(target = targets,
                      x = 0.2 + s * (0.6 * n - 0.2),
                      y = 0.6 * n - s * (0.6 * n - 0.3))

  edge_df <- merge(mantel_df, diag_pos, by = "var")
  edge_df <- merge(edge_df, t_pos, by = "target", suffixes = c("_v", "_t"))

  p_all <- p +
    ggnewscale::new_scale_color() +
    geom_segment(data = edge_df,
                 aes(x = x_t, y = y_t, xend = x_v, yend = y_v,
                     colour = p_cat, linewidth = r_bin, linetype = dir),
                 inherit.aes = FALSE, alpha = 0.9) +
    geom_point(data = t_pos, aes(x = x, y = y), shape = 21, size = 3, stroke = 1,
               fill = "white", inherit.aes = FALSE, show.legend = FALSE) +
    geom_text(data = t_pos, aes(x = x, y = y, label = target), inherit.aes = FALSE,
              nudge_y = -0.35, show.legend = FALSE) +
    scale_colour_manual(values = c(p_le_01 = "orange", p01_05 = "blue", p_gt_05 = "grey80"),
                        breaks = c("p_le_01", "p01_05", "p_gt_05"),
                        labels = c("p <= 0.01", "0.01 < p <= 0.05", "p > 0.05"),
                        name = "Mantel p", drop = FALSE,
                        guide = guide_legend(order = 3)) +
    scale_linewidth_manual(values = c(r_ge_0.5 = 1.6, r_0.25_0.5 = 1.0, r_lt_0.25 = 0.5),
                           breaks = c("r_ge_0.5", "r_0.25_0.5", "r_lt_0.25"),
                           labels = c("|r| >= 0.5", "0.25 <= |r| < 0.5", "|r| < 0.25"),
                           name = "Mantel r", drop = FALSE,
                           guide = guide_legend(order = 4)) +
    scale_linetype_manual(values = c(pos = "solid", neg = "dashed"),
                          breaks = c("pos", "neg"), labels = c("r > 0", "r < 0"),
                          name = "Direction", drop = FALSE,
                          guide = guide_legend(order = 5))
}

ggsave(paste0(out_prefix, ".pdf"), p_all, width = width_in, height = height_in)
ggsave(paste0(out_prefix, ".svg"), p_all, width = width_in, height = height_in,
       device = svglite::svglite)
message("Written ", out_prefix, ".pdf / .svg",
        if (!is.null(mantel_df)) paste0(" and ", out_prefix, "_mantel.csv") else "")
