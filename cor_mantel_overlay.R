# Install missing dependencies for reproducibility
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan", repos = "https://cloud.r-project.org")
if (!requireNamespace("ggnewscale", quietly = TRUE)) install.packages("ggnewscale", repos = "https://cloud.r-project.org")
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite", repos = "https://cloud.r-project.org")
library(ggplot2)
library(vegan)
library(ggnewscale)

# Read data and convert non-numeric variables to numeric codes for correlation/Mantel
input_path <- if (file.exists("demo.csv")) "demo.csv" else "Copy of data.csv"
df <- read.csv(input_path, check.names = FALSE)
vars <- c("Year","cultivar","Treatment","0-10cm AP","10-20cm AP","20-30cm AP","30-40cm AP","P-aMT","P-aMTR","P-aMCR","A-aMCR","PA","PT","P-aPCR","A-aPCR","PHI","PEP","BR","SP")
df_sel <- df[, vars, drop = FALSE]
df_num <- as.data.frame(lapply(df_sel, function(x) if (is.numeric(x)) x else as.numeric(as.factor(x))))

n <- ncol(df_num)
varnames <- colnames(df_num)

# Map p-values to significance stars
stars <- function(p) {
  if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
}

# Compute upper-triangular Pearson correlations and p-values
grid <- expand.grid(i = seq_len(n), j = seq_len(n))
grid <- subset(grid, j > i)
rmat <- cor(df_num, use = "pairwise.complete.obs", method = "pearson")
pmat <- matrix(NA_real_, n, n)
for (i in seq_len(n)) for (j in seq_len(n)) {
  if (i != j) {
    xi <- df_num[[i]]; xj <- df_num[[j]]; ok <- is.finite(xi) & is.finite(xj)
    pmat[i, j] <- if (sum(ok) >= 3) tryCatch(cor.test(xi[ok], xj[ok])$p.value, error = function(e) NA_real_) else NA_real_
  }
}

# Build plotting dataframe in long format
df_plot <- data.frame(
  row_i = grid$i,
  col_j = grid$j,
  row = factor(varnames[grid$i], levels = varnames),
  col = factor(varnames[grid$j], levels = varnames),
  r = mapply(function(a, b) rmat[a, b], grid$i, grid$j),
  p = mapply(function(a, b) pmat[a, b], grid$i, grid$j)
)
df_plot$label <- vapply(df_plot$p, stars, character(1))
df_plot$sig_cat <- ifelse(is.na(df_plot$p), "ns",
                          ifelse(df_plot$p < 0.001, "***",
                                 ifelse(df_plot$p < 0.01, "**",
                                        ifelse(df_plot$p < 0.05, "*", "ns"))))
df_plot$sig_cat <- factor(df_plot$sig_cat, levels = c("***","**","*","ns"))

# Coordinates for heatmap grid and overlay elements
df_plot$xi <- df_plot$col_j
df_plot$yi <- n - df_plot$row_i + 1
diag_pos <- data.frame(var = varnames, x = seq_len(n), y = n - seq_len(n) + 1)

# Correlation heatmap (upper triangle) + significance stars
p <- ggplot(df_plot, aes(x = xi, y = yi, fill = r)) +
  geom_tile(width = 1, height = 1, color = "grey80", linewidth = 0.3) +
  geom_text(aes(label = label, colour = sig_cat), size = 5, fontface = "bold", na.rm = TRUE, show.legend = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1), guide = guide_colorbar(order = 1)) +
  scale_colour_manual(
    values = c("***" = "black", "**" = "black", "*" = "black", "ns" = "black"),
    breaks = c("***","**","*"),
    labels = c("p<0.001","p<0.01","p<0.05"),
    guide = guide_legend(title = "Significance", override.aes = list(label = c("***","**","*")), order = 2)
  ) +
    geom_point(data = diag_pos, aes(x = x, y = y), shape = 16, size = 5, inherit.aes = FALSE, show.legend = FALSE) +
  scale_x_continuous(position = "top", breaks = seq_len(n), labels = varnames, expand = c(0, 0), limits = c(-0.2, n + 0.6)) +
  scale_y_continuous(position = "right", breaks = seq_len(n), labels = rev(varnames), expand = c(0, 0), limits = c(-0.4, n + 0.6)) +
  coord_fixed(clip = "off") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, margin = margin(b = 8)),
    axis.text.y.right = element_text(hjust = 0, margin = margin(l = 6)),
    plot.margin = margin(t = 20, r = 20, b = 40, l = 40)
  )

# Mantel test: relate distance matrices of targets (Yield/PUE/SB) to each explanatory variable
targets <- c("Yield","PUE","SB")
df_tar <- df[, targets, drop = FALSE]
df_tar_num <- as.data.frame(lapply(df_tar, function(x) if (is.numeric(x)) x else as.numeric(as.factor(x))))

# Numeric vars use absolute difference; categorical use 0/1 mismatch; return dist object
var_dist <- function(x) {
  if (is.numeric(x)) as.dist(outer(x, x, function(a, b) abs(a - b))) else as.dist(outer(as.factor(x), as.factor(x), function(a, b) as.numeric(a != b)))
}

mantel_results <- list()
set.seed(123)
for (t in targets) {
  xr <- df_tar_num[[t]]
  for (v in varnames) {
    yv <- df_num[[v]]
    ok <- is.finite(xr) & is.finite(yv)
    if (sum(ok) >= 3 && length(unique(xr[ok])) > 1 && length(unique(yv[ok])) > 1) {
      dx <- var_dist(xr[ok]); dy <- var_dist(yv[ok])
      mt <- tryCatch(vegan::mantel(dx, dy, method = "pearson", permutations = 999), error = function(e) NULL)
      if (!is.null(mt)) mantel_results[[length(mantel_results) + 1]] <- data.frame(target = t, var = v, r = unname(mt$statistic), p = unname(mt$signif))
    }
  }
}

mantel_df <- if (length(mantel_results)) do.call(rbind, mantel_results) else NULL
if (!is.null(mantel_df) && nrow(mantel_df) > 0) {
  # Bin p-values, r magnitude, and sign
  mantel_df$p_cat <- ifelse(mantel_df$p > 0.05, "pgt05", ifelse(mantel_df$p > 0.01, "p01_05", "p001_01"))
  mantel_df$r_bin <- ifelse(mantel_df$r >= 0.5, "r_ge_0.5", ifelse(mantel_df$r >= 0.25, "r_ge_0.25", "r_lt_0.25"))
  mantel_df$dir <- ifelse(mantel_df$r >= 0, "pos", "neg")

  # Layout of target nodes and variable nodes; edges connect them
  var_pos <- data.frame(var = varnames, x = seq_len(n), y = n - seq_len(n) + 1)
  t_pos <- rbind(
    data.frame(target = "Yield", x = 0.2, y = (n + 0.5) / 2),
    data.frame(target = "PUE",   x = 0.2, y = 1.0),
    data.frame(target = "SB",    x = (n + 0.5) / 2, y = 0.3)
  )
  edge_df <- merge(mantel_df, var_pos, by.x = "var", by.y = "var")
  edge_df <- merge(edge_df, t_pos, by.x = "target", by.y = "target", suffixes = c("_v", "_t"))

  # Overlay Mantel results: color = p bins, linewidth = r bins, linetype = sign
  p_all <- p +
    ggnewscale::new_scale_color() +
    geom_segment(data = edge_df, aes(x = x_t, y = y_t, xend = x_v, yend = y_v, colour = p_cat, linewidth = r_bin, linetype = dir), inherit.aes = FALSE, alpha = 0.9) +
    geom_point(data = t_pos, aes(x = x, y = y), shape = 21, size = 3, stroke = 1, fill = "white", inherit.aes = FALSE, show.legend = FALSE) +
    geom_text(data = t_pos, aes(x = x, y = y, label = target), inherit.aes = FALSE, hjust = 0.5, angle = 0, nudge_y = -0.35, show.legend = FALSE) +
    scale_colour_manual(values = c(p001_01 = "orange", p01_05 = "blue", pgt05 = "grey80"), breaks = c("p001_01","p01_05","pgt05"), labels = c("0.001 < p <= 0.01","0.01 < p <= 0.05","p > 0.05"), name = "Mantel p", guide = guide_legend(order = 3, override.aes = list(shape = NA, label = ""))) +
    scale_linewidth_manual(values = c(r_ge_0.5 = 1.6, r_ge_0.25 = 1.0, r_lt_0.25 = 0.7), breaks = c("r_ge_0.5","r_ge_0.25","r_lt_0.25"), labels = c("r >= 0.5","0.25 <= r < 0.5","r < 0.25"), name = "Mantel r", guide = guide_legend(order = 4, override.aes = list(shape = NA, label = ""))) + 
    scale_linetype_manual(values = c(pos = "solid", neg = "dashed"), breaks = c("pos","neg"), labels = c("r > 0","r < 0"), name = "Direction", guide = guide_legend(order = 5, override.aes = list(shape = NA, label = "")))

  # Export PDF and SVG for downstream use
  ggsave("cor_upper_heatmap_mantel.pdf", p_all, width = 10, height = 8)
  ggsave("cor_upper_heatmap_mantel.svg", p_all, width = 10, height = 8, device = svglite::svglite)
}
