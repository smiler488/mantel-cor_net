# Run the figure script on demo.csv and check its outputs.
#   Rscript tests/check_outputs.R
out <- file.path(tempdir(), "fig")
status <- system2("Rscript", c("cor_mantel_overlay.R", "demo.csv", out))
stopifnot(status == 0)
stopifnot(file.exists(paste0(out, ".pdf")), file.exists(paste0(out, ".svg")))
m <- read.csv(paste0(out, "_mantel.csv"))
stopifnot(nrow(m) == 3 * 19)                       # 3 targets x 19 variables
stopifnot(all(m$p > 0 & m$p <= 1), all(abs(m$r) <= 1))
# categorical variables use the mismatch distance; Treatment should relate to Yield
stopifnot(m$p[m$target == "Yield" & m$var == "Treatment"] <= 0.01)
cat("mantel-cor_net: outputs OK\n")
