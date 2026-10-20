# Dunn's post-hoc test following Kruskal–Wallis, with Benjamini–Hochberg FDR
# correction (manuscript Fig 5 / Table S7).
#
# Input: tab-separated file with columns "population" and "value"
#        (one row per individual).
# Usage: Rscript dunn_test.R load_per_group.tsv dunn_results.tsv

suppressMessages({
  library(dunn.test)
})

args <- commandArgs(trailingOnly = TRUE)
infile  <- ifelse(length(args) >= 1, args[1], "load_per_group.tsv")
outfile <- ifelse(length(args) >= 2, args[2], "dunn_results.tsv")

d <- read.delim(infile, header = TRUE)
stopifnot(all(c("population", "value") %in% names(d)))

kw <- kruskal.test(value ~ population, data = d)
cat(sprintf("Kruskal-Wallis chi2 = %.3f, df = %d, P = %.4g\n",
            kw$statistic, kw$parameter, kw$p.value))

dt <- dunn.test(d$value, d$population, method = "bh", kw = TRUE, list = TRUE)
res <- data.frame(
  comparison = dt$comparisons,
  Z          = dt$Z,
  P          = dt$P,
  P_adj_BH   = dt$P.adjusted
)
res$signif <- cut(res$P_adj_BH,
                  breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
                  labels = c("****", "***", "**", "*", "ns"))
write.table(res, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Written:", outfile, "\n")
