# Extract airway counts and run a simplified paired differential expression (base R only).
# Usage: Rscript prep.R <airway.RData> <sample_table.csv> <out.js>
args <- commandArgs(TRUE)
e <- new.env(); load(args[1], e)
x <- get(ls(e)[1], e)
counts <- attr(attr(attributes(x)$assays, "data"), "listData")[[1]]
ids <- attr(attributes(x)$rowRanges@partitioning, "NAMES")
st <- read.csv(args[2], row.names = 1)
runs <- rownames(st)
stopifnot(length(ids) == nrow(counts), ncol(counts) == length(runs))
colnames(counts) <- runs
keep <- rowSums(counts) >= 10
counts <- counts[keep, ]; ids <- ids[keep]
# median-of-ratios size factors
lg <- log(counts); ok <- is.finite(rowSums(lg))
sf <- exp(apply(lg[ok, ] - rowMeans(lg[ok, ]), 2, median))
norm <- sweep(counts, 2, sf, "/")
l <- log2(norm + 1)
trt <- st$dex == "trt"
# pair treated/untreated by cell line
o_t <- order(st$cell[trt]); o_u <- order(st$cell[!trt])
d <- l[, trt][, o_t] - l[, !trt][, o_u]
n <- ncol(d); lfc <- rowMeans(d); se <- apply(d, 1, sd) / sqrt(n)
tstat <- ifelse(se > 0, lfc / se, 0)
p <- 2 * pt(-abs(tstat), n - 1)
padj <- p.adjust(p, "BH")
sym <- c(ENSG00000152583="SPARCL1", ENSG00000179094="PER1", ENSG00000116584="ARHGEF2",
         ENSG00000189221="MAOA", ENSG00000120129="DUSP1", ENSG00000163884="KLF15",
         ENSG00000157514="TSC22D3", ENSG00000109906="ZBTB16", ENSG00000165995="CACNB2")
genes <- lapply(seq_along(ids), function(i) list(
  id = ids[i], sym = unname(sym[ids[i]]), bm = round(mean(norm[i, ]), 1),
  lfc = round(lfc[i], 3), p = signif(p[i], 3), q = signif(padj[i], 3),
  c = as.integer(round(norm[i, ]))))
j <- function(v) paste0("[", paste0('"', v, '"', collapse = ","), "]")
con <- file(args[3], "w")
writeLines(paste0("window.AIRWAY={samples:", j(runs), ",dex:", j(st$dex), ",cell:", j(st$cell),
  ",genes:["), con)
for (g in genes) {
  writeLines(paste0('{id:"', g$id, '",sym:', if (is.null(g$sym) || is.na(g$sym)) "null" else paste0('"', g$sym, '"'),
    ",bm:", g$bm, ",lfc:", g$lfc, ",p:", g$p, ",q:", g$q, ",c:[", paste(g$c, collapse = ","), "]},"), con)
}
writeLines("]};", con); close(con)
cat("genes:", length(genes), " sig(q<0.05):", sum(padj < 0.05, na.rm = TRUE), "\n")
print(head(data.frame(id = ids, lfc, padj)[order(padj), ], 8))
