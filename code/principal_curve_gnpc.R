# Adapter of the historical R princurve chat script. This is not Slingshot.
# Usage: Rscript code/principal_curve_gnpc.R embedding.csv NEW_OUTPUT_DIRECTORY
# Do not substitute this for the later Python IADRC curve backend.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Expected embedding.csv and new output directory")
if (!requireNamespace("princurve", quietly=TRUE)) stop("Install princurve")
if (dir.exists(args[2])) stop("Output directory already exists")
ph <- read.csv(args[1], colClasses=c(sample_id="character"))
stopifnot(all(c("sample_id", "age", "PHATE1", "PHATE2") %in% names(ph)))
if (anyDuplicated(ph$sample_id)) stop("Duplicate sample_id")
X <- as.matrix(ph[, c("PHATE1", "PHATE2")])
if (any(!is.finite(X)) || any(!is.finite(ph$age))) stop("Nonfinite data")
set.seed(1)
fit <- princurve::principal_curve(X)
if (!isTRUE(fit$converged)) stop("Principal curve did not converge")
raw <- fit$lambda
if (diff(range(raw)) <= 0) stop("Degenerate curve ordering")
score <- (raw - min(raw)) / diff(range(raw))
if (cor(score, ph$age, method="spearman") < 0) score <- 1 - score
# The recovered early script oriented scores negatively. This adapter writes
# the positive-with-age orientation used by the manuscript's downstream scores.
dir.create(args[2], recursive=TRUE)
write.csv(data.frame(sample_id=ph$sample_id, age=ph$age, proteomic_aging_score=score),
          file.path(args[2], "PRIVATE_scores.csv"), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(args[2], "session_info.txt"))
