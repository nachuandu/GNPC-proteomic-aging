# Portable adapter of the recovered preprocessing script; no cohort files bundled.
# Usage: Rscript code/preprocess_gnpc.R raw_matrix.rds NEW_OUTPUT_DIRECTORY
# Input must already contain the approved, QC-filtered control subset, with
# rows=samples and columns=assays, uniquely named. This script does not infer CN.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Expected matrix.rds and new output directory")
if (!requireNamespace("impute", quietly=TRUE)) stop("Install Bioconductor impute in the approved environment")
if (dir.exists(args[2])) stop("Output directory already exists")
x <- as.matrix(readRDS(args[1]))
storage.mode(x) <- "double"
if (is.null(rownames(x)) || is.null(colnames(x)) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) stop("Unique sample and assay names required")
x[x == -1] <- NA_real_
if (any(x < 0, na.rm=TRUE) || any(is.infinite(x))) stop("Unexpected non-RFU values")
keep <- colMeans(is.na(x)) <= 0.50
if (!any(keep)) stop("No assays remain")
x <- x[, keep, drop=FALSE]
# Recovered GNPC script used log2(RFU+1). This is distinct from the later
# IADRC adapter's log2(RFU); do not silently replace one with the other.
xlog <- log2(x + 1)
set.seed(1)
ximp <- t(impute::impute.knn(t(xlog), k=5)$data)
if (any(!is.finite(ximp))) stop("Imputation left nonfinite entries")
mu <- colMeans(ximp)
sdv <- apply(ximp, 2, sd)
z <- scale(ximp)
z[, !is.finite(sdv) | sdv == 0] <- 0
dir.create(args[2], recursive=TRUE)
write.csv(data.frame(sample_id=rownames(z), z, check.names=FALSE), file.path(args[2], "PRIVATE_standardized_matrix.csv"), row.names=FALSE)
write.csv(data.frame(feature=colnames(ximp), mean_log2=mu, sd_log2=sdv), file.path(args[2], "PRIVATE_scaling.csv"), row.names=FALSE)
saveRDS(ximp, file.path(args[2], "PRIVATE_imputed_log2.rds"))
writeLines(capture.output(sessionInfo()), file.path(args[2], "session_info.txt"))
