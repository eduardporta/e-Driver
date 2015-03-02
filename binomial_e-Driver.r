args <- commandArgs(trailingOnly = TRUE)

e.Driver <- function (MR, TM, LR, LP) {
	result <- binom.test (MR, TM, LR/LP, alternative = "greater")
	return (result$p.value)
}

mut.data <- read.table (file = args[1], sep = "\t", header = TRUE)

mut.data$p <- mapply (e.Driver, mut.data$MutationsRegion, mut.data$TotalMutations, mut.data$LengthRegion, mut.data$LengthProtein)

tissues <- unique (mut.data$Tissue)

for (t in 1:length (tissues)) {
	
	tissue.subset <- subset (mut.data, Tissue == tissues[t])
	tissue.subset$q <- p.adjust (tissue.subset$p, method = "fdr")
	
	outfile.name <- paste (args[2], tissues[t], "with_corrected_p_values.txt", sep = "_")
	write.table (tissue.subset, file = outfile.name, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}