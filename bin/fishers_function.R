kegg <- c()
id <- c()
stats <- c()

for (i in 1:nrow(b))
{
	numbers <- c(b[i,3], b[i,4], b[i,5], b[i,6])
	conting <- matrix(numbers, nrow = 2, ncol = 2, byrow = T)
	test = fisher.test(conting, alternative = "greater")
	stats[i] <- test$p
	kegg[i] <- b[i,1]
	id[i] <- formatC(b[i,2],width=5,format="d",flag="0")
}

output <- matrix(nrow = length(kegg), ncol = 4)
for (i in 1:length(stats))
{
	output[i,1] <- kegg[i]
	output[i,2] <- id[i]
	output[i,3] <- stats[i]
}

output[,4] = p.adjust(stats,method="fdr")

output[,3] = round(as.numeric(output[,3]),digits=5)
output[,4] = round(as.numeric(output[,4]),digits=5)

