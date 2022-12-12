range <- snakemake@params[[1]]

r <- as.vector(unlist(strsplit(range, ",")), "numeric")

thresh <- c()

for (i in 1:length(r))
{
  thresh[i] <- paste(r[i], " 0 ", r[i])
  
}

writeLines(thresh, snakemake@output[[1]])