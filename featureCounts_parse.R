# Transform featureCounts output into a simple count data frame

# load gplots package
suppressPackageStartupMessages(library(gplots))

# assign counts_file as the variable for featureCounts output from unix.
# "counts.txt" is hard-coded; please suit according to featureCounts output file name
counts_file <- "counts.txt"

#design.csv contains the sample and condition, useful for output table formatting.
#for example
# sample, condition
# UHR_1, UHR
# UHR_2, UHR
# UHR_3, UHR
# HBR_1, HBR ...etc...
design_file <- "design.csv"

# counts outpute file name
output_file <- "counts.csv"

# import counts file (output from featureCounts) 
table <- read.table(counts_file, header = T)

# subset counts file dataframe 
df <- table[ ,c(1, 7:length(names(table)))]

# import design_file (containing root names of each file)
design <- read.csv(design_file)
names(df) <- c("names", design[ ,"sample"])
 
# write the result into a stdout
write.csv(df, file=output_file, row.names = F, quote = F )

# inform the user
print(paste("# Output file :", output_file))

