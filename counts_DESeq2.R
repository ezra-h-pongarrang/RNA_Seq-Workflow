suppressPackageStartupMessages(library(DESeq2))

design_file <- "design.csv"
counts_file <- "counts.csv" #counts file from parsing output
output_file <- "results.csv"

# prepare colData
# for matrix input: a DataFrame or data.frame with at least a single column. 
# Rows of colData correspond to columns of countData
colData <- read.csv(design_file, stringsAsFactors = F)

# factorise colData$condition
colData$condition <- factor(colData$condition)

# make sure to relevel and that the first factor is the first entry in the file
colData$condition <- relevel(colData$condition, toString(colData$condition[1]))
sample_name <- colData$sample

# countData for matrix input: a matrix of non-negative integers
df <- read.csv(counts_file)

# housekeeping and making sure (later) countData has the matching column names 
#as colData
# round the supposedly countData
countData <- round(df[ ,sample_name])
otherCols <- df[!((names(df)%in%sample_name))]

# run DESeq
#
#
# DESeqDataSetFromMatrix(countData, colData, design).
dsd <- DESeqDataSetFromMatrix(countData, colData, design = ~condition) 
dse <- DESeq(dsd)
result <- results(dse)


# Table formatting
#
#
res_data <- cbind(otherCols, data.frame(result))
# add foldchange
res_data$foldChange = 2^res_data$log2FoldChange
# rename p adjusted column to FDR. Real padj will be made later
names(res_data)[names(res_data)=="padj"] <- "FDR"
# rename pvalue to P Value
names(res_data)[names(res_data)=="pvalue"] <- "PValue"

# Calculate real p adjusted
res_data$PAdj=p.adjust(res_data$PValue, method = "hochberg")

# FALSE DISCOVERY
# Sort the data by PValue to compute false discovery counts.
res_data = res_data[with(res_data, order(PValue, -foldChange)), ]
# Compute the false discovery counts on the sorted table.
# rowIndex*FD rate
res_data$falsePos = 1:nrow(res_data) * res_data$FDR

# Create the additional columns that we wish to present.
res_data$baseMeanA = 1
res_data$baseMeanB = 1

# Get the normalized counts
normed = counts(dse, normalized=TRUE)

# Round normalized counts to a single digit.
normed = round(normed, 1)

# Merge the two datasets by row names
final <- merge(res_data, normed, by=0)
# Sort again for output.
final = final[with(final, order(PValue, -foldChange)), ]

# Sample names for condition A
col_names_A = data.frame(split(colData, colData$condition)[1])[,1]

# Sample names for condition B
col_names_B = data.frame(split(colData, colData$condition)[2])[,1]

# Create the individual baseMean columns.
final$baseMeanA = rowMeans(final[, col_names_A])
final$baseMeanB = rowMeans(final[, col_names_B])

# round for the wellbeing of humanity
final$baseMean = round(final$baseMean, 1)
final$log2FoldChange = round(final$log2FoldChange, 1)
final$lfcSE = round(final$lfcSE, 2)
final$stat = round(final$stat, 2)
final$FDR = round(final$FDR, 4)
final$falsePos = round(final$falsePos, 0)
final$baseMeanA= round(final$baseMeanA, 1)
final$baseMeanB = round(final$baseMeanB, 1)
final$foldChange = round(final$foldChange, 3)

# Reformat these columns as string
final$PAdj = formatC(final$PAdj, format = "e", digits = 1)
final$PValue = formatC(final$PValue, format = "e", digits = 1)

# Rename the first column
colnames(final)[1] <- "name"

# Reorganize columns names to make more sense
new_col = c("name", names(otherCols), "baseMean","baseMeanA","baseMeanB","foldChange",
             "log2FoldChange","lfcSE","stat","PValue","PAdj", "FDR","falsePos",col_names_A, col_names_B)

# Slice the dataframe with new columns.
final = final[, new_col]

# write the result into a stdout
write.csv(final, file=output_file, quote = F, row.names = F)

print(paste("Output file:", output_file))
