# Load the library.
suppressPackageStartupMessages(library(gplots))

results_file <- "results.csv"
output_file <- "heatmap.pdf"
# Inform the user.
print("# Tool: Create Heatmap ")
print(paste("# Input: ", results_file))
print(paste("# Output: ", output_file))

#set minimun FDR
min_FDR <- 0.05

# Plot width
WIDTH = 12

# Plot height
HEIGHT = 13

# Set the margins
MARGINS = c(9, 12)

# Relative heights of the rows in the plot.
LHEI = c(1, 5)

# Import deseq2 output data
data <- read.csv(results_file,  header = T, as.is = T)

#subset data for values less than threshold set
data = subset(data, data$FDR <= min_FDR)

# The heatmap row names will be taken from the first column
row_names = data[ ,1]


# The code assumes that the normalized data matrix is listed to the right of the falsePos column.
idx = which(colnames(data) == "falsePos") + 1

# The normalized counts are on the right size (as mentioned in lin 35).
norm_counts = data[, idx : ncol(data)]

# Load the data from the second column on.
values = as.matrix(norm_counts)

# Adds a little noise to each element to avoid the
# clustering function failure on zero variance rows.
values = jitter(values, factor = 1, amount = 0.00001)

# Normalize each row to a z-score
zscores = NULL
for (i in 1 : nrow(values)) {
  row = values[i,]
  zrow = (row - mean(row)) / sd(row)
  zscores = rbind(zscores, zrow)
}

# Set the row names on the zscores.
row.names(zscores) = row_names

# Turn the data into a matrix for heatmap2.
zscores = as.matrix(zscores)

# Set the color palette.
col = greenred

# Create a PDF device
pdf(output_file, width = WIDTH, height = HEIGHT)

heatmap.2(zscores, col=col, density.info="none", Colv=NULL,
          dendrogram="row", trace="none", margins=MARGINS, lhei=LHEI)

#dev.off()
