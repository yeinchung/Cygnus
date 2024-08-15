# Load necessary library
library(ComplexHeatmap)

# Load the data
data <- read.csv("/Users/yeinchung/Downloads/all_BG.csv")

# Basic data exploration
hist(data$Number.MarkerPositive.Evs)
print(colnames(data))
print(unique(data$stage))

# Subset data based on stage
f1 <- data[data$stage == "F1", ][, 2:11]
t4 <- data[data$stage == "T4-1", ][, 2:11]
t4 <- data[data$stage == "T4-1", ][, 2:11]

# Check if the subsets are created correctly
print(dim(f1))
print(dim(t4))

# Define make_comb_mat if not already defined
if (!exists("make_comb_mat")) {
  make_comb_mat <- function(df) {
    comb_matrix <- as.matrix(df)
    return(comb_matrix)
  }
}

# Convert data to combination matrix format
m_f1 <- make_comb_mat(f1)
m_t4 <- make_comb_mat(t4)

# Check the dimensions and structure of the combination matrices
print(dim(m_f1))
print(dim(m_t4))
print(head(m_f1))
print(head(m_t4))

# Function to calculate deviation and p-values
calc_deviation_and_p <- function(comb_matrix, total_size) {
  results <- data.frame(
    Intersection = character(ncol(comb_matrix)),
    Deviation = numeric(ncol(comb_matrix)),
    Expected = numeric(ncol(comb_matrix)),
    Observed = numeric(ncol(comb_matrix)),
    P_Value = numeric(ncol(comb_matrix)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(ncol(comb_matrix))) {
    I <- comb_size(comb_matrix)[i]  # number of elements in the intersection I
    S_plus <- which(comb_matrix[, i] == 1)  # sets contained
    S_minus <- which(comb_matrix[, i] == 0)  # sets not contained in the intersection
    
    # calculate products for sets in S_plus and S_minus
    if (length(S_plus) > 0) {
      prod_plus <- prod(sapply(S_plus, function(j) set_size(comb_matrix)[j] / total_size))
    } else {
      prod_plus <- 1
    }
    
    if (length(S_minus) > 0) {
      prod_minus <- prod(sapply(S_minus, function(j) 1 - set_size(comb_matrix)[j] / total_size))
    } else {
      prod_minus <- 1
    }
    
    # calculate expected value and deviation
    expected_I <- prod_plus * prod_minus * total_size
    deviation <- I - expected_I
    
    # calculate p-value using binom.test
    result <- binom.test(x = I, n = total_size, p = prod_plus * prod_minus, alternative = "two.sided")
    p_value <- result$p.value
    
    results$Intersection[i] <- paste(i)
    results$Deviation[i] <- deviation
    results$P_Value[i] <- p_value
    results$Expected[i] <- expected_I
    results$Observed[i] <- I 
  }
  
  return(results)
}



# Calculate deviations and p-values for F1
result_f1 <- calc_deviation_and_p(m_f1, nrow(f1))
result_t4 <- calc_deviation_and_p(m_t4, nrow(t4))
print(result_f1)




# Visualize it in bar plot
library(ggplot2)

result_f1$Intersection <- factor(result_f1$Intersection, levels = result_f1$Intersection)

ggplot(result_f1, aes(x = Intersection, y = Deviation)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Deviation of Intersections F1") +
  xlab("Intersection") +
  ylab("Deviation")

ggplot(result_f1, aes(x = Intersection, y = -log10(P_Value))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("P-values of Intersections F1 (-log10 scale)") +
  xlab("Intersection") +
  ylab("-log10(P_Value)")

library(ggplot2)

ggplot(result_f1, aes(x = Intersection, y = -log10(P_Value), fill = Deviation)) +
  geom_bar(stat = "identity") +
  theme_minimal()+
  scale_fill_gradientn(limits = c(-max(result_f1$Deviation),max(result_f1$Deviation)), colors = rev(brewer.pal(9, 'RdYlBu'))) +  # Adjust colors as needed
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("P-values of Intersections F1 (-log10 scale) Colored by Deviation") +
  xlab("Intersection") +
  ylab("-log10(P_Value)") + 
  geom_hline(yintercept=-log(0.05/29), linetype="dashed", color = "red") 

library(ComplexHeatmap)

Heatmap(m_f1, name = "Combinations F1", 
        show_row_names = FALSE, 
        show_column_names = FALSE, 
        cluster_rows = F, 
        cluster_columns = F)

pheatmap(m_f1, cluster_rows = F, cluster_cols = F, color = c("lightgrey", "purple"), border_color = "lightgrey")

# 3-4 markers independent -> p -> how to interpret 
# Running simulation 
set.seed(123) 
num_rows <- 1000

# Create each column with the desired values
A <- sample(c(1, 0), num_rows, replace = T, prob = c(0.5, 0.5))
B <- sample(c(1, 0), num_rows, replace = T, prob = c(0.5, 0.5))
C <- sample(c(1, 0), num_rows, replace = T, prob = c(0.5, 0.5))

# Combine the columns into a data frame
exp <- data.frame(A, B, C)

# Randomly reorder columns
exp_reordered <- exp[, sample(ncol(exp))]

exp_reordered
m_exp <- make_comb_mat(exp)
result_exp <- calc_deviation_and_p(m_exp, nrow(exp))
print(result_exp)

ggplot(result_f1, aes(x = as.numeric(Intersection), y = -log10(P_Value), fill = Deviation)) +
  geom_bar(stat = "identity") +
  theme_minimal()+
  scale_fill_gradientn(colors = rev(brewer.pal(9, 'RdYlBu'))) +  # Adjust colors as needed
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("P-values of Intersections F1 (-log10 scale) Colored by Deviation") +
  xlab("Intersection") +
  ylab("-log10(P_Value)") + 
  geom_hline(yintercept=-log10(0.05/nrow(result_f1)), linetype="dashed", color = "red") 

pheatmap(m_f1, cluster_rows = F, cluster_cols = F, color = c("lightgrey", "orange"), border_color = "lightgrey")


ggplot(result_t4, aes(x = as.numeric(Intersection), y = -log10(P_Value), fill = Deviation)) +
  geom_bar(stat = "identity") +
  theme_minimal()+
  scale_fill_gradientn(limits = c(-35, 150), colors = rev(brewer.pal(9, 'RdYlBu'))) +  # Adjust colors as needed
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("P-values of Intersections T1 (-log10 scale) Colored by Deviation") +
  xlab("Intersection") +
  ylab("-log10(P_Value)") + 
  geom_hline(yintercept=-log10(0.05/nrow(result_t4)), linetype="dashed", color = "red") 

pheatmap(m_t4, cluster_rows = F, cluster_cols = F, color = c("white", "lightgrey"), border_color = "NA")

i = 1
for(i in 1:dim(result_t4)[1]){
  result_t4[["Set_names"]][i] <- comb_name(m_t4[, i])
}


# only show top 10 
result_t4_top100 <- result_t4[order(result_t4$P_Value), ][1:100, ]

result_t4_top100$P_Value <- as.character(-log10(result_t4_top100$P_Value))
result_t4_top100$P_Value <- factor(result_t4_top100$P_Value, levels=unique(result_t4_top100$P_Value))
ggplot(result_t4_top100, aes(x = Intersection, y = P_Value, fill = Deviation)) +
  geom_bar(stat = "identity") +
  theme_minimal()+
  scale_fill_gradientn(limits = c(-35, 150), colors = rev(brewer.pal(9, 'RdYlBu'))) +  # Adjust colors as needed
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("P-values of Intersections T1 (-log10 scale) Colored by Deviation") +
  xlab("Intersection") +
  ylab("-log10(P_Value)") + 
  geom_hline(yintercept=-log10(0.05/nrow(result_t4_top100)), linetype="dashed", color = "red") 





### reordering
# Load required libraries
# Load required libraries
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# Calculate the cutoff value
cutoff <- -log10(0.05 / nrow(result_t4))

# Filter the data frame to include only values above the cutoff
result_t4_filtered <- result_t4[-log10(result_t4$P_Value) > cutoff, ]

# Order the filtered data frame by -log10(P_Value) in decreasing order
result_t4_filtered <- result_t4_filtered[order(-log10(result_t4_filtered$P_Value), decreasing = TRUE), ]

# Exclude the highest value
result_t4_filtered <- result_t4_filtered[-1, ]

# Create the bar plot with ggplot2 and reverse x order
ggplot(result_t4_filtered, aes(x = reorder(as.factor(Intersection), log10(P_Value)), y = -log10(P_Value), fill = Deviation)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = rev(brewer.pal(9, 'RdYlBu'))) +  # Adjust colors as needed
  ggtitle("P-values of Intersections T1 (-log10 scale) Colored by Deviation") +
  xlab("Intersection") +
  ylab("-log10(P_Value)") + 
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "red") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

# Create the heatmap with pheatmap (if you want to subset this data too, modify accordingly)
pheatmap(m_t4, cluster_rows = FALSE, cluster_cols = FALSE, color = c("white", "lightgrey"), border_color = NA)

# subsetting to make the heatmap 
result_t4_filtered[["match"]] <- paste0("V", rownames(result_t4_filtered))
new_t4 <- as.data.frame(m_t4)
new_t4 <- new_t4[, result_t4_filtered$match]
pheatmap(new_t4, cluster_rows = FALSE, cluster_cols = FALSE, color = c("white", "lightgrey"), border_color = NA)


ggplot(result_t4_filtered, aes(x = reorder(as.factor(Intersection), log10(P_Value)), y = Observed/nrow(t4), group = 1)) +
  geom_area(fill = "red", alpha = 0.5)+
  theme_minimal()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        panel.grid = element_blank())

ggplot(result_t4_filtered, aes(x = reorder(as.factor(Intersection), log10(P_Value)), y = Observed/nrow(t4))) +
  geom_bar(stat = "identity", fill = "red", alpha = 0.3) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        panel.grid = element_blank())

ggplot(result_t4_filtered, aes(x = reorder(as.factor(Intersection), log10(P_Value)), y = (Observed-Expected)/nrow(t4))) +
  geom_bar(stat = "identity", fill = "red", alpha = 0.3) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),
        panel.grid = element_blank())




# make histogram 
sizes = as.data.frame(set_size(m_t4) )
# Load necessary libraries
library(ggplot2)

# Create the data frame
sizes <- data.frame(
  Gene = c("EpCAM", "MET", "SDC1", "EGFR", "SP.B", "CTSH", "PDL1", "ADAM10", "MUC1", "HER2"),
  Set_Size = c(3552, 4568, 993, 4969, 1989, 4827, 5247, 5191, 5188, 3663)
)

# Preserve the order of the Gene factor levels
sizes$Gene <- factor(sizes$Gene, levels = sizes$Gene)

# Create the bar plot
ggplot(sizes, aes(x = Gene, y = Set_Size)) +
  geom_bar(stat = "identity", fill = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Set Sizes of Various Genes") +
  xlab("Gene") +
  ylab("Set Size")


EnhancedVolcano(result_t4[-626,],
                lab = rownames(result_t4[-626,]),
                x = 'Deviation',
                y = 'P_Value', 
                pCutoff = 0.05 / nrow(result_t4), xlab = "Deviation", )
