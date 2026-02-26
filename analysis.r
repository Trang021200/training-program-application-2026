# ---------------------------------------------------------

# Melbourne Bioinformatics Training Program

# This exercise to assess your familiarity with R and git. Please follow
# the instructions on the README page and link to your repo in your application.
# If you do not link to your repo, your application will be automatically denied.

# Leave all code you used in this R script with comments as appropriate.
# Let us know if you have any questions!


# You can use the resources available on our training website for help:
# Intro to R: https://mbite.org/intro-to-r
# Version Control with Git: https://mbite.org/intro-to-git/

# ----------------------------------------------------------

# Load libraries -------------------
# You may use base R or tidyverse for this exercise
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
setwd("/Users/conchim/training-program-application-2026/")

# ex. library(tidyverse)

# Load data here ----------------------
# Load each file with a meaningful variable name.
GSE60450.meta <- read.csv("data/GSE60450_filtered_metadata.csv", check.names = FALSE,
                          stringsAsFactors = FALSE)

colnames(GSE60450.meta)[1] <- "id"

GSE60450.gene.expression <- read.csv("data/GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv",
                                          check.names = FALSE, stringsAsFactors = FALSE)
colnames(GSE60450.gene.expression)[1] <- "gene_id"

# Inspect the data -------------------------

# What are the dimensions of each data set? (How many rows/columns in each?)
# Keep the code here for each file.

## Expression data
dim(GSE60450.gene.expression)
#[1] 23735    14

## Metadata
dim(GSE60450.meta)
#[1] 12  4

##Prepare/combine the data for plotting ------------------------
##How can you combine this data into one data.frame?

#Pivot the expression data 
GSE60450.gene.expression.pivot <- GSE60450.gene.expression %>%
  pivot_longer(
    cols = starts_with("GSM"),
    names_to = "id",
    values_to = "expression_level"
  )

#Merge the metadata and gene expression by the column id 

GSE60450.merge <- merge(GSE60450.meta, GSE60450.gene.expression.pivot, by = "id",
                        all = T)

# Plot the data --------------------------
## Plot the expression by cell type

#Log transformation the the expression level 
GSE60450.merge <- GSE60450.merge %>% 
  mutate(log_expression_level = log2(expression_level + 1)) # add 1 to take care of 0 value

## Can use boxplot() or geom_boxplot() in ggplot2
#Reorder the characteristic 
char_order <- c("mammary gland, luminal cells, virgin",
                "mammary gland, basal cells, virgin",
                "mammary gland, luminal cells, 18.5 day pregnancy",
                "mammary gland, basal cells, 18.5 day pregnancy",
                "mammary gland, luminal cells, 2 day lactation",
                "mammary gland, basal cells, 2 day lactation")

Plot.expression.by.cell.type <- ggplot(GSE60450.merge, 
                                       aes(x=factor(characteristics, levels = char_order), 
                                           y=log_expression_level, fill = immunophenotype)) + 
  geom_boxplot() +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, max(GSE60450.merge$log_expression_level, na.rm = TRUE), by = 2)) + 
  theme_bw() +
  labs(
    x = "Cell type",
    y = "Log10 Gene Expression"
  )

Plot.expression.by.cell.type
## Save the plot
### Show code for saving the plot with ggsave() or a similar function
ggsave(
  "results/Expression_by_cell_type.png",
  plot = Plot.expression.by.cell.type,
  width = 8,
  height = 3,
  dpi = 300
)
