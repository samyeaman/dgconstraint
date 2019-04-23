#' Calculations for a Single Generation Functionin Matrix Form
#'
#' This function allows you to calculate the pairwise C-score using the hypergeometric approach, a p-value for 'all lineages' contrast using chi-square, and the estimates of the effective proportion of adaptive loci for a data set with a single generation.
#'
#' @param paper the data in csv that you want to analyze, in a folder named data-in
#' @param environment The environment in which the experiment occured
#' @param species Specify if the organism is "Sac" or "Ecoli_K12" or "Ecoli_O157-H7", or manually input the gene count of your species
#' @param population a list of populations in the data
#' @return a table with all the calculated infromation
#' @export
#' @examples
 #singlematrix("Author2018","YPD", "Sac", c("P1", "P2", "P3" ,"P4", "P5"))
#####################

singlematrix <- function(paper, environment, species, population, numGenes = NA){

# library(tidyverse)
# library(readr)
# library(devtools)
# library(dgconstraint)
# library(Hmisc)

geneNumbers <- read_csv(file.path(getwd(),"data-in/GeneDatabase.csv"))

data <- read_csv(file.path(getwd(), "data-in", paste0(paper, ".csv")))

if (species %in% geneNumbers$Species){
  numGenes <- filter(geneNumbers, Species == species)$NumGenes
}

if(is.na(numGenes)){
  prompt <- "Your species is unspecified or not in our database. How many genes does it have? \n"
  numGenes <- as.numeric(readline(prompt))
}


  data.1 <- data %>%
  arrange(Gene) %>%
  drop_na(Gene) %>%
  drop_na(population)
data.2 <- data.1 %>%
  select(population)

data.matrix <- as.matrix(data.2)
rownames(data.matrix) <-data.1$Gene

num_parallel <- data.frame(data.matrix, Count=rowSums(data.matrix, na.rm = FALSE, dims = 1), Genes = row.names(data.matrix))

genes_parallel <- num_parallel %>%
  as_tibble() %>%
  filter(Count > 1)


Non_genes_parallel <- num_parallel %>%
  as_tibble() %>%
  filter(Count == 1)


num_parallel_genes <- nrow(genes_parallel)
num_non_parallel_genes <- nrow(Non_genes_parallel)
total_genes <- num_non_parallel_genes + num_parallel_genes
parallel_genes <- paste0(genes_parallel$Genes, collapse=", ")

full_matrix <- rbind(data.matrix, array(0,c(numGenes-total_genes,ncol(data.matrix))))


c_hyper <- append(c_hyper, pairwise_c_hyper(full_matrix))
p_chisq <- append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200))
estimate <- append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T))

c_hyper[c_hyper <= 0] <- 0
c_hyper[c_hyper == "NaN"] <- 0

df <- tibble( paper = paper, environment = environment, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)

newdir <- file.path(getwd(), "data-out")
if (!file.exists(newdir)){
  dir.create(newdir, showWarnings = FALSE)
  cat(paste("\n\tCreating new directory: ", newdir), sep="")
}

filename <- file.path(getwd(), "data-out", paste(paper, "_Analysis.csv", sep=""))
write.csv(df, file=filename, row.names=FALSE)

}
