#' plot_manhattan
#'
#' This function takes the p-value output from GWASbyGLM and combines it with SNP positional information to generate a manhattan plot using the 'qqman' package
#'
#' @param mdpSNPinfo a data.frame with 3 columns containing snp name, chromosome, and basepair position. SNPs must be in the same order as the mdpN genotype matrix used for GWASbyGLM
#' @param GWASout a vector containing p-values generated with the GWASbyGLM function
#' @param error a value from 0 to 1 indicating what type-1 error threshold should be used for calculating a Bonferroni cutoff. Default set to 0.05.
#'
#' @return a manhattan plot with a horizontal line representing the Bonferroni correction cutoff.
#'
#' @examples
#' # Reading example data from the Github repository
#' GWASout <-read.csv(text=getURL("https://raw.githubusercontent.com/mtmcgowan/CROPS545_HW4_GLM_GWAS/master/example_GWASout"), header=T, sep = " ", row.names = 1)
#' mdpSNPinfo <-read.csv(text=getURL("https://raw.githubusercontent.com/mtmcgowan/CROPS545_HW4_GLM_GWAS/master/example_mdpSNPinfo"), header=T, sep = " ")
#'
#' # Generating the manhattan plot
#' plot_manhattan(mdpSNPinfo, GWASout)

plot_manhattan <- function(mdpSNPinfo, GWASout, error = 0.05) {
  # Generate the proper frame for the manhattan function
  manhattan_frame <- data.frame(mdpSNPinfo$Chromosome, mdpSNPinfo$Position, GWASout, mdpSNPinfo$SNP)
  names(manhattan_frame) <- c("CHR", "BP", "P", "SNP")

  # Calculate how many markers/tests were performed
  testct <- nrow(manhattan_frame)

  # Calculate the Bonferroni cutoff
  bonferroni <- -log10(error/testct)

  # Generate the manhattan plot
  manhattan(manhattan_frame, suggestiveline = FALSE, genomewideline = bonferroni)
}
