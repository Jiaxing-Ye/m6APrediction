#' @title Encode DNA sequences
#' @description Convert DNA sequences into categorical nucleotide-position features.
#' @param dna_strings A character vector of equal-length DNA sequences using uppercase A, T, C and G.
#' @return A data frame with one row per sequence and columns named nt_pos1, nt_pos2, and so on. Each column is a factor with levels A, T, C and G.
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' @title Predict multiple m6A sites
#' @description Predict m6A probabilities and classes for multiple sites using a fitted random forest model.
#' @param ml_fit A fitted randomForest classification model.
#' @param feature_df A data frame containing gc_content, RNA_type, RNA_region,
#'   exon_length, distance_to_junction, evolutionary_conservation and DNA_5mer.
#'   Length and distance inputs must already be log2-transformed.
#'   DNA_5mer must contain uppercase five-base DNA sequences with A at the centre.
#' @param positive_threshold A numeric classification threshold, defaulting to 0.5.
#'   Scores strictly greater than this threshold are classified as Positive;
#'   all other scores are classified as Negative.
#' @return The input data frame with two additional columns:
#'   predicted_m6A_prob (numeric) and predicted_m6A_status (character).
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction", mustWork = TRUE))
#' feature_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, feature_df[1:3, ], positive_threshold = 0.6)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' @title Predict a single m6A site
#' @description Predict the m6A probability and class for one site using a fitted random forest model.
#' @param ml_fit A fitted randomForest classification model.
#' @param gc_content A numeric GC-content value between 0 and 1.
#' @param RNA_type An RNA type: mRNA, lincRNA, lncRNA or pseudogene.
#' @param RNA_region An RNA region: CDS, intron, 3'UTR or 5'UTR.
#' @param exon_length A numeric exon length already transformed using log2.
#' @param distance_to_junction A numeric distance to the junction already transformed using log2.
#' @param evolutionary_conservation A numeric evolutionary conservation score.
#' @param DNA_5mer An uppercase five-base DNA sequence with A at the centre.
#' @param positive_threshold A numeric classification threshold, defaulting to 0.5.
#'   Scores strictly greater than this threshold are classified as Positive;
#'   all other scores are classified as Negative.
#' @return A named character vector containing predicted_m6A_prob and predicted_m6A_status.
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction", mustWork = TRUE))
#' prediction_single(
#'   ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
