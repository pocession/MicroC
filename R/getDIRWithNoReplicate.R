#' GetDIRWithNoReplicate
#'
#' @description
#' This function reads interaction counts and performs differential analysis between two conditions.
#' This function is specifically for experiments without biological replicate.
#' Hence the result has no statistical significance.
#' Please treat the log2 fold change, p value, or FDR as descriptive values.
#' 
#' @author Tsunghan Hsieh
#'
#' @param chr a character string specifying the interaction regions from which chromosome 
#' @param treat a character string specifying the name and path of the treatment file (.csv)
#' @param ctrl a character string specifying the name and path of the control file (.csv)
#' @param bcv a number specifying the hypothetical biological covaraiance, default is 0.4
#' See section 2.12 
#' https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#' @param output A character string specifying the name and path of output file
#' the default path is as same as the treatment file
#' the default output file name is the treatment file name prefixed with "DIR_"
#' @return a dataframe of the output data
#'
#' @export
#'
#' @importFrom here here
#' @importFrom edgeR DGEList
#' @importFrom assertthat assert_that
#' @importFrom dplyr select mutate left_join
#' @importFrom tidyr separate
#'
#' @examples
#' \dontrun{
#' df <- GetDIRWithNoReplicate(
#' chr = "chr2",
#' treat = here::here("./Results/processing/44112_A_bg_43615_mc6contact_map_extracted.csv"),
#' ctrl = here::here("./Results/processing/44111_A_ctrl_43614_mc5contact_map_extracted.csv"),
#' bcv = 0.4,
#' output = NULL
#' )
#' }

## Get differentiation regions from microC experiments

## TODO
## Check ./Python/subset_hic.ipynb to generate the dataframe
##

GetDIRWithNoReplicate <- function(chr, treat, ctrl, bcv, output) {
  # check input ----------------------------------------------------------------
  
  assertthat::assert_that(is.character(chr), 
                          msg = "The given argument is not a string.\n")
  
  assertthat::assert_that(is.character(treat), 
                          msg = "The given argument is not a string.\n")
  
  assertthat::assert_that(grepl("\\.csv$", treat), 
                          msg = "The given file to the input file is not a .csv file.\n")
  
  assertthat::assert_that(file.exists(here::here(treat)), 
                          msg = "The given file path to the input file does not exist.\n")
  
  assertthat::assert_that(is.character(ctrl), 
                          msg = "The given argument is not a string.\n")
  
  assertthat::assert_that(grepl("\\.csv$", ctrl), 
                          msg = "The given file to the experiment file is not a .csv file.\n")
  
  assertthat::assert_that(file.exists(here::here(ctrl)), 
                          msg = "The given file path to the expertiment file does not exist.\n")
  
  assertthat::assert_that(is.numeric(bcv), 
                          msg = "The given argument is not a number.\n")
  
  # check output ---------------------------------------------------------------
  
  if (is.null(output)) {
    # Use sub() to extract the desired part
    outputfname <- sub(".*/([^.]+)\\.csv", "\\1", treat)
    outputfname <- paste0("DIR_", outputfname, ".csv")
    outputdir <- dirname(treat)
    output <- here::here(outputdir, outputfname)
    remove(outputfname, outputdir)
  } else {
    assertthat::assert_that(is.character(output), 
                            msg = "The given argument is not a string.\n")
    
    assertthat::assert_that(grepl("\\.csv$", output), 
                            msg = "The given file to the output file is not a .csv file.\n")
    
    assertthat::assert_that(file.exists(here::here(dirname(output))), 
                            msg = "The given file path to the expertiment file does not exist.\n")
  }
  
  # Read and process interaction matrix  ---------------------------------------
  treat_df <- .generateInteractionMatrixCount(treat, "treat")
  ctrl_df <- .generateInteractionMatrixCount(ctrl, "ctrl")
  
  # Combine the dataframe ------------------------------------------------------
  combine_df <- treat_df |>
    dplyr::inner_join(ctrl_df, by = "ixid")
  
  # Filter data based on abundance ---------------------------------------------
  ## Filter the last 5% interactions
  combine_df <- combine_df |>
    dplyr::mutate(avgCPM = log2(counts_treat * counts_ctrl))
  
  # Calculate the 5th percentile (least 5% value) of the avg column
  least_5_percent_value <- quantile(combine_df$avgCPM, 0.05)
  
  combine_df <- combine_df |>
    dplyr::filter(avgCPM > least_5_percent_value) |>
    dplyr::select(-c(avgCPM))
  
  # Perform differential analysis ----------------------------------------------
  ## Generate counts df for differential analysis
  ## Note if the comparison is treat vs ctrl, 
  ## in edgeR it should be written as group = ctrl:treat
  counts_df <- combine_df |>
    dplyr::select(counts_treat, counts_ctrl)
  rownames(counts_df) <- combine_df$ixid
  dir_df <- edgeR::DGEList(counts=counts_df, group=2:1) # treat vs ctrl
  
  et <- edgeR::exactTest(dir_df, dispersion=bcv^2)
  et <- as.data.frame(et)
  et$index <- rownames(et)
  
  # Generate the result dataframe  ---------------------------------------------
  result <- data.frame(chr = chr, index = combine_df$ixid)
  
  result <- result |>
    tidyr::separate(index, c("region1", "region2"), remove = FALSE) |>
    dplyr::mutate(region1 = as.numeric(region1),
                  region2 = as.numeric(region2)) |>
    dplyr::left_join(et, by = "index") |>
    dplyr::select(-c(index))
  remove(et, counts_df, ctrl_df, treat_df, combine_df)
  
  write.csv(result, output)
  message("Done! The file is written into: ", output)
  return(result)
}


#' .generateInteractionMatrixCount
#' @description
#' Get the interaction count dataframe
#' @param file a character specifying the input file
#' The inputfile is a interaction matrix
#' @param sample a character specifying the sample type, i.e. treat or control
#' @return a dataframe containing the interaction index and count
#'  
.generateInteractionMatrixCount <- function(file, sample) {
  df <- read.csv(file)
  df <- df |>
    dplyr::mutate(ixid = paste0(region1, "-", region2)) |>
    dplyr::mutate(self = ifelse((region2-region1 == 0), TRUE, FALSE)) |>
    dplyr::filter(self != TRUE) |>
    dplyr::select(ixid, region1, region2, counts)
  
  ## Change the column names 
  for (i in 2:ncol(df)) {
    colnames(df)[i] <- paste0(colnames(df)[i], "_", sample)
  }
  
  return(df)
}