#' GetInteractionPosIndex
#'
#' @description
#' This function generates the position index for the interested genomic segment and write into a .txt file.
#' The position index is the mid-point of the genomic regions
#' For example, if there is a interaction between two genomic regions: 11000-16000 and 18000-23000
#' 
#' @author Tsunghan Hsieh
#'
#' @param chr a character string specifying the interested segment from which chromosome 
#' @param start a integer string specifying starting position of the interested segment
#' @param end a integer string specifying starting position of the interested segment
#' @param res a integer specifying the resolution 
#' The above parameters should be same as ./Python/subset_hic_data.py
#' @param output A character string specifying the name and path of output file
#' @return a dataframe of the output data
#'
#' @export
#'
#' @importFrom here here
#' @importFrom assertthat assert_that
#' @importFrom dplyr select mutate left_join
#' @importFrom tidyr separate
#'
#' @examples
#' \dontrun{
#' df <- GetInteractionPosIndex(
#' chr = "chr2",
#' start = 112735986,
#' end = 113204585,
#' res = 5000,
#' output = here::here("./Results/processing/start_position_index.txt")
#' )
#' }

GetInteractionPosIndex <- function(chr, start, end, res) {
  
  # check input ----------------------------------------------------------------
  
  assertthat::assert_that(is.character(chr), 
                          msg = "The given argument is not a string.\n")
  
  assertthat::assert_that(is.integer(start), 
                          msg = "The given argument is not a integer\n")
  
  assertthat::assert_that(is.integer(end), 
                          msg = "The given argument is not a integer\n")
  
  assertthat::assert_that(is.integer(res), 
                          msg = "The given argument is not a integer\n")
  
  assertthat::assert_that(file.exists(here::here(output)), 
                          msg = "The given file path to the output file does not exist.\n")
  
  fragment_num <- round((end - start)/res)
  fragment_index <- c()
  for (i in 1:(fragment_num-1)) {
    new_start <- start + 1 + (2*i-1)*round(res/2)
    fragment_index <- c(fragment_index, new_start)
  }
  
  # Handle the last index ------------------------------------------------------
  ## If the whole length of interested segment is not divisible by res
  ## Then the length between the last and last - 1 index is longer than the resolution
  last_interval <- round((end - fragment_index[(length(fragment_index)-1)])/2)
  fragment_index[length(fragment_index)] <- fragment_index[(length(fragment_index)-1)] + last_interval + 1
  
  # Write the final data -------------------------------------------------------
  data <- as.character(fragment_index)
  
  # Write the lines to the file
  writeLines(data, output)
}
