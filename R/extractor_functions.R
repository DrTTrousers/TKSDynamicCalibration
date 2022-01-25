#' Simulator AF summary data extractor
#'
#' This function extracts the axial force data from a single station knee simulator file
#' @param file a file path as a tab-delimited .txt.
#'             Defaults to an empty variable for the string
#' @param skiptorow a number which indicates the row to begin reading from
#'
#' @return
#' @export
#' @importFrom grDevices rainbow
#' @importFrom graphics frame
#' @importFrom stats lm predict
#' @importFrom utils read.delim
#' @examples
#' \dontrun{
#' simextract()
#' }
simextract <- function(file = simfile, skiptorow = 3){
  read.delim(file, skip=(skiptorow-1), stringsAsFactors = FALSE, sep = "\t") %>%
    select(Index, Stn.1.AF.Load) %>%
    mutate_if(is.character,as.numeric) %>%
    drop_na() %>%
    group_by(Index) %>%
    summarise(Mean_AF = mean(Stn.1.AF.Load))
}


#' Tekscan raw summary extractor
#'
#' This function extracts the raw data from a tekscan ascii file and creates an aggregate
#' @param file a file path as a comma separated .csv.
#'             Defaults to an empty variable for the string
#' @param headerrows rows to skip removing pre-amble. Defaults to 24 for an unsaved file
#'
#' @return
#' @export
#' @examples
#' \dontrun{
#' tekextract()
#' }
tekextract <- function(file = tekfile, headerrows = 24){
  read_csv(file,
           col_names = FALSE,
           col_types = cols(.default = "n"),
           skip = headerrows) %>%
    select_if(~!all(is.na(.))) %>%
    drop_na() %>%
    mutate(frame = rep(1:(length(X1)/26), each=26)) %>%
    group_by(frame) %>%
    summarise_all(sum) %>%
    mutate(Sum_Raw = rowSums(.[grep("X", names(.))], na.rm = TRUE)) %>%
    select(frame, Sum_Raw)
}


#' Tekscan full frame extractor
#'
#' This function extracts the raw data from a tekscan ascii file and creates a dataframe
#' retaining spatial information.
#' @param file a file path as a comma separated .csv.
#'             Defaults to an empty variable for the string
#' @param headerrows rows to skip removing pre-amble. Defaults to 24 for an unsaved file
#'
#' @return
#' @export
#' @examples
#' \dontrun{
#' sampleextract()
#' }
sampleextract <- function(file = samplefile, headerrows = 24){
  read_csv(file,
           col_names = FALSE,
           col_types = cols(`0` = col_number(),
                            X1 = col_number(),
                            X3 = col_number()),
           skip = headerrows) %>%
    drop_na() %>%
    mutate(frame = rep(1:(length(X1)/26), each=26))
}

