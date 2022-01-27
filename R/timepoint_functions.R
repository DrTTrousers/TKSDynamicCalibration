#' Full Simulator profile extractor
#'
#' As distinct from simextract(), this version extracts the whole simulator output rather than a summary
#' @param file simulator file paths, defaults to the default expected global environment variable for simfiles
#' @param skiptorow as with simextract() this variable controls the headspace skipping until relevant cycles
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' simext_full()
#' }
simext_full <- function(file = simfile, skiptorow = 3){
  read.delim(file, skip=(skiptorow-1),
             stringsAsFactors = FALSE,
             sep = "\t") %>%
    select(-X, -X.1) %>%
    mutate_if(is.character,as.numeric) %>%
    drop_na() %>%
    group_by(Index) %>%
    summarise(across(.fns = mean))
}


#' Results long format
#'
#' Takes the calibrated results, or optionally results list and reformats it to a long data structure to facilitate
#' plotting. Uses the "Position" variable as identifiers. Labels the variables in a standard format for other functions.
#' @param df results table or list of results tables to be produced.
#' @param s indicates which side to filter the results from
#'
#' @return
#' @export
#' @importFrom tidyr pivot_longer
#'
#' @examples
#' \dontrun{
#' results_pivot()
#' }
results_pivot <- function(df = results, s = "Medial"){
  res <- bind_rows(df, .id = "Position") %>% ungroup()

  names(res)[4:7] <- c("TF", "mP", "pP", "LA")

  res_long <-  res %>%
    filter(side == s) %>%
    group_by(Index, Position) %>%
    select(-frame, -side) %>%
    summarise(across(everything(), mean), .groups = "keep") %>%
    pivot_longer(!c(Index, Position), names_to = "x", values_to = "y")
}

#' Feature extracted bar plots
#'
#' This function creates barplots to compare features at different time instances in the tekscan recording
#' @param data the results values as a long variable format as set up using results_pivot()
#' @param variable the intended variable to be examined
#' @param timepoints a vector containing the index values to be plotted on the x axis
#' @param title graph title as a string in quotes to label the plot
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' featbarplot()
#' }
featbarplot <- function(data = "data",
                        variable = "LA, mP, pP, TF",
                        timepoints = "c() vector",
                        title = "Peaks in...."){

  target <- timepoints
  if(variable == "LA"){ yat <- "Loaded Area (mm^2)"}
  else if (variable == "pP"){ yat <- "Peak Pressure (MPa)"}
  else if (variable == "mP"){ yat <- "Mean Pressure (MPa)"}
  else if (variable == "TF"){ yat <- "Total Force (N)"}
  else {stop("No Variable Specified, Function Stopped")}

  data %>%
    filter(x == variable,
           Index %in% target) %>%
    ggplot(aes(as.factor(round(Index/128, 2)), y, fill = Position))+
    geom_col(width = 0.5,
             position = position_dodge(width = 0.6))+
    theme_bw()+
    xlab("Time(s)")+
    ylab(yat)+
    ggtitle(title)
}

#' Metric labels
#'
#' A sub function to organise the labels correctly. Not to be used alone
#' @param variable variable from plot
#' @param value value from plot
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' metric_labeller()
#' }
metric_labeller <- function(variable, value){

  metric_names <- list('LA' = "Loaded Area (mm^2)",
                       'mP' = "Mean Pressure (MPa)",
                       'pP' = "Peak Pressure (MPa)",
                       'TF' = "Total Force (N)")

  return(unlist(metric_names[value]))
}


#' Data Joining
#'
#' A function to combine the results with the simdata for use plotting the simulator peaks.
#' Currently only used for simulator peaks, but could be useful for regression analysis.
#' Only combines data where the simulator and recording have the same names.
#' @param r Calibrated Tekscan results from sumdataframe(), or list.
#' @param s Simulator data as extracted by simext_full()
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' data_join()
#' }
data_join <- function(r = results, s = simdata){

res <- bind_rows(r, .id = "Position") %>% ungroup()

names(res)[4:7] <- c("TF", "mP", "pP", "LA")

AA <- bind_rows(s, .id = "Recording")

joined <-  res %>%
  filter(Position == names(r)[1]) %>%
  dplyr::select(-frame, -Position) %>%
  group_by(Index, side) %>%
  summarise(across(everything(), mean)) %>%
  left_join(.,
            AA %>% filter(Recording == names(r)[1]),
            by = "Index")

joined
}

#' Profile peaks locator
#'
#' Similar to Syncplot, this function finds the peaks in the profile data and marks them for comparison
#' @param df Data containing the simulator output
#' @param plot a true/false switch which determines if the output of this function is graphical or numerical
#' @param s a string to filter the table by for each condyle
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' prof_peaks()
#' }
prof_peaks <- function(df = data_join(), s = "Medial", plot = TRUE){

varaxis <- df %>%
           filter(side == s) %>%
           ungroup() %>%
           select(Stn.1.AF.Load, Stn.1.FE.Pos, Index) %>%
           pivot_longer(!Index, names_to = "x", values_to = "y")

af <- find_peaks(filter(varaxis, x=="Stn.1.AF.Load")$y)-1
fe <- find_peaks(filter(varaxis, x=="Stn.1.FE.Pos")$y)-1

if(plot == FALSE){
  return(list(af, fe))
}else if(is.null(plot)){
  stop("Please indicate output")
}

v <- varaxis %>%
     ggplot(aes(Index/128, y))+
     geom_line()+
     facet_wrap(~x, scales = "free",
                strip.position = "top",
                nrow = 6) +
     theme_bw()+
     theme(strip.background = element_blank(),
           strip.placement = "outside")+
     ylab(NULL)+
   geom_vline(data = filter(varaxis, x =="Stn.1.AF.Load"), aes(xintercept = af[1]/128), linetype = "dashed")+
   geom_vline(data = filter(varaxis, x =="Stn.1.AF.Load"), aes(xintercept = af[2]/128), linetype = "dashed")+
   geom_vline(data = filter(varaxis, x =="Stn.1.AF.Load"), aes(xintercept = af[3]/128), linetype = "dashed")+
   geom_vline(data = filter(varaxis, x =="Stn.1.FE.Pos"), aes(xintercept = fe[1]/128), linetype = "dashed", colour = "red")+
   geom_vline(data = filter(varaxis, x =="Stn.1.FE.Pos"), aes(xintercept = fe[2]/128), linetype = "dashed", colour = "red")+
   geom_vline(data = filter(varaxis, x =="Stn.1.FE.Pos"), aes(xintercept = fe[3]/128), linetype = "dashed", colour = "red")+
   xlab("Time (s)")

v
}

#' Tekscan plotter
#'
#' This function finds the peaks in the tekscan recording data and marks them for comparison
#' @param df Data produced by results_pivot()
#' @param plot a true/false switch which determines if the output of this function is graphical or numerical
#' @param m1 an optional list vector of peak locations generated by prof_peaks which annotates peak positions
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' tks_peaks()
#' }
tks_peaks <- function(df = results_pivot(), plot = TRUE, m1 = NULL){

if(plot == FALSE){
tf <- find_peaks(filter(df, x=="TF")$y)-1
mp <- find_peaks(filter(df, x=="mP")$y)-1
pp <- find_peaks(filter(df, x=="pP")$y)-1
la <- find_peaks(filter(df, x=="LA")$y)-1
return(list(tf, mp, pp, la))
} else if(is.null(plot)){
  stop("Please indicate output")
}

P <- df %>% filter(x != "TF") %>%
  ggplot(aes(Index/128, y, colour = Position))+
  geom_point()+
  theme_bw()+
  facet_wrap(~x, scales = "free",
             labeller = metric_labeller,
             strip.position = "top",
             nrow = 4)+
  xlab("Time (s)")+
  ylab(NULL)+
  theme(strip.background = element_blank(),
        strip.placement = "outside")

if(!is.null(m1)){
  af <- m1[[1]]
  fe <- m1[[2]]

P<- P+
    geom_vline( aes(xintercept = af[1]/128), linetype = "dashed")+
    geom_vline( aes(xintercept = af[2]/128), linetype = "dashed")+
    geom_vline( aes(xintercept = af[3]/128), linetype = "dashed")+
    geom_vline( aes(xintercept = fe[1]/128), linetype = "dashed", colour = "red")+
    geom_vline( aes(xintercept = fe[2]/128), linetype = "dashed", colour = "red")+
    geom_vline( aes(xintercept = fe[3]/128), linetype = "dashed", colour = "red")
}

P
}


#' Tekscan peaks annotator
#'
#' A further set of annotations generated from tks_peaks can be optionally annotated by generating the plot
#' from tks_peaks and passing it to this function.
#' @param P a plot generated by tks_peaks()
#' @param m2 a list vector of peaks generated by tks_peaks(plot=NULL, m1=NULL)
#' @param res_long a dataframe of the results in a long format generated by results_pivot
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' tks_annotate()
#' }
tks_annotate <- function(P = tks_peaks(m1=NULL), m2 = tks_peaks(plot=FALSE, m1=NULL), res_long = rl){

  tf <- m2[[1]]
  mp <- m2[[2]]
  pp <- m2[[3]]
  la <- m2[[4]]

  P+
  geom_vline(data = filter(res_long, x=="LA"),
             aes(xintercept = la[1]/128), colour = "darkgreen",linetype = "dashed")+
  geom_vline(data = filter(res_long, x=="LA"),
             aes(xintercept = la[2]/128), colour = "darkgreen",linetype = "dotted")+
  geom_vline(data = filter(res_long, x=="LA"),
             aes(xintercept = la[3]/128), colour = "darkgreen",linetype = "dashed")+
  geom_vline(data = filter(res_long, x=="LA"),
             aes(xintercept = la[4]/128), colour = "darkgreen",linetype = "dashed")+
  geom_vline(data = filter(res_long, x=="mP"),
             aes(xintercept = mp[1]/128), colour = "purple", linetype = "dashed")+
  geom_vline(data = filter(res_long, x=="mP"),
             aes(xintercept = mp[2]/128), colour = "purple", linetype = "dashed")+
  geom_vline(data = filter(res_long, x=="mP"),
             aes(xintercept = mp[3]/128), colour = "purple", linetype = "dashed")+
  geom_vline(data = filter(res_long, x=="mP"),
             aes(xintercept = mp[4]/128), colour = "purple", linetype = "dashed")+
  geom_vline(data = filter(res_long, x=="pP"),
             aes(xintercept = pp[1]/128), colour = "Orange", linetype = "dashed")+
  geom_vline(data = filter(res_long, x=="pP"),
             aes(xintercept = pp[2]/128), colour = "Orange", linetype = "dashed")+
  geom_vline(data = filter(res_long, x=="pP"),
             aes(xintercept = pp[3]/128), colour = "Orange", linetype = "dashed")+
  geom_vline(data = filter(res_long, x=="pP"),
             aes(xintercept = pp[4]/128), colour = "Orange", linetype = "dashed")
}
