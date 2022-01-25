#' Peak finding helper function
#'
#' This function finds peaks in data and returns the numerical position in the vector
#' @param x values to be checked
#' @param m the range of values either side of each peak which must be different to identify a peak
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' find_peaks()
#' }
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

#' Index matching function
#'
#' This function uses the find_peaks function to match tekscan output and tekscan frames to the simulator index
#' This function has two outputs and therefore cannot be used directly but must be called alone then the outputs
#' collected from the global environment. This is used automatically in the autofit program, but can be run manually.
#' @param sim Output of the sim extractor function. Defaults to calling simextract()
#' @param tek Output of the TekScan extractor function. Defaults to calling tekextract()
#' @param numcycle a variable telling the operation how many cycles were recorded
#' @param peaknum identifies which peak to use based on its position in the output vector of find_peaks
#' @param mod a custom correction variable for minor adjustment of the position.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' indexer()
#' }
indexer <- function(sim = simextract(),
                    tek = tekextract(),
                    numcycle = 5,
                    peaknum = 1,
                    mod = 0){

  t_peak <- find_peaks(tek$Sum_Raw, m=15)[peaknum]
  s_peak <- find_peaks(sim$Mean_AF, m=15)[1] - mod

  frms <- length(tek$frame)
  stepval <- 128*numcycle/frms
  mind <- 127

  a <- (s_peak -(t_peak*stepval))

  if(s_peak%%stepval ==1){
    s_peak <- s_peak+1
  }

  if (a>=1){

    heads <- seq(a, s_peak, stepval)
  }else if(a<1){

    if(a%%stepval ==1){
      a <- a+stepval+1
    }
    b <- (mind+a)

    if(b%%stepval == 1){
      b <- b+1
    }
    heads <- c(seq(b, mind, stepval),
               seq(0, s_peak, stepval))
  }

  tails <- seq(s_peak, mind, stepval)

  rests <- seq(0, mind, stepval)

  Index <- vector(mode = "numeric", length = frms)
  Index[1:t_peak] <- heads

  y <- (t_peak + length(tails))-1
  Index[t_peak:y] <- tails

  ss <- 1:y
  Index[-ss] <- rep(rests)

  tek_indexed <- tek %>% mutate("Index" = Index)
  tek_indexed
}

#' Calibration model generator
#'
#' Performs a polynomial regression to generate an equation converting the raw TekScan values into Newtons.
#' required to calibrate the output of the TekScan
#' @param tek The synchronised TekScan summary output. One output of the indexer() function
#' @param sim The output of simextract()
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' calmodel(tek_sum, sim_sum)
#' }
calmodel <- function(tek = tek_sum, sim = sim_sum){
  lm(Mean_AF ~0+poly(Mean_Raw, 3, raw = T),
     data=merge(tek, sim, by="Index"))
}

#' Tidy tekscan data extractor
#'
#' Wrapper function using the indexer tool output to produce important summary data frames.
#' Has two possible outputs including a summary of the total value for each frame, or a each frame with its corresponding index values
#' @param x Output of the indexer function, regardless of if the peaks match
#' @param SumOrInd  string serving as a switch to produce either output
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' tidytek()
#' }
tidytek <- function(x = indexer(),
                    SumOrInd = "S"){
  t <- tolower(SumOrInd)
  if(t == "s"){
    x %>% group_by(Index) %>%
      summarise(Mean_Raw = mean(Sum_Raw), .groups = "keep")
  }else if (t == "i"){
    x %>% group_by(Index) %>%
      summarise(frame = frame, .groups = "keep") %>%
      mutate(frame = as.factor(frame))
  }else{
    message("Please choose either summary or index for output")
  }
}


#' Autofit function for calibration model
#'
#' This function iteratively controls the indexer parameters to optimise the calibration model
#' There are two loops, one which sets the peak number between 1 and the value of attempts,
#' and another which searches a small adjustment modifier from -10 to +10.
#' Stops when the calibration model fits to an R squared of >the target value.
#' Has three outputs to the global environment (using the tidytek and calmodel functions):
#' model
#' tek_sum
#' tek_ind
#' @param t Output of the tekextract() function to feed into the indexer() function
#' @param s Output of the simextract() function to feed into the indexer() function
#' @param attempts The number of different peaks the matching algorithm should iterate through
#' @param cycles The number of cycles recorded, assumed to be 5
#' @param target The target R squared value to meet or exceed before returning
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' autofit()
#' }
autofit <- function(t = tekextract(),
                    s = simextract(),
                    attempts = 6,
                    cycles = 5,
                    target = 0.97){

  minRsq <- 0
  peak <- 0
  fudge <- 0

  while(minRsq < target-0.1){
    peak <- peak+1
    tryCatch({
      tek_sum <- tidytek(indexer(sim = s, tek = t, peaknum = peak, mod = fudge, numcycle = cycles),"s")
      tek_ind <- tidytek(indexer(sim = s, tek = t, peaknum = peak, mod = fudge, numcycle = cycles),"i")

      model <- calmodel(tek_sum, s)

      minRsq <- summary(model)$adj.r.squared},

      error=function(e){})

    if(peak == attempts){
      stop(paste("Peak search tried", peak, "peaks without optimising R^2. Check Simulator file corresponds to TekScan file!"))
    }

    message(paste("Adjusted R squared is:", minRsq, "At peak number:", peak))
  }

  if(minRsq < target){
    fudge <- -10

    while(minRsq < target){

      tryCatch({
        tek_sum <- tidytek(indexer(sim = s, tek = t, peaknum = peak, mod = fudge, numcycle = cycles),"s")
        tek_ind <- tidytek(indexer(sim = s, tek = t, peaknum = peak, mod = fudge, numcycle = cycles),"i")

        model <- calmodel(tek_sum, s)

        minRsq <- summary(model)$adj.r.squared},

        error=function(e){})

      fudge <- fudge+1

      if(fudge == 10){
        stop(paste("Peak search tried", peak, "peaks and", fudge, "adjustments without optimising R^2. Check Simulator file corresponds to TekScan file!"))
      }
    }
  }

  message(paste("Adjusted R squared is:", minRsq, "At peak number:", peak, "Using a", fudge, "adjustment"))

  assign("model", model, envir = .GlobalEnv)
  assign("tek_sum", tek_sum, envir = .GlobalEnv)
  assign("tek_ind", tek_ind, envir = .GlobalEnv)
}
