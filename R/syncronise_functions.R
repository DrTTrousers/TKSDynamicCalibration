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
#' find_peaks()
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
#' collected from the global environment.
#' @param sim Output of the sim extractor function. Defaults to calling simextract()
#' @param tek Output of the TekScan extractor function. Defaults to calling tekextract()
#' @param peaknum identifies which peak to use based on its position in the output vector of find_peaks
#' @param modifier a custom correction variable for minor adjustment of the position
#'
#' @return
#' @export
#'
#' @examples
#' indexer()
indexer <- function(sim = simextract(),
                    tek = tekextract(),
                    peaknum = 1,
                    modifier = 0){
  offset <- (trunc(find_peaks(tek$Sum_Raw, m=15)[peaknum] - (find_peaks(sim$Mean_AF, m=15)[1])/2)-modifier)
  set <- tek %>%
    mutate(Index = c(seq((128-(offset*2)), 127, 2),
                     rep(seq(0,127,2),
                         times = 4),
                     seq(0, ((126-(offset*2))), 2)))
  set_mr <- set %>%
                    group_by(Index) %>%
                    summarise(Mean_Raw = mean(Sum_Raw), .groups = "keep")
  set_fr <- set %>%
                    group_by(Index) %>%
                    summarise(frame = frame, .groups = "keep") %>%
                    mutate(frame = as.factor(frame))

  assign("tek_sum", set_mr, envir = .GlobalEnv)
  assign("tek_ind", set_fr, envir = .GlobalEnv)
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
#' calmodel(tek_sum, sim_sum)
calmodel <- function(tek = tek_sum, sim = sim_sum){
  lm(Mean_AF ~0+poly(Mean_Raw, 3, raw = T),
     data=merge(tek, sim, by="Index"))
}

#' Autofitting function for calibration model
#'
#' This function iteratively controls the indexer parameters to optimise the calibration model
#' There are two loops, one which sets the peak number between 1 and 4,
#' and another which increments the modifier up to 20.
#' Stops when the calibration model fits to an R squared of >0.99.
#' Has three outputs to the global environment:
#' model
#' tek_sum
#' tek_ind
#' @param t Output of tekextract()
#' @param s Output of simextract()
#'
#' @return
#' @export
#'
#' @examples
#' autofit(tek_sum, )
autofit <- function(t = tekextract(),
                    s = simextract()){

  minRsq <- 0
  peak <- 0
  mod <- 0
  tex <- t
  simx <- s

while(minRsq < 0.90){
  peak <- peak+1
  tryCatch({
    indexer(sim = simx, tek = tex, peaknum = peak, modifier = mod)

    model <- calmodel(tek_sum, simx)

    minRsq <- summary(model)$adj.r.squared},

    error=function(e){})

  if(peak == 5){
    message("Peak search complete without optimising R^2")
    break
  }
  message(paste("Tuning complete. Adjusted R squared is:", minRsq, "At peak number:", peak))
}

while(minRsq < 0.99){
  mod <- mod+1
  tryCatch({
    indexer(sim = simx, tek = tex, peaknum = peak, modifier = mod)

    model <- calmodel(tek_sum, simx)

    minRsq <- summary(model)$adj.r.squared},

    error=function(e){})

  if(mod == 20){
    message("Modifier search complete without optimising R^2")
    break
  }
  message(paste("Tuning complete. Adjusted R squared is:", minRsq, "At peak number:", peak, "And a +", mod, "Modifier"))
}
  assign("model", model, envir = .GlobalEnv)
}
