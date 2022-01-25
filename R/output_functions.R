#' Pressure map calibration
#'
#' This function uses the model to calculate the MPa value of each sensel in the pressure map.
#' Used as a helper function for other outputs, but can be called to produce results as a dataframe.
#' @param file TekScan output to be calibrated. Defaults to the output of sampleextract()
#' @param modselect the model to be used for prediction. Defaults to the output of autofit() from the GlobalEnvironment
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' calibrate()
#' }
calibrate <- function(file = sampleextract(), modselect = model){
  prd <- file %>%
         select_if(negate(is.character)) %>%
         group_by(frame) %>%
         summarise_all(sum) %>%
         mutate(Mean_Raw = rowSums(.[2:45])) %>%
         select(-(contains("X"))) %>%
         mutate(Cal_Sum = predict(modselect, .))

  map <- merge(file, prd, by = "frame") %>%
         select_if(negate(is.character)) %>%
         mutate_at(vars(contains("X")),funs(./Mean_Raw)) %>%
         mutate_all(~replace(., is.nan(.), 0)) %>%
         mutate_at(vars(contains("X")),funs(.*Cal_Sum)) %>%
         mutate_at(vars(contains("X")),funs(./1.6129))

  return(map)
}

#' Animated Tekscan Plot
#'
#' This is the big one. Takes the TekScan dynamic recording, calibrates it according to model.
#' Produces as an output raster plots of the recording including summary statistics, labels and
#' centre-of-force markers.
#' Outputs a complex plotly object with several animated graphs.
#' Requires a large amount of memory to Knit to HTML.
#' @param file TekScan recording to be analysed. Defaults to the output of sampleextract()
#' @param calmodel Model for calibration of data. Defaults to the output of autofit()
#' @param knee right or left knee, as this dictates labelling
#' @param index_labels output from the Indexer function or autofit. Indicates gait progression on each frame.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' tekplot2()
#' }
tekplot2 <- function(file = sampleextract(),
                     calmodel = model,
                     knee = "rtkn",
                     index_labels = tek_ind ){ #Need to incorporate the indexing function to label the index step.
  l <- "Lateral"
  m <- "Medial"
  k <- quo(c(rep(m, times=(length(variable)/2)),
             rep(l, times=(length(variable)/2))))
  fac <- c(m, l)

  if (knee == "ltkn"){
    k <- quo(c(rep(l, times=(length(variable)/2)),
               rep(m, times=(length(variable)/2)))) #does not work as a reversed vector as a quosure?
  }
  if (knee == "ltkn"){
    fac <- c(l, m)
  }

  c <- calibrate(file, calmodel) %>%
    mutate_all(~replace(., is.nan(.), 0))%>%
    select(-Mean_Raw, -Cal_Sum) %>%
    mutate(rownum = rep(1:26, times=(length(X1)/26)))

  p1<- c %>% group_by(frame) %>%
    select(-rownum) %>%
    summarise(across(where(is.numeric), max)) %>% #each column max value per frame.
    melt(id.vars="frame") %>%
    ggplot(aes(as.numeric(variable), value, frame=frame))+
    geom_point()+
    theme_bw()+
    ylim(0, 10)

  c2 <- c %>% mutate(rownum = rep(1:26, times=(length(X1)/26))) %>%
    melt(id.vars=c("frame", "rownum")) %>%
    mutate(side = !! k)


  p2_l<- c2 %>%
    filter(side==l) %>%
    select(-side) %>%
    group_by(frame, rownum) %>%
    summarise(across(where(is.numeric), max)) %>% #each row max value per frame
    ggplot(aes(as.numeric(rownum), value, frame=frame))+
    geom_point()+
    coord_flip()+
    theme_bw()+
    ylim(0, 10)

  p2_m<- c2%>%
    filter(side==m) %>%
    select(-side) %>%
    group_by(frame, rownum) %>%
    summarise(across(where(is.numeric), max)) %>% #each row max value per frame
    ggplot(aes(as.numeric(rownum), value, frame=frame))+
    geom_point()+
    coord_flip()+
    theme_bw()+
    ylim(0, 10)

  x <- c2 %>% select(-rownum) %>% group_by(frame, variable, side) %>%
    summarise(across(where(is.numeric), max)) %>%
    group_by(frame, side) %>%
    slice_max(value, n=1) %>%
    slice_head(n=1) %>%  # Sepearate out a vector of the column locations for maxes per frame and side.
    select(-value)

  y <- c2 %>% select(-variable) %>% group_by(frame, rownum, side) %>%
    summarise(across(where(is.numeric), max)) %>%
    group_by(frame, side) %>%
    slice_max(value, n=1) %>%
    slice_head(n=1)%>%  # Sepearate out a vector of the row locations for maxes per frame and side.
    select(-value)

  cof <- merge(x, y, by=c("frame", "side")) %>%
    rename(xcof = variable,
           ycof = rownum)

  p3<- merge(c2, cof, by=c("frame", "side")) %>%
    mutate(side = factor(side, levels = fac)) %>%
    ggplot(aes(variable, rownum, fill=value, frame=frame))+
    geom_raster()+
    geom_point(aes(xcof, ycof),colour = "black", fill = "white")+
    facet_wrap("side", scales = "free")+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())+
    scale_fill_gradientn(colours = rev(rainbow(5)),
                         name = "MPa",
                         limits = c(0, 10))

  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )


  pcnt <- index_labels$Index/127*100

  pltly1 <- ggplotly(p1)
  pltly2l <- ggplotly(p2_l)
  pltly2m <- ggplotly(p2_m)
  pltly3 <- ggplotly(p3)
  filler <- plotly_empty()
  text <- plot_ly()%>%
    add_text(x = 1,
             y = 1,
             yshift = 1,
             text = paste("Gait<br> Progress:<br>", pcnt, "%"),
             frame = index_labels$frame,
             showlegend = FALSE,
             textfont = list(size = 11)) %>%
    layout(xaxis = ax, yaxis = ax)

  layout <- list(pltly2m, pltly3, pltly2l, filler, pltly1, text)
  if (knee == "ltkn"){
    layout <- list(pltly2l, pltly3, pltly2m, filler, pltly1, text)
  }

  subplot(layout, nrows = 2,
          heights = c(0.8, 0.2), widths = c(0.1, 0.8, 0.1),   margin = 0.02, #shareX = TRUE,
          shareY = TRUE,
          titleX = FALSE, titleY = FALSE) %>%
    animation_opts(frame = 750,
                   transition = 0,
                   redraw = TRUE,
                   easing = "quad",
                   mode = "next") %>% partial_bundle()

}


#' Syncronisation Plot
#'
#' Graphically and interactively illustrates the synchronisation of TekScan output with Simulator output.
#' Useful graphical tool to check autofit() has correctly optimised the variables.
#' Best combined with modplot() as two halves of the same analysis, using plotly.
#' Can be partial bundled.
#' @param x TekScan summary from indexer() or autofit()
#' @param title A string to label the plot
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' syncplot()
#' subplot(syncplot(), modplot(), nrows=2) %>% partial_bundle()
#' }
syncplot <- function(x = tek_sum, title = "Peak Synchronisation Check"){
  ggplotly(
    x %>% melt (id.vars = "Index") %>%
      ggplot(aes(Index, value, colour = variable))+
      geom_line(size = 1)+
      geom_vline(xintercept = find_peaks(sim_sum$Mean_AF, m=15)[1], linetype = "dotted", colour = "blue", size = 0.8)+
      geom_vline(xintercept = find_peaks(sim_sum$Mean_AF, m=15)[2], linetype = "dotted", colour = "blue", size = 0.8)+
      geom_vline(xintercept = find_peaks(sim_sum$Mean_AF, m=15)[3], linetype = "dotted", colour = "blue", size = 0.8)+
      scale_colour_brewer(palette = "Set1", name = "Profile")+
      theme_bw()+
      ylab("Raw Sum")+
      xlab("Index")+
      ggtitle(title)
  )
}

#' Model fit plot
#'
#' Graphically and interactively illustrates the fit of the calibration curve.
#' Useful graphical tool to check autofit() has correctly minimised R squared and fit the data
#' Best combined with syncplot() as two halves of the same analysis, using plotly.
#' Can be partial bundled.
#' @param tek TekScan summary from tidytek() or autofit()
#' @param sim Simulator summary from simextract()
#' @param modselect model output of calmodel() or autofit()
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' modplot()
#' subplot(syncplot(), modplot(), nrows=2) %>% partial_bundle()
#' }
modplot <- function(tek = tek_sum, sim = sim_sum, modselect = model){
  merged <- merge(tek, sim, by="Index") %>%
    mutate(mod_out = predict(modselect, .))
  ggplotly(
    ggplot(data=merged)+
      geom_point(aes(Mean_AF, Mean_Raw, colour = Index))+
      geom_line(aes(mod_out, Mean_Raw), size =1, colour = "red")+
      theme_bw()+
      ylab("Raw Sum")+
      xlab("Axial Force")+
      ggtitle("Model Plot")
  )
}


#' Output Table
#'
#' Summary of the calibrated data as a table. Suited for storage without rounding for future use
#' @param file TekScan recording. Defaults to the output of sampleextract()
#' @param knee right or left knee, dictates labelling
#' @param calibrationmodel model for calibration as in calibrate()
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' sumdataframe()
#' }
sumdataframe <- function(file = sampleextract(),
                         knee = "rtkn",
                         calibrationmodel = model){
  l <- "Lateral"
  m <- "Medial"
  k <- quo(c(rep(m, times=(length(variable)/2)),
             rep(l, times=(length(variable)/2))))

  if (knee == "ltkn"){
    k <- quo(c(rep(l, times=(length(variable)/2)),
               rep(m, times=(length(variable)/2))))
  }

  calibrate(file, calibrationmodel) %>%
    mutate(rownum = rep(1:26, times=(length(X1)/26))) %>%
    select(-Mean_Raw, -Cal_Sum) %>%
    mutate(frame = as.factor(frame)) %>%
    melt(id.vars =c("frame", "rownum")) %>%
    mutate(side = !! k) %>%
    group_by(frame, side) %>%
    summarise("TotalForce(N)" = (sum(value[value!=0])*1.6129),
              "AveragePressure(MPa)" = mean(value[value!=0]),
              "PeakPressure(MPa)" = max(value),
              "LoadedArea(mm^2)" = (length(value[value!=0])*1.6129)) %>%
    group_by(side) %>%
    mutate_all(~replace(., is.nan(.), 0)) %>%
    left_join(., tek_ind, by="frame")
}

#' Results table
#'
#' Summary of the calibrated data as a table. Suited for presentation with rounding and detailed headers
#' Includes slicers for interaction after knitting
#' Could be a wrapper function for sumdataframe?
#' @param df a dataframe to be tidied up for knit, defaults to sumdataframe() output.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' sumtab(sumdataframe())
#' }
sumtab <- function(df = sumdataframe()){

  datatable(df %>%
      mutate_if(is.numeric, round, 3) %>%
      rename("Tibial Side" = side,
             "Tekscan Frame" = frame),
    filter = "top",
    rownames = FALSE)
}
