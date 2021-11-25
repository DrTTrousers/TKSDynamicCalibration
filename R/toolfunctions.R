# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#Extractor functions

simextract <- function(file = simfile, skiptorow = 3){
  read.delim(file, skip=(skiptorow-1), stringsAsFactors = FALSE, sep = "\t") %>%
    select(Index, Stn.1.AF.Load) %>%
    mutate_if(is.character,as.numeric) %>%
    drop_na() %>%
    group_by(Index) %>%
    summarise(Mean_AF = mean(Stn.1.AF.Load))
}

tekextract <- function(file = tekfile, headerrows = headspace){
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
} #used as part of indexer

sampleextract <- function(file = samplefile, headerrows = headspace){
  read_csv(file,
           col_names = FALSE,
           col_types = cols(`0` = col_number(),
                            X1 = col_number(),
                            X3 = col_number()),
           skip = headerrows) %>%
    drop_na() %>%
    mutate(frame = rep(1:(length(X1)/26), each=26))
} #used all over the place

#indexing and peak matching

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

indexer <- function(sim = simextract(),
                    tek = tekextract(),
                    peaknum = 1, #identifies which peak to use based on its position in the output vector of find_peaks
                    modifier = 0){ #a custom correction variable for debugging, may not still be needed
  offset <- (trunc(find_peaks(tek$Sum_Raw, m=15)[peaknum] - (find_peaks(sim$Mean_AF, m=15)[1])/2)-modifier) #calculates an offset to line up the peaks in one data set with the other
  set <- tek %>%
    mutate(Index = c(seq((128-(offset*2)), 127, 2),
                     rep(seq(0,127,2),
                         times = 4),
                     seq(0, ((126-(offset*2))), 2)))
  tek_sum <<- set %>% #this creates a sequence of index values based off the peak matching so all data is identified.
    group_by(Index) %>%
    summarise(Mean_Raw = mean(Sum_Raw), .groups = "keep")
  tek_ind <<- set %>%
    group_by(Index) %>%
    summarise(frame = frame, .groups = "keep") %>%
    mutate(frame = as.factor(frame))
  #Outputs controlled by superassignment THIS IS A PROBLEM. Refactor outputs with IF statement and call twice, once for each assignment?
}

#model generation

calmodel <- function(tek = tek_sum, sim = sim_sum){
  lm(Mean_AF ~0+poly(Mean_Raw, 3, raw = T),
     data=merge(tek, sim, by="Index"))
}

#predicting functions

predictor <- function(file = sampleextract(samplefile), modin = model){
  file %>% # Reduce the overheads for repeatedly importing the same thing.
    select_if(negate(is.character)) %>%
    group_by(frame) %>%
    summarise_all(sum) %>%
    mutate(Mean_Raw = rowSums(.[2:45])) %>%
    select(-(contains("X"))) %>%
    mutate(Cal_Sum = predict(modin, .))
}

calibrate <- function(file = sampleextract(samplefile), modselect = model){
  merge(file, predictor(file, modin = modselect),
        by = "frame") %>%
    select_if(negate(is.character)) %>%
    mutate_at(vars(contains("X")),funs(./Mean_Raw)) %>%
    mutate_at(vars(contains("X")),funs(.*Cal_Sum)) %>%
    mutate_at(vars(contains("X")),funs(./1.6129))
}


#Plotting functions

tekplot2 <- function(file = sampleextract(samplefile),
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


  pltly1 <- ggplotly(p1)
  pltly2l <- ggplotly(p2_l)
  pltly2m <- ggplotly(p2_m)
  pltly3 <- ggplotly(p3)
  filler <- plotly_empty()
  text <- plot_ly()%>%
    add_text(x = 1,
             y = 1,
             yshift = 1,
             text = paste("Gait<br> Index:<br>", index_labels$Index),
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

#tabulating functions

sumdataframe <- function(file = sampleextract(samplefile),
                         tek = tekextract(samplefile),
                         sim = sim_sum,
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

sumtab <- function(file = sampleextract(samplefile),
                   calmodel = model,
                   knee = "rtkn"){
  l <- "Lateral"
  m <- "Medial"
  k <- quo(c(rep(m, times=(length(variable)/2)),
             rep(l, times=(length(variable)/2))))

  if (knee == "ltkn"){
    k <- quo(c(rep(l, times=(length(variable)/2)),
               rep(m, times=(length(variable)/2))))
  }
  library(DT)
  datatable(
    calibrate(file, calmodel) %>%
      mutate(rownum = rep(1:26, times=(length(X1)/26))) %>%
      select(-Mean_Raw, -Cal_Sum) %>%
      mutate(frame = as.factor(frame)) %>%
      melt(id.vars =c("frame", "rownum")) %>%
      mutate(side = !! k) %>%
      group_by(frame, side) %>%
      summarise("Total Force (N)" = (sum(value[value!=0]*1.6129)),
                "Average Pressure (MPa)" = mean(value[value!=0]),
                "Peak Pressure (MPa)" = max(value),
                "Loaded area (mm^2)" = (length(value[value!=0])*1.6129)) %>%
      mutate_all(~replace(., is.nan(.), 0)) %>%
      mutate_if(is.numeric, round, 3) %>%
      left_join(., tek_ind, by="frame") %>%
      rename("Tibial Side" = side,
             "Tekscan Frame" = frame),
    filter = "top",
    rownames = FALSE)
}
