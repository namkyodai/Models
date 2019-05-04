stackedPlot <- function(data, time=NULL, col=1:length(data),...) {
  if (is.null(time)) {
    time <- 1:length(data[[1]])
    plot(0,0, xlim = range(time), ylim = c(0,max(rowSums(data))), axes=FALSE,t="n" , ...)
    axis(2, ylim=c(0,100),col="black",las=1)  ## las=1 makes horizontal labels
    axis(1,c(seq(1,YearMax,by=1)),c(seq(1,YearMax,by=1)))
    box()
    
    for (i in length(data):1) {
      #the summup to the current collumn
      prep.data <- rowSums(data[1:i]);
      # The polygon must have his first and last point on the zero line
      prep.y <- c(0, prep.data,0)
      prep.x <- c(time[1], time, time[length(time)] )
      polygon(prep.x, prep.y, col=col[i], border = NA);
    }
  }
}