hist_logy <- function(x) {
  hd <- hist(x,plot=F);
  hd$counts <- log10(hd$counts);
  plot(hd, ylab='log10(Frequency)', axes=F);
  axis(1);
  cnts <- seq(floor(min(hd$counts)), floor(max(hd$counts)), length.out=4);
  axis(2,at=cnts,labels=10**cnts);
}

############################
## Multiple plot function
############################
##
## ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
## - cols:   Number of columns in layout
## - layout: A matrix specifying the layout. If present, 'cols' is ignored.
##
## If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
## then plot 1 will go in the upper left, 2 will go in the upper right, and
## 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## vertical x
wzGGverticalX <- theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))

########
## IO
########

wzPDF <- function(fn, dwidth=0, dheight=0) {
  pdf(fn, width = 3.5 + dwidth, height = 3.5 + dheight)
  par(mar=c(5,4,1,1))
}

## PDF is measured in "in" (inches). We should be consistent.
wzPNG <- function(fn, dwidth=0, dheight=0) {
  png(fn, width = 3.5 + dwidth, height = 3.5 + dheight, unit='in', res=300)
  par(mar=c(4,4.5,1,1),
      cex.axis=1.2, cex.lab=1.3)
}

###################
## Formatting
###################

## Exponential Label in Log Scale
## wzExpAxis in Base Graphics
## ... includes: cex.axis=,  padj=, las=,
wzBaseExpAxis <- function(at, ax='x', base=10, ...) {
  if (ax == 'x') ax <- 1;
  if (ax == 'y') ax <- 2;
  axis(ax, at, sapply(at, function(x) parse(text=paste0(base, '^', x))), ...)
}

## This function exists for controling the distance of the label to the axis
## Unfortunately, in R, there is no simple way to move label away from the axis
## This is particularly necessary when the axis label is exponential
## This should be used with xlab='' in the main plot function
## line= is what you need to move label away from axis
wzBaseAxisLabel <- function(lab, ax='x', line=3, ...) {
  if (ax == 'x') ax <- 1;
  if (ax == 'y') ax <- 2;
  if (!hasArg(cex)) cex = par('cex.lab');
  if (!hasArg(cex)) font = par('font.lab');
  mtext(side=ax, text=lab, line=line, cex=cex, font=font, ...)
}


#######################################
## Scatter plot with coloring density
#######################################

wzPlotDens <- function(x1,x2,pch=16,...) {
  df <- data.frame(x1,x2)
  x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  plot(x2~x1, data=df[order(df$dens),], pch=pch, col=col, ...)
}



###############################################
################## OBSOLETE ###################
###############################################

#########################
## 1-D density plot
#########################
plotdens <- function(x, normalize=FALSE, add=FALSE, ...) {
  d <- density(na.omit(x))
  if (normalize) {
    d$y <- d$y / max(d$y)
  }
  if (add) {
    lines(d, ...)
  } else {
    plot(d, ...)
  }
}

wzSmoothDensity <- function(x, y, xlim=c(-2,2), nrpoints=100, nbins=256, ...) {
    palette <- colorRampPalette(
        c("white","lightblue","blue","green","yellow","orange","red","darkred"),
        space = "Lab")

    smoothScatter(x, y, xlim = xlim,
        nrpoints=nrpoints,
        nbin=c(nbins,nbins),
        colramp=palette, col='blue', ...)
}

