library(RColorBrewer)
library(tidyverse)

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


# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#######################################
## 2D Scatter plot with coloring density
#######################################
## old name wzPlotDens
wzPlotDens2d <- function(x1,x2,pch=16,colorstops=c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"),...) {
  df <- data.frame(x1,x2)
  x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  cols <-  colorRampPalette(colorstops)(256)
  df$col <- cols[df$dens]
  plot(x2~x1, data=df[order(df$dens),], pch=pch, col=col, ...)
}


#########################
## 1-D density plot
#########################
## old plotdens
wzPlotDens1d <- function(x, normalize=FALSE, bw="nrd0", add=FALSE, ...) {
    d <- density(na.omit(x), bw=bw)
    if (normalize) {
        d$y <- d$y / max(d$y)
    }
    if (add) {
        lines(d, ...)
    } else {
        plot(d, ...)
    }
}

wzPlotDens1d.fromMatrixGG <- function(x) {
    xx <- melt(x)
    ggplot(xx) + geom_density(aes(value, color=Var2))
}

wzPlotDens1d.fromMatrix <- function(x, bw="nrd0", normalize=FALSE, ...) {
    for (i in seq_len(ncol(x))) {
        if (i==1) {
            add <- FALSE;
        } else {
            add <- TRUE;
        }
        wzPlotDens1d(x[,i], normalize=normalize, bw=bw, add=add, ...);
    }
}

## this is the better way
wzPlotDens1d.SE <- function(se, color_var = NULL) {
    ad = assays(se)[[1]]
    cd = as_tibble(colData(se))
    add = melt(ad, value.name="betas", varnames=c("Probe_ID", "Sample_ID"))
    if (is.null(color_var)) {
        ggplot(add) + geom_density(aes(betas, group = "Sample_ID"))
    } else {
        add[[color_var]] = cd[[color_var]][match(add$Sample_ID, cd[["Sample_ID"]])]
        ggplot(add) + geom_density(aes_string("betas", group = "Sample_ID", color=color_var))
    }
}

## old name: wzSmoothDensity
wzPlotDens2d.smooth <- function(x, y, xlim=c(-2,2), nrpoints=100, nbins=256, ...) {
    palette <- colorRampPalette(
        c("white","lightblue","blue","green","yellow","orange","red","darkred"),
        space = "Lab")

    smoothScatter(x, y, xlim = xlim,
        nrpoints=nrpoints,
        nbin=c(nbins,nbins),
        colramp=palette, col='blue', ...)
}

## try to limit ncol(mtx) to < 8
wzPlotDens2d.smoothPairs <- function(mtx) {
    palette <- colorRampPalette(
        c("white","lightblue","blue","green","yellow","orange","red","darkred"),
        space = "Lab")
    
    ## Correlation panel
    upper <- function(x, y){
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- round(cor(x, y, method='spearman', use = "na.or.complete"), digits=2)
        txt <- paste0("R = ", r)
        cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
    }
    ## Customize upper panel
    lower <- function(x, y){
        smoothScatter(x, y, add=TRUE, nrpoints = 0, colramp = palette, col='blue')
        abline(0,1,lty='dotted')
    }
    ## histogram (source: psych::pairs.panel)
    panel.hist = function(x, ...) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y,col="grey")
    }
    ## not used but could try
    panel.hist.density = function(x,...) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y,col="grey")
        tryd <- try( d <- density(x,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
        if(class(tryd) != "try-error") {
            d$y <- d$y/max(d$y)
            lines(d)
        }
    }
    ## Create the plots
    pairs(mtx,
        diag.panel = panel.hist.density,
        lower.panel = lower,
        upper.panel = upper)
}


wzGetColors <- function(grouping, palette.name='Paired', group_name=TRUE) {
    groups <- unique(grouping)
    if (palette.name == 'Set1') nclr = 9
    else if (palette.name == 'Dark2') nclr = 8
    else nclr = 12

    colors = colorRampPalette(brewer.pal(nclr, palette.name))(length(groups))
    if (group_name == TRUE) {
        setNames(colors, groups)
    } else {
        colors
    }
}

wzRepelLabel <- function(
    grouping, group2color,
    segment.square = FALSE,
    segment.inflect = FALSE,
    segment.curvature = 0.2,
    nudge_x = 0.1) {

    df <- data.frame(pos=sapply(split(seq_along(grouping), grouping), mean)) %>%
        rownames_to_column(var='grouping')

    library(ggrepel)
    ggplot(df) +
        geom_text_repel(aes(y=pos, x=1, label=grouping),
            color=group2color[df$grouping],
            segment.square = segment.square,
            segment.inflect = segment.inflect,
            segment.size = 0.6,
            segment.curvature = segment.curvature,
            force = 0.2, nudge_x = nudge_x,
            direction = "y", hjust = 0) +
        theme_minimal() +
        theme(axis.line.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.title.x = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            ) +
        xlim(1,2.3) + scale_y_reverse() + ylab('')
}

ggVerticalX <- function() {
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
ggPercentage <- function() {
    scale_y_continuous(labels=scales::percent)
}
## vertical x NOTE: use ggVerticalX above
wzGGverticalX <- theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))

# example from https://stackoverflow.com/questions/35324892/ggplot2-setting-geom-bar-baseline-to-1-instead-of-zero
## ggplot(dat, aes(x, ratio, fill = ifelse(ratio > 1,"GT1","LT1"))) +
##   geom_bar(stat="identity") +
##   scale_fill_manual(values=c("blue","red"), name="LT or GT 1") +
##   scale_y_continuous(trans = shift_trans(1))
shift_trans = function(d = 0) {
  scales::trans_new("shift", transform = function(x) x - d, inverse = function(x) x + d)
}


plotSesameSpeciesAUC <- function(sdf, title = NULL) {
    addr <- sesameDataGet(sprintf("%s.addressSpecies", sdfPlatform(sdf)))
    auc <- inferSpecies(sdf, return.auc = TRUE)
    species <- inferSpecies(sdf, return.species = TRUE)$scientificName
    detection <- as.integer(probeSuccessRate(sdf) * 100)
    df <- data.frame(auc=auc, species=names(auc))
    ggplot() +
        geom_segment(data=df[df$species!=addr$reference,],
            aes(x=auc,y=0,xend=auc,yend=1)) +
        geom_segment(data=df[df$species==addr$reference,],
            aes(x=auc,y=0,xend=auc,yend=1), size=2, color="red") +
        geom_segment(data=df[df$species==species,],
            aes(x=auc,y=0,xend=auc,yend=1), size=2, color="darkgreen") +
        xlim(c(0,1)) + ylab("") +
        xlab(sprintf("AUC (max: %1.2f; species: %s; detection: %1.0f%%)",
            max(auc), species, detection)) +
        theme_bw() +
        ggtitle(title) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank())
}
