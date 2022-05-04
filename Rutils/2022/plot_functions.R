
wzPlotDens2d.smooth <- function(x, y, xlim=c(-2,2), nrpoints=100, nbins=256, ...) {
    palette <- colorRampPalette(
        c("white","lightblue","blue","green","yellow","orange","red","darkred"),
        space = "Lab")

    smoothScatter(x, y, xlim = xlim,
        nrpoints=nrpoints,
        nbin=c(nbins,nbins),
        colramp=palette, col='blue', ...)
}

wzvenn <- function(..., dnames=NULL) {
  data <- list(...)
  if (is.null(dnames))
    dnames <- 1:length(data)
  names(data) <- dnames
  library(VennDiagram)
  g <- venn.diagram(data, filename = NULL)
  grid.newpage()
  pushViewport(viewport(unit(0.1,'npc'),unit(0.1,'npc'),unit(0.8,'npc'),unit(0.8,'npc'), just=c('left','bottom')))
  grid.draw(g)
}

