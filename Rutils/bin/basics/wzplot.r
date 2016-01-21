#!/usr/bin/env r

suppressMessages(library(docopt))

"Usage:
  wzplot.r barplot [-n NAMECOL] [-c PLOTCOL] [-t <INPUT>] -o <OUTPUT> [--main MAINTEXT] [--xlab XLABTEXT] [--ylab YLABTEXT] [--figwid FIGWIDTH] [--fighei FIGHEIGHT]

Options:
  -t <INPUT>          input table [default: stdin]
  -c PLOTCOL          plot column
  -n NAMECOL          name column
  -o <OUTPUT>         output figure [default: a]
  --figwid FIGWIDTH   figure width [default: 10]
  --fighei FIGHEIGHT  figure height [default: 8]
  --main MAINTEXT     main text
  --xlab XLABTEXT     xlab text
  --ylab YLABTEXT     ylab text
"-> doc

opt <- docopt(doc)

#suppressMessages(library(methylKit))
#cat(str(opt))

if (opt$barplot) {
    in.table <- read.table(opt$t, header=FALSE)

    pdf(opt$o, width=as.numeric(opt$figwid), height=as.numeric(opt$fighei))
    par(cex=1.5)
    barplot(in.table[,as.numeric(opt$c)], names.arg=in.table[,as.numeric(opt$n)])
    if (!is.null(opt$main)) {
        title(main=opt$main)
    }
    if (!is.null(opt$xlab)) {
        title(xlab=opt$xlab)
    }
    if (!is.null(opt$ylab)) {
        title(ylab=opt$ylab)
    }
    dev.off()
}
