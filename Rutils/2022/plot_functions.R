
wzGetColors <- function(grouping, palette_func = brewer.paired, group_name=TRUE) {
    ## some good color palettes:
    ## brewer.paired, brewer.set1, brewer.set2, brewer.dark2
    ## polychrome (n<=36), galsbey (n<=32),
    ## alphabet, alphabet2, brewer.set3, parula, turbo
    library(pals)
    groups <- unique(grouping)
    ## if (palette.name == 'Set1') nclr = 9
    ## else if (palette.name == 'Dark2') nclr = 8
    ## else nclr = 12
    colors = palette_func(length(groups))
    ## colors = colorRampPalette(brewer.pal(nclr, palette.name))(length(groups))
    if (group_name == TRUE) {
        if(length(colors) != length(groups)) {
            stop("the palette cannot provide so many colors")
        }
        setNames(colors, groups)
    } else {
        colors
    }
}

wzPlotDens1d.SE <- function(se, color_var = NULL) {
    ad = assays(se)[[1]]
    cd = as_tibble(colData(se))
    add = melt(ad, value.name="betas", varnames=c("Probe_ID", "Sample_ID"))
    if (is.null(color_var)) {
        ggplot(add) + geom_density(aes(betas, color = Sample_ID))
    } else {
        add[[color_var]] = cd[[color_var]][match(add$Sample_ID, cd[["Sample_ID"]])]
        ggplot(add) + geom_density(aes_string("betas", group = "Sample_ID", color=color_var))
    }
}

wzPlotDens2d <- function(x1,x2,pch=16,colorstops=c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"),...) {
    df <- data.frame(x1,x2)
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
    df$dens <- col2rgb(x)[1,] + 1L
    cols <-  colorRampPalette(colorstops)(256)
    df$col <- cols[df$dens]
    plot(x2~x1, data=df[order(df$dens),], pch=pch, col=col, ...)
}

wzPlotDens2d.smooth <- function(x, y, xlim=c(-2,2), stop.points=c("white","lightblue","blue","green","yellow","orange","red","darkred"), nrpoints=100, nbins=256, ...) {
    palette <- colorRampPalette(stop.points, space = "Lab")
    test = cor.test(x, y, method="spearman", use="na.or.complete")
    smoothScatter(x, y, xlim = xlim,
        nrpoints=nrpoints,
        nbin=c(nbins,nbins),
        colramp=palette, col='blue',
        main=sprintf("rho=%1.3f, p=%1.2f", test$estimate, test$p.value), ...)
}

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

#' @param label counts or percent
wzvenn <- function(data, proportional = TRUE, label = "counts") {
    ## data is a list of elements
    ## there are several packages for drawing venn / Euler diagrams
    ## - VennDiagram (maintained by Hanbo Chen and Paul Boutros)
    ## - venneuler
    ## - eVenn
    ## - venn
    ## - colorfulVennPlot
    ## - eulerr, based on venneuler, the most advanced, see https://cran.r-project.org/web/packages/eulerr/vignettes/introduction.html
    if (proportional) {
        library(eulerr) # recommended
        fit = euler(data)
        ## https://cran.r-project.org/web/packages/eulerr/eulerr.pdf
        plot(fit, labels = list(font = 4), quantities = list(type=label))
    } else {
        library(VennDiagram)
        g <- venn.diagram(data, filename = NULL)
        grid.newpage()
        pushViewport(viewport(unit(0.1,'npc'),unit(0.1,'npc'),unit(0.8,'npc'),unit(0.8,'npc'), just=c('left','bottom')))
        grid.draw(g)
    }
}

wzupset <- function(data, nsets=5, mb.ratio = c(0.7, 0.3)) {
    library(UpSetR)
    upset(fromList(data), order.by="freq", nsets=25, mb.ratio = mb.ratio)
}

ggVerticalX <- function() {
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

gg45X <- function() {
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

ggPercentage <- function() {
    scale_y_continuous(labels=scales::percent)
}

ggNoLegend <- function() {
    guides(fill="none")
    ## or theme(legend.position="none")
    ## or scale_fill_discrete(guide="none")
}

visualizeSEsimple <- function(se, rows=NULL, cols=NULL,
    name_base="a", legend_hpad=0.1, stop.points=turbo(10), #or parula(10)
    dmin=0, dmax=1, label.use.data = FALSE, show_legend=TRUE,
    show_row_names = FALSE, show_column_names = FALSE) {

    plt = WHeatmap(assay(se), cmp=CMPar(stop.points=stop.points, dmin=dmin, dmax=dmax),
        name=paste0(name_base, 'matrix'),
        xticklabels = show_column_names, xticklabels.n = ncol(se),
        yticklabels = show_row_names, yticklabels.n = nrow(se))
    last_matrix = paste0(name_base, 'matrix')
    
    ## row bars
    last_name = last_matrix
    for (i in seq_along(rows)) {
        bar = rows[i]
        if (paste0(bar,".colors") %in% names(metadata(se))) {
            bar.colors = metadata(se)[[paste0(bar,".colors")]]
        } else { bar.colors = NULL }

        if (label.use.data) { lud = (i==length(rows))
        } else { lud = FALSE }

        plt = plt + WColorBarH(colData(se)[,bar],
            TopOf(last_name, height=0.05),
            cmp = CMPar(label2color = bar.colors),
            name=paste0(name_base, "row", bar), label=bar,
            label.use.data=lud, label.pad=0.3,
            xticklabel.side='t')
        last_name = paste0(name_base, "row", bar)
    }

    ## column bars
    last_name = last_matrix
    for (i in seq_along(cols)) {
        bar = cols[i]
        if (paste0(bar,".colors") %in% names(metadata(se))) {
            bar.colors = metadata(se)[[paste0(bar,".colors")]]
        } else { bar.colors = NULL }
        if (label.use.data) { lud = (i==length(cols))
        } else { lud = FALSE }
        plt = plt + WColorBarV(rowData(se)[,bar],
            RightOf(last_name, width=0.05),
            cmp = CMPar(label2color = bar.colors),
            name=paste0(name_base, "col", bar), label=bar,
            label.side='left', label.pad=0.3,
            label.use.data=lud,
            yticklabel.side='l')
        last_name = paste0(name_base, "col", bar)
    }
    last_column_name = last_name

    ## legend
    if (show_legend) {
        for (bar in rows) {
            plt = plt + WLegendV(paste0(name_base, "row", bar),
                BottomLeftOf(last_column_name, just=c("left","top"), h.pad=legend_hpad),
                name=paste0(name_base, "legendrow", bar))
            last_name = paste0(name_base, "legendrow", bar)
        }
        for (bar in cols) {
            plt = plt + WLegendV(paste0(name_base, "col", bar),
                BottomLeftOf(last_column_name, just=c("left","top"), h.pad=legend_hpad),
                name=paste0(name_base, "legendcol", bar))
            last_name = paste0(name_base, "legendcol", bar)
        }
    }

    plt
}

## this is a fancier version of visualizeSE, it creates gaps between the column groups and row groups
## se = SummarizedExperiment(assays=list(betas=bt), colData=meta[match(colnames(bt), meta$Sample_ID),])
## metadata(se) = readExcelColors("~/samplesheets/2022/20220109_TCGA.MTAP.annoation.xlsx")
## visualizeSE(se, rows="rowbar1_name")
visualizeSE <- function(se, rows=NULL, cols=NULL,
    column_split=NA, column_split_nms=NULL, column_split_pad=0.05, column_cluster=FALSE,
    name_base="a", legend_hpad=0.3, stop.points=parula(20), show_row_names = FALSE, show_column_names = FALSE) {

    if (is.na(column_split)) {
        ses = list(main=se)
    } else {
        if (is.null(column_split_nms)) {
            column_split_nms = sort(unique(colData(se)[[column_split]]))
        }
        ses = lapply(column_split_nms, function(nm) {
            se[,colData(se)[[column_split]] == nm]
        })
        names(ses) = column_split_nms
    }

    if (column_cluster) {
        ses = lapply(ses, columnClusterSE)
    }

    for (i in seq_along(ses)) {
        column_nm = names(ses)[i]
        se = ses[[i]]
        if (i==1) {
            plt = WHeatmap(assay(se), cmp=CMPar(stop.points=stop.points, dmin=0, dmax=1), name=paste0(name_base, column_nm, 'matrix'), xticklabels = show_column_names, xticklabels.n = ncol(se), yticklabels = show_row_names, yticklabels.n = nrow(se))
        } else {
            plt = plt + WHeatmap(assay(se), cmp=CMPar(stop.points=stop.points, dmin=0, dmax=1),
                dm = RightOf(last_column, pad=column_split_pad), name=paste0(name_base, column_nm, 'matrix'))
        }
        last_column = paste0(name_base, column_nm, 'matrix')

        ## row bars
        last_name = paste0(name_base, column_nm, 'matrix')
        for (bar in rows) {
            if (paste0(bar,".colors") %in% names(metadata(se))) {
                bar.colors = metadata(se)[[paste0(bar,".colors")]]
            } else {
                bar.colors = NULL
            }
            ## only label the last column
            if (i == length(ses)) {
                label = bar
            } else {
                label = ""
            }
            plt = plt + WColorBarH(colData(se)[,bar], TopOf(last_name),
                cmp = CMPar(label2color = bar.colors),
                name=paste0(name_base, column_nm, bar), label=label)
            last_name = paste0(name_base, column_nm, bar)
        }

        plt = plt + WLabel(column_nm, TopOf(last_name))
    }

    last_name = last_column
    for (bar in rows) {
        for (i in seq_along(ses)) {
            plt = plt + WLegendV(paste0(name_base, names(ses)[i], bar),
                TopRightOf(last_name, just=c("left","top"), h.pad=legend_hpad),
                name=paste0(name_base, last_column, "legend", bar, i))
            last_name = paste0(name_base, last_column, "legend", bar, i)
        }
    }
    plt
}

## usage: p <- ggplot(USA.states, aes(x=state.region,y=Income))+geom_violin(trim = F)
## fancy_violin(p)
fancy_violin <- function(p, size=0.3) {
    mywidth <- .35 # bit of trial and error
    ## This is all you need for the fill:
    vl_fill <- data.frame(ggplot_build(p)$data) %>%
        mutate(xnew = x- mywidth*violinwidth, xend = x+ mywidth*violinwidth) 

    ## Bit convoluted for the outline, need to be rearranged: the order matters
    vl_poly <- vl_fill %>% 
        dplyr::select(xnew, xend, y, group) %>%
        pivot_longer(-c(y, group), names_to = "oldx", values_to = "x") %>% 
        arrange(y) %>%
        split(., .$oldx) %>%
        map(., function(x) {
            if(all(x$oldx == "xnew")) x <- arrange(x, desc(y))
            x
        }) %>%
        bind_rows()

    ggplot() +
        geom_polygon(data = vl_poly, aes(x, y, group = group), color= "black", size = size, fill = NA) +
        geom_segment(data = vl_fill, aes(x = xnew, xend = xend, y = y, yend = y, color = y))
}
