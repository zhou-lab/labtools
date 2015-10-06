hist_logy <- function(x) {
  hd <- hist(x,plot=F);
  hd$counts <- log10(hd$counts);
  plot(hd, ylab='log10(Frequency)', axes=F);
  axis(1);
  cnts <- seq(floor(min(hd$counts)), floor(max(hd$counts)), length.out=4);
  axis(2,at=cnts,labels=10**cnts);
}