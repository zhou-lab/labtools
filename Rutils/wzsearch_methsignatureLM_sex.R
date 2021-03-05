#!/usr/bin/env R
library(stringr)
OpenRead <- function(arg) {

  if (arg %in% c("-", "/dev/stdin")) {
    file("stdin", open="r")
  } else if (grepl("^/dev/fd/", arg)) {
    fifo(arg, open="r", blocking=TRUE)
  } else {
    file(arg, open = "r")
  }
}

args <- commandArgs(trailingOnly = TRUE)
## print(readLines(OpenRead(a[1])))
## read in the group, only the 2nd column is used
a <- read.table(OpenRead(args[1]))
groups <- setNames(a$V2, a$V1)
group_alphabet <- sort(unique(groups))

deque <- lapply(1:5, function(i) NULL)
deque_ind <- 1
deque_loc <- lapply(1:5, function(i) NULL)

cat('chrm','beg','end','range','mostU','gapU','baseline',group_alphabet,'\n', sep="\t")
con <- OpenRead(args[2])
while (TRUE) {
  line = readLines(con, n=1)
  if (length(line) == 0) { break; }
  fields <- scan(text=line, what = character(1), sep="\t", quiet=TRUE)
  betas <- as.numeric(fields[4:length(fields)])

  deque[[deque_ind]] = betas
  deque_loc[[deque_ind]] = fields[2]
  deque_ind <- deque_ind + 1
  if(deque_ind > 5) {deque_ind = 1;}
  b1 = do.call(c, lapply(deque, function(bt) if(is.null(bt)) {NA;} else {bt}))
  g1 = do.call(c, lapply(deque, function(bt) if(is.null(bt)) {NA;} else {groups}))
#cat(str(deque),'\n')
#cat(b1,'\n')
#cat(g1,'\n')
#cat("tgt", length(b1), length(g1),'\n')
  b2 = b1[!is.na(b1)]
  g2 = g1[!is.na(b1)]
#cat("tgt2", length(b2), length(g2),'\n')
#cat(table(g2),'\n')
#cat(names(table(g2)),'\n')
  tg <- g2 %in% names(which(table(g2) >=10)) # min data points
  b2 = b2[tg]
  g2 = g2[tg]

  if(length(unique(g2)) < 2) { next; }
  g2levels <- sort(unique(g2))
  g2 <- factor(g2, levels=g2levels)
  ft = lm(b2 ~ g2)$coef
  ft = ft[2:length(ft)]
  baseline = ft[1] + ft['g2Male']
  names(ft) <- str_replace(names(ft),'g2','')
  delta <- -sum(ft)/(length(ft)+1); ft <- ft+delta; ft[g2levels[1]] <- delta; # normalize

  sft <- sort(ft[group_alphabet])
  rg <- max(sft) - min(sft)

  cat(fields[1],deque_loc[[deque_ind]],fields[3],rg,names(sft)[1],sft[2]-sft[1],baseline,ft[group_alphabet],'\n', sep="\t")
}



