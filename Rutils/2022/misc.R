## redefine the factor level of factor_name using the order of val_name (a column) in df
reorder_df_factor <- function(df, factor_name, val_name, group_df = df, decreasing = TRUE) {
    group_val = sapply(split(group_df[[val_name]], group_df[[factor_name]]), median, na.rm=TRUE)
    df[[factor_name]] = factor(df[[factor_name]], levels=names(sort(group_val, decreasing = decreasing)))
    df
}

view <- function(data, autofilter=TRUE) {
    ## source: https://jeromyanglim.tumblr.com/post/33825729070/function-to-view-r-data-frame-in-spreadsheet
    ## adapted to using writexl
    ## data: data frame
    ## autofilter: whether to apply a filter to make sorting and filtering easier
    data = as.data.frame(data)
    open_command <- switch(Sys.info()[['sysname']],
        Windows= 'open',
        Linux  = 'xdg-open',
        Darwin = 'open')

    require(writexl)
    temp_file <- paste0(tempfile(), '.xlsx')
    wb <- writexl::write_xlsx(data, temp_file)
    
    ## require(XLConnect)
    ## wb <- loadWorkbook(temp_file, create = TRUE)
    ## createSheet(wb, name = "temp")
    ## writeWorksheet(wb, data, sheet = "temp", startRow = 1, startCol = 1)
    ## if (autofilter) setAutoFilter(wb, 'temp', aref('A1', dim(data)))
    ## saveWorkbook(wb, )
    system(paste(open_command, temp_file))
}

xtab_set <- function(A,B,AnotB = FALSE, BnotA = FALSE, AandB = FALSE){
    both    <-  union(A,B)
    inA     <-  both %in% A
    inB     <-  both %in% B

    if (AnotB) {
        both[!(both %in% B)]
    } else if (BnotA) {
        both[!(both %in% A)]
    } else if (AandB) {
        both[(both %in% A) & (both %in% B)]
    } else {
        table(inA,inB)
    }
}

splitString <- function (x, s, i) {
    l <- strsplit(x, s)
    vapply(l, function(x) x[[i]], character(1))
}

calculate_density <- function(data, eps) {
  densities = numeric(nrow(data))
  for (i in 1:nrow(data)) {
    point = data[i, ]
    distances = sqrt((data$x - point$x)^2 + (data$y - point$y)^2)
    densities[i] = sum(distances < eps) - 1  # Subtract 1 to exclude the point itself
  }
  return(densities)
}

include_neighbors <- function(df, cluster, target_cluster, radius) {
    target_points = df[(!is.na(cluster)) & cluster == target_cluster, c("x", "y")]
    within_radius = sapply(seq_len(nrow(df)), function(i) {
        distances = sqrt((target_points$x - df[i,"x"])^2 + (target_points$y - df[i,"y"])^2)
        any(distances <= radius)
    })
    cluster[within_radius & is.na(cluster)] = target_cluster
    cluster
}

subcluster <- function(df0, target_cluster, eps=6, minPts=5) {
    df1 = df0[target_cluster,]
    subc = dbscan(df1[,c("x","y")], eps=eps, minPts=minPts)$cluster
    ## increase eps for bigger clusters
    res = ifelse(target_cluster, "1", "0")
    res[target_cluster] = paste0("sub.", subc)
    res
}

normalizeNumVec <- function(x) {
    x <- x - median(x, na.rm=T)
    x <- x / sd(x)
    x
}

chunkvec_bychunknumber <- function(x, chunk_number) {
      unname(split(x, cut(seq_along(x), chunk_number, labels=F)))
}

chunkvec_bychunksize <- function(x, chunk_size) {
      num_chunks <- ceiling(length(x) / chunk_size)
      groups <- rep(1:num_chunks, each = chunk_size, length.out = length(x))
      unname(split(x, groups))
}

oversample <- function(df, col, nmax) {
    do.call(bind_rows, lapply(unique(df[[col]]), function(level) {
        subset_df = df %>% dplyr::filter(.data[[col]] == level)
        rows_to_add = nmax - nrow(subset_df)
        if (rows_to_add > 0) {
            subset_df = bind_rows(subset_df,
                subset_df %>% slice_sample(n = rows_to_add, replace = TRUE))
        }
        subset_df
    }))
}

