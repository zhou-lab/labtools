get_points_close_to_hull <- function(cluster_data, radius_threshold=0.1) {
    library(geometry)
    library(dplyr)

    cluster_data = cluster_data[,c("x","y")]

    ## Compute convex hull indices
    hull_indices <- chull(cluster_data$x, cluster_data$y)
    hull_points <- cluster_data[hull_indices, ]
    
    ## Initialize a vector to store minimum distances of points to the hull
    min_distances <- numeric(nrow(cluster_data))
    
    ## Calculate distance from each point to each segment of the hull
    for (i in 1:nrow(cluster_data)) {
        point <- cluster_data[i, ]
        distances <- sapply(1:nrow(hull_points), function(j) {
            next_index <- ifelse(j == nrow(hull_points), 1, j + 1)
            dist_to_segment(point, hull_points[j, ], hull_points[next_index, ])
        })
        min_distances[i] <- min(distances)
    }
    
    ## Select points where the minimum distance is less than the threshold
    min_distances <= radius_threshold
}

## Calculate distance from a point to a line segment
dist_to_segment <- function(point, line_start, line_end) {
    ## Line segment vector
    line_vec <- c(line_end$x - line_start$x, line_end$y - line_start$y)
    
    ## Point vector relative to line_start
    point_vec <- c(point$x - line_start$x, point$y - line_start$y)
    
    ## Project point_vec onto line_vec
    proj <- sum(point_vec * line_vec) / sum(line_vec * line_vec)
    proj_point <- line_start + proj * line_vec
    
    ## Check if projection is outside the line segment
    if (proj < 0) {
        closest_point <- line_start
    } else if (proj > 1) {
        closest_point <- line_end
    } else {
        closest_point <- proj_point
    }
    
    ## Return distance from the point to the closest point on the segment
    sqrt(sum((point - closest_point)^2))
}
