#' Find the elbow point in a vector of variances or eigenvalues.
#'
#' @param x A numeric vector of variances or eigenvalues.
#'
#' @returns The index of the elbow point in the vector.
#' @keywords internal
#'
find_elbow <- function(x) {
  # x = vector of variances or eigenvalues
  n <- length(x)
  all_points <- cbind(1:n, x)

  # line between first and last point
  line_vec <- all_points[n,] - all_points[1,]
  line_vec <- line_vec / sqrt(sum(line_vec^2))

  # distances
  distances <- sapply(1:n, function(i) {
    point_vec <- all_points[i,] - all_points[1,]
    crossprod(rbind(point_vec[1], point_vec[2]),
              rbind(-line_vec[2], line_vec[1]))
  })

  which.max(abs(distances))
}
