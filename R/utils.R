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


# We use this function from the RAFSIL package.
#' ELBOW DETECTION
#'
#' @param scores vector
elbow_detection <-function(scores) {
  # We included this function from uSORT package, version 1.6.0 .
  #=====================
  num_scores <- length(scores)
  if (num_scores < 2) {
    stop("Input scores must be a vector with length more than 1!")
  }
  scores <- data.frame(id = seq_len(num_scores), value = scores)
  sorted_scores <- scores[order(scores$value, decreasing = TRUE),
  ]
  xy_coordinates <- cbind(x = seq_len(num_scores), y = sorted_scores$value)
  start_point <- xy_coordinates[1, ]
  end_point <- xy_coordinates[num_scores, ]
  x1 <- start_point[1]
  x2 <- end_point[1]
  y1 <- start_point[2]
  y2 <- end_point[2]
  a <- y1 - y2
  b <- x2 - x1
  c <- x1 * y2 - x2 * y1
  dist_to_line <- abs(a * xy_coordinates[, "x"] + b * xy_coordinates[,
                                                                     "y"] + c)/sqrt(a^2 + b^2)
  best_point_id <- which.max(dist_to_line)
  score_cutoff <- xy_coordinates[best_point_id, "y"]
  select_ID <- scores$id[which(scores$value >= score_cutoff)]
  return(select_ID)
}


##RSV

#' RANDOM PROJECTION SVD
#'
#' @param A matrix
#' @param K rank
#' @return list with components U S and V
#' @importFrom pracma orth
#' @importFrom stats rnorm
fast.rsvd <- function( A, K ) {
  #============================

  M = dim(A)[1]
  N = dim(A)[2]
  P = min(2*K,N)
  X = matrix(rnorm(N*P),nrow=N,ncol=P)
  Y = A%*%X
  W1 = pracma::orth(Y)
  B = t(W1)%*%A
  res = svd(B,nu=min(dim(B)),nv=min(dim(B)))
  W2 = res$u
  tmp_S = res$d
  S = array(0,c(length(tmp_S),length(tmp_S)))
  diag(S) = tmp_S
  V = res$v
  U = W1%*%W2
  K = min(K,dim(U)[2])
  U = U[,1:K]
  S = S[1:K,1:K]
  V = V[,1:K]

  return(list(U=U,S=S,V=V))
}



#' PCA built on fast.rsvd
#'
#' @param X matrix
#' @param K number of principal components
#' @return projection of X on principal top principal components
#' @importFrom matrixStats colVars
fast.pca <- function(X, K = 50){
  #=====================

  #X = t(X)
  tmp_val = as.vector(colSums(X)/nrow(X))
  # X = X - t(apply(array(0,dim(X)),MARGIN=1,FUN=function(x) {x=tmp_val}))
  X <- sweep(X, 1, tmp_val, FUN = "-")
  res = fast.rsvd(X,K)
  U = res$U
  S = res$S
  K = min(dim(S)[2],K)
  diag_val = sqrt(diag(S[1:K,1:K]))
  diag_mat = array(0,c(length(diag_val),length(diag_val)))
  diag(diag_mat) = diag_val
  X = U[,1:K]%*%diag_mat
  normalization_val = sqrt(rowSums(X*X))
  X = X / apply(array(0,c(length(normalization_val),K)),MARGIN=2,FUN=function(x) {x=normalization_val})
  pcgeneres<-X
  varlist<-colVars(pcgeneres)
  ordered_varlist<-order(varlist,decreasing = TRUE)
  LM1<-pcgeneres[,ordered_varlist]
  varf<-colVars(LM1)

  num<-length(elbow_detection(varf))
  pcgene<-LM1[,c(1:num)]

  return(pcgene)
}
