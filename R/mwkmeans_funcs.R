#' Calculates the value of the objective function for optimizing cluster assignments.
#'
#' @param j
#' @param jPrime
#' @param x
#' @return
ObjectiveValue <- function(j, jPrime, x) {
  # calculate u
  dif <- jPrime - j
  u <- drop((t(dif) %*% (jPrime - x)) / (t(dif) %*% dif))
  # account for domain
  if (u < 0) {u = 0}
  if (u > 1) {u = 1}
  # find the objective
  A <- x - u*j - (1-u)*jPrime
  obj <- t(A) %*% A
  # return
  return(c(u, obj))
}

#' Calculates the Phi matrix that minimizes the objective for each data.
#' 
#' @param centers Current cluster cents
#' @param data Data matrix (features X samples)
MinAssignments <- function(centers, data) {
  # initialize phi as kxn matrix of all zeros
  k <- dim(centers)[2]
  n <- dim(data)[2]
  phi <- matrix(0L, k, n)
  
  # find the optimal assignment for each data point
  for (i in (1:n)) {
    # get the objective value for each pair of clusters
    objVals <- matrix(ncol = 4) 
    for (j in (1:(k - 1))) {
      for (jPrime in ((j + 1):k)) {
        obj <- ObjectiveValue(centers[,j], centers[,jPrime], data[,i])
        objVals <- rbind(objVals, c(j, jPrime, obj))
      }
    }
    # find the lowest score
    min <- which.min(objVals[,4])
    phi[objVals[min,1], i] <- objVals[min, 3]
    phi[objVals[min,2], i] <- 1 - objVals[min, 3]
  }
  return(phi)
}

#' Calculates the optimal centers given data and cluster assignments.
#' 
#' @param phi Matrix of phi values for each data.
#' @param data Data matrix (features X samples).
#' @return Optimal center matrix
MinClusters <- function(phi, data) {
  # identify new centers
  new.centers <- phi %*% t(phi)
  new.centers <- solve(new.centers)
  new.centers <- new.centers %*% phi
  new.centers <- new.centers %*% t(data)
  # return transpose
  return(t(new.centers))
}

#' Generates a multiway k-means given the data and starting center points.
#' 
#' @param data Matrix of data (features X samples)
#' @param centers Starting centers in matrix format (features X centers)
#' @return A list; 'phi' with the cluster assignments for each data and 'centers' with the calculated cluster centers.
#' @export
MWKMeans <- function(data, centers) {
  # find dimension variables
  d <- dim(data)[2]
  n <- dim(data)[1]
  k <- dim(centers)[2]
  # tracker variables
  phi <- MinAssignments(centers, data)
  newCenters <- MinClusters(phi, data)
  while (!all( abs(newCenters-centers) < 1e-3)) {
    centers <- newCenters
    phi <- MinAssignments(centers, data)
    newCenters <- MinClusters(phi, data)
  }
  colnames(phi) <- colnames(data)
  # return the cluster centers and the weight matrix
  return(list(phi, newCenters))
}
