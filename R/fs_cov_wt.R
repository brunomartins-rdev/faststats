#' Weighted or unweighted covariance and optional correlation
#'
#' Computes the covariance matrix of a numeric matrix `x`, optionally weighted
#' by a vector `wt`. Can also return the correlation matrix. Supports unbiased
#' or maximum likelihood estimation.
#'
#' @param x A numeric matrix or data frame. Each row is an observation, each
#' column a variable.
#' @param wt Optional numeric vector of weights for each row. Must be
#' non-negative and sum > 0.
#' @param cor Logical. If TRUE, also returns the correlation matrix.
#' @param center Logical or numeric vector specifying the column means to use
#' for centering. If TRUE (default), the weighted or unweighted mean is used.
#' If a numeric vector is provided, that is used for centering.
#' @param method Character string, either "unbiased" (default) or "ML".
#' "unbiased" divides by 1 - sum(wt^2), "ML" divides by sum(wt).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{cov}{The covariance matrix.}
#'   \item{cor}{The correlation matrix, if `cor = TRUE`.}
#'   \item{center}{The column means used for centering.}
#'   \item{n.obs}{Number of observations (rows in `x`).}
#'   \item{wt}{The normalized weights used, if `wt` was provided.}
#' }
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60), nrow = 20, ncol = 3)
#' wt <- runif(nrow(x))
#' res1 <- fs_cov_wt(x) # unweighted
#' res2 <- fs_cov_wt(x, wt, TRUE, method = "ML") # weighted with correlation
#'
#' @export
fs_cov_wt <- function(
    x,
    wt = NULL,
    cor = FALSE,
    center = TRUE,
    method = "unbiased") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x))) stop("'x' must contain finite values only")

  cov_wt_cpp(
    x,
    wt_ = wt,
    cor = cor,
    center_ = if (isTRUE(center)) NULL else center,
    method = method
  )
}
