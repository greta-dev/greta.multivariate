#' @name multivariate_probit
#' @title multivariate probit distribution
#'
#' @description greta probability distribution over a K-dimensional vector of
#'   binary variables.
#'
#' @param mean matrix or column vector (of length K) of unconstrained parameters
#'   giving the mean of the latent normal parameters
#' @param C K-dimensional correlation matrix
#' @param dim a scalar giving the number of rows in the resulting greta array
#'
#' @export
multivariate_probit <- function (mean, C, dim = 1)
  distrib('multivariate_probit', mean, C, dim)

distrib <- greta::.internals$nodes$constructors$distrib
tf_iprobit <- greta::.internals$tensors$tf_iprobit

# multivariate probit distribution
multivariate_probit_distribution <- R6Class (
  'multivariate_probit_distribution',
  inherit = .internals$nodes$node_classes$distribution_node,
  public = list(

    initialize = function (mean, C, dim) {

      # coerce to greta arrays
      mu <- as.greta_array(mean)
      R <- as.greta_array(C)

      # check dimensions of mean
      if (ncol(mean) != 1 | length(dim(mean)) != 2) {

        stop ('mean must be a 2D greta array with one column, but has dimensions ',
              paste(dim(mean), collapse = ' x '),
              call. = FALSE)

      }

      # check dimensions of C
      if (nrow(C) != ncol(C) | length(dim(C)) != 2) {

        stop ('C must be a square 2D greta array, but has dimensions ',
              paste(dim(C), collapse = ' x '),
              call. = FALSE)

      }

      # compare possible dimensions
      dim_mean <- nrow(mean)
      dim_C <- nrow(C)

      if (dim_mean != dim_C) {

        stop ('mean and C have different dimensions, ',
              dim_mean, ' vs ', dim_C,
              call. = FALSE)

      }

      if (dim_mean == 1) {

        stop ('the multivariate probit distribution is for vectors, ',
              'but the parameters were scalar',
              call. = FALSE)

      }

      # check dim is a positive scalar integer
      dim_old <- dim
      dim <- as.integer(dim)
      if (length(dim) > 1 || dim <= 0 || !is.finite(dim)) {

        stop ('dim must be a scalar positive integer, but was: ',
              capture.output(dput(dim_old)),
              call. = FALSE)

      }

      # coerce the parameter arguments to nodes and add as children and
      # parameters
      super$initialize('multivariate_probit',
                       dim = c(dim, length(mean)),
                       discrete = TRUE)
      self$add_parameter(mean, 'mean')
      self$add_parameter(C, 'C')

      # create latent variable and add as a parameter
      u = variable(0, 1, dim = self$dim)
      self$add_parameter(u, 'u')
      # this needs to be overwritted on distribution assignment to get the
      # correct dimensions!

    },

    tf_distrib = function (parameters) {

      mean <- parameters$mean
      C <- parameters$C
      u <- parameters$u

      # change this to use the cholesky factor if available!
      L <- tf$cholesky(C)

      N <- self$dim[1]
      K <- self$dim[2]

      stop ("not implemented")

      # return a tf function, taking the binary vector and returning the density
      # (including the log jacobian transform?)

      log_pdf <- function (x) {

        # find the elements of x that are 1, and those that are 0, count them
        # and split u too. Do this on the values matrix, rather than the tensor?

        # Add this as a parameter and ignore x? This will make it harder to
        # switch to discrete inference though.

        N_pos
        N_neg
        idx_neg
        idx_pos

        u_pos <- u[idx_pos]
        u_neg <- u[idx_neg]

        bound <- tf$zeros(self$dim)
        lj_pos <- tf$zeros(c(N_pos, K))
        lj_neg <- tf$zeros(c(N_neg, K))
        prev <- tf$zeros(N)

        for (k in seq_len(K)) {

          bound[, k] <- tf_iprobit(-(mu[, k] + prev) / L)
          z_pos <- tf_iprobit(bound + (1 - bound) * u_pos)
          z_neg <- tf_iprobit(bound * u_neg)
          z[, k] <- c(z_pos, z_neg)
          prev <- L[k, 1:k] %*% z[1:k, ]
          lj_pos[, k] <- tf$log(1 - bound_pos[, k])
          lj_neg[, k] <- tf$log(bound_neg[, k])

        }

      }

      list(log_pdf = log_pdf,
           cdf = NULL,
           log_cdf = NULL)

    },

    # no CDF for multivariate distributions
    tf_cdf_function = NULL,
    tf_log_cdf_function = NULL

  )
)

