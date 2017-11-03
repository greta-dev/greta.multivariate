#' @name conditional_bernoulli
#' @title conditional bernoulli distribution
#'
#' @description greta probability distribution over a K-dimensional vector of
#'   binary variables, arising from independent Bernoulli draws each conditioned
#'   on a single draw from another Bernoulli draw.
#'
#' @details A compound distribution, where elements of the bernoulli vector
#'   variable \emph{y} can only be 1 if a scalar latent bernoulli variable
#'   \emph{z} takes value 1, i.e.:
#'
#'   \deqn{ y_i ~ bernoulli(z * p_i)}
#'   \deqn{z ~ bernoulli(psi)}
#'   where
#'   \deqn{p_i = p(y_i = 1 | z = 1)}
#'   \deqn{psi = p(z = 1)}
#'
#'   \emph{p} and \emph{psi} are distinguishable provided there are multiple
#'   trials in each observation of \emph{y}. The density of this compound
#'   distribution can be calculated directly, explicitly integrating over the
#'   latent variable \emph{z}, as:
#'
#'   \deqn{psi * prod((p ^ y) * (1 - p) ^ (1 - y)) + max(y) * (1 - psi)}
#'
#'   This formulation underpins the ecological imperfect-detection model of
#'   MacKenzie et al. where \emph{y} and \emph{p} are vectors indicating whether
#'   a species was detected at each visit, and the probability of detection
#'   (which may vary between visits), and \emph{z} and \emph{psi} are scalars
#'   indicating whether the species was present (assumed to be the same at all
#'   visits) and the probability of being present.
#'
#' @references MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S.,
#'   Andrew Royle, J., & Langtimm, C. A. (2002). Estimating site occupancy rates
#'   when detection probabilities are less than one. \emph{Ecology}, 83(8),
#'   2248-2255.
#'
#' @param p matrix (of dimension \code{dim} x K) of (conditional) probabilities
#'   of success
#' @param psi scalar or column vector (of length \code{dim}) of probabilities
#'   for the latent bernoulli variable
#' @param dim a scalar giving the number of rows in the resulting greta array
#'
#' @importFrom R6 R6Class
#' @export
conditional_bernoulli <- function (p, psi, dim = 1)
  distrib('conditional_bernoulli', p, psi, dim)

# multivariate probit distribution
conditional_bernoulli_distribution <- R6Class (
  'conditional_bernoulli_distribution',
  inherit = .internals$nodes$node_classes$distribution_node,
  public = list(

    initialize = function (p, psi, dim) {

      # coerce to greta arrays
      p <- as.greta_array(p)
      psi <- as.greta_array(psi)

      # check dimensions of p
      if (ncol(p) < 2 | length(dim(p)) != 2) {

        stop ("p must be a 2D greta array with at least two columns, ",
              "but has dimensions ",
              paste(dim(p), collapse = " x "),
              call. = FALSE)

      }

      # check dimensions of psi
      if (ncol(psi) != 1 | length(dim(psi)) != 2) {

        stop ("psi must be a 2D greta array with one column, ",
              "but has dimensions ",
              paste(dim(psi), collapse = " x "),
              call. = FALSE)

      }

      # compare possible dimensions
      dim_p <- nrow(p)
      dim_psi <- nrow(psi)

      if (dim_p != dim_psi) {

        stop ("p and psi have different dimensions, ",
              dim_p, " vs ", dim_psi,
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
      super$initialize('conditional_bernoulli',
                       dim = c(dim, ncol(p)),
                       discrete = TRUE)
      self$add_parameter(p, "p")
      self$add_parameter(psi, "psi")

    },

    tf_distrib = function (parameters) {

      p <- parameters$p
      psi <- parameters$psi

      # return a tf function, taking the binary vector and returning the density

      log_prob <- function (x) {

        # for each row, were all elements 0?
        none <- tf_as_float(tf_rowsums(x, 1L) == 0)

        one <- fl(1)

        # log conditional probability
        # cp <- (p ^ x) * (one - p) ^ (one - x)
        # log_cp <- log(cp)
        log_cp <- x * log(p) + (one - x) * log(one - p)

        # log probability
        # log(rowProds(cp) * psi + nd * (1 - psi))
        prob <- exp(tf_rowsums(log_cp, 1L) + log(psi)) + none * (one - psi)
        log(prob)

      }

      list(log_prob = log_prob,
           cdf = NULL,
           log_cdf = NULL)

    },

    # no CDF for multivariate distributions
    tf_cdf_function = NULL,
    tf_log_cdf_function = NULL

  )
)
