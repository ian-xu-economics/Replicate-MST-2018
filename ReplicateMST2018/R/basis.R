#' Create an MTRBasis list
#'
#' This function creates an MTRBasis list containing three components: a, b, and ib.
#'
#' @param a A list of functions for the 'a' component.
#' @param b A list of functions for the 'b' component.
#' @param ib A list of functions for the 'ib' component.
#' @return A list containing the components a, b, and ib.
#' @export
MTRBasis = function(a, b, ib) {
  stopifnot(length(b) == length(ib))
  list(a = a, b = b, ib = ib)
}

#' Create a Bernstein basis
#'
#' This function generates a Bernstein basis of degree K.
#'
#' @param K The degree of the Bernstein basis.
#' @return An MTRBasis list for the Bernstein basis.
#' @export
bernstein_basis = function(K) {
  a = list(function(z) 1)
  b = lapply(0:K, function(k) function(u) bernstein_polynomial(u, k, K))
  ib = lapply(0:K, function(k) function(lower_pt, upper_pt) integrate_bernstein_polynomial(lower_pt, upper_pt, k, K))

  MTRBasis(a, b, ib)
}

#' kth Bernstein basis polynomial of degree K
#'
#' This function calculates the kth Bernstein basis polynomial of degree K at a given point u.
#'
#' @param u The point at which to evaluate the polynomial.
#' @param k The index of the Bernstein basis polynomial.
#' @param K The degree of the Bernstein basis.
#' @return The value of the kth Bernstein basis polynomial of degree K at point u.
#' @export
bernstein_polynomial = function(u, k, K) {
  choose(K, k) * u^k * (1 - u)^(K - k)
}

#' Integrate kth Bernstein basis polynomial of degree K
#'
#' This function calculates the integral of the kth Bernstein basis polynomial of degree K between lower_pt and upper_pt.
#'
#' @param lower_pt The lower limit of integration.
#' @param upper_pt The upper limit of integration.
#' @param k The index of the Bernstein basis polynomial.
#' @param K The degree of the Bernstein basis.
#' @return The integral value.
integrate_bernstein_polynomial = function(lower_pt, upper_pt, k, K) {
  sum(sapply(k:K, function(i) {
    (-1)^(i - k) * choose(K, i) * choose(i, k) * (upper_pt^(i + 1) - lower_pt^(i + 1)) / (i + 1)
  }))
}

#' Weighted Kth degree Bernstein basis polynomial
#'
#' This function calculates a weighted sum of Bernstein basis polynomials of degree K at a given point u.
#'
#' @param u The point(s) at which to evaluate the polynomial.
#' @param K The degree of the Bernstein basis.
#' @param weights A vector of weights for each basis polynomial.
#' @return The weighted sum of the Bernstein basis polynomials.
#' @export
weighted_bernstein_basis = function(u, K, weights) {

  basis_values = bernstein_basis(K)$b
  stopifnot(length(weights) == length(basis_values))

  results = NULL
  for(j in u){
    results = results %>%
      c(sum(sapply(1:length(basis_values), function(i) weights[i] * basis_values[[i]](j))))
  }

  return(results)

}

#' Create a constant spline basis
#'
#' This function generates a constant spline basis based on the given knots.
#'
#' @param knots A vector of knot points, must include 0 and 1.
#' @return An MTRBasis list for the constant spline basis.
#' @export
constantspline_basis = function(knots) {
  stopifnot(all(knots >= 0) && all(knots <= 1))
  stopifnot(any(knots == 0) && any(knots == 1))

  knots = unique(sort(knots))

  a = list(function(z) 1)

  in_partition = function(u, k) {
    if (k != tail(knots, n = 1)) {
      return(as.integer(knots[k - 1] <= u & u < knots[k]))
    } else {
      return(as.integer(knots[k - 1] <= u & u <= knots[k]))
    }
  }

  b = lapply(2:length(knots), function(k) function(u) in_partition(u, k))
  ib = lapply(2:length(knots), function(k) function(u, v) max(0, (min(v, knots[k]) - max(u, knots[k - 1]))))

  MTRBasis(a, b, ib)
}
