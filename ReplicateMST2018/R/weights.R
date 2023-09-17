#' Compute Average Weights for Target Parameters
#'
#' This function computes the average weights for different target parameters.
#'
#' @param tp A list containing the name of the target parameter and its integration limits.
#' @param dgp An optional DGP object containing the data generating process.
#' @return A data frame containing the average weights.
#' @export
compute_average_weights = function(tp, dgp = NULL) {

  if (tp$name == "LATE(u1, u2)") {
    u1 = tp$int_limits[1]
    u2 = tp$int_limits[2]
    results = data.frame(u = c(0, u1, u2, 1),
                          `average weight for d = 1` = rep(0, 4),
                          `average weight for d = 0` = rep(0, 4))
    results[2, 2] = 1 / (u2 - u1)
    results[2, 3] = -1 / (u2 - u1)
  } else {
    stop("Weight computation is unsupported for this target parameter.")
  }

  return(results)
}

#' Compute Average Weights for IV-like Parameters
#'
#' This function computes the average weights for IV-like parameters.
#'
#' @param ivlike A list containing the name of the IV-like parameter and its function 's'.
#' @param dgp A DGP object containing the data generating process.
#' @return A data frame containing the average weights.
#' @export
compute_average_weights_ivlike = function(ivlike, dgp) {

  results = data.frame(u = unique(c(0, dgp$pscore)))
  order = order(dgp$pscore)

  s = ivlike$s

  if(ivlike$name == "Saturated"){
    terms = function(d) {
      s_values = sapply(1:length(dgp$suppZ), function(i) s(d, dgp$suppZ[i]))
      return(s_values * dgp$densZ)
    }
  } else{
    terms = function(d){
      s(d, dgp$suppZ) * dgp$densZ
    }
  }

  # d = 1 weights
  d1terms = terms(1)
  summands = d1terms[rev(order)]  # order by decreasing pscore
  summands = c(0, summands)
  results$`average weight for d = 1` = rev(cumsum(summands))

  # d = 0 weights
  d0terms = terms(0)
  summands = d0terms[order]  # order by increasing pscore
  summands = c(0, summands)
  results$`average weight for d = 0` = cumsum(summands)

  return(results)
}
