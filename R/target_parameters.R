#' Create a Target Parameter Structure
#'
#' @param name Name of the target parameter.
#' @param int_limits Function for integration limits.
#' @param int_constant Function for integration constant.
#' @param legendtitle Optional legend title.
#' @return A list representing the target parameter.
#' @export
TargetParameter = function(name, int_limits, int_constant, legendtitle=NULL) {
  if (is.null(legendtitle)) {
    legendtitle = name
  }
  return(list(name=name, int_limits=int_limits, int_constant=int_constant, legendtitle=legendtitle))
}

#' Evaluate Target Parameter
#'
#' @param tp Target parameter object.
#' @param mtrs List of MTR objects.
#' @param dgp Data Generating Process object.
#' @return A numeric value representing the evaluated target parameter.
#' @export
eval_tp = function(tp, mtrs, dgp) {

  gamma_star = compute_gamma_star(tp, list(mtrs[[1]]$basis, mtrs[[2]]$basis), dgp)
  total = 0
  for (l in seq_along(mtrs)) {
    for (d in 0:1) {
      total = total + sum(gamma_star[[l]][[d + 1]] * mtrs[[d + 1]]$theta)
    }
  }
  return(total)
}

#' Compute Gamma Star for Each Model
#'
#' @param tp Target parameter object.
#' @param bases List of basis functions for each model.
#' @param dgp Data Generating Process object.
#' @return A list of Gamma_star matrices for each model.
#' @export
compute_gamma_star = function(tp, bases, dgp) {
  lapply(seq_along(bases), function(l) {
    lapply(0:1, function(d) {
      compute_single_gamma_star(tp, bases[[l]], d, l, dgp)
    })
  })
}

#' Compute Single Gamma Star
#'
#' @param tp Target parameter object.
#' @param basis Basis functions for a specific model.
#' @param d Treatment indicator (0 or 1).
#' @param l Model index.
#' @param dgp Data Generating Process object.
#' @return A Gamma_star matrix for the specific basis.
#' @export
compute_single_gamma_star = function(tp, basis, d, l, dgp) {
  gamma_star = matrix(0, nrow=length(basis$a), ncol=length(basis$b))
  for (z in 1:length(dgp$suppZ)) {
    il = tp$int_limits(dgp$suppZ[z])
    gamma_star = gamma_star + outer(
      sapply(basis$a, function(aj) aj(dgp$suppZ[z])),
      sapply(basis$ib, function(ibk) ibk(il[1], il[2])),
      function(a, b) a * b * tp$int_constant(l, d, dgp$suppZ[z])
    )
  }
  return(gamma_star)
}

#' Create LATE Target Parameter
#'
#' @param dgp Data Generating Process object.
#' @param u1 Lower limit for LATE.
#' @param u2 Upper limit for LATE.
#' @param l Model index (default is 1).
#' @return A TargetParameter object for LATE.
#' @export
late = function(dgp, u1, u2, l = 1){
  stopifnot(u1 <= u2)

  name = "LATE(u1, u2)"
  int_limits = function(z){
    return(c(u1,u2))
  }
  int_constant = function(l, d, z){
    return((l==1)*(2*d-1)*dgp$find_density(z)/(u2-u1))
  }
  legendtitle = glue("LATE({u1}, {u2})")

  return(TargetParameter(name, int_limits, int_constant, legendtitle))
}

#' Create LATE Target Parameter
#'
#' @param dgp Data Generating Process object.
#' @param l Model index (default is 1).
#' @return A TargetParameter object for ATT.
#' @export
att = function(dgp, l = 1){
  prd1 = dgp$pscore * dgp$densZ

  name = "ATT"

  int_limits = function(z){
    return(c(0, dgp$find_pscore(0)))
  }
  int_constant = function(l,d,z){
    return((l==1)*(2*d-1)/prd1*dgp$find_density(z))
  }


  return(TargetParameter(name, int_limits, int_constant, name))

}
