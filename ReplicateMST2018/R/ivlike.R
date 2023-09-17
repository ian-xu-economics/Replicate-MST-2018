#' Compute IV Slope
#'
#' @param dgp Data Generating Process object.
#' @return A list containing the name and s function.
#' @export
ivslope = function(dgp) {

  name = "IV Slope"

  expZ = sum(dgp$suppZ * dgp$densZ)
  expD = sum(dgp$pscore * dgp$densZ)
  expDZ = sum(dgp$pscore * dgp$densZ * dgp$suppZ)

  covDZ = expDZ - expD * expZ

  s = function(d, z) {
    return((z - expZ) / covDZ)
  }

  ivlike = list(name = name, s = s, extra_info = NULL)

  return(ivlike)
}

#' Compute OLS Slope
#'
#' @param dgp Data Generating Process object.
#' @return A list containing the name and s function.
#' @export
olsslope = function(dgp) {

  name = "OLS Slope"

  prd1 = sum(dgp$pscore * dgp$densZ)

  s = function(d, z) {
    return((d - prd1) / (prd1 * (1 - prd1)))
  }

  IVLike = list(name = name, s = s, extra_info = NULL)
  return(IVLike)
}

#' Compute IV Slope Indicator
#'
#' @param dgp Data Generating Process object.
#' @param support Support set for the indicator.
#' @return A list containing the name, s function, support set, and legend titles.
#' @export
ivslope_indicator = function(dgp, support) {

  support = unique(support)
  support = sort(support)

  indices = match(support, dgp$suppZ)
  stopifnot(!is.na(indices))  # ensure support is in dgp.suppZ

  name = glue("IV Slope for 1(Z == z) for z in {[toString(support)]}.",
               .open = "[",
               .close = "]")

  expZind = function(i) {
    return(dgp$densZ[i])
  }

  expDZind = function(i) {
    return(dgp$pscore[i] * dgp$densZ[i])
  }

  expD = sum(dgp$pscore * dgp$densZ)

  covDZind = function(i) {
    return(expDZind(i) - expD * expZind(i))
  }

  s = lapply(indices, function(i) {
    return(function(d, z) {
      return(((z == dgp$suppZ[i]) - expZind(i)) / covDZind(i))
    })
  })

  legendtitle = c()
  for (z in support) {
    legendtitle = append(legendtitle, glue("IV Slope (1[Z = {z}])"))
  }

  ivlike = list(name = name, s = s, support = support, legendtitle = legendtitle)

  return(ivlike)
}

#' Create Saturated IV-like List
#'
#' @param dgp Data Generating Process object.
#' @return A list containing the name, s functions, and legend titles.
#' @export
make_slist = function(dgp) {
  name = "Saturated"

  combinations = expand.grid(d_bar = 0:1, z_bar = 1:length(dgp$suppZ))

  s = mapply(function(d_bar, z_bar) {
    function(d, z) {
      return((d == d_bar) * (z == dgp$suppZ[z_bar]))
    }
  }, combinations$d_bar, combinations$z_bar, SIMPLIFY = FALSE)

  d_string = c("(1 - D)", "D")
  z_string = glue("1[Z = {as.character(dgp$suppZ)}]")

  # Create legend titles
  legendtitle = c()
  for (z in z_string) {
    for (d in d_string) {
      legendtitle = append(legendtitle, paste0(d, z))
    }
  }

  ivlike = list(name = name, s = s, legendtitle = legendtitle)

  return(ivlike)
}

#' Compute Beta_s
#'
#' @param dgp Data Generating Process object.
#' @param slist List of s functions.
#' @param param Optional parameter for slist.
#' @return A numeric vector of beta_s values.
#' @export
compute_beta_s = function(dgp, slist, param = NA) {
  gamma_s = compute_gamma_s(list(dgp$mtrs[[1]]$basis,
                                  dgp$mtrs[[2]]$basis),
                             dgp,
                             slist = slist,
                             param = param)

  # Initialize beta_s
  beta_s = rep(NA, dim(gamma_s[[1]])[1])

  # Compute beta_s
  for (s in 1:length(beta_s)) {
    beta_s[s] = sum(gamma_s[[1]][s,,] * dgp$mtrs[[1]]$theta +
                       gamma_s[[2]][s,,] * dgp$mtrs[[2]]$theta)
  }

  return(beta_s)
}

#' Compute gamma_s for Each Model
#'
#' @param bases List of basis functions for each model.
#' @param dgp Data Generating Process object.
#' @param slist List of s functions.
#' @param param Optional parameter for slist.
#' @return A list of gamma_s arrays for each model.
#' @export
compute_gamma_s = function(bases, dgp, slist, param = NA) {
  # Placeholder for future functionality
  # d takes on value 0 and 1
  return(lapply(0:1, function(d) {
    compute_gamma_s_for_basis(bases[[d + 1]], d, dgp, slist = slist, param = param)
  }))
}

#' Compute gamma_s for a given basis
#'
#' @param basis Basis functions for a specific model.
#' @param d Treatment indicator (0 or 1).
#' @param dgp Data Generating Process object.
#' @param slist List of s functions.
#' @param param Optional parameter for slist.
#' @return A gamma_s array for the specific basis.
#' @export
compute_gamma_s_for_basis = function(basis, d, dgp, slist, param = NA) {

  if(slist == "ivslope"){
    slist = list(ivslope(dgp)$s)
  } else if(slist == "olsslope"){
    slist = list(olsslope(dgp)$s)
  } else if(slist == "ivslopeind"){
    slist = ivslope_indicator(dgp, support = param)$s
  } else if(slist == "saturated"){
    slist = make_slist(dgp)$s
  }

  # Initialize gamma_s
  gamma_s = array(0, dim = c(length(slist), length(basis$a), length(basis$b)))

  # Loop through each row of suppZ
  for (i in 1:length(dgp$suppZ)){
    z = dgp$suppZ[i]
    intlb = (1 - d) * dgp$pscore[i]
    intub = d * dgp$pscore[i] + (1 - d) * 1

    for(s in 1:length(slist)){
      for(aj in 1:length(basis$a)){
        for(ibk in 1:length(basis$ib)){
          gamma_s[s,aj,ibk] = gamma_s[s,aj,ibk] + basis$a[[aj]](z) * basis$ib[[ibk]](intlb, intub) * slist[[s]](d,z) * dgp$densZ[i]
        }
      }
    }
  }

  return(gamma_s)
}
