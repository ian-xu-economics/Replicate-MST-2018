#' DGP R6 Class
#'
#' An R6 class for representing a Data Generating Process (DGP).
#'
#' @field suppZ Support of Z.
#' @field densZ Density of Z.
#' @field pscore Propensity score.
#' @field mtrs List of Marginal Treatment Response (MTR) functions.
DGP = R6::R6Class("DGP",
              public = list(
                suppZ = NULL,
                densZ = NULL,
                pscore = NULL,
                mtrs = NULL,

                #' Initialize DGP object
                #'
                #' @param suppZ Support of Z.
                #' @param densZ Density of Z.
                #' @param pscore Propensity score.
                #' @param mtrs List of MTR functions.
                #' @return A new DGP object.
                initialize = function(suppZ, densZ, pscore, mtrs) {
                  # Check that dimensions of suppZ and lengths of densZ and pscore match
                  if (length(suppZ) != length(densZ)) {
                    stop("Dimensions of suppZ and length of densZ do not match.")
                  }

                  if (length(densZ) != length(pscore)) {
                    stop("Lengths of densZ and pscore do not match.")
                  }

                  # Check that densZ values are non-negative
                  if (any(densZ < 0)) {
                    stop("densZ contains negative values.")
                  }

                  # Check that densZ values sum to 1
                  if (sum(densZ) != 1) {
                    stop("densZ values do not sum to 1.")
                  }

                  self$suppZ = suppZ
                  self$densZ = densZ
                  self$pscore = pscore
                  self$mtrs = mtrs
                },

                #' Find Propensity Score for Given Z
                #'
                #' @param z Value of Z.
                #' @return Propensity score corresponding to Z.
                find_pscore = function(z) {
                  idx = which(self$suppZ == z)
                  if (length(idx) != 1) {
                    stop("Expected to find one and only one match for z.")
                  }
                  return(self$pscore[idx])
                },

                #' Find Density for Given Z
                #'
                #' @param z Value of Z.
                #' @return Density corresponding to Z.
                find_density = function(z) {
                  idx = which(self$suppZ == z)
                  if (length(idx) != 1) {
                    stop("Expected to find one and only one match for z.")
                  }
                  return(self$densZ[idx])
                }
              )
            )

#' Create DGP Object for MST2018 Model
#'
#' This function creates a DGP object based on the MST2018 model.
#'
#' @return A new DGP object.
#' @export
dgp_MST2018 = function() {
  suppZ = 0:2
  densZ = c(0.5, 0.4, 0.1)
  pscoreZ = c(0.35, 0.6, 0.7)

  basis0 = bernstein_basis(2)
  theta0 = c(0.6, 0.4, 0.3)
  mtr0 = MTR(basis0, theta0)

  basis1 = bernstein_basis(2)
  theta1 = c(0.75, 0.5, 0.25)
  mtr1 = MTR(basis1, theta1)

  return(DGP$new(suppZ, densZ, pscoreZ, list(mtr0, mtr1)))
}

#' Create DGP Object for MST2018 Model
#'
#' This function creates a DGP object based on inputs.
#' @field suppZ Support of Z.
#' @field densZ Density of Z.
#' @field pscore Propensity score.
#' @field mtrs List of Marginal Treatment Response (MTR) functions.
#' @return A new DGP object.
#' @export
create_dgp = function(suppZ, densZ, pscoreZ, basis0, theta0, basis1, theta1){

  mtr0 = MTR(basis0, theta0)
  mtr1 = MTR(basis1, theta1)

  return(DGP$new(suppZ, densZ, pscoreZ, list(mtr0, mtr1)))

}
