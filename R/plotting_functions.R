#' Get Non-Zero Weights from Weight Data Frame
#'
#' @param weight_df Data frame containing weights.
#' @param d Treatment indicator (0 or 1).
#' @return A data frame containing non-zero weights.
#' @export
get_non_zero_weights = function(weight_df, d){

  stopifnot(d == 0 | d == 1)

  # Find column for d = 0 case
  column_d0 = which(str_detect(colnames(weight_df), "0") == TRUE)

  # Find column for d = 1 case
  column_d1 = which(str_detect(colnames(weight_df), "1") == TRUE)

  num_rows = nrow(weight_df)

  new_weight_df = NULL

  if(weight_df$u[num_rows] != 1){
    weight_df = weight_df %>%
      rbind(0)
  }

  new_weight_df = data.frame(
    u_start = head(weight_df$u, -1),
    u_end = weight_df$u[-1]
    )

  if(new_weight_df$u_end[nrow(new_weight_df)] == 0){
    new_weight_df$u_end[nrow(new_weight_df)] = 1
  }

  if(d == 0){
    new_weight_df$avg_weight = head(weight_df[,column_d0], -1)
  } else{
    new_weight_df$avg_weight = head(weight_df[,column_d1], -1)
  }

  # Define a small tolerance value
  tolerance = 1e-12

  # Return non-zero weights
  return(new_weight_df[abs(new_weight_df$avg_weight) > tolerance, ])

}

#' Get Weights from Weight Data Frame
#'
#' @param weight_df Data frame containing weights.
#' @param d Treatment indicator (0 or 1).
#' @return A data frame containing weights.
#' @export
get_weights = function(weight_df, d){

  stopifnot(d == 0 | d == 1)

  # Find column for d = 0 case
  column_d0 = which(str_detect(colnames(weight_df), "0") == TRUE)

  # Find column for d = 1 case
  column_d1 = which(str_detect(colnames(weight_df), "1") == TRUE)

  num_rows = nrow(weight_df)

  new_weight_df = NULL

  if(weight_df$u[num_rows] != 1){
    weight_df = weight_df %>%
      rbind(0)
  }

  new_weight_df = data.frame(
    u_start = head(weight_df$u, -1),
    u_end = weight_df$u[-1]
  )

  if(new_weight_df$u_end[nrow(new_weight_df)] == 0){
    new_weight_df$u_end[nrow(new_weight_df)] = 1
  }

  if(d == 0){
    new_weight_df$avg_weight = head(weight_df[,column_d0], -1)
  } else{
    new_weight_df$avg_weight = head(weight_df[,column_d1], -1)
  }

  return(new_weight_df)

}

#' Create a ggplot based on various conditions and data
#'
#' This function generates a ggplot object based on the provided data, conditions, and parameters.
#' The function supports different basis types and allows for upper and lower bounds.
#'
#' @param dgp Data Generating Process object containing necessary data and parameters.
#' @param d The treatment indicator (0 or 1).
#' @param upper Logical indicating whether to use upper bounds.
#' @param weights_tp List of weights for the treatment parameters.
#' @param weights_ivlike List of weights for IV-like estimators.
#' @param bounds Object containing upper and lower bound solutions.
#' @param basis The type of basis to use ("constant_spline", "weighted_bernstein", or "dgp_mtr").
#' @param title The title for the plot.
#'
#' @return A ggplot object.
#' @export
create_plot = function(dgp, d, upper, weights_tp, weights_ivlike, bounds, basis, title){

  legend_names = NULL

  final_plot = ggplot()

  if(basis == "constant_spline"){
    constant_spline_weight = c(0, 1, 0.35, 0.9, dgp$pscore) %>%
      unique() %>%
      sort() %>%
      data.frame(u = .) %>%
      filter(u != 1)

    if(upper == TRUE){
      constant_spline_weight = constant_spline_weight %>%
        cbind(d0 = bounds$upper_solution_d0) %>%
        cbind(d1 = bounds$upper_solution_d1)

      final_plot = final_plot +
        geom_segment(data = constant_spline_weight %>%
                       get_weights(d = d),
                     aes(x = u_start,
                         xend = u_end,
                         y = avg_weight,
                         yend = avg_weight,
                         color = "Maximizing MTR"),
                     linewidth = 1)

      legend_names = c(legend_names, "Maximizing MTR")
    } else{
      constant_spline_weight = constant_spline_weight %>%
        cbind(d0 = bounds$lower_solution_d0) %>%
        cbind(d1 = bounds$lower_solution_d1)

      final_plot = final_plot +
        geom_segment(data = constant_spline_weight %>%
                       get_weights(d = d),
                     aes(x = u_start,
                         xend = u_end,
                         y = avg_weight,
                         yend = avg_weight,
                         color = "Minimizing MTR"),
                     linewidth = 1)

      legend_names = c(legend_names, "Minimizing MTR")
    }
  } else if(basis == "weighted_bernstein"){
    if(upper == TRUE){
      if(d == 1){
        final_plot = final_plot +
          stat_function(fun = weighted_bernstein_basis,
                        args = list(K = length(bounds$upper_solution_d1) - 1,
                                    weights = bounds$upper_solution_d1),
                        aes(color = "Maximizing MTR"))
      } else{
        final_plot = final_plot +
          stat_function(fun = weighted_bernstein_basis,
                        args = list(K = length(bounds$upper_solution_d0) - 1,
                                    weights = bounds$upper_solution_d0),
                        aes(color = "Maximizing MTR"))
      }
      legend_names = c(legend_names, "Maximizing MTR")
    } else{
      if(d == 1){
        final_plot = final_plot +
          stat_function(fun = weighted_bernstein_basis,
                        args = list(K = length(bounds$lower_solution_d1) - 1,
                                    weights = bounds$lower_solution_d1),
                        aes(color = "Minimizing MTR"))
      } else{
        final_plot = final_plot +
          stat_function(fun = weighted_bernstein_basis,
                        args = list(K = length(bounds$lower_solution_d0) - 1,
                                    weights = bounds$lower_solution_d0),
                        aes(color = "Minimizing MTR"))
      }
      legend_names = c(legend_names, "Minimizing MTR")
    }
  } else if(basis == "dgp_mtr"){
    if(d == 1){
      final_plot = final_plot +
        stat_function(fun = weighted_bernstein_basis,
                      args = list(K = length(dgp$mtrs[[2]]$theta) - 1,
                                  weights = dgp$mtrs[[2]]$theta),
                      aes(color = "DGP MTR"))
    } else{
      final_plot = final_plot +
        stat_function(fun = weighted_bernstein_basis,
                      args = list(K = length(dgp$mtrs[[1]]$theta) - 1,
                                  weights = dgp$mtrs[[1]]$theta),
                      aes(color = "DGP MTR"))
    }

    legend_names = c(legend_names, "DGP MTR")
  }

  for(i in 1:length(weights_tp)){
    weight = compute_average_weights(weights_tp[[i]],
                                     dgp)

    if(weights_tp[[i]]$name == "LATE(u1, u2)"){
      legendTitle = glue("LATE({weights_tp[[i]]$int_limits[1]}, {weights_tp[[i]]$int_limits[2]})")
    }

    final_plot = final_plot +
      geom_segment(data = weight %>%
                     get_non_zero_weights(d = d),
                   aes(x = u_start,
                       xend = u_end,
                       y = avg_weight,
                       yend = avg_weight,
                       color = legendTitle),
                   linewidth = 1)

    legend_names = c(legend_names, legendTitle)
  }

  weights_ivlike_df = NULL

  for(i in 1:length(weights_ivlike)){

    if(str_detect(weights_ivlike[[i]]$name, "IV Slope for") == TRUE){
      for(j in 1:length(weights_ivlike[[i]]$s)){
        weights_ivlike_df = weights_ivlike_df %>%
          rbind(compute_average_weights_ivlike(list(name = weights_ivlike[[i]]$legendtitle[j], s = weights_ivlike[[i]]$s[[j]]), dgp) %>%
                  get_non_zero_weights(d = d) %>%
                  mutate(name = weights_ivlike[[i]]$legendtitle[j]))

        legend_names = c(legend_names, weights_ivlike[[i]]$legendtitle[j])
      }
    } else if(weights_ivlike[[i]]$name == "Saturated"){
      if(d == 0){
        for(j in seq(1, length(weights_ivlike[[i]]$s), 2)){
          weights_ivlike_df = weights_ivlike_df %>%
            rbind(compute_average_weights_ivlike(list(name = weights_ivlike[[i]]$legendtitle[j],
                                                      s = weights_ivlike[[i]]$s[[j]]),
                                                 dgp) %>%
                    get_non_zero_weights(d = d) %>%
                    mutate(name = weights_ivlike[[i]]$legendtitle[j]))

          legend_names = c(legend_names, weights_ivlike[[i]]$legendtitle[j])
        }
      } else{
        for(j in seq(2, length(weights_ivlike[[i]]$s), 2)){
          weights_ivlike_df = weights_ivlike_df %>%
            rbind(compute_average_weights_ivlike(list(name = weights_ivlike[[i]]$legendtitle[j],
                                                      s = weights_ivlike[[i]]$s[[j]]),
                                                 dgp) %>%
                    get_non_zero_weights(d = d) %>%
                    mutate(name = weights_ivlike[[i]]$legendtitle[j]))

          legend_names = c(legend_names, weights_ivlike[[i]]$legendtitle[j])
        }
      }
    } else{
      weights_ivlike_df = weights_ivlike_df %>%
        rbind(compute_average_weights_ivlike(weights_ivlike[[i]], dgp) %>%
                get_non_zero_weights(d = d) %>%
                mutate(name = weights_ivlike[[i]]$name))

      legend_names = c(legend_names, weights_ivlike[[i]]$name)
    }

  }

  final_plot = final_plot +
    geom_segment(data = weights_ivlike_df,
                 aes(x = u_start,
                     xend = u_end,
                     y = avg_weight,
                     yend = avg_weight,
                     color = weights_ivlike_df$name),
                 linewidth = 1) +
    scale_y_continuous(name = NULL) +
    scale_x_continuous(name = "u",
                       limits = c(0,1)) +
    scale_color_brewer(name = NULL,
                       palette = "Set1",
                       breaks = legend_names) +
    theme(legend.text=element_text(size = 14),
          legend.position = "bottom") +
    ggtitle(glue("{title} (D = {d})"))

  if(length(legend_names) > 3){
    final_plot = final_plot +
      guides(color = guide_legend(nrow = 2, byrow = TRUE))
  }

  return(final_plot)

}
