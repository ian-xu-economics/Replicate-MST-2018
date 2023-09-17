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
