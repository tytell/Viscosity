require(tidyverse)
require(hdf5r)

write_hdf5 <- function(df, filename, mode='w') {
  file.h5 <- H5File$new(filename, mode)
  
  df <-
    df %>%
    mutate(across(where(is.factor), ~fct_explicit_na(.x, na_level='(Missing)')))
  
  cols <- colnames(df)
  h5attr(file.h5, 'colnames') <- cols
  for (c1 in cols) {
    file.h5[[c1]] <- df[[c1]]
  }
  
  file.h5$close_all()
}

read_hdf5 <- function(filename) {
  file.h5 <- H5File$new(filename, 'r')
  
  cols <- h5attr(file.h5, 'colnames')
  df <- map_dfc(cols, ~ file.h5[[.x]]$read())
  colnames(df) <- cols
  
  file.h5$close_all()
  
  df %>%
    mutate(across(where(is.factor), ~ if_else(.x != '(Missing)', .x, NA_integer_) %>%
                    fct_drop('(Missing)')))
}

write_hdf5_old <- function(df, filename, mode='w') {
  file.h5 <- H5File$new(filename, mode)
  
  df <-
    df %>%
    mutate(across(where(is.factor), ~fct_explicit_na(.x, na_level='(Missing)')))
  
  file.h5[['data']] <- df
  
  file.h5$close_all()
}

read_hdf5_old <- function(filename) {
  file.h5 <- H5File$new(filename, 'r')
  dset <- file.h5[['data']]
  
  df <- as_tibble(dset$read())
  
  file.h5$close_all()
  
  df %>%
    mutate(across(where(is.factor), ~ if_else(.x != '(Missing)', .x, NA_integer_) %>%
                    fct_drop('(Missing)')))
}
