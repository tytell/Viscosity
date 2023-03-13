require(tidyverse)

fuzzy_time_match <- function(t1, t2) {
  tmatch <- rep_along(t2, NA)
  
  for (i in seq_along(t1)) {
    if ((i > 1) & (i < length(t1))) {
      dt = (t1[i] - t1[i-1])/2 + (t1[i+1] - t1[i])/2
    }
    else if (i == 1) {
      dt = (t1[i+1] - t1[i])/2
    }
    else {
      dt = (t1[i] - t1[i-1])/2
    }
    
    isclose = abs((t2 - t1[i])) < dt
    if (any(isclose, na.rm = TRUE)) {
      k <- which(isclose)
      
      if (length(k) > 1) {
        j <- which.min(abs(t2[k] - t1[i]))
        k <- k[j]
      }
      tmatch[k] <- t1[i]
    }
  }
  tmatch
}

fuzzy_time_match_group <- function(df, key, df2) {
  df2 <- filter(df2, trial == key$trial)
  
  df %>%
    mutate(tmatch = fuzzy_time_match(df2$t, df$t))
}