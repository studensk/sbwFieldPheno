#' Predict days to development given a set of daily or hourly temperature data
#' @param weather user-defined weather dataset; requires columns:
#' Temp (temperature in degrees Celsius), Date (YYYY-MM-DD), Hour (integer, 0-23)
#' @param sbwcolony colony identifier; one of 'AB', 'IPQL', 'NB', 'NWT', 'ON', 'QC'
#' @param period observation interval of weather data; can be one of 'hour' or 'day'
#' @param stage larval development stage to be predicted; one of 'L3', 'L4', 'L5',
#' 'L6', 'Pupa' or 'All'
#' @param ecdf if TRUE, return empirical CDF (density object) of predicted days
#' @param n.post number of posterior samples to simulate from (max 100)
#' @param individuals number of individuals in the population per posterior sample
#' @returns Matrix or array of Julian dates; dimensions are `individuals*n.post*stages`
#' @export

dev_days <- function(weather,
                     sbwcolony = 'NB',
                     period = 'hour',
                     stage = 'Pupa',
                     ecdf = FALSE,
                     n.post = 100,
                     individuals = 100) {

  n.post <- pmin(n.post, 100)

  if (ecdf & stage == 'All') {
    warning('Cannot track individuals across stages with ecdf = TRUE')
  }

  div <- ifelse(period == 'hour', 24, 1)

  stages <- paste0('L', 2:6)
  stages.end <- c(stages[2:5], 'Pupa')

  if (stage == 'All') {
    ind.ind <- 1:5
  }
  else {
    ind.ind <- which(stages.end == stage)
  }

  weather.clean <- weather |>
    dplyr::mutate(Date = as.Date(Date),
           Temp = pmax(Temp, 0),
           DateTime = Date + lubridate::hours(Hour),
           Year = lubridate::year(Date),
           JDay = lubridate::yday(Date)) |>
    dplyr::arrange(DateTime)
  all.dates <- weather.clean |>
    dplyr::select(DateTime) |>
    dplyr::reframe(DateTime = seq(min(DateTime), max(DateTime), by = period)) |>
    dplyr::left_join(weather.clean, by = dplyr::join_by(DateTime)) |>
    subset(is.na(Temp)) |>
    unique() |>
    nrow()
  if (all.dates > 0) {
    warning(paste0('Dataset contains ', all.dates,
                   ' missing site-', period,
                   ifelse(all.dates > 1, 's', '')))
  }

  info.sub <- subset(curve_info, colony == sbwcolony) |>
    dplyr::arrange(run.index, stage)
  inds <- sort(info.sub$index)
  run.inds <- sort(sample(unique(info.sub$run.index),
                          n.post, replace = FALSE))

  stage.mat <- array(dim = c(individuals,
                             length(run.inds),
                             length(ind.ind)))
  nweather <- nrow(weather.clean)
  isplines <- matrix(nrow = nweather, ncol = length(stages))

  for (i in seq(run.inds)) {
    sub <- subset(info.sub, run.index == run.inds[i])
    new.inds <- sub$index

    for (j in seq(new.inds)) {
      ind <- new.inds[[j]]
      isplines[,j] <- with(subset(curves, index == ind), {
        spl <- splines::interpSpline(temp, rate/div)
        predict(spl, weather.clean$Temp)$y
      })
    }

    stage.mat[,i,] <- t(replicate(individuals, {
      ind <- 0
      inds <- vector(length = length(stages))
      mult <- exp(stats::rnorm(nrow(sub), 0, sub$s_eps))

      for (s in seq(stages)) {
        if (stage == stages[s]) {break}
        ind <- match(TRUE,
                     cumsum(isplines[,s]*(seq(nweather)>ind)*mult[s]) > 1)

        inds[s] <- weather.clean$JDay[ind]
      }
      return(inds[ind.ind])
    }))
  }
  if (ecdf) {
    den <- apply(stage.mat, 3, function(x) {
      d <- stats::density(x)
      d$y <- d$y/sum(d$y)
      return(d)
    })
    if (length(den) == 1) {
      return(den[[1]])
    }
    else {
      names(den) <- stages.end
      return(den)
    }
  }
  else {
    return(drop(stage.mat))
  }
}

