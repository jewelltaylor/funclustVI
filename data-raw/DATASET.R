## code to prepare `DATASET` dataset goes here


get_real_data <- function() {
  
  res_temp = as.vector(weathercan::weather_dl(station_ids = 1776, start = "1994-01-01", end = "2003-12-31", interval = "day")$max_temp)
  tor_temp = as.vector(weathercan::weather_dl(station_ids = 5051, start = "1994-01-01", end = "2003-12-31", interval = "day")$max_temp)

  cl = matrix(0, 18, 31)
  count = 1
  for (i in seq(1, 3285, 365)) {
    upperbound1 = i + 364
    res_temp_curr = tidyr::replace_na(res_temp[i:upperbound1], 0)[187:217]
    cl[count, ] = res_temp_curr
    
    count = count + 1
  }
  
  for (j in seq(1, 3285, 365)) {
    upperbound2 = j + 364
    tor_temp_curr = tidyr::replace_na(tor_temp[j:upperbound2], 0)[187:217]
    cl[count, ] = tor_temp_curr
    count = count + 1
  }
  
  return(cl)
}

dataset = get_real_data()

usethis::use_data(dataset, overwrite = TRUE)


