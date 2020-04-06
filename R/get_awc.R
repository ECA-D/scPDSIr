

get_awc <- function(lat, lon)
{
  l <- length(lon)
  awc <- array(dim=l)
  
  lon_index <- which.min(abs(lon - awc_lon))
  
  # lat index
  lat_index <- which.min(abs(lat - awc_lat))
  # lon index
  lon_index <- which.min(abs(lon - awc_lon))
  indat <- awc_dat[lon_index, lat_index]
  if(is.na(indat))
  {
    new_coord <- array(dim=c(8,2))
    new_coord[1,] <- c(1,-1)
    new_coord[2,] <- c(1,0)
    new_coord[3,] <- c(1,1)
    new_coord[4,] <- c(0,-1)
    new_coord[5,] <- c(0,1)
    new_coord[6,] <- c(-1,-1)
    new_coord[7,] <- c(-1,0)
    new_coord[8,] <- c(-1,1)
    
    for(c in 1:8) # if it doesn't find one then it is a ship, so leave coords as is...
    {
      lon_index_new <- lon_index+new_coord[c,1]
      if(lon_index_new <= 0) lon_index_new <- lon_index
      
      lat_index_new <- lat_index+new_coord[c,2]
      if(lat_index_new <= 0) lat_index_new <- lat_index
      
      check <- awc_dat[lon_index_new, lat_index_new]
      if(is.na(check) == 0)
      {
        lon_index <- lon_index_new
        lat_index <- lat_index_new
        break
      }
    }
  }
  awc <- awc_dat[lon_index, lat_index]
  return(awc)
}