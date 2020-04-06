library(gridclimind)
library(ncdf4)

in.dir <- "input/"
author.data <- list(ERA ="5")
out.dir <- "scPDSI_output/"

############ scPDSI ############

scPDSI.file <- "scPDSI_testfile.nc"
out.file <- sprintf("%s%s", out.dir, scPDSI.file)
if(file.exists(out.file)){
  file.remove(out.file)
}
input.files <- c(paste0(in.dir,"pet_0.1deg_regular_1980-2010_sub.nc"),
                 paste0(in.dir,"rr_0.1deg_regular_1980-2010_sub.nc"))

create.scPDSI.from.files(input.files, out.file, input.files[1], author.data,  parallel=FALSE, start=1980, end=2010, cal_start=1980, cal_end=2010)


##### test output #######

in.nc_1= nc_open("scPDSI_output/scPDSI_testfile_master.nc")
lats <- ncvar_get( in.nc_1, "latitude")   # coordinate variable
nlat <- length(lats)
lons <- ncvar_get( in.nc_1, "longitude")   # coordinate variable
nlon <- length(lons)
tm <- ncvar_get( in.nc_1, "time")/24
nt <- length(tm)
indat1 <- as.matrix(ncvar_get( in.nc_1, "scPDSI", start=c(1,1,1), count=c(nlon,nlat,nt)) )

in.nc_2= nc_open("scPDSI_output/scPDSI_testfile.nc")
indat2 <- as.matrix(ncvar_get( in.nc_2, "scPDSI", start=c(1,1,1), count=c(nlon,nlat,nt)) )

if(all.equal.numeric(indat1, indat2, tolerance=0.0001)==TRUE)
{
  print("scPSDI output matches master file")
}else{
  print("Warning: scPDSI output differs from master file")
}
nc_close(in.nc_1)
nc_close(in.nc_2)
