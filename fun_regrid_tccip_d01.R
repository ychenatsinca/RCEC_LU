

# this funtion try to resample the input array to a reference WRF domain/grid


fun_regrid <- function(nav.lon,nav.lat,array.in,type) {
print("regriding")		       
#default GFS grid
#nav_lon	<- c( seq(-180+(0.234375/2),180-(0.234375/2), 0.234375) )
#nav_lat <- c( seq(-90+(0.234375/2),  90-(0.234375/2), 0.234375) ) 
#str(array_in)
library(raster)
nx <- length(nav.lon)
ny <- length(nav.lat)
#ix1=1250
#ix2=1300
#iy1=450
#iy2=500
ix1=1
ix2=nx
iy1=1
iy2=ny
old.grid  <- raster(t(array.in[,ny:1]),
		    xmn=nav.lon[ix1],xmx=nav.lon[ix2], 
		    ymn=nav.lat[iy1],ymx=nav.lat[iy2],
		    crs= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") ) 
#
new.grid <- raster(ncol=249, nrow=249)
extent(new.grid) <- extent(108.4, 133.6 ,12.56, 35.49 )

if (type ==1) {
new.grid <- resample(old.grid, new.grid, method='bilinear')
}else{
new.grid <- resample(old.grid, new.grid, method='bilinear')
}
#
print(old.grid)
print(new.grid)
#image(new_grid) 
#par(mfrow=c(1,2))
#library(fields)
#image(old.grid)
#image(new.grid)
return(as.matrix(new.grid))
}
