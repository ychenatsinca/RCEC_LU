

#
source("src_function_ncdf4.R")
source("fun_regrid_tccip_d02.R")

#load tccip d02 geo data
d02.fname = "../geo_em.d02.nc"
d02.geo <- fun_read_nc(arg1=d02.fname)

mask_arr_a  <- d02.geo$LANDMASK
nav_lon_a   <- d02.geo$CLONG
nav_lat_a   <- d02.geo$CLAT

nx_a=220 
ny_a=280
nz_a=24



work.year = 1982 

bn.5km <- fun_read_nc(arg1="BN_Ratio_map_tw_5km.nc")
#regrid 5km bn.ratio to 2km 
bn.2km <- fun_regrid(nav.lon = bn.5km$lon, 
	             nav.lat = bn.5km$lat,
		     array.in = bn.5km$BN_Ratio,
		     type=1)
#
bn.2km <- t(bn.2km)[,ny_a:1]


# load 7pft land cover data @ 500m
#/work/ychen/ycmeet/SPOT_Classification/finalmap/TWD_1997_TM_Taiwan/7_class/7classfinal
dir.name <- c("/work/ychen/ycmeet/SPOT_Classification/finalmap/TWD_1997_TM_Taiwan/7_class/7classfinal/")
lc.data <- fun_read_nc(arg1=paste(dir.name,work.year,"_7PFT.nc",sep=""))

# calculate the 7pft fraction in 2km resolution grid
nx_b=747
ny_b=900

lon_b  <- array(NA,dim=c(nx_b*ny_b))
lat_b  <- array(NA,dim=c(nx_b*ny_b))
arr_b  <- array(NA,dim=c(nx_b*ny_b))



#convert RCEC-PFT NO to WRF-PFT(USGS) NO
lc.data$LC[lc.data$LC==13 ] <- 16 #water
lc.data$LC[lc.data$LC==1 ]  <- 13 #forest
lc.data$LC[lc.data$LC==2 ]  <- 5  #grass
lc.data$LC[lc.data$LC==8 ]  <- 1  #built
lc.data$LC[lc.data$LC==26]  <- 19 #bare
lc.data$LC[lc.data$LC==4 ]  <- 2  #dry crop
lc.data$LC[lc.data$LC==52]  <- 3  #irr crop



# for irregular grid spacing a data.frame to store lon,lat information by nx*ny elements
for (j in 1:ny_b) {
   for (i in 1:nx_b) {
        id <- i + (j-1)*nx_b
        lon_b[id] <- lc.data$nav_lon[i,j]
        lat_b[id] <- lc.data$nav_lat[i,j]
        arr_b[id] <- lc.data$LC[i,j]
   }
}

# combine lon,lat 
lonlat_b <- cbind(lon_b,lat_b)


lcfrac <- array(0, dim=c(nx_a,ny_a,nz_a) )
lcdom  <- array(0, dim=c(nx_a,ny_a))
tmp.arr <- array(0, dim=c(nx_a, ny_a))
# copy the original data, we only update the land cover data over the Taiwan
lcfrac <- d02.geo$LANDUSEF
lcdom  <- d02.geo$LU_INDEX
# set resolution at 2km
dx=0.02
dy=0.02 

ld.go <- TRUE
if(ld.go) { 
for (i in 1:(nx_a)) {
    for (j in 1:(ny_a) ) {
#for (i in 100:130 ) {
#    for (j in 1:210) {
                      # we only do the land points over Taiwan
		        if ( (d02.geo$LANDMASK[i,j] == 1 ) &
			     (nav_lon_a[i,j] >= 120.0) & (nav_lon_a[i,j] <= 122.05) &
                             (nav_lat_a[i,j] >= 21.75) & (nav_lat_a[i,j] <=   25.5) )  {
                                # initialized lcfrac[i,j,] = 0 at the point
			        lcfrac[i,j, ] = 0	
				# copy lon_lat_b information to lalo
			        lalo      <- lonlat_b
			        arr_lalo  <- arr_b
			        # find the pixels in the each TCCIP geo_grid (lon_a, lat_a)
			        tmp_id <- which(lalo[,1] >=  nav_lon_a[i,j]-dx/2 & lalo[,1] < nav_lon_a[i,j]+dx/2 &
					        lalo[,2] >=  nav_lat_a[i,j]-dy/2 & lalo[,2] < nav_lat_a[i,j]+dy/2 )
		                grid_lc <- arr_lalo[tmp_id]	
				grid_cont <- length(grid_lc)
                                #rm na values
			        grid_lc <- grid_lc[!is.na(grid_lc)]	
			        #if ( grid_cont >= 16 ) {
		                    # calculate the fraction of land cover       
			            for ( k in 1:nz_a) { 
                                            #print(grid_lc)		          	            
					    if ( (!is.null(grid_lc[grid_lc == k]) ) & (length(grid_lc[grid_lc == k]) > 1 )) {
				   	      lcfrac[i,j,k] <- as.numeric(length(grid_lc[grid_lc == k ])) / as.numeric(grid_cont)
			                    #  print(lcfrac[i,j,k])
					    } 
			            }# end k loop
	      		#print(paste("dominate type:",aaa))
				#print(paste("land point:","lon1:",nav_lon_a[i,j], " lon2:", nav_lon_a[i+1,j+1], sep="" ) )
					  #print(paste("           ","lat1:",nav_lat_a[i,j], " lat2:", nav_lat_a[i+1,j+1], sep="" ) )
		                #} # end if of grid_cont 
			} # end if of mask
                }# end of j
     print(paste("i=",i,sep="")) 
}# end of i
}# end ld.go

# use 2km BN ratio to allocate Forests Cover Fraction to differnt type of forests  
# to Evergreen Broadleaf(13) BN > 0.55
# to Evergreen Needleaf(14) if BN <= 0.45 
# to Mix (15) if 45 < BN <= 0.55

for (i in 1:nx_a) {
    for (j in 1:ny_a) {
       # do the region over Taiwan 
       if ( (d02.geo$LANDMASK[i,j] == 1 )  &
	    (nav_lon_a[i,j] >= 120.0) & (nav_lon_a[i,j] <= 122.05)  &
            (nav_lat_a[i,j] >= 21.75) & (nav_lat_a[i,j] <=   25.5)  ) {

	     tmp.frac <- lcfrac[i,j,13]
      #  if ( bn.2km[i,j] <= 0.3){
#	    	lcfrac[i,j,13] <- tmp.frac*.2
#		lcfrac[i,j,14] <- tmp.frac*.7
#		lcfrac[i,j,15] <- tmp.frac*.1
#	}
 #       if ( (bn.2km[i,j] > 0.5) & (bn.2km[i,j] <= 0.7) ) {
	     lcfrac[i,j,13] <- tmp.frac*(bn.2km[i,j])
             lcfrac[i,j,14] <- tmp.frac*(1.-bn.2km[i,j])*0.3
	     lcfrac[i,j,15] <- tmp.frac*(1.-bn.2km[i,j])*0.7
#	}
#	if ( bn.2km[i,j] > 0.7) {
#	   	lcfrac[i,j,13] <- tmp.frac*.7
#		lcfrac[i,j,14] <- tmp.frac*.2
#		lcfrac[i,j,15] <- tmp.frac*.1
#	}	
       # add a loop for checking total area
             tmp.sum <- sum(lcfrac[i,j,],na.rm=T) 
             tmp.rsl <- tmp.sum - 1.0  
	     if ( abs(tmp.rsl)  > 1E-3)  {
              print( paste("warning: total cover fraction error!", "Total sum:", tmp.sum,
	  	       "lat:", nav_lat_a[i,j], "lon:",nav_lon_a[i,j],sep=""))  
               # allocare the residule to bare soilfor the condiction of the tmp.sum  below 1
               if ( tmp.rsl < 0) {
	          lcfrac[i,j,19] <- lcfrac[i,j,19] - tmp.rsl
	       }
	       # do nothing for residule > 0, keep the ratio

	      }# end if   
       }# end if 
       
       #find dominate land cover type	   
       tmp.lc <- lcfrac[i,j,]
       lcdom[i,j] <- which( tmp.lc  == max(tmp.lc, na.rm=T))[1]  
	
    }#end for j 	    
}#end for i 	


#generate the geo.d02 file


ld.create.nc <- TRUE
if (ld.create.nc) {

      

  #====== setup a NCDF file "define dimensions and  variables"  
  #====== define the  dimension ================================ 
  dim1 <- ncdim_def("Time","", 1, unlim=TRUE) 
  dim2 <- ncdim_def("DateStrLen","",  1:19, create_dimvar=FALSE) 
  dim3 <- ncdim_def("west_east","",  1:220, create_dimvar=FALSE) 
  dim4 <- ncdim_def("south_north","",1:280, create_dimvar=FALSE)
  dim5 <- ncdim_def("south_north_stag","",1:281, create_dimvar=FALSE)
  dim6 <- ncdim_def("west_east_stag","",1:221,   create_dimvar=FALSE) 
  dim7 <- ncdim_def("land_cat","",1:24,     create_dimvar=FALSE)  
  dim8 <- ncdim_def("soil_cat","",1:16,     create_dimvar=FALSE) 
  dim9 <- ncdim_def("month","",1:12,        create_dimvar=FALSE)  
  #============================================================= 
  mv <- 1e+10 
  #      ncvar_ def( name, units, dim, missval, longname=name, prec="float",
  #      shuffle=FALSE, compression=NA, chunksizes=NA, verbose=FALSE )

  var0 <- ncvar_def(name="Times",units="",  dim=list(dim2,dim1), ,prec="char") 
  var1 <-  ncvar_def(name="XLAT_M", units="degrees latitude", dim=list(dim3,dim4,dim1),  missval=mv ,prec="float") 
  var2 <-  ncvar_def(name="XLONG_M",units="degrees longitude", dim=list(dim3,dim4,dim1), missval=mv, prec="float")
  var3 <-  ncvar_def(name="XLAT_V",units="degrees latitude",  dim=list(dim3,dim5,dim1),  missval=mv, prec="float")
  var4 <-  ncvar_def(name="XLONG_V",units="degrees longitude", dim=list(dim3,dim5,dim1), missval=mv, prec="float")
  var5 <-  ncvar_def(name="XLAT_U",units="degrees latitude",  dim=list(dim6,dim4,dim1),  missval=mv, prec="float")
  var6 <-  ncvar_def(name="XLONG_U",units="degrees longitude", dim=list(dim6,dim4,dim1), missval=mv, prec="float")
  var7 <-  ncvar_def(name="CLAT",units="degrees latitude",  dim=list(dim3,dim4,dim1),    missval=mv, prec="float")
  var8 <-  ncvar_def(name="CLONG",units="degrees longitude", dim=list(dim3,dim4,dim1),   missval=mv, prec="float")
  var9 <-  ncvar_def(name="MAPFAC_M",units="none",dim=list(dim3,dim4,dim1),    missval=mv, prec="float")
  var10 <- ncvar_def(name="MAPFAC_V",units="none",dim=list(dim3,dim5,dim1),   missval=mv, prec="float")
  var11 <- ncvar_def(name="MAPFAC_U",units="none",dim=list(dim6,dim4,dim1),   missval=mv, prec="float")
  var12 <- ncvar_def(name="MAPFAC_MX",units="none",dim=list(dim3,dim4,dim1),  missval=mv, prec="float")
  var13 <- ncvar_def(name="MAPFAC_VX",units="none",dim=list(dim3,dim5,dim1),  missval=mv, prec="float")
  var14 <- ncvar_def(name="MAPFAC_UX",units="none",dim=list(dim6,dim4,dim1),  missval=mv, prec="float")
  var15 <- ncvar_def(name="MAPFAC_MY",units="none",dim=list(dim3,dim4,dim1), missval=mv, prec="float")
  var16 <- ncvar_def(name="MAPFAC_VY",units="none",dim=list(dim3,dim5,dim1), missval=mv, prec="float")
  var17 <- ncvar_def(name="MAPFAC_UY",units="none",dim=list(dim6,dim4,dim1),  missval=mv, prec="float")
  var18 <- ncvar_def(name="E",units="-", dim=list(dim3,dim4,dim1),          missval=mv, prec="float")
  var19 <- ncvar_def(name="F",units="-",  dim=list(dim3,dim4,dim1), 	    missval=mv, prec="float")
  var20 <- ncvar_def(name="SINALPHA",units="none", dim=list(dim3,dim4,dim1),  missval=mv, prec="float")
  var21 <- ncvar_def(name="COSALPHA",units="none", dim=list(dim3,dim4,dim1),  missval=mv, prec="float")
  var22 <- ncvar_def(name="LANDMASK",units="none", dim=list(dim3,dim4,dim1),  missval=mv, prec="float")
  var23 <- ncvar_def(name="XLAT_C",units="degrees latitude",dim=list(dim6,dim5,dim1),  missval=mv, prec="float")
  var24 <- ncvar_def(name="XLONG_C",units="degrees longitude",dim=list(dim6,dim5,dim1),missval=mv, prec="float")
  var25 <- ncvar_def(name="SINALPHA_U",units="none",dim=list(dim6,dim4,dim1),          missval=mv, prec="float")
  var26 <- ncvar_def(name="COSALPHA_U",units="none",dim=list(dim6,dim4,dim1), 	    missval=mv, prec="float")
  var27 <- ncvar_def(name="SINALPHA_V",units="none",dim=list(dim3,dim5,dim1), 	    missval=mv, prec="float")
  var28 <- ncvar_def(name="COSALPHA_V",units="none",dim=list(dim3,dim5,dim1),       missval=mv, prec="float")
  var29 <- ncvar_def(name="LANDUSEF",units="category",dim=list(dim3,dim4,dim7,dim1),  missval=mv, prec="float")
  var30 <- ncvar_def(name="LU_INDEX",units="category",dim=list(dim3,dim4,dim1),       missval=mv, prec="float")
  var31 <- ncvar_def(name="HGT_M",units="meters MSL", dim=list(dim3,dim4,dim1),       missval=mv, prec="float")
  var32 <- ncvar_def(name="SOILTEMP",units="Kelvin",  dim=list(dim3,dim4,dim1),       missval=mv, prec="float")
  var33 <- ncvar_def(name="SOILCTOP",units="category",dim=list(dim3,dim4,dim8,dim1),  missval=mv, prec="float")
  var34 <- ncvar_def(name="SCT_DOM",units="category",  dim=list(dim3,dim4,dim1),      missval=mv, prec="float")
  var35 <- ncvar_def(name="SOILCBOT",units="category", dim=list(dim3,dim4,dim8,dim1), missval=mv, prec="float")
  var36 <- ncvar_def(name="SCB_DOM",units="category", dim=list(dim3,dim4,dim1),       missval=mv, prec="float")
  var37 <- ncvar_def(name="ALBEDO12M",units="percent",  dim=list(dim3,dim4,dim9,dim1),missval=mv, prec="float")
  var38 <- ncvar_def(name="GREENFRAC",units="fraction", dim=list(dim3,dim4,dim9,dim1),missval=mv, prec="float")
  var39 <- ncvar_def(name="LAI12M",units="Kelvin",  dim=list(dim3,dim4,dim9,dim1),    missval=mv, prec="float")
  var40 <- ncvar_def(name="SNOALB",units="percent",  dim=list(dim3,dim4,dim1),        missval=mv, prec="float")
  var41 <- ncvar_def(name="SLOPECAT",units="catelogry", dim=list(dim3,dim4,dim1),     missval=mv, prec="float")
  var42 <- ncvar_def(name="CON",units="", dim=list(dim3,dim4,dim1),  missval=mv, prec="float")
  var43 <- ncvar_def(name="VAR",units="m", dim=list(dim3,dim4,dim1), missval=mv, prec="float")
  var44 <- ncvar_def(name="OA1",units="", dim=list(dim3,dim4,dim1),  missval=mv, prec="float")
  var45 <- ncvar_def(name="OA2",units="",  dim=list(dim3,dim4,dim1), missval=mv, prec="float")
  var46 <- ncvar_def(name="OA3",units="",  dim=list(dim3,dim4,dim1), missval=mv, prec="float")
  var47 <- ncvar_def(name="OA4",units="",  dim=list(dim3,dim4,dim1), missval=mv, prec="float")
  var48 <- ncvar_def(name="OL1",units="fraction",  dim=list(dim3,dim4,dim1),  missval=mv, prec="float")
  var49 <- ncvar_def(name="OL2",units="fraction",  dim=list(dim3,dim4,dim1),  missval=mv, prec="float")
  var50 <- ncvar_def(name="OL3",units="fraction",  dim=list(dim3,dim4,dim1),  missval=mv, prec="float")
  var51 <- ncvar_def(name="OL4",units="fraction",  dim=list(dim3,dim4,dim1),  missval=mv, prec="float")
  var52 <- ncvar_def(name="VAR_SSO",units="meters2 MSL", dim=list(dim3,dim4,dim1),   missval=mv, prec="float")
  var53 <- ncvar_def(name="LAKE_DEPTH",units="meters MSL", dim=list(dim3,dim4,dim1), missval=mv, prec="float")
  #==== assign file name ===
  out_filename <- paste("d02.geo_RCEC_",work.year,".nc",sep="") 
  #==== creat the nc file based on definition ==========================  
  names<-list(var0,var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,
	      var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,
	      var21,var22,var23,var24,var25,var26,var27,var28,var29,var30,
	      var31,var32,var33,var34,var35,var36,var37,var38,var39,var40,
	      var41,var42,var43,var44,var45,var46,var47,var48,var49,var50,
	      var51,var52,var53) 
  #setting the input_information from "out_filename"   
  output_nc <- nc_create( out_filename, names) 
  #===== set attributes of variables in nc files for various types ===== 
  print(paste("Writing dataset to new nc file: ", out_filename, sep="")) 
  # load old data keep original dimension  
  org.data <- fun_read_nc(arg1=d02.fname,arg2=FALSE)

  tt <- array(paste(work.year,"-01-01_00:00:00",sep=""), dim=c(1,19)) 
  # put new data into file  
  ncvar_put(output_nc, var0,  tt  )
  ncvar_put(output_nc, var1,   org.data$XLAT_M   )  
  ncvar_put(output_nc, var2,   org.data$XLONG_M  )
  ncvar_put(output_nc, var3,   org.data$XLAT_V   )
  ncvar_put(output_nc, var4,   org.data$XLONG_V  )
  ncvar_put(output_nc, var5,   org.data$XLAT_U   )
  ncvar_put(output_nc, var6,   org.data$XLONG_U  )  
  ncvar_put(output_nc, var7,   org.data$CLAT     )
  ncvar_put(output_nc, var8,   org.data$CLONG    )
  ncvar_put(output_nc, var9,   org.data$MAPFAC_M ) 
  ncvar_put(output_nc, var10,  org.data$MAPFAC_V ) 
  ncvar_put(output_nc, var11,  org.data$MAPFAC_U )
  ncvar_put(output_nc, var12,  org.data$MAPFAC_MX) 
  ncvar_put(output_nc, var13,  org.data$MAPFAC_VX) 
  ncvar_put(output_nc, var14,  org.data$MAPFAC_UX) 
  ncvar_put(output_nc, var15,  org.data$MAPFAC_MY) 
  ncvar_put(output_nc, var16,  org.data$MAPFAC_VY) 
  ncvar_put(output_nc, var17,  org.data$MAPFAC_UY) 
  ncvar_put(output_nc, var18,  org.data$E )
  ncvar_put(output_nc, var19,  org.data$F )
  ncvar_put(output_nc, var20,  org.data$SINALPHA )
  ncvar_put(output_nc, var21,  org.data$COSALPHA )
  ncvar_put(output_nc, var22,  org.data$LANDMASK ) 
  ncvar_put(output_nc, var23,  org.data$XLAT_C   )
  ncvar_put(output_nc, var24,  org.data$XLONG_C  )
  ncvar_put(output_nc, var25,  org.data$SINALPHA_U) 
  ncvar_put(output_nc, var26,  org.data$COSALPHA_U) 
  ncvar_put(output_nc, var27,  org.data$SINALPHA_V)
  ncvar_put(output_nc, var28,  org.data$COSALPHA_V) 
  # replace the land cover datat from reconstruction 
  #org.data$LANDUSEF[,,,1] <- lcfrac
  ncvar_put(output_nc, var29,  lcfrac)
  #org.data$LU_INDEX[,,1] <- lcdom
  ncvar_put(output_nc, var30,  lcdom)
  ncvar_put(output_nc, var31,  org.data$HGT_M    )
  ncvar_put(output_nc, var32,  org.data$SOILTEMP )
  ncvar_put(output_nc, var33,  org.data$SOILCTOP )
  ncvar_put(output_nc, var34,  org.data$SCT_DOM  )
  ncvar_put(output_nc, var35,  org.data$SOILCBOT )
  ncvar_put(output_nc, var36,  org.data$SCB_DOM  )
  ncvar_put(output_nc, var37,  org.data$ALBEDO12M)
  ncvar_put(output_nc, var38,  org.data$GREENFRAC) 
  ncvar_put(output_nc, var39,  org.data$LAI12M )
  ncvar_put(output_nc, var40,  org.data$SNOALB )
  ncvar_put(output_nc, var41,  org.data$SLOPECAT) 
  ncvar_put(output_nc, var42,  org.data$CON) 
  ncvar_put(output_nc, var43,  org.data$VAR) 
  ncvar_put(output_nc, var44,  org.data$OA1) 
  ncvar_put(output_nc, var45,  org.data$OA2) 
  ncvar_put(output_nc, var46,  org.data$OA3) 
  ncvar_put(output_nc, var47,  org.data$OA4) 
  ncvar_put(output_nc, var48,  org.data$OL1) 
  ncvar_put(output_nc, var49,  org.data$OL2) 
  ncvar_put(output_nc, var50,  org.data$OL3)
  ncvar_put(output_nc, var51,  org.data$OL4) 
  ncvar_put(output_nc, var52,  org.data$VAR_SSO   )
  ncvar_put(output_nc, var53,  org.data$LAKE_DEPTH) 
  # add atribute to the variables
  att.fid <- rep(104, 53)
  att.mod <- c("XY ","XY ","XY ","XY ","XY ","XY ","XY ","XY ","XY ","XY ","XY ","XY ","XY ",
	       "XY ","XY ","XY ","XY ","XY ","XY ","XY ","XY ","XY ","XY ",
	       "XY ","XY ","XY ", "XY ","XY ","XYZ","XY ","XY ","XY ","XYZ",
	       "XY ","XYZ","XY ","XYZ","XYZ","XYZ","XY ","XY ","XY ","XY ","XY ","XY ","XY ",
	       "XY ","XY ","XY ","XY ","XY ","XY ","XY ")	
  att.des <- c("Latitude on mass grid","Longitude on mass grid","Latitude on V grid","Longitude on V grid",
	       "Latitude on U grid","Longitude on U grid","Computational latitude on mass grid","Computational longitude on mass grid",
	       "Mapfactor on mass grid","Mapfactor on V grid","Mapfactor on U grid","Mapfactor (x-dir) on mass grid",
	       "Mapfactor (x-dir) on V grid","Mapfactor (x-dir) on U grid","Mapfactor (y-dir) on mass grid",
	       "Mapfactor (y-dir) on V grid","Mapfactor (y-dir) on U grid","Coriolis E parameter","Coriolis F parameter",
	       "Sine of rotation angle","Cosine of rotation angle","Landmask : 1=land, 0=water",
	       "Latitude at grid cell corners","Longitude at grid cell corners",
	       "Sine of rotation angle on U grid","Cosine of rotation angle on U grid",
	       "Sine of rotation angle on V grid","Cosine of rotation angle on V grid",
	       "24-category USGS landuse","Dominant category","GMTED2010 30-arc-second topography height","Annual mean deep soil temperature",
	       "16-category top-layer soil type","Dominant category","16-category bottom-layer soil type",
	       "Dominant category","Monthly surface albedo","MODIS FPAR","MODIS LAI","Maximum snow albedo",
	       "Dominant category","orographic convexity","stdev of subgrid-scale orographic height",
	       "orographic asymmetry","orographic asymmetry","orographic asymmetry","orographic asymmetry","effective orographic length",
	       "effective orographic length","effective orographic length","effective orographic length",
	       "Variance of Subgrid Scale Orography","Topography height")
  att.sta <- c("M","M","V","V","U","U","M","M","M","V","U","M","V","U","M","V","U","M","M","M","M","M",
	       "CORNER","CORNER","U","U","V","V","M","M","M","M","M","M","M","M","M","M","M","M","M","M",
               "M","M","M","M","M","M","M","M","M","M","M")
  att.srx <- rep(1,53);  att.sry <- rep(1,53)
  for (i in 1:53) {
       ncatt_put( nc=output_nc, varid=get(paste("var",i,sep="")), attname="FieldType", attval=att.fid[i], prec=NA)
       ncatt_put( nc=output_nc, varid=get(paste("var",i,sep="")), attname="MemoryOrder", attval=att.mod[i], prec="text")
       ncatt_put( nc=output_nc, varid=get(paste("var",i,sep="")), attname="description", attval=att.des[i], prec="text")
       ncatt_put( nc=output_nc, varid=get(paste("var",i,sep="")), attname="stagger", attval=att.sta[i], prec="text")
       ncatt_put( nc=output_nc, varid=get(paste("var",i,sep="")), attname="sr_x", attval=att.srx[i], prec=NA)
       ncatt_put( nc=output_nc, varid=get(paste("var",i,sep="")), attname="sr_y", attval=att.sry[i], prec=NA)
  }
              
  # put the global attribution to the file
  ncatt_put(output_nc, 0, "File description:", 
	    paste("The land cover over Taiwan was generated by a regional map reconstruction",
		 "made by  (Chen et al. 2019.) https://www.nature.com/articles/s41598-019-40063-1",sep="")
		 )
  # close the nc file 
  nc_close(output_nc) 
# ==== end of create nc ====
} # end of ld.creat.nc

