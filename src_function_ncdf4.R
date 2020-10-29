#==============================================================
#Function: fun_read_nc
#
#First date: 2014-05-18 
#Revised: 2017-03-03
#
#Purpose: get the variable for any nc file 
#         and output a list for every varibles in nc file 
#==============================================================
fun_read_nc<- function(arg1,arg2="TRUE") {
  #load  ncdf library
  #
  library(ncdf4)
  #arg1: filepath fo the nc file from James multilayer output
  print(paste("arg1: for reading file path ;", arg1)) 
  # open the read in file and copy the variables to the dataframe for analysis
  input_nc <- nc_open(arg1)
  #str(input_nc)
  # creat a result list for storage varibales as results list for output
  result<-list() 
  
  for (i in 1:input_nc$ndims ) {
 #   print(i)
    # store each  variable with respect of the dim_name 
    result[[input_nc$dim[[i]]$name]] <- input_nc$dim[[i]]$vals
  }
  
  for (i in 1:length(input_nc$var) ) {
    # store each variable with respect of the var_name 
    print(paste("creat var:",as.character(input_nc$var[[i]]$name),sep=""))	  
    try( result[[as.character(input_nc$var[[i]]$name)]] <- ncvar_get(input_nc, as.character(input_nc$var[[i]]$name), collapse_degen=arg2))
  }
 
 
  nc_close(input_nc)
  # show result structure 
  print(str(result))
  # export the datatable		
  return(result)
} #end of fun_read_nc FUNCTION

#=========================================================
# Creating point forcing file for ORCHIDEE 
#
# Author: Yiying CHEN
#
#First Date: 2014-05-18
#Revised : 2017-03-03
#Purpose:preparing the forcing data for the nc file 
#        & doing the unit convertion 
#=========================================================
fun_make_forcing <- function (input_list1, input_list2) 
{
  print(names(input_list1) )
  print(input_list2)
  #===load ncdf4 library ====
  library(ncdf4)
  #======= Start loading the variables in the list1 ===============
  swdown   <- input_list1[[1]]    #short wave rad.   [W/m2]
  lwdown   <- input_list1[[2]]    #long wave rad.    [W/m2]
  rainfall <- input_list1[[3]]    #rainfall rate     [kg/m2/s]  if (mm/dt)  data*(1/dt) for convertion 
  snowfall <- input_list1[[4]]    #snowfall rate     [kg/m2/s] 
    
  psurf    <- input_list1[[5]]    #air pressure      [Pa]       if (kPa)  data*1000.0
  tair     <- input_list1[[6]]    #air temperature   [K]        if (oC)   data+273.15
  qair     <- input_list1[[7]]    #specific humidity [kg/kg]
  wind     <- input_list1[[8]]    #wind speed        [m/s]
  dt       <- input_list1[[9]]    #delta T           [s]
  #=== calculated the specific humidity from relative humidity and air temperarute 
  #=== saturate water vapor pressure calculation
  #esat <- 611.2*exp((17.67 * tair)/(tair+243.5)) # Bolton (1980) tair in (oC); qstat in(Pa)
  #eair <- esat* (rhair/100.)
  #qair<- ( (0.622*eair) / ((psurf*1000.) - (0.378*eair)) )
  
  #=== unit conservation
  psurf    <- psurf    #   *1000.
  tair     <- tair     #   + 273.15
  rainfall <- rainfall #   * (1./dt) 

  #=== get the tatal datatable lenth 
  n_tsteps       <- length(input_list1[[1]])
  
  #=== inital variables
  steps_since   <- NULL
  seconds_since <- as.double(NA, length=n_tsteps)
  
  #=== create time step indexs with seconds
  for (i in 1:n_tsteps) {
    steps_since[i]   <- 1.0*(i-1)
    seconds_since[i] <- dt*(i-1)  # for half-hourly
  }
  #=== set some constant and genrate time index 
  missing_r = 1.0e+20
  #=== arrange a datatable(forcing) for the forcing data
  forcing <- data.frame (
    swdown=swdown,
    lwdown=lwdown,
    rainfall=rainfall,
    snowfall=snowfall,
    qair=qair,
    psurf=psurf,
    tair=tair,
    wind=wind,
    tsteps_since=steps_since,
    seconds_since=seconds_since
    )	
  
  lon_arr           <- input_list2[[1]]
  lat_arr           <- input_list2[[2]]
  output_fname      <- input_list2[[3]]
  time_start_origin <- input_list2[[4]]
  
  #=== setup a NCDF file "define dimensions and  variables"
  
  #latitude/latitude/timestep 
  dimX <- ncdim_def("lon", "",  1:length(lon_arr), create_dimvar=FALSE )
  dimY <- ncdim_def('lat' ,"",  1:length(lat_arr), create_dimvar=FALSE)
  dimT <- ncdim_def("tstep", "", 1:n_tsteps,  create_dimvar=FALSE)
  
  mv <- 1e+20
  d1 <- ncvar_def( "lon", "degrees_east",  dimX,longname="Longitude", mv, prec="double")
  d2 <- ncvar_def( "lat", "degrees_north", dimY,longname="Latitude", mv, prec="double")
  d3 <- ncvar_def( "tstep",  paste("seconds since ",time_start_origin ,sep=""), dimT, mv, prec="double")
  
  v1 <- ncvar_def( "SWdown", "W/m^2", list(dimX,dimY,dimT),longname="Surface incident SW radiation ", mv, prec="double")
  v2 <- ncvar_def( "LWdown", "W/m^2", list(dimX,dimY,dimT),longname="Surface incident LW radiation ", mv, prec="double")
  v3 <- ncvar_def( "Rainf", "kg/m^2/s", list(dimX,dimY,dimT),longname="Rainfall rate" ,mv, prec="double")
  v4 <- ncvar_def( "Snowf", "kg/m^2/s", list(dimX,dimY,dimT),longname="Snowfall rate" ,mv, prec="double")
  v5 <- ncvar_def( "Qair", "kg/kg", list(dimX,dimY,dimT),longname="Specific humidity at 2m" ,mv, prec="double")
  v6 <- ncvar_def( "PSurf", "Pa", list(dimX,dimY,dimT),longname="Surfce pressure at 2m" ,mv, prec="double")
  v7 <- ncvar_def( "Tair", "K", list(dimX,dimY,dimT),longname="Temperature at 2m" ,mv, prec="double")
  v8 <- ncvar_def( "Wind", "m/s", list(dimX,dimY,dimT),longname="Surf wind speed at 10m" ,mv, prec="double")
  
  #====
  #creat the nc file based on definition 
  names<-list(d1,d2,d3,v1,v2,v3,v4,v5,v6,v7,v8)
  new_nc <- nc_create(paste(output_fname, sep=""), names)
  #====
  
  #set attributes of variables in nc files for various types
  #longitude/latitude/time
  #att.put.ncdf( new_nc, d3, "calendar","noleap" , prec="text" )
  ncatt_put( new_nc, d3, "calendar","gregorian" , prec="text" )
  ncatt_put( new_nc, d3, "origin",time_start_origin , prec="text" )
  #land sea mask
  
  #short wave radiation
  ncatt_put( new_nc, v1, "axis", "TYX" , prec="text" )
  ncatt_put( new_nc, v1, "associate", 'time (nav_lat nav_lon)'   , prec="text" )
  #long wave radiation
  ncatt_put( new_nc, v2, "axis", "TYX" , prec="text" )
  ncatt_put( new_nc, v2, "associate", 'time (nav_lat nav_lon)'   , prec="text" )
  #rainfall
  ncatt_put( new_nc, v3, "axis", "TYX" , prec="text" )
  ncatt_put( new_nc, v3, "associate", 'time (nav_lat nav_lon)'   , prec="text" )
  #snowfall
  ncatt_put( new_nc, v4, "axis", "TYX" , prec="text" )
  ncatt_put( new_nc, v4, "associate", 'time (nav_lat nav_lon)'   , prec="text" )
  #specific humidity
  ncatt_put( new_nc, v5, "axis", "TYX" , prec="text" )
  ncatt_put( new_nc, v5, "associate", 'time (nav_lat nav_lon)'   , prec="text" )
  #air pressure
  ncatt_put( new_nc, v6, "axis", "TYX" , prec="text" )
  ncatt_put( new_nc, v6, "associate", 'time (nav_lat nav_lon)'   , prec="text" )
  #air temperature
  ncatt_put( new_nc, v7, "axis", "TYX" , prec="text" )
  ncatt_put( new_nc, v7, "associate", 'time (nav_lat nav_lon)'   , prec="text" )
  #wind 
  ncatt_put( new_nc, v8, "axis", "TYX" , prec="text" )
  ncatt_put( new_nc, v8, "associate", 'time (nav_lat nav_lon)'   , prec="text" )
  
  #put coordinate & time  
  ncvar_put(new_nc, d1,lon_arr)
  ncvar_put(new_nc, d2,lat_arr)
  ncvar_put(new_nc, d3, forcing$seconds_since)
  
  #put forcing data
  ncvar_put(new_nc, v1, forcing$swdown   )
  ncvar_put(new_nc, v2, forcing$lwdown   )
  ncvar_put(new_nc, v3, forcing$rainfall )
  ncvar_put(new_nc, v4, forcing$snowfall )
  ncvar_put(new_nc, v5, forcing$qair     )
  ncvar_put(new_nc, v6, forcing$psurf    )
  ncvar_put(new_nc, v7, forcing$tair     )
  ncvar_put(new_nc, v8, forcing$wind     )

  #close the nc file
  nc_close(new_nc)
  print(paste("fun_make_focring() creates the forcing file for ORCHIDEE: ", input_list2[[3]],sep='')) 
} #End Fun_make_forcing
