##  Version 3.0
##  Author: Jasper Van doninck
##  Contact: jasper.vandoninck@utu.fi, vandoninck.jasper@gmail.com
##
##  Description:
##    Process Landsat Collection 1 Level 1 image to surface reflectance
##    Includes BRDF+topographic correction
##    WARNING: required ESPA software to be installed!
##      For instructions:
##        * ESPA product formatter: https://github.com/USGS-EROS/espa-product-formatter
##        * LEDAPS:                 https://github.com/USGS-EROS/espa-surface-reflectance/tree/master/ledaps
##
##  Usage:
##    surfaceReflectance(fid,inDir,outDir,tempDirRoot,
##                       removeTmpFiles=TRUE, ptstart=proc.time(), 
##                       bs=100, 
##                       brdf=NULL, outSZA=30,
##                       topo=NULL, demFile=NULL)
##
##  Arguments:
##    fid:            Character. Landsat file ID
##    inDir:          Character. Directory of Landsat L1 .tar.gz file.
##    outDir:         Character. Directory of output L2 files.
##    tempDirRoot:    Character. Path to directory in which temporary directory will be created.
##    removeTmpFiles: Logical. Defines whether temporary files are to be removed.
##    ptstart:        Object of class "proc.time". Used to track progress.
##    bs:             Numeric. Defines number of image lines read into memory simultaneously
##    brdf:           NULL or character (one of "Roy"/"Flood"/"UTU"). See details.
##    outSZA:         Numeric. Solar zenith angle (in degrees) of BRDF-normalized image. Only required if argument "brdf" is not NULL.
##    topo:           Null or character ("C"). See details.
##    demFile:        Character. Filename of digital elevation model covering extent of Landsat image. Only required if parameter "topo" is not NULL.
##
##  Details:
##    Parameter "brdf" defines the parameters of the BRDF correction to be applied. 
##    Possible values: 
##      "Roy": set of global parameters obtained from MODIS  (Roy, D. P., H. K. Zhang, J. Ju, J. L. Gomez-Dans, P. E. Lewis, C. B. Schaaf, Q. Sun, J. Li, H. Huang, and V. Kovalskyy. 2016. "A General Method to Normalize Landsat Reflectance Data to Nadir BRDF Adjusted Reflectance." Remote Sensing of Environment 176 (April): 255-71. https://doi.org/10.1016/j.rse.2016.01.023.)
##      "Flood": set of parameters derived for Australia (Flood, N., T. Danaher, T. Gill, and S. Gillingham. 2013. "An Operational Scheme for Deriving Standardised Surface Reflectance from Landsat TM/ETM+ and SPOT HRG Imagery for Eastern Australia." Remote Sensing 5 (1): 83-109. https://doi.org/10.3390/rs5010083.)
##      "UTU": set of parameters derived for Amazonian forests (Van doninck, J., and H. Tuomisto. 2017. "Evaluation of Directional Normalization Methods for Landsat TM/ETM+ over Primary Amazonian Lowland Forests." International Journal of Applied Earth Observation and Geoinformation 58 (June): 249-63. https://doi.org/10.1016/j.jag.2017.01.017.)
##
##    Parameter "topo" defines topographic correction. Only C-correction is implemented.
##
##  Value:
##    Numeric. 
##      0:  Processing completed correctly.
##      1:  Error in unpacking .tar.gz file or while running ESPA software.
##      3:  Processing skipped. Insufficient cloud-free pixels in image.
##      5:  Processing skipped. Incorrect input argument(s).
##
##  Used packages: 
##    raster
##    rgdal
##    RStoolbox
#install.packages(c("raster", "rgdal", "RStoolbox")) 

surfaceReflectance <- function(fid,inDir,outDir,tempDirRoot,
                               removeTmpFiles=TRUE, 
                               ptstart=proc.time(), 
                               bs=100,
                               brdf=NULL,
                               outSZA=30,
                               topo=NULL,
                               demFile=NULL){
  
  ## Load libraries and functions
  require(raster)
  require(rgdal)

  #Kvol - RossThick kernel for volume scattering (Roujean et al., 1992)
  Kvol <- function(sunza, satza, relaz, degrees=FALSE){
    if (degrees){
      sunza <- sunza*pi/180
      satza <- satza*pi/180
      relaz <- relaz*pi/180
    }
    ph_an <- acos((cos(sunza)*cos(satza))+(sin(sunza)*sin(satza)*cos(relaz)))
    result <- (((pi/2-ph_an)*cos(ph_an)+sin(ph_an))/(cos(sunza)+cos(satza)))-(pi/4)
    return(result)
  }
  
  #Kgeo - LiSparse-R kernel for geometric-optical surface scattering (Wanner et al., 1995)
  Kgeo <- function(sunza, satza, relaz, degrees=FALSE){
    if (degrees){
      sunza <- sunza*pi/180
      satza <- satza*pi/180
      relaz <- relaz*pi/180
    }
    #br <- 1
    hb <- 2
    ph_an <- acos((cos(sunza)*cos(satza))+(sin(sunza)*sin(satza)*cos(relaz)))
    D <- sqrt((tan(sunza))^2+(tan(satza))^2-2*tan(sunza)*tan(satza)*cos(relaz))
    cost <- (hb*(sqrt(D^2+(tan(sunza)*tan(satza)*sin(relaz))^2)/((1/cos(sunza))+(1/cos(satza)))))
    cost[cost > 1] <- 1
    cost[cost < -1] <- -1
    t <- acos(cost)
    O=(1/pi)*(t-sin(t)*cos(t))*((1/cos(sunza))+(1/cos(satza)))
    result <- O-(1/cos(sunza))-(1/cos(satza))+0.5*(1+cos(ph_an))*(1/cos(sunza))*(1/cos(satza))
    return(result)
  }
  
  #Some other helper functions
  nearest.vec <- function(x, vec){
    smallCandidate <- findInterval(x, vec, all.inside=TRUE)
    largeCandidate <- smallCandidate + 1
    #nudge is TRUE if large candidate is nearer, FALSE otherwise
    nudge <- 2 * x > vec[smallCandidate] + vec[largeCandidate]
    return(smallCandidate + nudge)
  }
  timeStamp <- function(pt.start){
    pt <- proc.time()-pt.start
    pt <- pt["elapsed"]/(60*60)
    hh <- floor(pt)
    pt <- (pt-hh)*60
    mm <- floor(pt)
    ss <- floor((pt-mm)*60)
    if (hh<10) hhtxt <- paste0("0",as.character(hh)) else hhtxt <- as.character(hh)
    if (mm<10) mmtxt <- paste0("0",as.character(mm)) else mmtxt <- as.character(mm)  
    if (ss<10) sstxt <- paste0("0",as.character(ss)) else sstxt <- as.character(ss)
    return(paste0("[",hhtxt,":",mmtxt,":",sstxt,"]"))
  }
  
  ##  Define constants
  # Sensor
  sensor <- substr(fid,1,4)
  # Landsat bands
  sbands <- c(1,2,3,4,5,7)

  # Normalized angles
  if(!is.null(brdf)){
    Kgeo.norm <- Kgeo(outSZA, 0,0, degrees=TRUE)
    Kvol.norm <- Kvol(outSZA, 0,0, degrees=TRUE)
  }

  ##  Check input parameters
  if (!is.null(topo)){
    # Check if parameter "topo" is accepted value
    if (!isTRUE(topo=="C")){
      cat("[ERROR] Parameter \"topo\" must one of (NULL, \"C\").\n", sep="")
      return(5)
    }
    # Check if DEM raster exists (in case topographic correction performed)
    if(!is.character(demFile)){
      cat("[ERROR] Parameter \"demFile\" is not a character variable.\n", sep="")
      return(5)
    }
    if(isFALSE(file.exists(demFile))){
      cat("[ERROR] File \"",demFile ,"\" does not exist.\n", sep="")
      return(5)
    }
  }
    # Check BRDF model parameters 
  fpData <- switch(tolower(brdf),
                   "roy"   = c(0.0372, 0.0580, 0.0574, 0.1535, 0.1154, 0.0639,
                               0.0079, 0.0178, 0.0227, 0.0330, 0.0453, 0.0387)/
                                  c(0.0774, 0.1306, 0.1690, 0.3093, 0.3430, 0.2658),
                   "flood" = c(0.93125413991, 0.687401438519, 0.645033011917, 0.704036740665, 0.360201003097, 0.290061903555,
                               0.260953557124, 0.213872135374, 0.180032152925, 0.093518142066, 0.162796996525, 0.147723009593),  
                   "utu"   = c(1.2933, 0.5550, 0.9937, 0.3489, 1.0129, 1.2742,
                               0.1902, 0.2456, 0.0000, 0.2755, 0.2348, 0.2426),
                    NA)
  if(isTRUE(sum(is.na(fpData)))){
    cat("[ERROR] Parameter \"brdf\" must be one of (\"Roy\", \"Flood\", \"UTU\").\n", sep="")
    return(5)
  } else fp <- matrix(data= fpData, ncol=2, dimnames=list(sbands, c("vol", "geo")))
 
  # Create/set directories
  tempDir <- paste0(tempDirRoot,"/",fid)
  if (!file.exists(tempDir)) dir.create(tempDir, recursive=TRUE)
  if (!file.exists(outDir)) dir.create(outDir, recursive=TRUE)
  wdir <- tempDir
  cwdir <- getwd()
  setwd(wdir)
  
  # Set temporal directory for "raster" package 
  rasterTmpDir <- file.path(wdir,"rasterTmpDir")
  if (!file.exists(rasterTmpDir)) dir.create(rasterTmpDir)
  options(rasterTmpDir=rasterTmpDir)
  
  # Check if file is available untar .tar.gz file
  cat(timeStamp(ptstart)," ",fid,": Untar/LEDAPS/CFMASK/angles\n", sep="")
  if (!file.exists(file.path(inDir,paste0(fid,".tar.gz")))){
    cat("[ERROR] ",fid,": File not found\n", sep="")
    setwd(cwdir)
    unlink(wdir, recursive=TRUE)
    return(1)
  }
  untar(file.path(inDir,paste0(fid,".tar.gz")), exdir=wdir)
  #Convert LPGS to ESPA format:
  if (system(paste0("convert_lpgs_to_espa --mtl ",fid,"_MTL.txt --xml ",fid,".xml --del_src_files"), intern=FALSE, ignore.stdout=TRUE)){
    cat("[ERROR] ",fid,": Error converting LPGS to ESPA\n", sep="")
    setwd(cwdir)
    unlink(wdir, recursive=TRUE)
    return(1)
  }
  # Apply LEDAPS
  if (system(paste0("do_ledaps.py --xml ",fid,".xml"), intern=FALSE, ignore.stdout=TRUE, ignore.stderr=TRUE)){
    ledaps.fail <- 1
    nTry <- 2
    while (ledaps.fail & nTry <5){
      cat("[WARNING] ",fid,": Error running LEDAPS, attempt ",nTry,"\n", sep="")
      ledaps.fail <- system(paste0("do_ledaps.py --xml ",fid,".xml"), intern=FALSE, ignore.stdout=TRUE, ignore.stderr=TRUE)
      nTry <- nTry+1
    }
    if (ledaps.fail){
      cat("[ERROR] ",fid,": Error running LEDAPS\n", sep="")
      setwd(cwdir)
      unlink(wdir, recursive=TRUE)
      return(1)
    }
  }
  # CFMASK
  if (system(paste0("cfmask --xml ",fid,".xml"), intern=FALSE, ignore.stdout=TRUE, ignore.stderr=TRUE)){
    cat("[ERROR] ",fid,": Error running CFMASK\n", sep="")
    setwd(cwdir)
    unlink(wdir, recursive=TRUE)
    return(1)
  }
  ##  Landsat angles   --> ONLY IF !is.null(brdf)
  if (system(paste0("landsat_angles ",fid,"_ANG.txt"), intern=FALSE, ignore.stdout=TRUE, ignore.stderr=TRUE)){
    cat("[ERROR] ",fid,": Error running Landsat angles\n", sep="")
    setwd(cwdir)
    unlink(wdir, recursive=TRUE)
    return(1)
  }
  
  # Clean up intermediate files (DN + TOA)
  if(removeTmpFiles) dummy <- file.remove(list.files(wdir, pattern=paste0(fid,"_B"), full.names=TRUE))
  if(removeTmpFiles) dummy <- file.remove(list.files(wdir, pattern=paste0(fid,"_toa"), full.names=TRUE))
  
  # Create empty raster template
  nullR <- raster(file.path(wdir, paste0(fid, "_sr_fill_qa.img")), vals=FALSE)
  nullR <- raster(nrows=nrow(nullR), ncols=ncol(nullR), ext=extent(nullR), crs=crs(nullR), resolution=res(nullR))
  nullR_ext <- extent(nullR)
  nullR_ext_p <- extent(t(project(t(as.matrix(nullR_ext)), projection(nullR), inv=TRUE)))
  
  # Create masks from QA bands and check number of unmasked pixels
  cat(timeStamp(ptstart)," ",fid,": Create mask\n",sep="")
  maskR<- overlay(raster(file.path(wdir, paste0(fid, "_sr_fill_qa.img"))),
                  raster(file.path(wdir, paste0(fid, "_sr_cloud_qa.img"))),
                  raster(file.path(wdir, paste0(fid, "_sr_cloud_shadow_qa.img"))),
                  raster(file.path(wdir, paste0(fid, "_sr_adjacent_cloud_qa.img"))),
                  raster(file.path(wdir, paste0(fid, "_sr_band1.img"))),
                  raster(file.path(wdir, paste0(fid, "_sr_band2.img"))),
                  raster(file.path(wdir, paste0(fid, "_sr_band3.img"))),
                  raster(file.path(wdir, paste0(fid, "_sr_band4.img"))),
                  raster(file.path(wdir, paste0(fid, "_sr_band5.img"))),
                  raster(file.path(wdir, paste0(fid, "_sr_band7.img"))),
                  raster(file.path(wdir, paste0(fid, "_cfmask.img"))),
                  fun=function(fi,cl,cs,ac,a,b,c,d,e,f,g)
                  {return(fi==0 & cl==0 & cs==0 & ac==0 & a!=2000 & b!=2000 & c!=2000 & d!=2000 & e!=2000 & f!=2000 & (!is.na(g) & g<2) )},
                  filename=file.path(wdir,paste0(fid,"_mask.grd")), overwrite=TRUE)
  if (cellStats(maskR,stat='sum')<10000){ # Exit if too few unmasked pixels
    cat(timeStamp(ptstart)," ",fid,": <10000 pixels, abort\n",sep="")
    setwd(cwdir)
    unlink(wdir, recursive=TRUE)
    return(2)
  }
  
  ##  BRDF correction
  if(is.null(brdf)){
    #No BRDF correction
    rho_brdf <- stack(as.list(file.path(tempDir, paste0(fid,"_sr_band",sbands,".img"))))
    rho_brdf <- mask(rho_brdf,maskR, maskvalue=0,updatevalue=NA, 
                     file.path(rasterTmpDir,paste0(fid,"_rho_brdf.grd")), overwrite=TRUE, datatype="INT2S")
  } else {
    #Perform BRDF correction
    cat(timeStamp(ptstart)," ",fid,": BRDF normalisation.\n",sep="")

    rho_brdf <- brick(nullR, nl=length(sbands))
    rho_brdf <- writeStart(rho_brdf, file.path(rasterTmpDir,paste0(fid,"_rho_brdf.grd")), overwrite=TRUE, datatype="INT2S")
    for(row_in in seq(1,nrow(nullR),bs)){
      #Read row(s) from rasters
      sun.zen_row <- getValues(raster(file.path(wdir,"angle_solar_B04.img"),band=1),row_in,bs)*0.01*pi/180
      sun.azi_row <- getValues(raster(file.path(wdir,"angle_solar_B04.img"),band=2),row_in,bs)*0.01*pi/180
      
      sat.zen_row <- getValues(raster(file.path(wdir,"angle_sensor_B04.img"),band=1),row_in,bs)*0.01*pi/180
      sat.azi_row <- getValues(raster(file.path(wdir,"angle_sensor_B04.img"),band=2),row_in,bs)*0.01*pi/180
      
      Kgeo.init.h_row <-Kgeo(sun.zen_row, sat.zen_row, sun.azi_row-sat.azi_row)
      Kvol.init.h_row <-Kvol(sun.zen_row, sat.zen_row, sun.azi_row-sat.azi_row)

      mask_row <- getValues(maskR,row_in,bs)

      rho_std_row <- matrix(data=NA, ncol=length(sbands), nrow=length(mask_row))
      for (b in 1:6){
        rho_dir <- getValues(raster(file.path(wdir, paste0(fid,"_sr_band",sbands[b],".img"))),row_in,bs)
        gamma <- (1+fp[b,"vol"]*Kvol.norm+fp[b,"geo"]*Kgeo.norm)/(1+fp[b,"vol"]*Kvol.init.h_row+fp[b,"geo"]*Kgeo.init.h_row)
        rho_norm <- round(gamma*rho_dir)
        
        rho_norm[mask_row==0] <- NA
        rho_norm[rho_norm<0] <- 0
        rho_norm[rho_norm>10000] <- 10000
        
        rho_std_row[,b] <- rho_norm
      }
      rho_brdf <- writeValues(rho_brdf,rho_std_row,row_in)
    }
    rho_brdf <- writeStop(rho_brdf)

  }
  
  ##  Topographic correction
  if(is.null(topo)){
    rho_topo <- rho_brdf
  } else {
    #Combined BRDF/topograpgic correction is not implemented in this version, 
    #use "topCor" in RStoolbox package for now, only "C" implemented. 
    cat(timeStamp(ptstart)," ",fid,": Topographic correction\n",sep="") 
    library(RStoolbox)
    ##  Topographic correction (to be move )
    dem.crop <- crop(raster(demFile), nullR_ext_p)
    dem.elev <- projectRaster(dem.crop, nullR, method="bilinear",
                              filename=file.path(rasterTmpDir,paste0(fid,"_dem_elev.grd")), overwrite=TRUE)
    rho_topo <- topCor(rho_brdf, dem.elev, file.path(tempDir,paste0(fid,"_MTL.txt")), method=topo)
  }
  
  rho_out <- writeRaster(rho_topo, file.path(outDir, paste0(fid,".tif")), overwrite=TRUE, datatype="INT2S")
  
  #Clean up
  cat(timeStamp(ptstart)," ",fid,": Clean up \n", sep="")
  setwd(cwdir)
  options(rasterTmpDir=tempdir())
  if(removeTmpFiles) unlink(wdir, recursive=TRUE)
  rm(list=ls())
  dummy <- gc()
  return(0)
  
}

