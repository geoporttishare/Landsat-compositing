##  Author:   Jasper Van doninck
##  Contact:  jasper.vandoninck(at)utu.fi, vandoninck.jasper(at)gmail.com
##
##  Description:
##    Reproject set of Landat L2 surface reflectance images to common extent and composite using predefined criterion 
##
##  Usage:
##    OGIIR_composite(inBricks, inDates, outRef, filename, tempDirRoot=tempdir(), int="bilinear", comp="medoid")
##
##  Arguments:
##    inBricks:     List of RasterBrick objects. Each list element represent a multilayer Landsat TM/ETM+ image
##    inDates:      Numeric or character vector. Each element corresponds to the acquisition date of the respective images is "inBricks". Typically in format YYYYDDD or YYYYMMDD, but can be free-form. If unknown, can be set to 1:length(inBricks).
##    outRef:       Raster* object. Reference raster indicating extent/resolution of output composite image
##    filename:     Character. Filename of the output Raster*
##    tempRootDir:  Character. Path to directory in which temporary directory will be created.
##    int:          Character. One of c("bilinear","ngb") for bilinear and nearest neigbour interpolation, respectively
##    comp:         Character. Indicating compositing criterion. One of c("medoid", "median", "maxNDIV) for multidimensional median, band-wise median of maximum NDVI compositing, respectively.
##
##  Value:
##    Raster* object
##

OGIIR_composite <- function(inBricks, inDates, outRef, filename, 
                            tempDirRoot=tempdir(),
                            int="bilinear",
                            comp="medoid"){
  
  library(raster)
  require(abind)
  require(cluster)
  
    # #For parallel processing 
    #require(parallel)
    #cl <- makeCluster(2) #Set number of cores here. For details, see "parallel" package
    #clusterExport(cl, "pam")
  
  ##  Create tempdir
  tmpDir <- tempfile("ogiir_comp_", tmpdir=tempDirRoot)
  dir.create(tmpDir, recursive=TRUE)
  
  ##  Step 1: Reproject all images to common extent
  projFun <- function(n, inBricks, outRef, tmpDir, int){
    require(raster)
    tryVal <- try(projBrick <- projectRaster(inBricks[[n]], outRef, method=int,
                                             filename=file.path(tmpDir, paste0("proj_",n,".grd")), overwrite=TRUE), 
                  silent=TRUE)
    if (class(tryVal)=="try-error") return(1) else return(0)
  }
  projFailed <- sapply(1:length(inBricks), projFun, inBricks, outRef, tmpDir, int)
    # #For parallel processing
    #projFailed <- parSapply(cl, 1:length(inBricks), projFun, inBricks, outRef, tmpDir, int)
    
  #(Do something with rasters failed to reproject)
  inDates <- inDates[!projFailed]

  ##  Step 2: pixel-based compositing
    #Define compositing function
  compFun <- switch(tolower(comp),
                    "medoid" = function(x,d){  #Function for medoid (multidimensional median) compositing
                        #Remove duplicated observations (in overlap area scenes same track)
                      dup <- duplicated(d)
                      x <- x[!dup,,drop=FALSE]
                      d <- d[!dup]
                        #Check number of unmasked observation
                      valid <- !is.na(rowSums(x))
                      if (sum(valid) == 0) return(list(r=numeric(6+1)+NA, n=0, d=NA)) #No unmasked observations, return NAs
                      if (sum(valid) == 1) return(list(r=x[valid,,drop=FALSE],n=0,d=d[valid])) #Single unmasked observation, return that observation
                      x <- x[valid,,drop=FALSE]
                      d <- d[valid]
                      if (sum(valid) == 2) {
                        #Two unmasked observations. No median compositing possible, maxNDVI backup
                        ndviMax <- which.max((x[,4]-x[,3])/(x[,4]+x[,3]))
                        return(list(r=x[ndviMax,], n=2, d=d[ndviMax]))
                      }
                        #Three or more unmasked observation
                      med.obj <- pam(x, k=1)
                      return(list(r=x[med.obj$id.med,],n=sum(valid),d=d[med.obj$id.med]))
                    },
                    "median" = function(x,d){  #Function for band-wise median compositing
                      
                      dup <- duplicated(d)
                      x <- x[!dup,,drop=FALSE]
                      valid <- !is.na(rowSums(x))
                      
                      if (sum(valid) == 0) return(list(r=numeric(6+1)+NA, n=0, d=NA)) #No unmasked observations, return NAs
                      return(list(r=apply(x[valid,,drop=FALSE],2,quantile,probs=0.5), n=sum(valid),d=NA))
                    },
                    "maxndvi"= function(x,d){  #Function for maximum NDVI compositing
                      valid <- !is.na(rowSums(x))
                      if (sum(valid) == 0) return(list(r=numeric(6+1)+NA, n=0, d=NA)) #No unmasked observations, return NAs
                      x <- x[valid,,drop=FALSE]
                      d <- d[valid]
                      ndviMax <- which.max((x[,4]-x[,3])/(x[,4]+x[,3]))
                      return(list(r=x[ndviMax,], n=2, d=d[ndviMax]))
                    }  
                    )
  
  fl.proj <- list.files(tmpDir, full.names=TRUE)   #Add error message here if length(fl.proj)==0
  brickList <- lapply(fl.proj, brick)
  
  #Surface reflectance
  if(!dir.exists(dirname(filename))) dir.create(dirname(filename), recursive=TRUE)
  out.comp <- brick(outRef, nl=6)
  out.comp <- writeStart(out.comp, filename=filename, overwrite=TRUE)
  
  bs <- 50
  ##  Pixel-based compositing will be done by blocks of image rows, to avoid that all data has to be read into memory simultaneously.
  ##  Parameter "bs" defines this block size. setting this to a larger number will reduce the number of read/write commands and speed up processing time, provided that enough memory is available.
  for (startRow in seq(1,nrow(outRef),bs)){
    rowValues <- abind(lapply(brickList, getValues, row=startRow, nrows=bs), along=3)
    npix <- nrow(rowValues)
    
    list.img.band <- lapply(seq_len(npix), function(i) t(unname(rowValues[i,,])))
    compReturn <- lapply(list.img.band, compFun, inDates)
      #Parallel
      #compReturn <- parLapply(cl, list.img.band, compFun, imgDates)
    
    compReturn <- parLapply(cl, list.img.band, compFun, inDates)
    outValues <- do.call(rbind,lapply(compReturn, function(x){x$r}))
    out.comp <- writeValues(out.comp, outValues, startRow)
  }
  
  # Close output raster
  out.comp <- writeStop(out.comp)
  
  ##  Clean up tempdir
  unlink(tmpDir, recursive=TRUE)

    #For parallel processing: stop cluster
    #stopCluster(cl)
}


