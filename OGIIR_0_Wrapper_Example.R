##  Author: Jasper Van doninck
##
##  Description: This is an example of a simple wrapper for browsing, processing and compositing Landsat imagery.


##  Set EarthExplorer username and password
  #username
userName <- "my.EE.username"
  #password
require(getPass)
pwd <- getPass()

##  Set search parameters
  #Min/max longitude and latitude
xMin <- 22
xMax <- 22.5
yMin <- 60
yMax <- 60.5
  #Start and end date of search window
startDate <- 20110501
endDate <- 20110701
  #Maximum scene cloud cover
cloudCover <- 60

##  Set directories
  #dowloaded Level-1 images
l1Dir <- "C/OGIIR_example/L1"
  #topographically and BRDF normalized surface reflectance images
srDir <- "C:/OGIIR_example/SR"
  #output composite images
compDir <- "C:/OGIIR_example/comp"
  #directory in which temporary directories will be created
tempDirRoot <- "C:/OGIIR_example"

# Step 1: Query and download all available Landsat images
source(file.path("directory_of_R_files", "OGIIR_1_browseEE.R"))
fl <- browseEE(userName, pwd, xMin, xMax, yMin, yMax, 
               searchFrom=startDate, searchTo=endDate,
               landCC=cloudCover,
               download=l1Dir)

# Step 2: Process all downloaded images to surface reflectance, apply mask, normalize for BRDF effects (skipped topographic correction)
source(file.path("directory_of_R_files", "OGIIR_2_srTopoBrdf.R"))
sapply(fl, surfaceReflectance,
       inDir=l1Dir,outDir=srDir,tempDirRoot=tempDirRoot,
       brdf="Roy", outSZA=30,
       topo=NULL, demFile=NULL)

#(Steps 1-2 could be skipped when downloading L2 data directly, but this would exclude )

# Step 3: Reprojection of images to common extent and pixel-based compositing
source(file.path("directory_of_R_files", "OGIIR_3_composite.R"))
  #File IDs of surface reflectance images
fIDs <- substr(list.files(srDir, pattern="_sr_band1.tif"), 1,48)
  #Dates of image acquisition
inDates <- substr(list.files(srDir, pattern="_sr_band1.tif"), 18,25)
  #Reference output extent-dimensions
library(raster)
outRef <- raster(xmn=xMin, xmx=xMax, ymn=yMin, ymx=yMax, crs=CRS("+proj=longlat +datum=WGS84"), nrows=1800, ncols=1800)
  #Load all Landat images into list of raster bricks
srBricks <- sapply(fIDs, function(x){stack(as.list(list.files(srDir, pattern=x, full.names=TRUE)))})
  #Apply bilinear interpolation and multidimensional median compositing
LandsatComposite <- OGIIR_composite(srBricks, inDates, outRef, file.path(compDir,"LandsatComposite.tif"), 
                            tempDirRoot=tempDirRoot, int="bilinear", comp="medoid")
