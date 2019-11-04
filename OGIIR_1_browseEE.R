##  Version 2.2
##  Author:   Jasper Van doninck
##  Contact:  jasper.vandoninck(at)utu.fi, vandoninck.jasper(at)gmail.com
##
##  Description:
##    JSON search and download Landat TM/ETM+ images from USGS EarthExplorer
##      Provide search bounding box, dates, maximum cloud cover and output directory
##      Downloads files to output directory and returns IDs of images matching the criteria
##      WARNING: Requires USGS machine-to-machine access permission.
##
##  Usage:
##    browseEE(userName, pwd, xMin, xMax, yMin, yMax, searchFrom, searchTo,
##             searchMonths=NULL,
##             dataSets=c("LANDSAT_TM_C1", "LANDSAT_ETM_C1"),
##             maxResults=5000,
##             collectionCategory="T1", dataTypeLevel1="L1TP",
##             landCC=NULL, sceneCC=NULL,
##             download=NULL,
##             forceOverwrite=FALSE, setVerbose=TRUE)
##
##  Arguments:
##    userName:           Character. USGS EarthExplorer user name
##    pwd:                Character. USGS EarthExplorer password
##    xMin:               Numeric. Minimum longitude search bounding box
##    xMax:               Numeric. Maximum longitude search bounding box
##    yMin:               Numeric. Minimum latitude search bounding box
##    yMax:               Numeric. Maximum latitude search bounding box
##    searchFrom:         Numeric or character. Search start date. Format YYYYMMDD.
##    searchTo:           Numeric or character. Search end date. Format YYYYMMDD.
##    searchMonths:       Integer vector or NULL. Specify months (1-12) to be retrieved.
##    dataSets            Character. Datasets to be searched. Only "Landsat_TM_C1" and "LANDSAT_ETM_C1" implemented.
##    maxResults:         Integer. Maximum number of search results returned (maximum 5000)
##    collectionCategory: Character. Collection category to be searched. Only "T1" (Tier 1) implemented.
##    dataTypeLevel1:     Character. Data type level to be searched. Only "L1TP" implemented.
##    landCC:             Numeric. Maximum land cloud cover.
##    sceneCC:            Numeric. Maximum scene cloud cover. 
##    download:           Character or NULL. If character, specify directory in which Landsat image are downloaded. If NULL, no images are downloaded.
##    forceOverwrite:     Logical. Specify whether download has to overwrite or skip existing files.
##    setVerbose:         Logical. 
##
##  Value:
##    Character vector. Vecotr of Landsat scene IDs
##

browseEE <- function(userName, pwd, xMin, xMax, yMin, yMax, searchFrom, searchTo,
                        searchMonths=NULL,
                        dataSets=c("LANDSAT_TM_C1", "LANDSAT_ETM_C1"),
                        maxResults=5000,
                        collectionCategory="T1",
                        dataTypeLevel1="L1TP",
                        landCC=NULL,
                        sceneCC=NULL,
                        download=NULL,
                        forceOverwrite=FALSE,
                        setVerbose=TRUE){
  
  if (setVerbose){
    cat("*******************************\n")
    cat("**  EarthExplorer Inventory  **\n")
    cat("*******************************\n")
  }
  ##  Load libraries
  require(httr)
  require(jsonlite)
  
  serviceURL <- 'https://earthexplorer.usgs.gov/inventory/json/v/1.3.0/'
  
  ##  Request API key
  if (setVerbose) cat("USGS login \n")
  req <- POST(paste0(serviceURL,'login'), 
              body = URLencode(paste0('jsonRequest={"username":"', userName, '","password":"', pwd, '","authType":"EROS","catalogId":"EE"}')),
              content_type("application/x-www-form-urlencoded; charset=UTF-8"))
  stop_for_status(req, "connect to server.")
  warn_for_status(req)
  responseContent <- content(req)
  if (!is.null(responseContent$errorCode)){
    if (setVerbose) cat(paste0("Authentication error: ",responseContent$error,"\n"))
    return(NULL)
  }
  apiKey <- responseContent$data
  
  ##  Query field ID's
  fieldIDs <- matrix(data=c(19873, 19884,
                            19879, 19887,
                            19881, 19893,
                            19874, 19883,
                            19880, 19890,
                            19876, 19886),
                     ncol=2, byrow=TRUE)
  colnames(fieldIDs) <- c("LANDSAT_TM_C1", "LANDSAT_ETM_C1")
  rownames(fieldIDs) <- c("wrsPath", "wrsRow", "landCC", "sceneCC", "collectionCategory", "dataTypeLevel1")
  
  ##  Format search start and search end string
  startDate <- paste0(substr(searchFrom,1,4),"-",substr(searchFrom,5,6),"-",substr(searchFrom,7,8),"T00:00:00")
  endDate <- paste0(substr(searchTo,1,4),"-",substr(searchTo,5,6),"-",substr(searchTo,7,8),"T23:59:59")
  
  ##  Format seach monhts string
  if (!is.null(searchMonths)) searchMonths <- paste0('"months":','[',paste(searchMonths, collapse=","),'],')
  
  allProductIDs <- NULL
  for (dataSet in dataSets){
    if (setVerbose) cat(paste0("Searching ",dataSet),"\n")
    
    ##  Format additional criteria string
    if (is.null(collectionCategory) & is.null(dataTypeLevel1) & is.null(landCC) & is.null(sceneCC)){
      additionalCriteria <- NULL
    } else {
      additionalCriteria <- '"additionalCriteria":{"filterType":"and","childFilters":['
      filt <- NULL
      #Collection category
      if(!is.null(collectionCategory)) filt <- c(filt,paste0('{"filterType":"value","fieldId":"',fieldIDs["collectionCategory",dataSet],'","value":','"',collectionCategory,'"',',"operand":"="}'))
      #Data Type Level-1
      if(!is.null(dataTypeLevel1)) filt <- c(filt,paste0('{"filterType":"value","fieldId":"',fieldIDs["dataTypeLevel1",dataSet],'","value":','"',dataTypeLevel1,'"',',"operand":"="}'))
      #Scene cloud cover
      if(!is.null(sceneCC)) filt <- c(filt,paste0('{"filterType":"value","fieldId":"',fieldIDs["sceneCC",dataSet],'","value":','"',sceneCC,'"',',"operand":"<"}'))
      #Land cloud cover
      if(!is.null(landCC)) filt <- c(filt,paste0('{"filterType":"value","fieldId":"',fieldIDs["landCC",dataSet],'","value":','"',landCC,'"',',"operand":"<"}'))
      #Combine criteria
      additionalCriteria <- paste0(additionalCriteria,paste(filt,collapse=","),']},')
    }
    
    ##  Query
    requestURL <- GET(paste0(serviceURL,
                             'search?jsonRequest={',
                             '"datasetName":"',dataSet,'",',
                             '"lowerLeft":{"latitude":"',yMin,'","longitude":"',xMin,'"},',
                             '"upperRight":{"latitude":"',yMax,'","longitude":"',xMax,'"},',
                             '"startDate":"',startDate,'","endDate":"',endDate,'",',
                             searchMonths,
                             additionalCriteria,
                             '"includeUnknownCloudCover":false,',
                             '"maxResults":"',maxResults,'",',
                             '"sortOrder":"ASC",',
                             '"apiKey":"',apiKey,'",',
                             '"node":"EE"}'))
    responseContent <- fromJSON(content(requestURL, "text", encoding="UTF-8"))
    
    ##  Check number of hits and download
    if (responseContent$data$totalHits > 0){
      sceneIDs <- responseContent$data$results$entityId
      displayIDs <- responseContent$data$results$displayId
      allProductIDs <- c(allProductIDs, displayIDs)
      
      ##  Download data
      if(!is.null(download)){
        if (!file.exists(download)) dir.create(download, recursive=TRUE)
        #Remove ID's of previously downloaded scenes
        
        if(!forceOverwrite){
          locs <- !(file.size(file.path(download,paste0(displayIDs,".tar.gz")))>0)
          locs[is.na(locs)] <- TRUE
        } else {
          locs <- !logical(length(displayIDs))
        }
  
        if (sum(locs)>0){
          if (setVerbose) cat(paste0("Downloading ",sum(locs)," files:\n"))
          sceneIDs <- sceneIDs[locs]
          displayIDs <- displayIDs[locs]
          for (s in 1:length(sceneIDs)){
            if (setVerbose) cat(paste0("  ",displayIDs[s]," "))
            requestURL <- GET(paste0(serviceURL,
                                     'download?jsonRequest={',
                                     '"datasetName":"',dataSet,'",',
                                     '"products":["STANDARD"],',
                                     '"entityIds":["',sceneIDs[s],'"],',
                                     '"apiKey":"',apiKey,'","node":"EE"}'))
            responseContent <- fromJSON(content(requestURL, "text", encoding="UTF-8"))
            if (!is.null(responseContent$errorCode)){
              if (setVerbose) cat("    Error retrieving download url, skipping\n")
              #Do something to remove file from fileslist
              #return(NULL)
            } else {
              if (setVerbose) cat(" ... ")
              downURL <- responseContent$data
                #correct file naming in case extra text is added at back filename 
              downURL_corr <- file.path(dirname(downURL), substr(basename(downURL),1,47))

              fn.out <- file.path(download,paste0(displayIDs[s],".tar.gz"))
              ptd <- proc.time()
              ret <- download.file(downURL_corr, fn.out, "wget", quiet = TRUE)
              ptd <- proc.time()-ptd
              if (ret>0 & setVerbose) 
                cat("    Error downloading file, skipping\n") 
              #Do something to remove file from fileslist
              if (ret==0 & setVerbose) 
                cat(" Completed: ",round(unname(ptd[3]),0)," sec (",round((file.size(fn.out)/1048576)/unname(ptd[3]),2),"MB/sec) \n", sep="")
            }
          } #for s
        } else {
          if (setVerbose) cat("All files downloaded\n")
        }
      } #if
    } #if
  } #for dataSet
  
  ##  Logout
  requestURL <- paste0(serviceURL,'logout?jsonRequest={"apiKey":"',apiKey,'"}')
  responseContent <- fromJSON(requestURL)
  
  ##  Return code for successfull completion
  if (setVerbose) cat("**  EarthExplorer Inventory Completed  **\n")
  return(allProductIDs)
  
}


