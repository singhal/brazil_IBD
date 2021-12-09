require(alphahull)
require(rgeos)
require(rgdal)
require(maps)
require(mapdata)
require(maptools)
require(rJava)
library(readxl)
require(dismo)
require(raster)
require(sp)
data(wrld_simpl)

make_polygon <- function(coords) {
  polygon = NA
  
  #if enough points for alpha hulls
  if (nrow(coords) > 3) {
    polygon <-try(getAlphaHullRange(coords, percent=1, partCount= partcount, buff=50000, coordHeaders=c('origLong','origLat'), verbose=TRUE), silent=TRUE)
    while ('try-error' %in% class(q)) {
      polygon <-try(getAlphaHullRange(coords, percent=1, partCount= partcount, buff=50000, coordHeaders=c('origLong','origLat'), verbose=TRUE), silent=TRUE)
    }
    
    #if 3 points, make minimum convex hull
  } else if (nrow(coords) == 3) {
    pt <- SpatialPoints(coords, proj4string=CRS('+proj=longlat +datum=WGS84'))
    pt <- spTransform(pt, CRS("+init=epsg:3395"))
    b <- gBuffer(gConvexHull(pt), width=50000)
    b <- spTransform(b, CRS('+proj=longlat +datum=WGS84'))
    polygon <- list(b, NA, 3)
    
    #if 2 points, make tube polygon
  } else if (nrow(coords) == 2){
    
    polygon <- list(linkTwoPoints(SpatialPoints(coords, proj4string=CRS('+proj=longlat +datum=WGS84'))), NA, 2)
    
    # if 1 point, make buffered point
  } else if (nrow(coords) == 1) {
    
    pt <- SpatialPoints(coords, proj4string=CRS('+proj=longlat +datum=WGS84'))
    pt <- spTransform(pt, CRS("+init=epsg:3395"))
    b <- gBuffer(pt, width=50000)
    b <- spTransform(b, CRS('+proj=longlat +datum=WGS84'))
    polygon <- list(b, NA, 1)
  }
  return(polygon)	
}


make_clip_polygon <- function(polygon) {
  # clip polygon to the coast of sa
  clip_polygon = polygon
  clip_polygon[[1]] = gIntersection(polygon[[1]], sa2)
  return(clip_polygon)
}


make_buffer_polygon <- function(polygon) {
  # put a buffer around the polygon 
  temp <- spTransform(polygon[[1]], CRS("+init=epsg:3395"))
  temp <- gBuffer(temp, width=buffdist)
  temp = Filter(function(f){f@ringDir == 1}, temp@polygons[[1]]@Polygons)
  temp = SpatialPolygons(list(Polygons(temp,ID=1)), proj4string=CRS('+init=epsg:3395'))
  
  buffer_polygon <- spTransform(temp, CRS('+proj=longlat +datum=WGS84'))
  # tests if the buffer geometry is valid
  if (gIsValid(buffer_polygon) == FALSE) {
    buffer_polygon = gBuffer(polygon, width=0)
  }
  
  buffer_polygon = gIntersection(buffer_polygon, sa2)
  return(buffer_polygon)
}


get_env_data <- function(rast_dir) {
  #load rasters
  setwd(rast_dir)
  env <- stack(list.files(pattern='.tif$'))
  names(env) <- gsub('CHELSA_bio10_', '', names(env))
  return(env)
}

run_maxent <- function(env, coords, bg) {
  occ = coords
  
  # don't want to run maxent if all the coordinates
  # are in the same grid - no information there
  if (nrow(occ) < 10) {
    # number unique cells
    num_cells = length(unique(cellFromXY(env[[1]], occ)))
    if (num_cells == 1) {
      # will halt things in the next step
      occ = occ[1,]
    }
  }
  
  # only run MaxEnt if you have a certain number of points
  if (nrow(occ) >= min_pts) {
    occ <- SpatialPointsDataFrame(coords=occ[,c('origLong','origLat')],data=occ,proj4string=CRS('+proj=longlat +datum=WGS84'))
    
    # figure out what occurrence points are associated with missing data
    x <- raster::extract(env, occ)
    occ <- occ[which(complete.cases(x)),]
    
    if (!missing(bg)) {
      bgocc <- SpatialPointsDataFrame(coords=bg[,c('origLong','origLat')], data=bg, proj4string=CRS('+proj=longlat +datum=WGS84'))
      
      # figure out what occurrence points are associated with missing data
      x <- raster::extract(env, bgocc)
      bgocc <- bgocc[which(complete.cases(x)),]
    }
    
    print("run maxent 1")
    if (missing(bg)) {
      # actually run maxent!
      xm <- maxent(x=env,p=occ,args=c("randomseed=true"), silent = F)
    } else {
      xm <- maxent(x=env,p=occ, a= bgocc, args=c("randomseed=true"))		
    }
    
    # thins variables so no overfitting	
    perm <- xm@results
    # gets permutation significance results
    perm <- perm[grep('permutation', rownames(perm)),]
    names(perm) <- gsub('.permutation.importance', '', names(perm))
    # selects only those that are significant above 5
    topVar <- names(which(perm > 5))
    # cannot build a model with a single variable
    # so if there is only one top variable, then add the second
    if (length(topVar) == 1) { 
      topVar <- c(topVar, names(sort(perm, decreasing=TRUE))[2])
    }
    
    print("run maxent 2")
    if (missing(bg)) {
      xm <- maxent(x=env[[topVar]], p=occ, args=c("randomseed=true"))
    } else {
      xm <- maxent(x=env[[topVar]],p=occ, a= bgocc, args=c("randomseed=true"))
    }
    print("predict")
    px <- predict(xm, env[[topVar]], progress='text')
    
    #least training presence threshold
    presenceProbs <- raster::extract(px,occ)
    thresh <- quantile(presenceProbs,threshold)
    # converts 0 - 1 habitat suitability to binary presence / absence
    px[which(values(px) < thresh)] <- NA
    px[which(values(px) >= thresh)] <- 1
    
    print("binary to poly")
    master <- binaryToPolygon(px, projection = proj4string(occ), buffer=10000, occ, dropNoOcc=TRUE)
    results  = list(master, xm)
    # results = list(px, xm)
    return(results)
  } else {
    # returns NA if too few points
    # or if all points in the same grid cell
    return(list(NA, NA))
  }
}


clip_hull_enm <- function(clip_polygon, buffer_polygon, enm) {
  # only clip if actually an ENM
  # if no ENM, the range becomes one of the simple polygons defiend by occurrences
  if (!class(enm) %in% c('SpatialPolygons','SpatialPolygonsDataFrame')) {
    res = clip_polygon[[1]]
  } else {
    res = gIntersection(buffer_polygon, enm)
  }
  return(res)
}

write_range <- function(range, sp, outdir) {
  # save our new ranges!
  sp_name = gsub("\\. ", "_", sp)
  sp_name = gsub(" ", "_", sp_name)
  filename = paste(outdir, sp_name, sep="")
  if (class(range) == 'SpatialPolygons') {
    IDs <- sapply(slot(range, "polygons"), function(x) slot(x, "ID"))
    df <- data.frame(rep(0, length(IDs)), row.names=IDs)
    writeSpatialShape(SpatialPolygonsDataFrame(range, data=df), filename)
  } else {
    writeSpatialShape(range, filename)
  }
}

# location of pascal scripts
dir = '/Users/sonal/Dropbox/scripts/eco_IBD/pascal_scripts/'
# location of raster files
rast_dir = '/Users/sonal/Desktop/wc2.1_2.5m_bio/'

# source some of pascal's helper scripts
source(paste(dir, 'ah2sp.R', sep=""))
source(paste(dir, 'ENM2range.R', sep=""))

sa = readOGR("~/Dropbox/brazil_gene_flow/data/ne_10m_land/ne_10m_land.shp")
sa2 = crop(sa, extent(-85, -34, -56, 13))

env1 = get_env_data(rast_dir)
env = crop(env1, extent(-85, -34, -56, 13))

d = read.table("~/Dropbox/brazil_gene_flow/data/distributions/Bothrops_moojeni.txt", stringsAsFactors = F, header = T)

# coords = as.data.frame(unique(d[,c("Long", "Lat")]))
# coords = coords %>% mutate_all(as.numeric)
# names(coords) = c("origLong", "origLat")
# # remove some obvious outliers & errors
# coords2 = coords[coords$origLong < -35, ]
# coords2 = coords2[coords2$origLong != -63.1414, ] 

d = readxl::read_xlsx("~/Dropbox/brazil_gene_flow/data/distributions/Table Sxx. Micrablepharus atticolus distribution.xlsx", skip = 1)

coords = as.data.frame(unique(d[,c("Longitude do Local", "Latitude do Local")]))
coords2 = coords %>% mutate_all(as.numeric)
names(coords2) = c("origLong", "origLat")

ll = read.csv("~/Dropbox/brazil_gene_flow/data/brazil_samples_v7.csv",
              stringsAsFactors = F)
ll2 =  ll[complete.cases(ll$LAT), ]
# ll3 = ll2[ll2$lineage == "Bothrops_moojeni", c("LON", "LAT")]
ll3 = ll2[ll2$lineage == "Micrablepharus_atticolus", c("LON", "LAT")]
names(ll3) = c("origLong", "origLat")
coords3 = unique(rbind(coords2, ll3))

# buffer distance to use in ENM clipping
buffdist <- 150000
# minimum number of points to use for ENM
min_pts = 3
# threshold for ENM
# anything in the 97% of the distribution 
# for suitability will be treated as part of the range
threshold <- 0.01

partcount <- 1
# get alpha hull / MCP / tube / buffered point polygon
polygon = make_polygon(coords3)
# check that it is a polygon
class(polygon[[1]])
# clip polygon to sa's coast
clip_polygon = make_clip_polygon(polygon)
buffer_polygon = make_buffer_polygon(polygon)
# get the range map based on enm
results = run_maxent(env, coords3)

# clip enm to the hull
range = clip_hull_enm(clip_polygon, buffer_polygon, results[[1]])
range2 <- spTransform(range, CRS('+proj=longlat +datum=WGS84'))
write_range(range2, "Micrablepharus_atticolus", "~/Desktop/distributions/")