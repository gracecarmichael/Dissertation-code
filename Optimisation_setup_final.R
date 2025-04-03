################################################################################
# R code used to set up optimisation models 
# author: Grace Carmichael
################################################################################

# read in packages needed for reading in shapefiles and plotting 
library(sf)
library(tmap)

#-------------------------------------------------------------------------------
# Get ward data and set up ward plots 
#-------------------------------------------------------------------------------

# read in ward location dataset
wards <- sf::st_read(dsn = "MDBWard2016.gdb")
# remove Western Cape wards
wards <- wards[which(wards$ProvinceName!="Western Cape"),]

# convert pick-up point coordinates to sf object
pup_sf <- st_as_sf(pup_loc, coords = c("Longitude", "Latitude"), crs = 4326)

# setup plot color palette
palette("Okabe-Ito")
pal <- palette()
pal[2] <- "#F5C35C"
pal[3] <- "#9DC3F7"
pal[4] <- "#D3E6A9"
pal[6] <- "#D3F2EA"
pal[7] <- "#C8C2F2"
pal[8] <- "#FEB1A3"
pal[9] <- "#EEC9E5"

# plot wards and overlay pick-up points
tmap_mode("plot")
tmap_options(check.and.fix = TRUE)
tm_shape(wards) + tm_fill("ProvinceName", n = 8, palette = pal[-1], title = "Province") + tm_borders() + tm_shape(pup_sf) + tm_dots(col = "firebrick", size=0.07) 


# store clinic coordinates as sf object
DHIS_geo <- DHIS_geo[!is.na(DHIS_geo$coordinates),]
clinic_sf <- st_as_sf(DHIS_geo, coords = c("longitude", "latitude"), crs = 4326)

# plot wards with clinics overlayed
tm_shape(wards) + tm_fill("ProvinceName", n = 8, palette = pal[-1], title = "Province") + tm_borders() + tm_shape(clinic_sf[which(clinic_sf$Province!="wc Western Cape Province"),]) + tm_dots(col = "dodgerblue4", size=0.07) 


# read in population data for wards
pop_by_ward <- box_read_rds(file_id = 1809633143383)
pop_by_ward <- as.data.frame(pop_by_ward)

# order by ward ID
pop_by_ward <- pop_by_ward[order(pop_by_ward$WardID),]


#-------------------------------------------------------------------------------
# create variables needed for optimisation 
#-------------------------------------------------------------------------------

# create dataset with clinics and which ward they fall into
pnts <- st_transform(clinic_sf, 2263)
pnts$Ward <- rep(0, nrow(pnts))

# get logical matrix with clinics on rows and wards on columns, indicating 
# whether a pup is inside a ward (true) or not (false)
int <- st_intersects(pnts, st_transform(wards,2263), sparse = FALSE)
sum(int)

# populate Ward column with the ward that pup i is located in
for(i in 1:nrow(int)){
  j <- which(int[i,] == TRUE)
  pnts$Ward[i] <- wards[j,]$WardID
}

# save as file
write.csv(pnts, file = "pnts.csv")
pnts <- read.csv("pnts.csv",header=TRUE,row.names=1, check.names=FALSE)


## create h (average total headcounts by ward) - use average monthly headcounts 
# over last 5 years

# create dataframe with clinic names and average monthly headcount over 5 years
h_clinic <- data.frame(name=clinics_all$organisationunitname, 
                       headcount=apply(clinics_all[,73:133], 1, mean, 
                                       na.rm=TRUE))
# create variable (h) to store headcounts by ward
h <- rep(0, nrow(wards))
# populate h with headcounts
for(i in 1:nrow(wards)){
  #only execute if there is a clinic in the ward
  if(length(which(pnts$Ward == wards$WardID[i]))!=0){  
    # get names of clinics in ward i
    clinics_in_ward <- pnts[which(pnts$Ward == wards$WardID[i]),4] 
    # sum headcounts of all clinics in ward i
    h[i] <- sum(h_clinic[h_clinic$name %in% clinics_in_ward,2]) 
  }
  # if no clinic in ward i assign headcount of 0
  else{h[i]<-0} 
}

# add ward IDs as entry names
names(h) <- wards$WardID

# save as file
write.table(h,file="h.txt", col.names = TRUE)
h <- read.table("h.txt",header=TRUE,row.names=1, check.names=FALSE)

# order h to be in alphabetic order of ward IDs
h <- h[order(rownames(h)),]
# if headcounts are NA replace with 0
h[is.na(h)] <- 0


## create N (matrix of all wards rows and columns - entry is 1 if a pup in row 
# ward can serve column ward - ie. if travel time is less than D) 
# set value of D (max travel time, will vary this later)
D <- 60
# read in travel times dataset
traveltimes <- read.csv("Wards_OD_Traveltime.csv") # this is as pairwise observations - want as a matrix

# store row ward IDs
rows <- sub(" -.*", "", traveltimes$Name)
# store column ward IDs
columns <- sub(".*- ", "", traveltimes$Name)
#row_wards <- row_wards[which(row_wards %in% wards_no_WC)]
#row_wards <- row_wards[which(column_wards %in% wards_no_WC)]
row_wards <- rows[which(rows %in% wards_no_WC & columns %in% wards_no_WC)]
column_wards <- columns[which(rows %in% wards_no_WC & columns %in% wards_no_WC)]
# store travel times between row wards and column wards
times <- traveltimes$Total_Minutes
times <- times[which(rows %in% wards_no_WC & columns %in% wards_no_WC)]
# convert into matrix
tt_mat <- as.matrix(xtabs(times ~ row_wards + column_wards))
# if travel time between i and j <= D then TRUE (1), else FALSE (0)
N <- tt_mat<=D

# save as file
write.table(N,file="N.txt", col.names = TRUE)
N <- read.table("N.txt",header=TRUE,row.names=1, check.names=FALSE)

write.table(tt_mat,file="tt_mat.txt")
tt_mat <- read.table("tt_mat.txt",header=TRUE,row.names=1, check.names=FALSE)
