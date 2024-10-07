##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
##          This code has been cleaned after the study was finalized       ##
## If you notice errors in the code, contact Pyry Toivonen (pntoiv@utu.fi) ##
##        NOTE: code includes heavy analyses that take processing time     ##
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##

#-#-#-# Code chapters #-#-#-#-#-#_#-#-#
## 1 Settings and starting libraries -#
## 2 Hop event classification        -#
## 3 Creating variables              -#
## 4 Modelling (GAM)                 -#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

## With the data given in the repository, chapters 2 and 3 can be skipped
## Chapter 2 can be performed with data given in https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study4656276921

## GPS data is produced by Saaristoluonnon hoito- ja suojeluyhdistys SLHSY ry and published with their permission


##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#### 1 Settings and starting libraries ####
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)

setwd("//utuhome.utu.fi/pntoiv/CoastRacc/Submission/DataRepository") # Make sure to change this to your working directory


##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#### 2 Hop event classification  ####
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##


Saaret_AOI <- st_read("Islands.gpkg") # This dataset is obtained from the GitHub data repository

Saaret_AOI <- Saaret_AOI %>%
  rename(SaarID = ID)
Saaret_buffer <- st_buffer(Saaret_AOI, 5)

CoastRaccDogs <- read.csv("CoastRacc data.csv") # GPS data from MoveBank
CoastRaccDogs$timestamp <- as.POSIXct(CoastRaccDogs$timestamp, tz="UTC", format="%Y-%m-%d %H:%M:%S")

RaccSF <- st_as_sf(CoastRaccDogs, coords=c("location.long","location.lat"), crs=4326, remove=F) %>%
  st_transform(3067) %>%
  mutate(X = st_coordinates(.)[,1],
         Y = st_coordinates(.)[,2]) %>%
  rename(ID = "individual.local.identifier")

RaccData <- st_join(RaccSF, Saaret_buffer, join=st_intersects, largest=T)
RaccData$SaarID <- ifelse(is.na(RaccData$SaarID), 0, RaccData$SaarID)


DataSplit <- split(RaccData, f=RaccData$ID)

TrackList <- list()
for (i in unique(RaccData$ID)) {
  TrackList[[i]] <- make_track(DataSplit[[i]],
                               .x=X,
                               .y=Y,
                               .t=timestamp,
                               ID=ID,
                               SaarID=SaarID,
                               crs = st_crs(3067))
  
  TrackList[[i]] <- TrackList[[i]] %>%
    mutate(Difference = SaarID - lag(SaarID)) # Calculates the difference between SaarID (detects if it changes between two subsequent rows)
}


TrackDF <- do.call(rbind, TrackList)

TrackDF %>%
  filter(SaarID == 0) %>% # locations in water
  count()

# Creates new column based on detected island hops with following values:
# Value 0 = individual remained in the same island,
# Value 1 = individual left island into water (>5 meters from the shore)
# Value 2 = individual arrived at new island or same island after 3 or more consecutive steps of being in water
# Value 3 = left island but came back within two next steps


# Code for hop event classification based on if else clauses (?case_when)

TrackList_IslandHop <- list()

for (i in unique(RaccData$ID)) {
  
  TrackList_IslandHop[[i]] <- TrackList[[i]] %>%
    mutate(IslandHop = case_when(Difference != 0 & SaarID == 0 ~ 1,
                                 Difference != 0 & SaarID != 0 ~ 2,
                                 is.na(Difference) ~ NA,
                                 .default = 0)) %>%
    mutate(IslandHop = case_when(Difference != 0 & Difference == -lag(Difference) ~ 3,
                                 Difference != 0 & Difference == -lead(Difference) ~ 3,
                                 is.na(Difference) ~ NA,
                                 .default = IslandHop)) %>%
    mutate(IslandHop = case_when(IslandHop == 2 & lag(SaarID) == 0 & lag(SaarID, 3) == SaarID ~ 3,
                                 IslandHop == 2 & lead(SaarID) == 0 & lead(SaarID, 2) == lag(SaarID) ~ 3,
                                 .default = IslandHop))
  
}



TrackDF <- do.call(rbind, TrackList_IslandHop)

# Some statistics
HopStats <- TrackDF %>%
  group_by(ID, IslandHop) %>%
  count() %>%
  pivot_wider(id_cols=ID, names_from=IslandHop, values_from=n, values_fill=0) %>%
  rename(Left = "1", Explore = "3", Hop = "2", Indiffer = "0") %>%
  mutate(Total = Left+Explore+Hop+Indiffer) %>%
  select(ID, Explore, Left, Hop, Indiffer, Total) 



# Creating data-frame including switching events between islands (also includes events where individuals left shore and came back to the same island)
HopDates1 <- TrackDF %>%
  group_by(ID) %>%
  filter(IslandHop %in% c(2,NA)) %>% # we also include the first date (has NA class)
  pull(t_, name=ID)

Hoppers <- data.frame(ID = names(HopDates1),
                      Dates = HopDates1)

RacStepList <- split(TrackDF, f=TrackDF$ID)
HopperList <- split(Hoppers, f=Hoppers$ID)
EventList1 <- list()

for (i in unique(TrackDF$ID)) {
  
  for (j in 1:nrow(HopperList[[i]])) {
    paste(j)
    EventList1[[i]][[j]] <- RacStepList[[i]] %>%
      filter(t_ >= HopperList[[i]][["Dates"]][j] & t_ <= HopperList[[i]][["Dates"]][j+1]) %>%
      mutate(EventID = paste0(i,"_",j))
  }
  
}

# This creates a flat data.frame including all the events
HopEvents1 <- map_dfr(names(EventList1), function(id) {
  bind_rows(EventList1[[id]], .id = "DataFrame") %>%
    mutate(ID = id)
})

##-#-#-#-#-#-#-#-#-#-#-#-#-#-##
#### 3 Creating variables  ####
##-#-#-#-#-#-#-#-#-#-#-#-#-#-##

EventData <- HopEvents1
TrackData <- TrackDF

EventData$t_ <- as.POSIXct(EventData$t_, "%Y-%m-%d %H:%M:%S", tz="UTC")

TrackData$t_ <- as.POSIXct(TrackData$t_, "%Y-%m-%d %H:%M:%S", tz="UTC")

print(TrackData %>% group_by(ID) %>% count(), n=30)


TrackData %>%
  group_by(IslandHop) %>%
  count()


# Function for mode
mode_value <- function(x) {
  x_nonzero <- x[x != 0]
  ux <- unique(x_nonzero)
  ux[which.max(tabulate(match(x_nonzero, ux)))]
}

ModelData <- EventData %>%
  group_by(ID, EventID) %>%
  summarise(TimeBefore = difftime(last(t_), first(t_), unit="hours"), # excludes the step duration of the last step
            IslandBefore = mode_value(SaarID), # The most visited island before hopping to new island
            IslandAfter = last(SaarID), # Island where the individual hopped to
            EventDate = last(t_), # The reported date of hop event
            EventTimeDiff = difftime(last(t_), tail(t_,2)[1], unit="hours"), # Time difference between locations before and after event
            EventMonth = month(label=F, last(t_)), # The reported month of hop event
            EventX1 = tail(x_,2)[1],
            EventY1 = tail(y_,2)[1],
            EventX2 = last(x_),
            EventY2 = last(y_)) %>% 
  ungroup()

Saaret_AOI$SaarArea <- st_area(Saaret_AOI)

ModelData <- ModelData %>%
  left_join(st_drop_geometry(Saaret_AOI)
            %>% rename(IslandBefore = ID, AreaBefore = SaarArea), join_by(IslandBefore)) %>%
  left_join(st_drop_geometry(Saaret_AOI)
            %>% rename(IslandAfter = ID, AreaAfter = SaarArea), join_by(IslandAfter))


Distances <- rep(NA, nrow(ModelData))

for (i in 1:nrow(ModelData)) {
  if (ModelData$IslandBefore[i] != 0) {
    Distances[i] <- st_distance(filter(Saaret_AOI, ID == ModelData$IslandBefore[i]), filter(Saaret_AOI, ID == ModelData$IslandAfter[i]))
  }
}

ModelData$IslandDistance <- Distances

ModelData %>%
  filter(IslandDistance == 0)# 6 individuals with above mentioned events. These might need to to be removed for scewing data

ModelData %>% 
  group_by(ID) %>%
  count() 

length(unique(ModelData$ID)) # Total 24 individuals with island hops (1054 events)

hist(ModelData$IslandDistance)


DF <- st_as_sf(select(ModelData, EventX1, EventY1, ID, EventID, IslandBefore, IslandAfter), coords=c("EventX1","EventY1"), crs=3067, remove=F)

DistMatrix <- st_distance(DF, Saaret_AOI)

nrow(DistMatrix)
colnames(DistMatrix) <- Saaret_AOI$ID
rownames(DistMatrix) <- DF$EventID

min(DistMatrix)

MinDist <- rep(NA, nrow(DistMatrix))

MeanDist5 <- rep(NA, nrow(DistMatrix))

MeanDist10 <- rep(NA, nrow(DistMatrix))

MeanDistAll <- rep(NA, nrow(DistMatrix))


sort_nonzero <- function(x) {
  x_nonzero <- x[x != 0]
  sort(x_nonzero)
}


for (i in 1:nrow(DistMatrix)) {
  
  MinDist[i] <- sort_nonzero(as.vector(DistMatrix[i,]))[1]
  
  MeanDist5[i] <- mean(sort_nonzero(as.vector(DistMatrix[i,]))[1:5], na.rm=T) 
  
  MeanDist10[i] <- mean(sort_nonzero(as.vector(DistMatrix[i,]))[1:10], na.rm=T) 
  
  MeanDistAll[i] <- mean(sort_nonzero(as.vector(DistMatrix[i,])), na.rm=T) # laskee keskiarvon ignooraamalla nollan
  
}

ModelData <- ModelData %>%
  mutate(MinDistIslands = MinDist,
         MeanDist_k5 = MeanDist5,
         MeanDist_k10 = MeanDist10,
         MeanDist_all = MeanDistAll)
# Dist unit = meter


# Adding the judas variable (days since last judas date)
JudasTable <- read.csv("PartnerRemovals.csv", sep=";")

JudasTable$JudasDate <- as.POSIXct(JudasTable$JudasDate, tz="UTC", format="%Y-%m-%d")
JudasTable$ID <- as.factor(JudasTable$ID)
ModelData$ID <- as.factor(ModelData$ID)
ModelData$Date <- date(ModelData$EventDate)



ModelData <- ModelData %>%
  left_join(JudasTable, join_by(ID,  closest(Date >= JudasDate)))

ModelData <- ModelData %>%
  mutate(SinceJudas = difftime(Date, JudasDate, units="days"))

ModelData5 <- ModelData %>%
  left_join(rename(JudasTable, JudasBefore = JudasDate), join_by(ID,  closest(Date <= JudasBefore)))

ModelData <- ModelData5 %>%
  mutate(BeforeJudas = difftime(date(EventDate), JudasBefore, units="days"))


Icedata <- nc_open("FMI_IceConc.nc") # https://doi.org/10.48670/moi-00132

time <- ncvar_get(Icedata, "time")
time_units <- ncatt_get(Icedata, "time", "units")
time_obs <- as.POSIXct(time, origin = "1981-01-01 00:00:00", tz="UTC")

IceRaster_C <- rast("FMI_IceConc.nc", subds="ice_concentration")
IceRaster_T <- rast("FMI_IceThick.nc", subds="ice_thickness")


IceC.values <- rep(NA, nrow(ModelData))
IceT.values <- rep(NA, nrow(ModelData))

for (i in 1:nrow(ModelData)) {
  xvec <- as.vector(unlist(c(ModelData[i,9], ModelData[i,11])))
  yvec <- as.vector(unlist(c(ModelData[i,10], ModelData[i,12])))
  
  ExtVec <- c(min(xvec), max(xvec), min(yvec), max(yvec))
  
  e <- ext(ExtVec)
  
  w <- matrix(1, nrow=5, ncol=5)
  r_filled <- focal(IceRaster_C[[time(IceRaster_C) == ModelData[i,]$Date]], w=w, fun="mean", na.policy="only", na.rm=T)
  r_filled2 <- focal(IceRaster_T[[time(IceRaster_T) == ModelData[i,]$Date]], w=w, fun="mean", na.policy="only", na.rm=T)
  
  r <- crop(r_filled, project(x=e, from="epsg:3067", to = crs(IceRaster_C)))
  
  IceC.values[i] <- global(r, fun="mean", na.rm=T)$mean
  
  r2 <- crop(r_filled2, project(x=e, from="epsg:3067", to = crs(IceRaster_T)))
  
  IceT.values[i] <- global(r2, fun="mean", na.rm=T)$mean
  
}


ModelData$IceCover <- IceC.values
ModelData$IceThickness <- IceT.values

ModelData <- ModelData %>%
  mutate(IceCover = if_else(is.na(IceCover), 0, IceCover),
         IceThickness = if_else(is.na(IceThickness), 0, IceThickness))



# Find out the last coordinates before hop from the most visited island
N0_EventX1 <- rep(NA,length(unique(EventData$EventID)))
N0_EventY1 <- rep(NA,length(unique(EventData$EventID)))

for (i in 1:length(unique(EventData$EventID))) {
  
  j <- unique(EventData$EventID)[i]
  
  df <- filter(EventData, EventID == j)
  
  N0_EventX1[i] <- last(filter(slice_tail(df,n=10), SaarID == mode_value(df$SaarID))$x_)
  N0_EventY1[i] <- last(filter(slice_tail(df,n=10), SaarID == mode_value(df$SaarID))$y_)
}

Addition <- data.frame(EventID = unique(EventData$EventID),
                       N0_EventX1 = N0_EventX1,
                       N0_EventY1 = N0_EventY1)

ModelData <- ModelData %>%
  left_join(Addition, join_by("EventID"))

ModelData <- ModelData %>%
  filter(IslandDistance != 0)


## Ice cover again

IceRaster_C <- rast("FMI_IceConc.nc", subds="ice_concentration")
IceRaster_T <- rast("FMI_IceThick.nc", subds="ice_thickness")
ModelData$Date <- as.POSIXct(ModelData$Date, tz="UTC", format="%Y-%m-%d")


IceC.values <- rep(NA, nrow(ModelData))
IceT.values <- rep(NA, nrow(ModelData))

for (i in 1:nrow(ModelData)) {
  xvec <- as.vector(unlist(c(ModelData[i,27], ModelData[i,11])))
  yvec <- as.vector(unlist(c(ModelData[i,28], ModelData[i,12])))
  
  ExtVec <- c(min(xvec), max(xvec), min(yvec), max(yvec))
  
  e <- ext(ExtVec)
  
  w <- matrix(1, nrow=5, ncol=5)
  r_filled <- focal(IceRaster_C[[time(IceRaster_C) == ModelData[i,]$Date]], w=w, fun="mean", na.policy="only", na.rm=T)
  r_filled2 <- focal(IceRaster_T[[time(IceRaster_T) == ModelData[i,]$Date]], w=w, fun="mean", na.policy="only", na.rm=T)
  
  r <- crop(r_filled, project(x=e, from="epsg:3067", to = crs(IceRaster_C)))
  
  IceC.values[i] <- global(r, fun="mean", na.rm=T)$mean
  
  r2 <- crop(r_filled2, project(x=e, from="epsg:3067", to = crs(IceRaster_T)))
  
  IceT.values[i] <- global(r2, fun="mean", na.rm=T)$mean
  
}


ModelData$IceCover <- IceC.values
ModelData$IceThickness <- IceT.values

ModelData <- ModelData %>%
  mutate(IceCover = if_else(is.na(IceCover), 0, IceCover),
         IceThickness = if_else(is.na(IceThickness), 0, IceThickness))

Sex <- read.csv("Sex.csv", sep=";")

# Lets to some unit transformations
ModelData2 <- ModelData %>%
  mutate(AreaBefore_km2 = as.vector(AreaBefore/1000000),
         AreaAfter_km2 = as.vector(AreaAfter/1000000),
         IslandDistance_km = IslandDistance/1000,
         MinDistIslands = MinDistIslands/1000,
         MeanDist_k5 = MeanDist_k5/1000,
         MeanDist_k10 = MeanDist_k10/1000,
         MeanDist_all = MeanDist_all/1000,
         TimeBefore = as.vector(TimeBefore),
         Judas = as.factor(case_when(is.na(SinceJudas) ~ "Before",
                                     .default = "After")),
         Ice = as.factor(case_when(IceCover < 25 ~ "Water",
                                   .default = "Ice"))) %>%
  filter(IslandDistance > 0) %>%
  left_join(Sex, join_by(ID))

ModelData2$ID <- as.factor(ModelData2$ID)
ModelData2$EventMonth <- as.factor(ModelData2$EventMonth)
ModelData2$Sex <- as.factor(ModelData2$Sex)
ModelData2$IslandAfter <- as.factor(ModelData2$IslandAfter)
ModelData2$IslandBefore <- as.factor(ModelData2$IslandBefore)


Saaret_AOI$ID <- as.factor(Saaret_AOI$ID)

StartPoints <- sf::st_as_sf(ModelData2, coords=c("N0_EventX1","N0_EventY1"), crs=3067, remove=F)
EndPoints <- sf::st_as_sf(ModelData2, coords=c("EventX2","EventY2"), crs=3067, remove=F)



## Least cost path analysis ##
library(leastcostpath)


## For loop
SteppingStonesList <- rep(NA, nrow(StartPoints))
WaterTravelList <- rep(NA, nrow(StartPoints))
LeastCostLength <- rep(NA, nrow(StartPoints))
IslandNumber <- rep(NA, nrow(StartPoints))
IslandDensity <- rep(NA, nrow(StartPoints))
SS_PlotList <- list()

for (i in 1:nrow(StartPoints)) {
  
  BBOX_Ext <- ext(st_bbox(rbind(StartPoints[i,],EndPoints[i,])))
  
  new_extent <- c(BBOX_Ext[1] - 1000, BBOX_Ext[2] + 1000, 
                  BBOX_Ext[3] - 1000, BBOX_Ext[4] + 1000) 
  
  NewExt <- ext(new_extent)
  
  Cropped_Islands <- st_crop(Saaret_AOI, NewExt)
  
  r <- terra::rast(NewExt, resolution=5, crs=crs(Cropped_Islands))
  
  r2 <- terra::rasterize(Cropped_Islands, r, field = 5, background = 1) # Polygons get weight of 10 while sea gets 1
  
  cs1 <- create_cs(r2, neighbours=16)
  
  path <- create_lcp(cs1, StartPoints[i,], EndPoints[i,])
  
  path_sf <- st_transform(path, st_crs(Cropped_Islands))
  
  SteppingStones <- st_intersection(Cropped_Islands, path_sf)
  
  WaterTravel <- st_difference(path_sf, st_union(Cropped_Islands))
  
  WaterTravel$Length <- st_length(WaterTravel)
  
  LeastCostLength[i] <- st_length(path_sf)
  
  
  if ("Length" %in% colnames(WaterTravel) && nrow(WaterTravel) > 0) {
    WaterTravelList[i] <- pull(WaterTravel, Length)
  } else {
    # Handle the case where Length column doesn't exist or WaterTravel is empty
    WaterTravelList[i] <- NA
  }
  
  SteppingStones <- SteppingStones %>%
    filter(!ID %in% c(StartPoints[i,]$IslandBefore, StartPoints[i,]$IslandAfter))
  
  SteppingStonesList[i] <- length(unique(pull(SteppingStones, ID)))
  
  IslandNumber[i] <- length(unique(pull(Cropped_Islands, ID)))
  
  IslandFreq <- freq(r2)
  
  IslandDensity[i] <- IslandFreq[IslandFreq$value==5,3]/sum(IslandFreq[,3])
  
  SS_PlotList[[i]] <- ggplot() +
    geom_sf(data=Cropped_Islands) +
    geom_sf(data=path_sf, colour="blue") +
    geom_sf(data=rbind(StartPoints[i,],EndPoints[i,]), color="red")
  
}

# With plots started 12:19, ended somewhere before 14:23
# Without plots started 10:50, ended 11:35
# With plots started 14:52, ended (ennuste: 17:00)

ModelData2$SteppingStones <- SteppingStonesList
ModelData2$IslandDensity <- IslandDensity
ModelData2$Islands_within <- IslandNumber
ModelData2$WaterTravelDist <- WaterTravelList
ModelData2$LeastCostLength <- LeastCostLength

ModData.G <- ModelData2 %>%
  dplyr::select(-JudasDate,-JudasBefore,-BeforeJudas) %>%
  mutate(SinceJudas_inv = case_when(is.na(SinceJudas) ~ 0,
                                    .default = 1/SinceJudas)) %>%
  mutate(SinceJudas_inv = case_when(is.infinite(SinceJudas_inv) ~ 1,
                                    .default = SinceJudas_inv)) %>%
  mutate(SinceJudas = case_when(is.na(SinceJudas) ~ -1,
                                .default = SinceJudas))


#write.csv(ModelData.G, "ModelData.csv", row.names=F)



##-#-#-#-#-#-#-#-#-#-#-#-##
#### 4 Modelling (GAM) ####
##-#-#-#-#-#-#-#-#-#-#-#-##


library(gamlss)
library(tidyverse)
library(ggeffects)

ModData.G <- read.csv("ModelData.csv") %>%
  mutate(Judas = as.factor(Judas),
         EventMonth = as.factor(EventMonth),
         Sex = as.factor(Sex),
         ID = as.factor(ID),
         IslandBefore = as.factor(IslandBefore),
         IslandAfter = as.factor(IslandAfter)) %>%
  mutate(Islands_within = Islands_within - 2) %>%
  mutate(Islands_within = case_when(Islands_within < 0 ~ 0,
                                    .default=Islands_within)) %>%
  mutate(WaterTravelKM = WaterTravelDist/1000,
         LeastCostKM = LeastCostLength/1000,
         DistDifference = WaterTravelKM - IslandDistance_km,
         PerDifference = DistDifference/IslandDistance_km,
         SteppingStoneDensity = SteppingStones/LeastCostKM) 


# Removing too influential outliers
ModData.G <- ModData.G %>% filter(!EventID %in% c("R23_53","R23_57","R23_39","R23_24","R23_58"))

# Removing male from Bolax (R4) because it already has partner in the data set and move together
ModData.G <- ModData.G %>% filter(!ID %in% c("R4"))

ModData.G$Judas <- relevel(ModData.G$Judas, "Before")


NullModel <- gamlss(log(IslandDistance_km) ~ 1,
                    data = ModData.G,
                    method = RS(),
                    family=NO,
                    control = gamlss.control(n.cyc = 500))

summary(NullModel) # Deviance 1165.339 with RE, without RE 3283.639


GAMLSS_LOGNO_5 <- gamlss(log(IslandDistance_km) ~
                           IceCover +
                           pb(Islands_within) +
                           pb(SteppingStones) +
                           Sex + 
                           log(AreaBefore_km2) +
                           log(AreaAfter_km2) +
                           Judas * EventMonth +
                           random(IslandBefore) +
                           random(IslandAfter),
                         data = ModData.G,
                         method = RS(),
                         family=NO,
                         control = gamlss.control(n.cyc = 500))

term.plot(GAMLSS_LOGNO_5, partial.resid=F)

plot(GAMLSS_LOGNO_5)
wp(GAMLSS_LOGNO_5, ylim.all=5)

summary(GAMLSS_LOGNO_5) 

# Plotting and other extra code (different model runs etc.) have been removed from this code for simplicity