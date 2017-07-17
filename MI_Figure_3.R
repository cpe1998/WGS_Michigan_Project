####################################################################################
# 
# Spatial distribution of the different Mycobacterium bovis clades
# 
# This code generates Figure 3 of paper in (1).
# Figure caption: Spatial distribution of Mycobacterium bovis clades identified by 
# Bayesian phylogenetic analysis. Colours of points and polygons correspond to the 
# colours in the phylogeny presented in Figure 2. Each polygon represents the minimum 
# convex polygon of the sampled locations of the isolates of each clade.
#
# (1) Manuscript "Implications for disease management at the wildlife/livestock 
# interface: using whole-genome sequencing to study the role of elk in bovine Tuberculosis 
# transmission in Michigan, USA" by L.C.M. Salvador, D.J. O’Brien, M.K. Cosgrove, T.P. Stuber, 
# A. Schooley, J. Crispell, S. Church, Y.T., Grohn, S. Robbe-Austerman, R.R. Kao
# 
# @developed by lcmsalvador, July 2017
#
# Input files:
# 1. Data/Michigan_sequence_traits_with_clades.csv  
# 2. Maps/Counties/Counties.dbf
# 3. Maps/Counties/Counties.shp
# 4. Maps/Counties/Counties.shx
#  
# Output file: michigan_clades_convex_hull.pdf
# 
####################################################################################

# packages
install.packages("maptools")
install.packages("rgeos")
library("maptools")
library("rgeos")


 # read traits file with clade ids for each isolate
clades <- read.csv("Data/Michigan_sequence_traits_with_clades.csv")
clade1 <- subset(clades, Clade==1)
clade2 <- subset(clades, Clade==2)
clade3 <- subset(clades, Clade==3)
clade4 <- subset(clades, Clade==4)


# read in shapefiles; here we use the specialized readShape* functions,
# but readShapeSpatial would produce identical output in all three cases
counties.mp <- readShapeSpatial("Maps/Counties/counties")


# specifying projection information is not strictly necessary for
# plotting, but does yield nicer default axis labels and aspect ratio in
# the case of geographic data
proj4string(counties.mp) <- "+proj=longlat +datum=WGS84"

# subset counties of interest:
# names <- c("Presque Isle", "Montmorency", "Otsego", "Oscoda", "Alcona", "Alpena", "Cheboygan", "Crawford")
names <- c("Presque Isle", "Montmorency", "Otsego", "Oscoda", "Alcona", "Alpena", "Cheboygan", "Crawford")

# only extracts the counties of interest
counties.mp <- counties.mp[is.element(counties.mp$NAME, names),]

# create spatial points data frame for each clade
spdf1 <- SpatialPointsDataFrame(coords = clade1[,c(6,7)], data = clade1,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

spdf2 <- SpatialPointsDataFrame(coords = clade2[,c(6,7)], data = clade2,
                                proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

spdf3 <- SpatialPointsDataFrame(coords = clade3[,c(6,7)], data = clade3,
                                proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

spdf4 <- SpatialPointsDataFrame(coords = clade4[,c(6,7)], data = clade4,
                                proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))


# creates a convex hull for the isolates of each clade
ch1 <- gConvexHull(spdf1)
ch2 <- gConvexHull(spdf2)
ch3 <- gConvexHull(spdf3)
ch4 <- gConvexHull(spdf4)


pdf(file="michigan_clades_convex_hull.pdf", width=7, height=7)
# generates simple map with the convex hull of each clade
plot(counties.mp, axes=TRUE, border="black", cex.axis=1.3)
points(clades[clades$Clade==1,]$Latitude, clades[clades$Clade==1,]$Longitude, pch=19, cex=0.7, col="royalblue")
points(clades[clades$Clade==2,]$Latitude, clades[clades$Clade==2,]$Longitude, pch=19, cex=0.7, col="coral4")
points(clades[clades$Clade==3,]$Latitude, clades[clades$Clade==3,]$Longitude, pch=19, cex=0.7, col="orange")
points(clades[clades$Clade==4,]$Latitude, clades[clades$Clade==4,]$Longitude, pch=19, cex=0.7, col="red")
plot(ch1,add=TRUE, border="royalblue", lwd=2)
plot(ch2,add=TRUE, border="coral4", lwd=2)
plot(ch3,add=TRUE, border="orange", lwd=2)
plot(ch4,add=TRUE, border="red", lwd=2)
legend("topright", legend=c("Clade 1", "Clade 2", "Clade 3", "Clade 4"), col=c("royalblue", "coral 4", "orange", "red"),lwd=c(1.6,1.6, 1.6, 1.6))
dev.off()

