######################################################################################
# 
# Spatial distribution of the different Mycobacterium bovis clades
# 
# This code generates Figure 3 of paper in (1).
# Figure caption: Spatial analysis of Mycobacterium bovis isolates.  
# A. Spatial analysis of distribution of M. bovis clades identified by Bayesian 
# phylogenetic analysis.  Each polygon represents the minimum convex polygon of the 
# sampled locations of the isolates of each clade. 
# B. Comparison of spatial distances between estimated and permuted clades. For every 
# pair of clades being compared we have randomly selected 1000 isolates from each. For 
# each random pair of isolates we calculated the spatial distance between them. This 
# analysis was repeated with random (permuted) clade assignments. 
# C. K-means analysis with 3 clusters (represented by symbols) versus 3 clades 
# (represented by colors). D. Optimal number of clusters estimated by within group sum 
# of squares (distances between individuals within each cluster). The optimal number of
# clusters will be the number after which within cluster differences become minimal; 
# here this occurs after ~ 3 clusters.
#
# (1) Manuscript "Implications for disease management at the wildlife-livestock 
# interface: using whole-genome sequencing to study the role of elk in Mycobacterium bovis 
# transmission in Michigan, USA" by L.C.M. Salvador, D.J. O’Brien, M.K. Cosgrove, T.P. Stuber, 
# A. Schooley, J. Crispell, S. Church, Y.T., Grohn, S. Robbe-Austerman, R.R. Kao
# 
# @developed by Joseph Crispell, November 2018
#
# Input file:
# Data/Michigan_sequence_traits_with_clades.csv  
#  
# Output file: 
# MI_Figure3.pdf
# 
####################################################################################
#### Load libraries ####

library(contoureR)

#### Load data ####

# Read in the data 
# Note: Modified from the original file; Data about farms id and spatial positions 
# were deleted

file <- "Data/MI_Elk_Data_134isolates_traits_withClades_v3.csv"
table <- read.table(file, header=TRUE, sep=",")

#### Examine spatial clustering ####

# Open the pdf
pdf("MI_Figure3.pdf")

# Set the plotting dimensions
par(mfrow=c(2,2))

# Define colours for the clusters 
colours <- c(rgb(1,0,0, 0.75), rgb(0,1,0, 0.75), rgb(0,0,1, 0.75), rgb(0,0,0, 0.75))
table$Colour <- colours[table$Clade]

# Plot the locations and colour by clade
plotDataWithPolygons(table, colours)

# Compare random samples of the data with real data and permuted data
differences <- calculateDistancesBetweenRandomSamples(table, nComparison=1000)
differencesPermutedClades <- calculateDistancesBetweenRandomSamples(table, nComparison=1000, permute=TRUE)
plotDifferences(differences, differencesPermutedClades)

# Kmeans clustering
runKMeansClustering(table, colours)

# Reset the plotting window dimensions
par(mfrow=c(1,1))

# Close the PDF
dev.off()


#### FUNCTIONS ####

plotDataWithPolygons <- function(table, colours){
  
  # Plot the points
  plot(table$POINT_X, table$POINT_Y, col=table$Colour, pch=19, bty="n", xlab="X", ylab="Y", 
       main="", las=1)
  
  # Add a polygon for each clade
  for(clade in unique(table$Clade)){
    
    # Get the subset of the data for the current clade
    subset <- table[table$Clade == clade, ]
    
    # Calculate the convex hull
    hull <- getConvexHull(x=subset$POINT_X, y=subset$POINT_Y)
    
    # Plot the convex hull
    polygon(x=subset[hull, "POINT_X"], y=subset[hull, "POINT_Y"], border=colours[clade])
  }
  
  # Add a legend
  legend("topright", legend=c("Clade:", unique(table$Clade)), text.col=c("black", colours[unique(table$Clade)]), bty="n")
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  # Add a plot label
  mtext("A", side=3, line=1, at=axisLimits[1], cex=2.5)
}

runKMeansClustering <- function(table, colours){
  # Code mostly taken from: http://gsp.humboldt.edu/OLM/R/03_02_ClusterAnalysis.html
  
  # Create a matrix to store the coordinate data
  matrix <- as.matrix(data.frame(x=table$POINT_X, y=table$POINT_Y, clade=table$Clade))
  
  # Run a kmeans analysis for 4 clusters - the number of clades we have
  clusterInfo <- kmeans(matrix[, c("x", "y")], centers=3)
  
  # Plot the kmeans clusters and 
  plot(matrix, col=colours[matrix[, "clade"]], pch=c(21, 22, 23, 24)[clusterInfo$cluster], 
       bg=colours[matrix[, "clade"]], bty="n", las=1)
  
  # Add legend
  legend("topright", legend=c("Clade", unique(table$Clade)), 
         text.col=c("black", colours[unique(table$Clade)]), bty="n")
  legend("top", legend=c("K-means cluster", sort(unique(clusterInfo$cluster))), 
         pch=c(19, c(21, 22, 23, 24)[sort(unique(clusterInfo$cluster))]), bty="n",
         col=c("white", "black", "black", "black", "black"))
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  # Add a plot label
  mtext("C", side=3, line=1, at=axisLimits[1], cex=2.5)
  
  # Calculate the variance of each column
  variance <- apply(matrix, 2, var)
  
  # Find the sum of squares for 1 cluster
  withinClusterSumOfSquares = (nrow(matrix)-1)*sum(variance)
  
  # Find the sum of squares for 2 to 15 clusters
  for (i in 2:30) {
    clusterInfo <- kmeans(matrix, centers=i)
    withinClusterSumOfSquares[i] <- sum(clusterInfo$withinss)
  }
  
  # Plot the result
  plot(1:30, withinClusterSumOfSquares, type="o", xlab="Number of Clusters", ylab="Within groups sum of squares",
       bty="n", las=1)
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  # Add a plot label
  mtext("D", side=3, line=1, at=axisLimits[1], cex=2.5)
  
}

plotDifferences <- function(differences, differencesPermutedClades){
  
  # Get the Y limits
  yLim <- c(0, max(c(max(differences), max(differencesPermutedClades))))
  
  # Build an initial empty plot
  plot(x=NULL, y=NULL, xlim=c(0.75, 3.5), ylim=yLim, bty="n", xaxt="n", xlab="",
       ylab="Spatial difference (m)", las=1,
       main="")
  
  # Add an X axis
  axis(side=1, at=1:3, labels=colnames(differences))
  
  # Plot the 95% bounds of the distributions before and after permutation
  for(i in 1:3){
    
    # Get the column
    column <- colnames(differences)[i]
    
    # Plot the points from BEFORE permutating the clades
    bounds <- quantile(differences[, column], probs=c(0.025, 0.975))
    points(x=c(i-0.15, i-0.15), y=c(bounds[1], bounds[2]), type="l", col="red", lwd=3)
    points(x=c(i-0.1, i-0.2), y=c(bounds[1], bounds[1]), type="l", col="red", lwd=3)
    points(x=c(i-0.1, i-0.2), y=c(bounds[2], bounds[2]), type="l", col="red", lwd=3)
    points(x=i-0.15, y=median(differences[, column]), pch=19)
    
    # Plot the points from AFTER permutating the clades
    bounds <- quantile(differencesPermutedClades[, column], probs=c(0.025, 0.975))
    points(x=c(i+0.15, i+0.15), y=c(bounds[1], bounds[2]), type="l", col="blue", lwd=3)
    points(x=c(i+0.1, i+0.2), y=c(bounds[1], bounds[1]), type="l", col="blue", lwd=3)
    points(x=c(i+0.1, i+0.2), y=c(bounds[2], bounds[2]), type="l", col="blue", lwd=3)
    points(x=i+0.15, y=median(differencesPermutedClades[, column]), pch=19)
  }
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  # Add a plot label
  mtext("B", side=3, line=1, at=axisLimits[1], cex=2.5)
  
  # Add a legend
  legend("top", legend=c("Before permutation", "After permutation"), text.col=c("red", "blue"), bty="n")
}

calculateDistancesBetweenRandomSamples <- function(table, nComparisons=100, permute=FALSE){
  
  # Permute the Clade IDs if requested
  if(permute){
    table$Clade <- sample(table$Clade)
  }
  
  # Initialise a data frame to store the distances between random samples from each clade
  differences <- data.frame(Diff_1_2=NA, Diff_1_3=NA, Diff_2_3=NA)
  
  # Compare each clade once
  for(i in unique(table$Clade)){
    for(j in unique(table$Clade)){
      
      # Skipping when i >= j means all clades will only be compared once
      if(i >= j){
        next
      }
      
      # Subset the data for Clade i and j
      subsetI <- table[table$Clade == i, ]
      subsetJ <- table[table$Clade == j, ]
      
      # Randomly pick Clade i coordinates
      randomI <- sample(seq_len(nrow(subsetI)), replace=TRUE, size=nComparisons)
      
      # Randomly pick Clade j coordinates
      randomJ <- sample(seq_len(nrow(subsetJ)), replace=TRUE, size=nComparisons)
      
      # Calculate spatial distance between random coordinates
      for(n in seq_len(nComparisons)){
        
        differences[n, paste0("Diff_", i, "_", j)] <- 
          euclideanDistance(x1=subsetI[randomI[n], "POINT_X"], y1=subsetI[randomI[n], "POINT_Y"],
                            x2=subsetJ[randomJ[n], "POINT_X"], y2=subsetJ[randomJ[n], "POINT_Y"])
      }
    }
  }
  
  return(differences)
}

euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt(sum((x1 - x2)^2 + (y1 - y2)^2)))
}

