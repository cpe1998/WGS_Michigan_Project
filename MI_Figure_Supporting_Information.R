###################################################################################
# Code that generates Supporting Information Figures S1 to S5 of manuscript in (1).
#
# (1) Manuscript "Implications for disease management at the wildlife/livestock 
# interface: using whole-genome sequencing to study the role of elk in bovine Tuberculosis 
# transmission in Michigan, USA" by L.C.M. Salvador, D.J. O’Brien, M.K. Cosgrove, T.P. Stuber, 
# A. Schooley, J. Crispell, S. Church, Y.T., Grohn, S. Robbe-Austerman, R.R. Kao
# 
# @developed by lcmsalvador, July 2017
#
############################################################################
## Figure S1 ##
#
# Caption: Bayesian dating permutation test for Mycobacterium bovis 
# isolates sampled from white-tailed deer, cattle and elk during 
# 1999 and 2013 in Michigan USA. The y axis corresponds to estimates 
# of the substitution rate (subst/site/year) and the x axis indicates 
# different replicate datasets. The error bars show the 95% highest 
# posterior density intervals. The real estimate (in red) was obtained 
# from the original sampling times and it is compared with the estimates 
# from 20 data sets with randomized tip labels (in black).
#
# Input files for RandomCluster:
# 1. Cluster/MI_Sequences.xml
# 2. Cluster/MI_Sequences.clusters.csv
#  
# Output files: MI_SequencesRep"i".xml, 1<=i<=20
# 
# Input files for PlotDRT:
# 1. MI_Sequences.log
# 2. MI_SequencesRep"i".log, 1<=i<=20
#
###########################################################################
# Packages
install.packages("TipDatingBeast")
library("TipDatingBeast")

# Date Randomization test
# creats 20 new files with randomization of sample times 
#
# Randomize sampling dates from file Cluster/MI_Sequences.xml and it creates 
# 20 new files MI_Sequences"i".xml, with 1<=i<=20 and the sampling time randomized
# It loads cluster file 'MI_Sequences.cluster.csv'
RandomCluster("MI_Sequences", reps = 20, loadCluster = T, writeTrees = F)

# all the new sequences MI_SequencesRep"i".xml will need to be run with BEAST

#error plot with the original data and the randomization runs
# it reads the BEAST output: MI_Sequences.log and MI_SequencesRep"i".log files
# it outputs file "Fig_DRT_MI_Sequences.pdf
PlotDRT("MI_Sequences", reps = 20, burnin=0.1)





############################################################################
## Figure S2 ##
#
# Caption: Posterior probabilities (PP) of the phylogenetic tree internal 
# nodes from an ancestral state host reconstruction using discrete traits. 
# Every node with PP>0.75 is considered to have strong support.
#
# Input: posterior probability values of tree internal nodes
# Output file: histogram_internal_nodes.png
###########################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# histogram with probability of internal nodes     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

internal.nodes <- c(0.88, 0.98, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99, 0.97, 0.6, 0.94, 0.57, 0.98, 1, 1, 1, 0.99, 0.98, 1, 1, 1, 0.97, 0.97, 1, 0.78, 0.83, 0.96, 1, 1, 1, 1, 1, 0.98, 0.57, 0.53, 0.99, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

png(file=paste0("histogram_internal_nodes.png"), width=12,height=8.5,units="in",res=300)
plot(hist(internal.nodes, breaks=16), main="")
abline(v = 0.75, col = "red", lty = 3)
dev.off()




############################################################################
## Figure S3 ##
#
# Caption: The estimated evolutionary rate of the sampled Mycobacterium bovis
# population.  The sampled posterior distribution of the evolutionary rate was
# estimated by BEAST 2 analyses using a strict clock and constant population 
# size. The mean evolutionary rate was estimated to be 0.41 substitutions per
# genome per year.
#
# Input file: Data/MI_Sequences.log
# Output file: clock_rate_distribution.png
###########################################################################

table <- {}
nSites <- 698
data <- read.table("Data/MI_Sequences.log", header=TRUE, skip=2)

# Convert the per site value into a per genome nSites 
table$GenomeClockRate <- data$clockRate.MI_Sequences * nSites
table$GenomeClockRate[1] <- NA
mean <- mean(table$GenomeClockRate, na.rm=TRUE) 
lower <- min(table$GenomeClockRate, na.rm=TRUE) 
upper <- max(table$GenomeClockRate, na.rm=TRUE) 

png(file="clock_rate_distribution.png", width=7, height=7, units="in", res=300)
h<-hist(table$GenomeClockRate, breaks = 100, xlab="Mean clock rate/ Genome/ Year")
cuts <- cut(h$breaks, c(-Inf, lower, upper,Inf))
plot(h, col=c("red", "white", "red")[cuts], xlab="Mean clock rate/ Genome/ Year")
abline(v= mean, col="red")
dev.off()




############################################################################
## Figure S4 ##
#
# Caption: The minimum pairwise genetic distance between isolates sampled 
# from cattle (‘Cattle’), white tailed deer (‘WTD’) and elk (‘Elk’). The 
# genetic distance was computed as the number of sites that differ between
# each pair of sequences (no repetition) using the R package ‘ape’. ‘N pairs’
# represent the number of times each comparison between host species had the
# lowest genetic distance (some isolates had minimum pairwise genetic 
# distances with more than one isolate that could be from the same or different
# host species).
#
# Input files:
# 1. Data/MI_sequences.nexus
# 2. Data/MI_Traits.csv
# 
# Output file: boxplot_min_genetic_distance.png
# 
###########################################################################
# code for Figures S4 and S5

install.packages("ape")
install.packages("reshape")
install.packages("geosphere")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("grid")
install.packages("gridExtra")

library("ape")
library("reshape")
library("geosphere")
library("ggplot2")
library("cowplot")
library("grid")
library("gridExtra")


# read alignments
data <- read.nexus.data("Data/MI_sequences.nexus")

# read traits file
traits <- read.csv("Data/MI_Traits.csv")

# compute genetic distance between isolates
############################################
genetic.distance <- dist.dna(as.DNAbin(data), model="raw", as.matrix=TRUE) # This is simply the proportion or the number of sites that differ between each pair of sequences. 

genetic.distance.melt <- melt(genetic.distance)
colnames(genetic.distance.melt) <- c("isolate1", "isolate2", "genetic.distance")

# eliminate duplicate interactions
aux1 <- genetic.distance.melt[, c("isolate2", "isolate1", "genetic.distance")]
colnames(aux1) <- c("isolate1", "isolate2", "genetic.distance")
genetic.distance.melt <- unique(rbind(genetic.distance.melt[,c(1:3)], aux1))

# assign species to each isolate (species is on the traits file)
genetic.distance.melt$species1 <- traits$Species[match(genetic.distance.melt$isolate1,traits$Isolate)] 
genetic.distance.melt$species2 <- traits$Species[match(genetic.distance.melt$isolate2,traits$Isolate)] 

matrix <- genetic.distance.melt

# eliminate comparisons between the same isolates
matrix <- matrix[!(matrix$isolate1 == matrix$isolate2),]

# create group interaction
dd <- which((matrix$species1 == "Deer" & matrix$species2 == "Deer"), arr.ind=TRUE)
matrix$interaction[dd] <- "Deer-Deer"
dc <- which((matrix$species1 == "Deer" & matrix$species2 == "Cattle") | (matrix$species2 == "Deer" & matrix$species1 == "Cattle"), arr.ind=TRUE)
matrix$interaction[dc] <- "Deer-Cattle"
ed <- which((matrix$species1 == "Deer" & matrix$species2 == "Elk") | (matrix$species2 == "Deer" & matrix$species1 == "Elk"), arr.ind=TRUE)
matrix$interaction[ed] <- "Deer-Elk"
ce <- which((matrix$species1 == "Cattle" & matrix$species2 == "Elk") | (matrix$species2 == "Cattle" & matrix$species1 == "Elk"), arr.ind=TRUE)
matrix$interaction[ce] <- "Cattle-Elk"
cc <- which((matrix$species1 == "Cattle" & matrix$species2 == "Cattle"), arr.ind=TRUE)
matrix$interaction[cc] <- "Cattle-Cattle"
ee <- which((matrix$species1 == "Elk" & matrix$species2 == "Elk"), arr.ind=TRUE)
matrix$interaction[ee] <- "Elk-Elk"


N <- table(matrix$interaction)
N <- melt(N)

# keep the number of occurrences of each interaction type
a <- gsub("Deer-Deer", paste0("Deer-Deer\n", as.character(N[N$Var.1=="Deer-Deer",]$value)), matrix$interaction)
b <- gsub("Deer-Elk", paste0("Deer-Elk\n", as.character(N[N$Var.1=="Deer-Elk",]$value)), a)
c <- gsub("Deer-Cattle", paste0("Deer-Cattle\n", as.character(N[N$Var.1=="Deer-Cattle",]$value)), b)
d <- gsub("Elk-Elk", paste0("Elk-Elk\n", as.character(N[N$Var.1=="Elk-Elk",]$value)), c)
e <- gsub("Cattle-Cattle", paste0("Cattle-Cattle\n", as.character(N[N$Var.1=="Cattle-Cattle",]$value)), d)
matrix$interaction2 <- gsub("Cattle-Elk", paste0("Cattle-Elk\n", as.character(N[N$Var.1=="Cattle-Elk",]$value)), e)



# minimum genetic distance
min.gdist <-  unique(aggregate(genetic.distance~isolate1,matrix, min))

a <- unique(merge(min.gdist, matrix))
a <- unique(a[, c("isolate1", "genetic.distance", "interaction" )])


N<-melt(table(a$interaction))
a$interaction2 <- a$interaction

a1 <- gsub("Deer-Deer", paste0("Deer-Deer\n", as.character(N[N$Var.1=="Deer-Deer",]$value)), a$interaction)
b1 <- gsub("Deer-Elk", paste0("Deer-Elk\n", as.character(N[N$Var.1=="Deer-Elk",]$value)), a1)
c1 <- gsub("Deer-Cattle", paste0("Deer-Cattle\n", as.character(N[N$Var.1=="Deer-Cattle",]$value)), b1)
d1 <- gsub("Elk-Elk", paste0("Elk-Elk\n", as.character(N[N$Var.1=="Elk-Elk",]$value)), c1)
a$interaction2 <- gsub("Cattle-Cattle", paste0("Cattle-Cattle\n", as.character(N[N$Var.1=="Cattle-Cattle",]$value)), d1)
#a$interaction2 <- gsub("Cattle-Elk", paste0("Cattle-Elk\n", as.character(N[N$Var.1=="Cattle-Elk",]$value)), e1)



png(file=paste0("boxplot_min_genetic_distance.png"), width=12,height=8.5,units="in",res=300)


p<-ggplot(a, aes(x=interaction2, y=genetic.distance)) + 
  geom_boxplot() +
  labs(y="minimum pairwise genetic distance", x="")+  
  geom_boxplot(fill='#999999')+
  theme( 
    axis.title.x = element_text(size=26),
    axis.title.y = element_text(size=26),
    axis.text.x = element_text(size=22),
    axis.text.y = element_text(size=20)
  )

grid.newpage()
footnote <- "N pairs"
g <- arrangeGrob(p, bottom = textGrob(footnote, x = 0.01, hjust = -0.1, vjust=-3.1, gp = gpar(fontsize = 22)))
grid.draw(g)
dev.off()





##################################################################################
## Figure S5 ##
#
# The minimum pairwise spatial distance between isolates sampled from cattle
# (‘Cattle’), white-tailed-deer (‘WTD’) and elk (‘Elk’). The spatial distance
# was computed as the shortest distance between two points (i.e., the 
# great-circle-distance), according to the haversine method using the R package
# ‘geosphere’. ‘N pairs’ represent the number of times each comparison between 
# host species had the shortest spatial distance.
#
# Input files:
# 1. Data/MI_sequences.nexus
# 2. Data/MI_Traits.csv
# 
# Output file: boxplot_min_spatial_distance.png
#################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# pairwise minimum spatial distance (separate by type)            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

min.sdist <-  unique(aggregate(spatial.distance~isolate1,matrix, min))

a <- unique(merge(min.sdist, matrix))
a <- unique(a[, c("isolate1", "spatial.distance", "interaction" )])

N<-melt(table(a$interaction))
a$interaction2 <- a$interaction

a1 <- gsub("Deer-Deer", paste0("Deer-Deer\n", as.character(N[N$Var.1=="Deer-Deer",]$value)), a$interaction)
b1 <- gsub("Deer-Elk", paste0("Deer-Elk\n", as.character(N[N$Var.1=="Deer-Elk",]$value)), a1)
c1 <- gsub("Deer-Cattle", paste0("Deer-Cattle\n", as.character(N[N$Var.1=="Deer-Cattle",]$value)), b1)
d1 <- gsub("Elk-Elk", paste0("Elk-Elk\n", as.character(N[N$Var.1=="Elk-Elk",]$value)), c1)
a$interaction2 <- gsub("Cattle-Cattle", paste0("Cattle-Cattle\n", as.character(N[N$Var.1=="Cattle-Cattle",]$value)), d1)
#a$interaction2 <- gsub("Cattle-Elk", paste0("Cattle-Elk\n", as.character(N[N$Var.1=="Cattle-Elk",]$value)), e1)


png(file=paste0("boxplot_min_spatial_distance_interactions.png"), width=12,height=8.5,units="in",res=300)


p<-ggplot(a, aes(x=interaction2, y=spatial.distance)) + 
  geom_boxplot() +
  labs(y="minimum pairwise spatial distance", x="")+  
  geom_boxplot(fill='#999999')+
  theme( 
    axis.title.x = element_text(size=26),
    axis.title.y = element_text(size=26),
    axis.text.x = element_text(size=22),
    axis.text.y = element_text(size=20)
  )

grid.newpage()
footnote <- "N pairs"
g <- arrangeGrob(p, bottom = textGrob(footnote, x = 0.01, hjust = -0.1, vjust=-3.1, gp = gpar(fontsize = 22)))
grid.draw(g)

dev.off()


############################################################################