###################################################################################
# Code that generates Supporting Information Figures S2 to S7 of manuscript in (1).
#
# (1) Manuscript "Disease management at the wildlife-livestock 
# interface: using whole-genome sequencing to study the role of elk in Mycobacterium bovis 
# transmission in Michigan, USA" by L.C.M. Salvador, D.J. O’Brien, M.K. Cosgrove, T.P. Stuber, 
# A. Schooley, J. Crispell, S. Church, Y.T., Grohn, S. Robbe-Austerman, R.R. Kao
# 
# @developed by lcmsalvador, July 2017
# @updated by lcmsalvador, November 2018

############################################################################
## Figure S2 ##
#
# Caption: Temporal signal in the sampled data. Bayesian dating permutation 
# test for Mycobacterium bovis isolates sampled from deer, cattle and elk 
# during 1996 and 2013 in Michigan USA. The y axis corresponds to estimates 
# of the estimated evolutionary rate (subst/site/year) and the x axis indicates 
# different replicate datasets. The error bars show the 95% highest posterior 
# density intervals. The real estimate (in red) was obtained from the original 
# sampling times and it is compared with the estimates from 10 data sets with 
# randomized tip labels (in black).
#
# Input files for RandomCluster:
# 1. MI_Elk_134isolates_Traits_withClades.csv
# 2. MI_Elk_extra_isolates_HKY_relExp_extskyline_DTA.xml
#  
# Output files: MI_SequencesRep"i".xml, 1<=i<=20
# FigureS2.pdf
# FigureS2.png
# 
# 
###########################################################################
###############################
# Temporal signal
# randomization of year labels
###############################

# number of isolates
nisol <- 134

# number of iterations
n <- 10

# Rate indicators of the original run
cattle_deer <- 0.996
cattle_elk <- 0.391
deer_elk <- 0.989

# read traits file
traits <- read.csv('Data/MI_Elk_134isolates_Traits_withClades.csv', header=TRUE, sep=",", stringsAsFactors = default.stringsAsFactors())

# retrieving year
id <- as.character(traits$ID)

aux <- strsplit(id, "_")
year <- c(1:length(id))
for (i in 1:length(id)){
  year[i] <-  aux[[i]][1] 
}

indexes <- c(1:n)
for (j in 1:length(indexes)){
  
  text <- readLines('Data/MI_Elk_extra_isolates_HKY_relExp_extskyline_DTA.xml', n=-1)
  
  # delete commas after sequence id in taxa id and traitset
  ref <- which(grepl("trait id", text))
  lines <- text[c((ref+1):(ref+nisol))]
  lines.index <- c((ref+1):(ref+nisol))
  
  rand.year <- sample(year,nisol, replace=FALSE)
  print(table(rand.year))
  
  for (l in (ref+1):(ref+length(lines)-1)){
    a <- strsplit(lines[l-ref], "=")
    b <- paste0(as.character(rand.year[l-ref]), ",")
    text[l] <-paste0(a[[1]][1],'=', b)
  }
  a <- strsplit(lines[nisol], "=")
  b <- paste0(as.character(rand.year[nisol]), ",")
  text[ref+length(lines)] <-paste0(a[[1]][1],'=', substr(b,1,nchar(b)-1))
  
  dir.create('temporal_random_runs/')
  dir.create(paste0('temporal_random_runs/random_run_',as.character(indexes[j])))
  writeLines(text, paste0('temporal_random_runs/random_run_',as.character(indexes[j]),'/random_run.xml'))
}


#################################
# these will be run by BEAST.....
#################################


########################################################################
# Input files (from BEAST output): 
# temporal_random_runs/temp_logs/beast_Rep", i, ".log. 1<=i<=10 
########################################################################

random.traits <- {}
for (i in 0:n){
  print(i)
  data <- read.table(paste0("temporal_random_runs/temp_logs/beast_Rep", i, ".log"), header=TRUE, skip=2)
  data$i <- i
  random.traits <- rbind(random.traits, data)
}
clockrate <- random.traits[, c("ucedMean.MI_Elk_Data_134isolates391snps", "i")]
names(clockrate) <- c("rate", "i")

df2 <- data_summary(clockrate, varname="rate", 
                    groupnames=c("i"))

df2$data[1] <- "original"
df2$data[2:11] <- "randomized" 

############################################################
# Function to calculate the mean and the standard deviation
# for each group
# data: data frame
# varname: the name of a column containing the variable to 
# be summarized
# groupnames : vector of column names to be used as
# grouping variables
############################################################

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

pdf(file="FigureS2.pdf", width=7, height=5)
p<- ggplot(df2, aes(x=i, y=rate)) + 
  geom_pointrange(aes(ymin=rate-sd, ymax=rate+sd, color=data)) +
  scale_fill_grey() + 
  theme_classic() +
  scale_color_manual(values=c('red','black')) +
  xlab("datasets") + 
  ylab("estimated evolutionary rate") +
  theme(axis.text.x=element_blank())
print(p)
dev.off()


png(file="FigureS2.png", width=7, height=5, units="in", res=300)
p<- ggplot(df2, aes(x=i, y=rate)) + 
  geom_pointrange(aes(ymin=rate-sd, ymax=rate+sd, color=data)) +
  scale_fill_grey() + 
  theme_classic() +
  scale_color_manual(values=c('red','black')) +
  xlab("datasets") + 
  ylab("estimated evolutionary rate") +
  theme(axis.text.x=element_blank())
print(p)
dev.off()



############################################################################
## Figure S3 #
#
# Caption: Model selection with Path Sampling. Bayes Factors (BF) were 
# calculated for seven different models run two times (for comparison) and 
# model selection was based on the highest average BF value of the two runs. 
# See Table S4 for a description of each model.
#
# Input: MLE and BF values for Path Sampling Analysis
# Figures: Figure_S3.png
#          
###########################################################################


values <- read.csv("data/MLE_PS_runs.csv")
toplot <- cbind(values$BF1, values$BF2)
names(toplot) <- c("run 1", "run 2")
toplot <- toplot[-1,]
toplot <- toplot[-8,]

toplot.melt <- melt(toplot)
names(toplot.melt) <- c("model", "run", "value")
bf <- toplot.melt
colnames(bf) <- c("run 1", "run 2")
rownames(bf) <- c("M1", "M2", "M3", "M4", "M5","M6","M7")
bf <- t(bf)
bf <- as.table(bf)


png(file="Figure_S3.png", width=5, height=5, units="in", res=300)
barplot(bf, ylab="Bayes Factor", col=c("darkgrey", "lightgrey"), ylim=c(0, 200), beside=TRUE, legend=rownames(bf))
dev.off()



############################################################################
## Figure S4 ##
#
# Caption: The estimated evolutionary rate of the sampled Mycobacterium 
# bovis population.  The sampled posterior distribution of the evolutionary
# rate was estimated by BEAST 2 analyses using an uncorrelated relaxed 
# exponential clock and an extended skyline demographic model. The mean 
# evolutionary rate was estimated to be 0.37 substitutions per genome per 
# year.
#
# Input file: Data/Beast_runs_combined.log
# Output file: Figure_S4.png
###########################################################################

nSites <- 391
mean <- 9.4872E-4 *nSites
lower <- 6.2573E-4 * nSites
upper <- 1.306E-3 * nSites

# Read in the data
table <- read.table("Data/Beast_runs_combined.log", header=TRUE, skip=2)
# Convert the per site value into a per genome nSites <-
table$GenomeClockRate <- table$ucedMean.MI_Elk_Data_134isolates391snps * nSites

png(file="FigureS4.png", width=5, height=5, units="in", res=300)
h<-hist(table$GenomeClockRate, breaks = 100, xlab="Estimated evolutionary rate \n (substitution/genome/year)")

cuts <- cut(h$breaks, c(-Inf, lower, upper,Inf))
plot(h, col=c("red", "white", "red")[cuts], xlab="Estimated evolutionary rate \n (substitution/genome/year)", main="")

abline( v= mean, col="red")
dev.off()


############################################################################
## Figures S5 and S6 ##
#
# Caption S5: The minimum pairwise genetic distance between isolates sampled 
# from cattle deer and elk. The genetic distance was computed as the number 
# of sites that differ between each pair of sequences using the R package 
# ‘ape’. For each isolate, we computed the genetic distance to all the other  
# isolates and recorded the host species interaction that had the minimum 
# distance. ‘N pairs’ represent the number of times each comparison between
# host species had the lowest genetic distance, however, some isolates had 
# minimum pairwise genetic distances with more than one isolate that could 
# be from the same or different host species. Interactions between cattle  
# and elk never had the minimum pairwise genetic distance between them.
#
# Caption S6: The minimum pairwise spatial distance between isolates sampled
# from cattle, deer and elk. The spatial distance was computed as the shortest
# distance between two points (i.e., the great-circle-distance), according to
# the haversine method using the R package ‘geosphere’. For each isolate, we 
# computed the spatial distance to all the other isolates and recorded the host
# species interaction that had the minimum distance. ‘N pairs’ represent the 
# number of times each comparison between host species had the shortest spatial
# distance. Interactions between cattle and elk never had the minimum pairwise
# spatial distance between them.
#
# Input files:
# 1. Data/MI_Elk_Data_134isolates391snps.nexus
# 2. Data/MI_Elk_134isolates_Traits_withClades.csv
# 
# Output files: Figure_S5.png, Figure_S6.png
# 
###########################################################################
library("ape")
library("reshape")
library("geosphere")
library("ggplot2")
library("cowplot")
library("grid")
library("gridExtra")


# compute genetic distances between isolates
# read alignments
data <- read.nexus.data("Data/MI_Elk_Data_134isolates391snps_07102018.nexus")

# read traits file
traits <- read.csv("Data/MI_Elk_Data_134isolates_traits.csv")

# compute genetic distance between isolates
############################################
genetic.distance <- dist.dna(as.DNAbin(data), model="raw", as.matrix=TRUE) 
# This is simply the proportion or the number of sites that differ between 
# each pair of sequences. 

genetic.distance.melt <- melt(genetic.distance)
colnames(genetic.distance.melt) <- c("isolate1", "isolate2", "genetic.distance")

# eliminate duplicate interactions
aux1 <- genetic.distance.melt[, c("isolate2", "isolate1", "genetic.distance")]
colnames(aux1) <- c("isolate1", "isolate2", "genetic.distance")
genetic.distance.melt <- unique(rbind(genetic.distance.melt[,c(1:3)], aux1))

# assign species to each isolate (species is on the traits file)
genetic.distance.melt$species1 <- traits$SPECIES[match(genetic.distance.melt$isolate1,traits$ID)] 
genetic.distance.melt$species2 <- traits$SPECIES[match(genetic.distance.melt$isolate2,traits$ID)] 


# compute spatial distance between isolates
############################################
# compute spatial distance between isolates
# distm: Distance matrix of a set of points, or between two sets of points
# The shortest distance between two points (i.e., the great-circle-distance or as the crow flies), according to the haversine method. This method assumes a spherical earth, ignoring ellipsoidal effects.
#distHaversine(p1, p2, r=6378137 - radius of earth)

spatial.distance <- distm(cbind(traits$POINT_X,traits$POINT_Y), fun=distHaversine)
rownames(spatial.distance) <- traits$ID
colnames(spatial.distance) <- traits$ID

spatial.distance.melt <- melt(spatial.distance)
colnames(spatial.distance.melt) <- c("isolate1", "isolate2", "spatial.distance")

matrix <- merge(genetic.distance.melt, spatial.distance.melt)
#17956

# eliminate comparisons between the same isolates
matrix <- matrix[!(matrix$isolate1 == matrix$isolate2),]
# 17822

# create group interaction
dd <- which((matrix$species1 == "DEER" & matrix$species2 == "DEER"), arr.ind=TRUE)
matrix$interaction[dd] <- "deer-deer"
dc <- which((matrix$species1 == "DEER" & matrix$species2 == "CATTLE") | (matrix$species2 == "DEER" & matrix$species1 == "CATTLE"), arr.ind=TRUE)
matrix$interaction[dc] <- "deer-cattle"
ed <- which((matrix$species1 == "DEER" & matrix$species2 == "ELK") | (matrix$species2 == "DEER" & matrix$species1 == "ELK"), arr.ind=TRUE)
matrix$interaction[ed] <- "deer-elk"
ce <- which((matrix$species1 == "CATTLE" & matrix$species2 == "ELK") | (matrix$species2 == "CATTLE" & matrix$species1 == "ELK"), arr.ind=TRUE)
matrix$interaction[ce] <- "cattle-elk"
cc <- which((matrix$species1 == "CATTLE" & matrix$species2 == "CATTLE"), arr.ind=TRUE)
matrix$interaction[cc] <- "cattle-cattle"
ee <- which((matrix$species1 == "ELK" & matrix$species2 == "ELK"), arr.ind=TRUE)
matrix$interaction[ee] <- "elk-elk"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# boxplot - pairwise minimum genetic distance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# minimum genetic distance
#min.gdist <-  unique(aggregate(genetic.distance~isolate1+isolate2,matrix, min))

min.gdist <-  unique(aggregate(genetic.distance~isolate1,matrix, min))
# 134
a <- unique(merge(min.gdist, matrix))
# 274
a <- unique(a[, c("isolate1", "isolate2", "genetic.distance", "interaction" )])

id <- which(table(a$isolate1)>1)

N<-melt(table(a$interaction))
a$interaction2 <- a$interaction

a1 <- gsub("deer-deer", paste0("deer-deer\n", as.character(N[N$Var.1=="deer-deer",]$value)), a$interaction)
b1 <- gsub("deer-elk", paste0("deer-elk\n", as.character(N[N$Var.1=="deer-elk",]$value)), a1)
c1 <- gsub("deer-cattle", paste0("deer-cattle\n", as.character(N[N$Var.1=="deer-cattle",]$value)), b1)
d1 <- gsub("elk-elk", paste0("elk-elk\n", as.character(N[N$Var.1=="elk-elk",]$value)), c1)
e1 <- gsub("cattle-cattle", paste0("cattle-cattle\n", as.character(N[N$Var.1=="cattle-cattle",]$value)), d1)
a$interaction2 <- gsub("cattle-elk", paste0("cattle-elk\n", as.character(N[N$Var.1=="cattle-elk",]$value)), e1)


png(file="Figure_S5.png", width=12,height=8.5,units="in",res=300)
n <- dim(a)[1]
a[n+1, c(1:4)] <- NA
a[n+1, 5] <- "cattle-elk\n0"

p<-ggplot(a, aes(x=interaction2, y=genetic.distance)) + 
  geom_boxplot() +
  labs(y="minimum pairwise genetic distance", x="")+  
  geom_boxplot(fill='#999999')+
  ylim(0,0.050) +
  theme( 
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16)
  )
grid.newpage()
footnote <- "N pairs"
g <- arrangeGrob(p, bottom = textGrob(footnote, x = 0.01, hjust = -0.1, vjust=-4.0, gp = gpar(fontsize = 16)))
grid.draw(g)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# pairwise minimum spatial distance (separate by type)            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

min.sdist <-  unique(aggregate(spatial.distance~isolate1,matrix, min))
a <- unique(merge(min.sdist, matrix))
a <- unique(a[, c("isolate1", "spatial.distance", "interaction" )])

N<-melt(table(a$interaction))
a$interaction2 <- a$interaction

a1 <- gsub("deer-deer", paste0("deer-deer\n", as.character(N[N$Var.1=="deer-deer",]$value)), a$interaction)
b1 <- gsub("deer-elk", paste0("deer-elk\n", as.character(N[N$Var.1=="deer-elk",]$value)), a1)
c1 <- gsub("deer-cattle", paste0("deer-cattle\n", as.character(N[N$Var.1=="deer-cattle",]$value)), b1)
d1 <- gsub("elk-elk", paste0("elk-elk\n", as.character(N[N$Var.1=="elk-elk",]$value)), c1)
e1 <- gsub("cattle-cattle", paste0("cattle-cattle\n", as.character(N[N$Var.1=="cattle-cattle",]$value)), d1)
a$interaction2 <- gsub("cattle-elk", paste0("cattle-elk\n", as.character(N[N$Var.1=="cattle-elk",]$value)), e1)

png(file="Figure_S6.png", width=12,height=8.5,units="in",res=300)

n <- dim(a)[1]
a[n+1, c(1:3)] <- NA
a[n+1, 4] <- "cattle-elk\n0"

p<-ggplot(a, aes(x=interaction2, y=spatial.distance)) + 
  geom_boxplot() +
  labs(y="minimum pairwise spatial distance", x="")+
  ylim(0, 18000) +
  geom_boxplot(fill='#999999')+
  theme( 
    axis.title.x = element_text(size=26),
    axis.title.y = element_text(size=26),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20)
  )

grid.newpage()
footnote <- "N pairs"
g <- arrangeGrob(p, bottom = textGrob(footnote, x = 0.01, hjust = -0.1, vjust=-3.1, gp = gpar(fontsize = 22)))
grid.draw(g)

dev.off()


############################################################################
## Figure S7 ##
#
# Caption S7: Extra “elk” analyses. Effect of adding 1 (A) and 2 (B) more 
# elk to the population by comparing the estimated posterior support of 
# direct host species transition between permuted and observed data. Through
# simulations we have extended our dataset with 1 and 2 extra elk to see their
# effect in the dynamics of the disease. We focused on the clades where we 
# have elk and cattle isolates (clades 1-2) and randomly chose 1 (A) and then 
# 2 (B) isolates that were extracted from the deer dataset and assumed added 
# them to the original elk dataset. We repeated this analysis 10 times and 
# computed the probability of support for pathogen transition for each 
# simulation via Discrete Ancestral Trait Mapping performed in BEAST v2 
# (‘Permuted data’). We compared it to the original analysis 
# (with 5 elk, ‘Observed Data’). The estimated posterior mean probability 
# of each host species interaction is the posterior probability that a 
# particular transition rate is positive. If this probability is high, the
# data strongly support a model in which there is direct pathogen transition 
# between that particular pair of host species. These results show that two 
# more elk in the system are not enough to generate pathogen interactions 
# between cattle and elk.   
#
# Input files:
# 1. Data/MI_Elk_134isolates_Traits_withClades.csv
# 2. Data/MI_Elk_extra_isolates_HKY_relExp_extskyline_DTA.xml
# 
# Output files: Figure_S7A.png
#               Figure_S7B.png
# 
###########################################################################
library("doBy")
library("ggplot2")
library("reshape")


#####################################################
# shuffle randomly elk labels of traits at tip nodes
#####################################################
# 117 deer
# 5 elk
# 12 cattle
nisol <- 134

# read traits file
traits <- read.csv('Data/MI_Elk_134isolates_Traits_withClades.csv', header=TRUE, sep=",", stringsAsFactors = default.stringsAsFactors())
species <- traits[, "SPECIES"]
clade123 <- which(traits$Clade!=4)
deer <- which(species[clade123]=="DEER")

######################
#### add 1 elk #######
indexes <- c(1:20)
for (j in 1:length(indexes)){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #  eliminate lines with the sequences to del seq.del
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  text <- readLines('Data/MI_Elk_extra_isolates_HKY_relExp_extskyline_DTA.xml', n=-1)
  
  # delete commas after sequence id in taxa id and traitset
  ref <- which(grepl("traitSet id", text))
  lines <- text[c((ref+1):(ref+nisol))]
  lines.index <- c((ref+1):(ref+nisol))
  #text <- text[-(lines.index)]
  
  # select random deer and replace it by one elk
  r <- sample(deer, 1)
  print(r)
  rand.species <- species
  rand.species[r] <- "ELK"
  #rand.species <- sample(species,nisol-1, replace=FALSE)
  print(table(rand.species))
  
  for (l in (ref+1):(ref+length(lines)-1)){
    a <- strsplit(lines[l-ref], "=")
    b <- paste0(as.character(rand.species[l-ref]), ",")
    text[l] <-paste0(a[[1]][1],'=', b)
  }
  a <- strsplit(lines[nisol], "=")
  b <- paste0(as.character(rand.species[nisol]), ",")
  text[ref+length(lines)] <-paste0(a[[1]][1],'=', substr(b,1,nchar(b)-1))
  
  dir.create('elk_random_runs/')
  dir.create(paste0('elk_random_runs/random_run_',as.character(indexes[j])))
  writeLines(text, paste0('elk_random_runs/random_run_',as.character(indexes[j]),'/elk_random_run.xml'))
}

######################
##### add 2 elk ######
######################

# 117 deer
# 5 elk
# 12 cattle
nisol <- 134

# read traits file
traits <- read.csv('Data/Data/MI_Elk_134isolates_Traits_withClades.csv', header=TRUE, sep=",", stringsAsFactors = default.stringsAsFactors())
species <- traits[, "SPECIES"]
clade123 <- which(traits$Clade!=4)
deer <- which(species[clade123]=="DEER")


indexes <- c(1:20)
for (j in 1:length(indexes)){
  
  text <- readLines('Data/MI_Elk_extra_isolates_HKY_relExp_extskyline_DTA.xml', n=-1)
  
  # delete commas after sequence id in taxa id and traitset
  ref <- which(grepl("traitSet id", text))
  lines <- text[c((ref+1):(ref+nisol))]
  lines.index <- c((ref+1):(ref+nisol))
  #text <- text[-(lines.index)]
  
  # select random deer and replace it by one elk
  r <- sample(deer, 2)
  print(r)
  rand.species <- species
  rand.species[r] <- "ELK"
  #rand.species <- sample(species,nisol-1, replace=FALSE)
  print(table(rand.species))
  
  for (l in (ref+1):(ref+length(lines)-1)){
    a <- strsplit(lines[l-ref], "=")
    b <- paste0(as.character(rand.species[l-ref]), ",")
    text[l] <-paste0(a[[1]][1],'=', b)
  }
  a <- strsplit(lines[nisol], "=")
  b <- paste0(as.character(rand.species[nisol]), ",")
  text[ref+length(lines)] <-paste0(a[[1]][1],'=', substr(b,1,nchar(b)-1))
  
  dir.create('elk_random_runs_2/')
  dir.create(paste0('elk_random_runs_2/random_run_',as.character(indexes[j])))
  writeLines(text, paste0('elk_random_runs_2/random_run_',as.character(indexes[j]),'/elk_random_run.xml'))
}


#################################
# these will be run by BEAST.....
#################################


########################################################################
# Input files (from BEAST output): 
# elk_random_runs_1/random_run_',i,'/elk_random_run.xml. 1<=i<=10
# elk_random_runs_2/random_run_',i,'/elk_random_run.xml  1<=i<=10
########################################################################
n<-10

# 1 extra elk

n<-10

random.traits <- {}
# i <- 1 has different structure (it has to be ignored)
for (i in 1:n){
  print(i)
  data <- read.table(paste0("elk_random_runs/logs_elk1/beast_", i, ".log"), header=TRUE, skip=2)
  data$i <- i
  random.traits <- rbind(random.traits, data)
}

rateIndicators <- summaryBy(rateIndicator.Host1+rateIndicator.Host2+rateIndicator.Host3~i, data=random.traits, FUN=mean)
names(rateIndicators) <- c("i","Cattle-Deer", "Cattle-Elk", "Deer-Elk")

rateIndicators.melt <- melt(rateIndicators, id=c("i"))
names(rateIndicators.melt) <- c("i", "interaction", "value")
# values for original run
point <- data.frame(variable = c("Cattle-Deer", "Cattle-Elk", "Deer-Elk"), value = c(0.996,0.391,0.989))
names(point) <- c("interaction", "value")

png(file="Figure_S7A.png", width=7, height=5, units="in", res=300)
ggplot(rateIndicators.melt, aes(x=interaction, y=value, fill = interaction)) + 
  geom_boxplot()+scale_fill_grey() + theme_classic() +
  annotate("text", x = point$interaction, y = point$value-0.01, label = "*", color="darkred", size=12)+
  labs(x = "Symmetric transition between host species", y="Estimated posterior probability", size=12)
dev.off()


# 2 extra elk
random.traits <- {}
for (i in 1:n){
  print(i)
  data <- read.table(paste0("elk_random_runs_2/logs_elk2/beast_", i, ".log"), header=TRUE, skip=2)
  data$i <- i
  random.traits <- rbind(random.traits, data)
}

rateIndicators <- summaryBy(rateIndicator.Host1+rateIndicator.Host2+rateIndicator.Host3~i, data=random.traits, FUN=mean)
names(rateIndicators) <- c("i","Cattle-Deer", "Cattle-Elk", "Deer-Elk")

rateIndicators.melt <- melt(rateIndicators, id=c("i"))
names(rateIndicators.melt) <- c("i", "interaction", "value")
# values for original run
point <- data.frame(variable = c("Cattle-Deer", "Cattle-Elk", "Deer-Elk"), value = c(0.996,0.391,0.989))
names(point) <- c("interaction", "value")

png(file="Figure_S7B.png", width=7, height=5, units="in", res=300)
ggplot(rateIndicators.melt, aes(x=interaction, y=value, fill = interaction)) + 
  geom_boxplot()+scale_fill_grey() + theme_classic() +
  annotate("text", x = point$interaction, y = point$value-0.01, label = "*", color="darkred", size=12)+
  labs(x = "Symmetric transition between host species", y="Estimated posterior probability", size=12)
dev.off()

