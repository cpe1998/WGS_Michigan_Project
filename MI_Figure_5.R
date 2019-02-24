##############################################################################
# 
# This code parses the XML file to randomly shuffle the host species 
# associated to each isolate
# 
# It generates Figure 5 of paper in (1).
#
# Figure caption: Comparison of the estimated posterior support of direct host
# species transition between permuted and observed data. The estimated 
# posterior mean probability of each host species interaction is the posterior 
# probability that a particular transition rate is positive. If this probability
# is high, the data strongly support a model in which there is direct pathogen
# transition between that particular pair of host species. The posterior means
# were estimated via a Discrete Ancestral Trait Mapping performed in BEAST v2.
# The ‘Permuted data’ correspond to the posterior means of 10 BEAST runs of 
# each interaction after permuting the host species labels each time. The 
# ‘Observed data’ correspond to the posterior mean of each interaction using 
# the observed data. 
#
# (1) Manuscript "Disease management at the wildlife-livestock 
# interface: using whole-genome sequencing to study the role of elk in Mycobacterium bovis 
# transmission in Michigan, USA" by L.C.M. Salvador, D.J. O’Brien, M.K. Cosgrove, T.P. Stuber, 
# A. Schooley, J. Crispell, S. Church, Y.T., Grohn, S. Robbe-Austerman, R.R. Kao
# 
# @developed by lcmsalvador, July 2017
# @updated by lcmsalvador, November 2018
#
# Input files:
# 1. Data/MI_Elk_Data_134isolates_Traits_withClades.csv
# 2. Data/MI_Elk_134isolates_HKY_relExp_extskyline_DTA.xml
#  
# Output directory: 
# random_runs/
# 
# output files:
# random_runs/random_run_i.xml, 1<=i<=n
# 
# figures: Figure_5.pdf
#          Figure_5.png
####################################################################################

install.packages("doBy")
install.packages("ggplot2")
installed.packages("reshape")
library("doBy")
library("ggplot2")
library("reshape")

##################################
# RANDOMIZATION OF ALL TIP LABELS
##################################

# number of isolates
nisol <- 134

# number of iterations
n <- 10

# read traits file
traits <- read.csv('Data/MI_Elk_134isolates_Traits_withClades.csv', header=TRUE, sep=",", stringsAsFactors = default.stringsAsFactors())
species <- traits[, "SPECIES"]

indexes <- c(1:n)
for (j in 1:length(indexes)){

  text <- readLines('Data/MI_Elk_extra_isolates_HKY_relExp_extskyline_DTA.xml', n=-1)
  
  # delete commas after sequence id in taxa id and traitset
  ref <- which(grepl("traitSet id", text))
  lines <- text[c((ref+1):(ref+nisol))]
  lines.index <- c((ref+1):(ref+nisol))

  rand.species <- sample(species,nisol, replace=FALSE)
  print(table(rand.species))
  
  for (l in (ref+1):(ref+length(lines)-1)){
    a <- strsplit(lines[l-ref], "=")
    b <- paste0(as.character(rand.species[l-ref]), ",")
    text[l] <-paste0(a[[1]][1],'=', b)
  }
  a <- strsplit(lines[nisol], "=")
  b <- paste0(as.character(rand.species[nisol]), ",")
  text[ref+length(lines)] <-paste0(a[[1]][1],'=', substr(b,1,nchar(b)-1))
  
  dir.create('random_runs/')
  dir.create(paste0('random_runs/random_run_',as.character(indexes[j])))
  writeLines(text, paste0('random_runs/random_run_',as.character(indexes[j]),'/random_run.xml'))
}

#################################
# these will be run by BEAST.....
#################################


####################################################
# After BEAST analyses on the files generated above
# Code to generate Figure 5
####################################################

# rateIndicator values for original run:
cattle_deer <- 0.996
cattle_elk <- 0.391
deer_elk <- 0.989

random.traits <- {}
for (i in 1:n){
  print(i)
  data <- read.table(paste0("random_runs/random_run_", i, "/beast_", i, ".log"), header=TRUE, skip=2)
  data$i <- i
  random.traits <- rbind(random.traits, data)
}

# compute average of rateIndicators 
rateIndicators <- summaryBy(rateIndicator.Host1+rateIndicator.Host2+rateIndicator.Host3~i, data=random.traits, FUN=mean)

# states are attributed alphabetically
names(rateIndicators) <- c("i","Cattle-Deer", "Cattle-Elk", "Deer-Elk")

rateIndicators.melt <- melt(rateIndicators, id=c("i"))
names(rateIndicators.melt) <- c("i", "interaction", "value")

# values for original run
point <- data.frame(variable = c("Cattle-Deer", "Cattle-Elk", "Deer-Elk"), value = c(cattle_deer, cattle_elk, deer_elk))
names(point) <- c("interaction", "value")

# Create png file
png(file="Figure5.png", width=7, height=5, units="in", res=300)
ggplot(rateIndicators.melt, aes(x=interaction, y=value, fill = interaction)) + 
  geom_boxplot()+scale_fill_grey() + theme_classic() +
  annotate("text", x = point$interaction, y = point$value-0.01, label = "*", color="darkred", size=12)+
  labs(x = "Symmetric transition between host species", y="Estimated posterior probability", size=12)
dev.off()

# Create pdf file
pdf(file="Figure5.pdf", width=7, height=5)
ggplot(rateIndicators.melt, aes(x=interaction, y=value, fill = interaction)) + 
  geom_boxplot()+scale_fill_grey() + theme_classic() +
  annotate("text", x = point$interaction, y = point$value-0.01, label = "*", color="darkred", size=12)+
  labs(x = "Symmetric transition between host species", y="Estimated posterior probability", size=12)
dev.off()


