###########################################################################
# 
# This code randomly shuffles the host species associated to each isolate
# 
# It generates Figure 5 of paper in (1).
#
# Figure caption: Comparison of the estimated posterior support of direct host
# species transition between permuted and observed data. The estimated 
# posterior mean probability of each host species interaction is the 
# posterior probability that a particular transition rate is positive. 
# If this probability is high, the data strongly support a model in which
# there is direct pathogen transition between that particular pair of host
# species. The posterior means were estimated via a Discrete Ancestral Trait
# Mapping performed in BEAST v2. The ‘Permuted data’ correspond to the posterior
# means of 20 BEAST runs of each interaction after permuting the host species 
# labels each time. The ‘Observed data’ correspond to the posterior mean of 
# each interaction using the observed data. The label ‘Deer’ corresponds to the 
# host species ‘white-tailed deer’.
#
# (1) Manuscript "Implications for disease management at the wildlife/livestock 
# interface: using whole-genome sequencing to study the role of elk in bovine Tuberculosis 
# transmission in Michigan, USA" by L.C.M. Salvador, D.J. O’Brien, M.K. Cosgrove, T.P. Stuber, 
# A. Schooley, J. Crispell, S. Church, Y.T., Grohn, S. Robbe-Austerman, R.R. Kao
# 
# @developed by lcmsalvador, July 2017
#
# Input files:
# 1. Data/MI_Traits.csv
# 2. Data/MI_Sequences.xml
#  
# Output directory: 
# random_runs/
# 
# output files:
# random_runs/random_run_i.xml, 1<=i<=20
# 
# figure: prob_transition_species_randomization.png
####################################################################################

# read traits file
traits <- read.csv('Data/MI_Traits.csv', header=TRUE, sep=",", stringsAsFactors = default.stringsAsFactors())
species <- traits[, "Species"]

indexes <- c(1:20)

# loop to generate new 20 data sets with host species labels shuffled
for (j in 1:length(indexes)){
#  eliminate lines with the sequences to del seq.del
text <- readLines('Data/MI_Sequences.xml', n=-1)

# delete commas after sequence id in taxa id and traitset
ref <- which(grepl("traitSet id", text))
lines <- text[c((ref+1):(ref+53))]
lines.index <- c((ref+1):(ref+53))

# randomly shuffles host species labels
rand.species <- sample(species,53, replace=FALSE)
  
# identifies the label in the string (after '+') and it replaces it by the new one
for (l in (ref+1):(ref+length(lines)-1)){
  a <- strsplit(lines[l-ref], "=")
  b <- paste0(as.character(rand.species[l-ref]), ",")
  text[l] <-paste0(a[[1]][1],'=', b)
  }
  a <- strsplit(lines[53], "=")
  b <- paste0(as.character(rand.species[53]), ",")
  text[ref+length(lines)] <-paste0(a[[1]][1],'=', substr(b,1,nchar(b)-1))

  
  # creates directory with output (20 new .xml files with shuffled labels )
  dir.create('random_runs/')
  dir.create(paste0('random_runs/random_run_',as.character(indexes[j])))
  writeLines(text, paste0('random_runs/random_run_',as.character(indexes[j]),'/random_run.xml'))
}

#################################
# these will be run by BEAST.....
#################################


######################################################
# Input files (BEAST output): 
# random_runs/random_run"j"/random_run.log, 1<=j<=20 
######################################################


random.traits <- {}
# i <- 1 has different structure (it has to be ignored)
for (i in 1:20){
  print(i)
  data <- read.table(paste0("random_runs/random_run", i, "/random_run.log"), header=TRUE, skip=2)
  data$i <- i
  random.traits <- rbind(random.traits, data)
}

rateIndicators <- summaryBy(rateIndicator.species1+rateIndicator.species2+rateIndicator.species3~i, data=random.traits, FUN=mean)
names(rateIndicators) <- c("i","Cattle-Deer", "Cattle-Elk", "Deer-Elk")

rateIndicators.melt <- melt(rateIndicators, id=c("i"))

# values for original run
point <- data.frame(variable = c("Cattle-Deer", "Cattle-Elk", "Deer-Elk"), value = c(0.9691081,0.6030079,0.9871876))

png(file="prob_transition_species_randomization.png", width=5, height=5, units="in", res=300)
ggplot(rateIndicators.melt, aes(x=variable, y=value, fill = variable)) + 
  geom_boxplot()+
  annotate("text", x = point$variable, y = point$value-0.01, label = "*", color="red", size=6)+
  labs(x = "Symmetric transition between host species", y="Estimated posterior probability")
dev.off()


