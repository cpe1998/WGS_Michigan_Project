##############################################################################
# 
# This code parses the XML file to randomly downsample the number of isolates
# 
# It generates Figure 6 of paper in (1).
# 
# Figure caption: Comparison of the estimated posterior support of direct host species 
# transition between subsampled and observed data. The estimated posterior mean
# probability of each interaction is the posterior probability that a particular
# transition rate is positive. If this probability is high, then the data strongly 
# support a model in which there is direct pathogen transition between that 
# particular pair of host species. The posterior means were estimated via a 
# Discrete Ancestral Trait Mapping performed in BEAST v2. The ‘Subsampled data’
# correspond to three subsets of 20 files where the different isolates found in 
# each species were randomly chosen to be part of the new data set. Subsample A
# corresponds to isolates sampled from five elk (‘Elk’), five randomly chosen cattle 
# (‘Cattle’), and five randomly chosen white-tailed deer (WTD; ‘Deer’); Subsample B
# corresponds to isolates sampled from five elk, nine cattle, and nine randomly 
# chosen WTD; and Subsample C corresponds to isolates sampled from five elk, nine 
# cattle, and twenty four randomly chosen WTD. The ‘All data’ correspond to the 
# posterior mean of each host species interaction output by one DATM analysis using 
# all of the observed data, which consists of five elk, nine cattle, and thirty-nine 
# WTD. 
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
# Output directories:
# downsampling/downsamplingA/downsampling_run_i/, 1<=i<=20
# downsampling/downsamplingB/downsampling_run_i/, 1<=i<=20
# downsampling/downsamplingC/downsampling_run_i/, 1<=i<=20
#
# output files:
# downsampling/downsamplingA/downsampling_run_i/downsampling_run.xml', 1<=i<=20
# downsampling/downsamplingB/downsampling_run_i/downsampling_run.xml', 1<=i<=20
# downsampling/downsamplingC/downsampling_run_i/downsampling_run.xml', 1<=i<=20
#
# Figure: prob_transition_species_downsampling.png
####################################################################################


# number of data for each species from a total of elk=5, deer=39, cattle=9
nelk <- 5
ndeer <- 5
ncattle <-5

# read traits file
traits <- read.csv('Data/MI_Traits.csv', header=TRUE, sep=",", stringsAsFactors = default.stringsAsFactors())
df <- traits[, c("traits", "Species")]
names(df) <- c('Isolate', 'Species')
deer <- df[df$Species=='Deer',]
cattle <- df[df$Species=='Cattle',]
total.deer <- dim(deer)[1]
total.cattle <- dim(cattle)[1]

# loop to create 20 new data sets with only 5 isolates from each host species
# 5 elk, 5 random cattle, 5 random deer
indexes <- c(1:20)

for (j in 1:length(indexes)){
deer.to.delete <- as.character(sample(deer$Isolate, total.deer-ndeer))
cattle.to.delete <- as.character(sample(cattle$Isolate, total.cattle-ncattle))

# sequence of isolates to delete from deer and cattle
seq.del <- c(deer.to.delete, cattle.to.delete)

#  eliminate lines with the sequences to del seq.del
text <- readLines('Data/MI_Sequences.xml', n=-1)
lines <- c()
for(i in 1:length(seq.del)){
  lines <- c(which(grepl(seq.del[i], text)), lines)
}
text <- text[-lines]

# delete commas after sequence id in taxa id and traitset
comma3 <- which(grepl("taxa id", text))
aux <- which(grepl(",", text[comma3-1]))

if (length(aux) != 0){
  text[comma3-1] <- substr(text[comma3-1], 1, nchar(text[comma3-1])-1)
}

comma4 <- which(grepl("</traitSet", text))
aux <- which(grepl(",", text[comma4-1]))

if (length(aux) != 0){
  text[comma4-1] <- substr(text[comma4-1], 1, nchar(text[comma4-1])-1)
}

# create downsampling output directories
dir.create('downsampling/downsamplingA')
dir.create(paste0('downsampling/downsamplingA/downsampling_run_',as.character(indexes[j])))
writeLines(text, paste0('downsampling/downsamplingA/downsampling_run_',as.character(indexes[j]),'/downsampling_run.xml'))
}



#################
# subsampling B #
#################

# number of data for each species from a total of elk=5, deer=39, cattle=9
nelk <- 5
ndeer <- 9
ncattle <-9

# read traits file
traits <- read.csv('Data/MI_Traits.csv', header=TRUE, sep=",", stringsAsFactors = default.stringsAsFactors())
df <- traits[, c("traits", "Species")]
names(df) <- c('Isolate', 'Species')
deer <- df[df$Species=='Deer',]
cattle <- df[df$Species=='Cattle',]
total.deer <- dim(deer)[1]
total.cattle <- dim(cattle)[1]

indexes <- c(1:20)
for (j in 1:length(indexes)){
  deer.to.delete <- as.character(sample(deer$Isolate, total.deer-ndeer))
  cattle.to.delete <- as.character(sample(cattle$Isolate, total.cattle-ncattle))
  
  seq.del <- c(deer.to.delete, cattle.to.delete)
  
  #  eliminate lines with the sequences to del seq.del
  text <- readLines('Data/MI_Sequences.xml', n=-1)
  lines <- c()
  
  for(i in 1:length(seq.del)){
    lines <- c(which(grepl(seq.del[i], text)), lines)
  }
  text <- text[-lines]
  
  # delete commas after sequence id in taxa id and traitset
  comma3 <- which(grepl("taxa id", text))
  aux <- which(grepl(",", text[comma3-1]))
  
  if (length(aux) != 0){
    text[comma3-1] <- substr(text[comma3-1], 1, nchar(text[comma3-1])-1)
  }
  
  comma4 <- which(grepl("</traitSet", text))
  aux <- which(grepl(",", text[comma4-1]))
  
  if (length(aux) != 0){
    text[comma4-1] <- substr(text[comma4-1], 1, nchar(text[comma4-1])-1)
  }
  
  dir.create('downsampling/downsamplingB')
  dir.create(paste0('downsampling/downsamplingB/downsampling_run_',as.character(indexes[j])))
  writeLines(text, paste0('downsampling/downsamplingB/downsampling_run_',as.character(indexes[j]),'/downsampling_run.xml'))
}



#################
# Subsampling C #
#################

# number of data for each species availale from elk=5, deer=39, cattle=9
nelk <- 5
ndeer <- 24
ncattle <-9

# read traits file
traits <- read.csv('data/MI_Traits.csv', header=TRUE, sep=",", stringsAsFactors = default.stringsAsFactors())
df <- traits[, c("traits", "Species")]
names(df) <- c('Isolate', 'Species')
deer <- df[df$Species=='Deer',]
cattle <- df[df$Species=='Cattle',]
total.deer <- dim(deer)[1]
total.cattle <- dim(cattle)[1]

indexes <- c(1:20)

for (j in 1:length(indexes)){
  deer.to.delete <- as.character(sample(deer$Isolate, total.deer-ndeer))
  cattle.to.delete <- as.character(sample(cattle$Isolate, total.cattle-ncattle))
  seq.del <- c(deer.to.delete, cattle.to.delete)
  
  
  #  eliminate lines with the sequences to del seq.del
  text <- readLines('Data/MI_sequences.xml', n=-1)
  lines <- c()
  for(i in 1:length(seq.del)){
    lines <- c(which(grepl(seq.del[i], text)), lines)
  }
  text <- text[-lines]
  
  # delete commas after sequence id in taxa id and traitset
  comma3 <- which(grepl("taxa id", text))
  aux <- which(grepl(",", text[comma3-1]))
  
  if (length(aux) != 0){
    text[comma3-1] <- substr(text[comma3-1], 1, nchar(text[comma3-1])-1)
  }
  
  comma4 <- which(grepl("</traitSet", text))
  aux <- which(grepl(",", text[comma4-1]))
  
  if (length(aux) != 0){
    text[comma4-1] <- substr(text[comma4-1], 1, nchar(text[comma4-1])-1)
  }
  
  dir.create('downsampling/downsamplingC')
  dir.create(paste0('downsampling/downsamplingC/downsampling_run_',as.character(indexes[j])))
  writeLines(text, paste0('downsampling/downsamplingC/downsampling_run_',as.character(indexes[j]),'/downsampling_run.xml'))
}



#################################
# these will be run by BEAST.....
#################################


########################################################################
# Input files (BEAST output): 
# downsampling/A/downsampling_run_",i,"/downsampling_run.log". 1<=i<=20 
# downsampling/B/downsampling_run_",i,"/downsampling_run.log". 1<=i<=20 
# downsampling/C/downsampling_run_",i,"/downsampling_run.log". 1<=i<=20 
########################################################################


#########################
# Plot Figure 6
########################


subtable <- {}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# read downsampling files 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
for (i in 1:20){
  print(i)
  table555 <- read.table(paste0("downsampling/A/downsampling_run_",i,"/downsampling_run.log"), header=TRUE, skip=2)
  table555$SACountFBD <- NULL
  table555$data <- "subsampling A"
  table555$i <- i
  
  table599 <- read.table(paste0("downsampling/B/downsampling_run_",i,"/downsampling_run.log"), header=TRUE, skip=2)
  table599$SACountFBD <- NULL
  table599$data <- "subsampling B"
  table599$i <- i
  
  table5924 <- read.table(paste0("downsampling/C/downsampling_run_",i,"/downsampling_run.log"), header=TRUE, skip=2)
  table5924$SACountFBD <- NULL
  table5924$data <- "subsampling C"
  table5924$i <- i
  
  subtable <- rbind(subtable, table555, table599, table5924)
}

rateIndicators <- summaryBy(rateIndicator.species2+rateIndicator.species1+rateIndicator.species3~(i+data), data=subtable, FUN=mean)

names(rateIndicators) <- c("i","data","Cattle-Elk", "Cattle-Deer", "Deer-Elk")
rateIndicators.melt <- melt(rateIndicators, id=c("data","i"))
rateIndicators.melt <- rateIndicators.melt[, c("value", "variable", "i", "data" )]
names(rateIndicators.melt) <- c("mean", "interaction", "i", "data")

point <- data.frame(variable = c("Cattle-Elk", "Cattle-Deer", "Deer-Elk"), value = c(0.9691081,0.6030079,0.9871876))
cbPalette <- c("#000000","#808080", "#C0C0C0")


png(file="prob_transition_species_downsampling.png", width=7, height=4.5, units="in", res=300)
ggplot(rateIndicators.melt, aes(x=interaction, y=mean, cond=interaction, fill = data, colour=interaction)) + 
  geom_boxplot()+
  scale_x_discrete(name = "Symmetric transition between host species") +
  scale_y_continuous(name = "Estimated posterior probability") +
  scale_fill_brewer(palette = "Reds")+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=8))+
  scale_colour_manual(values=cbPalette) +
  annotate("text", x = as.numeric(point$variable)+0.50, y = point$value, label = "*", color="red", size=6)
dev.off()





