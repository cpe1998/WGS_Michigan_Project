##############################################################################
# 
# This code parses the XML file to randomly downsample the number of isolates
# 
# It generates Figure 6 of paper in (1).
#
# Figure caption: Comparison of the estimated posterior support of direct host species transition between subsampled and observed data. The estimated posterior mean probability of each interaction is the posterior probability that a particular transition rate is positive. If this probability is high, then the data strongly support a model in which there is direct pathogen transition between that particular pair of host species. The posterior means were estimated via a Discrete Ancestral Trait Mapping performed in BEAST v2. The ‘Subsampled data’ correspond to three subsets of 10 files where the different isolates found in each species were randomly chosen to be part of the new data set. Subsample A corresponds to isolates sampled from five elk (‘Elk’), five randomly chosen cattle (‘Cattle’), and five randomly chosen deer (‘Deer’); Subsample B corresponds to isolates sampled from five elk, nine cattle, and nine randomly chosen deer; and Subsample C corresponds to isolates sampled from five elk, nine cattle, and twenty four randomly chosen deer. The ‘All data’ correspond to the posterior mean of each host species interaction output by one DATM analysis using all of the observed data, which consists of five elk, twelve cattle, and 117 deer. 
#
# (1) Manuscript "Implications for disease management at the wildlife-livestock 
# interface: using whole-genome sequencing to study the role of elk in Mycobacterium bovis 
# transmission in Michigan, USA" by L.C.M. Salvador, D.J. O’Brien, M.K. Cosgrove, T.P. Stuber, 
# A. Schooley, J. Crispell, S. Church, Y.T., Grohn, S. Robbe-Austerman, R.R. Kao
# 
# @developed by lcmsalvador, July 2017
# @updated by lcmsalvador, November 2018
#
# Input files:
# 1. Data/MI_Elk_Data_134isolates_Traits_withClades_OnlineVersion.csv
# 2. Data/MI_Elk_134isolates_HKY_relExp_extskyline_DTA.xml
#  
# # Output directories:
# downsampling/downsamplingA/downsampling_run_i/, 1<=i<=10
# downsampling/downsamplingB/downsampling_run_i/, 1<=i<=10
# downsampling/downsamplingC/downsampling_run_i/, 1<=i<=10
# downsampling/downsamplingD/downsampling_run_i/, 1<=i<=10
#
# output files:
# downsampling/downsamplingA/downsampling_run_i/downsampling_run.xml', 1<=i<=10
# downsampling/downsamplingB/downsampling_run_i/downsampling_run.xml', 1<=i<=10
# downsampling/downsamplingC/downsampling_run_i/downsampling_run.xml', 1<=i<=10
# downsampling/downsamplingC/downsampling_run_i/downsampling_run.xml', 1<=i<=10
#
# Figures: MI_Figure_6.png
#          MI_Figure_6.pdf    
####################################################################################
######################################################################
# number of isolates
nisol <- 134

# number of iterations
n <- 10

# Downsampling A
# number of data for each species (elk=5, deer=117, cattle=12)
nelk <- 5
ndeer <- 5
ncattle <-5

# read traits file
traits <- read.csv('Data/MI_Elk_134isolates_Traits_withClades_OnlineVersion.csv', header=TRUE, sep=",", stringsAsFactors = default.stringsAsFactors())
df <- traits[, c("ID", "SPECIES")]
names(df) <- c('Isolate', 'Species')
deer <- df[df$Species=='DEER',]
cattle <- df[df$Species=='CATTLE',]
total.deer <- dim(deer)[1]
total.cattle <- dim(cattle)[1]

indexes <- c(1:n)
for (j in 1:length(indexes)){
deer.to.delete <- as.character(sample(deer$Isolate, total.deer-ndeer))
cattle.to.delete <- as.character(sample(cattle$Isolate, total.cattle-ncattle))

seq.del <- c(deer.to.delete, cattle.to.delete)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  eliminate lines with the sequences to del seq.del
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

text <- readLines('Data/MI_Elk_extra_isolates_HKY_relExp_extskyline_DTA.xml', n=-1)
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

dir.create('downsampling/downsamplingA')
dir.create(paste0('downsampling/downsamplingA/downsampling_run_',as.character(indexes[j])))
writeLines(text, paste0('downsampling/downsamplingA/downsampling_run_',as.character(indexes[j]),'/downsampling_run.xml'))
}


################
# downsamplingB
################

nelk <- 5
ndeer <- 12
ncattle <-12

# number of data for each species (elk=5, deer=39, cattle=9)
df <- traits[, c("ID", "SPECIES")]
names(df) <- c('Isolate', 'Species')
deer <- df[df$Species=='DEER',]
cattle <- df[df$Species=='CATTLE',]
total.deer <- dim(deer)[1]
total.cattle <- dim(cattle)[1]

indexes <- c(1:10)
for (j in 1:length(indexes)){
  deer.to.delete <- as.character(sample(deer$Isolate, total.deer-ndeer))
  cattle.to.delete <- as.character(sample(cattle$Isolate, total.cattle-ncattle))
  
  seq.del <- c(deer.to.delete, cattle.to.delete)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #  eliminate lines with the sequences to del seq.del
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  text <- readLines('Data/MI_Elk_extra_isolates_HKY_relExp_extskyline_DTA.xml', n=-1)
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
# Downsampling C
#################
# number of data for each species (elk=5, deer=117, cattle=12)
nelk <- 5
ndeer <- 36
ncattle <-12
df <- traits[, c("ID", "SPECIES")]
names(df) <- c('Isolate', 'Species')
deer <- df[df$Species=='DEER',]
cattle <- df[df$Species=='CATTLE',]
total.deer <- dim(deer)[1]
total.cattle <- dim(cattle)[1]

indexes <- c(1:10)
for (j in 1:length(indexes)){
  deer.to.delete <- as.character(sample(deer$Isolate, total.deer-ndeer))
  cattle.to.delete <- as.character(sample(cattle$Isolate, total.cattle-ncattle))
  
  seq.del <- c(deer.to.delete, cattle.to.delete)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #  eliminate lines with the sequences to del seq.del
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  text <- readLines('Data/MI_Elk_extra_isolates_HKY_relExp_extskyline_DTA.xml', n=-1)
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


#################
# Downsampling D
#################

# number of data for each species (elk=5, deer=117, cattle=12)
nelk <- 5
ndeer <- 76
ncattle <-12

df <- traits[, c("ID", "SPECIES")]
names(df) <- c('Isolate', 'Species')
deer <- df[df$Species=='DEER',]
cattle <- df[df$Species=='CATTLE',]
total.deer <- dim(deer)[1]
total.cattle <- dim(cattle)[1]

indexes <- c(1:10)
for (j in 1:length(indexes)){
  deer.to.delete <- as.character(sample(deer$Isolate, total.deer-ndeer))
  cattle.to.delete <- as.character(sample(cattle$Isolate, total.cattle-ncattle))
  
  seq.del <- c(deer.to.delete, cattle.to.delete)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #  eliminate lines with the sequences to del seq.del
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  text <- readLines('Data/MI_Elk_extra_isolates_HKY_relExp_extskyline_DTA.xml', n=-1)
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
  
  dir.create('downsampling/downsamplingD')
  dir.create(paste0('downsampling/downsamplingD/downsampling_run_',as.character(indexes[j])))
  writeLines(text, paste0('downsampling/downsamplingD/downsampling_run_',as.character(indexes[j]),'/downsampling_run.xml'))
}


#################################
# these will be run by BEAST.....
#################################


########################################################################
# Input files (from BEAST output): 
# downsampling/downsamplingA/downsampling_run_",i,"/beast.log. 1<=i<=10 
# downsampling/downsamplingB/downsampling_run_",i,"/beast.log. 1<=i<=10 
# downsampling/downsamplingC/downsampling_run_",i,"/beast.log. 1<=i<=10 
# downsampling/downsamplingD/downsampling_run_",i,"/beast.log. 1<=i<=10 
########################################################################

#########################
# Plot Figure 6
########################
# Rate indicators of the original run
cattle_deer <- 0.996
cattle_elk <- 0.391
deer_elk <- 0.989

# Number of iterations
n <- 10
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# read downsampling files 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
subtable <- {}
for (i in 1:n){
  print(i)
  table555 <- read.table(paste0("downsampling/downsamplingA/downsampling_run_",i,"/beast.log"), header=TRUE, skip=2)
  aux1 <- table555[, c("rateIndicator.Host1", "rateIndicator.Host2", "rateIndicator.Host3")]
  aux1$data <- "subsampling A"
  aux1$i <- i
 
  table599 <- read.table(paste0("downsampling/downsamplingB/downsampling_run_",i,"/beast.log"), header=TRUE, skip=2)
  aux2 <- table599[, c("rateIndicator.Host1", "rateIndicator.Host2", "rateIndicator.Host3")]
  aux2$data <- "subsampling B"
  aux2$i <- i
  
  table5924 <- read.table(paste0("downsampling/downsamplingC/downsampling_run_",i,"/beast_", i, ".log"), header=TRUE, skip=2)
  aux3 <- table5924[, c("rateIndicator.Host1", "rateIndicator.Host2", "rateIndicator.Host3")]
  aux3$data <- "subsampling C"
  aux3$i <- i
  
  table5976 <- read.table(paste0("downsampling/downsamplingD/downsampling_run_",i,"/beast_", i, ".log"), header=TRUE, skip=2)
  aux4 <- table5976[, c("rateIndicator.Host1", "rateIndicator.Host2", "rateIndicator.Host3")]
  aux4$data <- "subsampling D"
  aux4$i <- i
  
  subtable <- rbind(subtable, aux1, aux2, aux3, aux4)
}

rateIndicators <- summaryBy(rateIndicator.Host2+rateIndicator.Host1+rateIndicator.Host3~(i+data), data=subtable, FUN=mean)

names(rateIndicators) <- c("i","data","Cattle-Elk", "Cattle-Deer", "Deer-Elk")
rateIndicators.melt <- melt(rateIndicators, id=c("data","i"))
rateIndicators.melt <- rateIndicators.melt[, c("value", "variable", "i", "data" )]
names(rateIndicators.melt) <- c("mean", "interaction", "i", "data")

point <- data.frame(variable = c("Cattle-Elk", "Cattle-Deer", "Deer-Elk"), value = c(cattle_deer,cattle_elk,deer_elk))

cbPalette <- c("#000000","#808080", "#C0C0C0")


# generage png file
png(file="Figure6.png", width=7, height=4.5, units="in", res=300)
ggplot(rateIndicators.melt, aes(x=interaction, y=mean, fill = data)) + 
  geom_boxplot()+
  scale_x_discrete(name = "Symmetric transition between host species") +
  scale_y_continuous(name = "Estimated posterior probability") +
  scale_fill_brewer(palette = "Reds")+
  theme_bw(
    #legend.title = element_blank(),
    legend.text=element_text(size=10))+
  scale_colour_manual(values=cbPalette) +
  annotate("text", x = as.numeric(point$variable)+0.50, y = point$value, label = "*", color="darkred", size=12)
dev.off()


# generage pdf file
pdf(file="Figure6.pdf", width=7, height=4.5)
ggplot(rateIndicators.melt, aes(x=interaction, y=mean, fill = data)) + 
  geom_boxplot()+
  scale_x_discrete(name = "Symmetric transition between host species") +
  scale_y_continuous(name = "Estimated posterior probability") +
  scale_fill_brewer(palette = "Reds")+
  theme_classic()+
  scale_colour_manual(values=cbPalette) +
  annotate("text", x = as.numeric(point$variable)+0.50, y = point$value, label = "*", color="darkred", size=12)
dev.off()





