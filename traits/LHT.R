## a script for organising trait data from 3 sources
# amniote, HBW and a Nabholz paper
# use ggplot to visually assess data for outliers
library(ggplot2)
library(reshape2)

# read in the data

# read in amniote data with added rows
amniote = read.csv("path/to/Amniote_accepted_names.csv", stringsAsFactors = F)
amniote[amniote == (-999)] = NA
benoit = read.csv("path/to/benoit.csv", stringsAsFactors = F)
hbw = read.csv(file = "path/to/HBW.csv", stringsAsFactors = F)
species = read.csv("path/to/taxon_names.csv", stringsAsFactors = F)

### MASS ###

# keep only mass and species names columns
amniote_mass = amniote[,c('species','adult_body_mass_g')]
# change mass column to specify the data
colnames(amniote_mass) [2] = "body_mass_am"

# use the 'species' dataset to merge and subset data
x = merge(species, amniote_mass, by = "matched_name", all.x = T)
# some species match 'species' name column not 'accepted name' column in amniote database
colnames(amniote_mass)[1] = 'species'
colnames(amniote_mass)[2] = 'mass2'
mass = merge(x, amniote_mass, by = 'species', all.x = T)

# move 'mass2' values to mass amniote column
mass$body_mass_am[is.na(mass$body_mass_am)] = mass$mass2[is.na(mass$body_mass_am)]
mass = mass[,c(1:4)]

# keep only mass and species names columns
ben_mass = benoit[,c('matched_name', 'Mass')]
# change mass column to specify the data
colnames(ben_mass) [2] = "body_mass_ben"

mass2 = merge(mass, ben_mass, by = 'Tips', all.x = T)

# if amniote not availiable use benoit
mass2$body_mass = mass2$body_mass_am
mass2$body_mass[is.na(mass2$body_mass)] = mass2$body_mass_ben[is.na(mass2$body_mass)]

# keep only matched name and mass columns
mass = mass2[,c(1:3,6)]

# plot the data and visually check for outliers
ggplot(mass2, aes(x=matched_name, y=body_mass)) + geom_point() 
ggplot(mass2, aes(x=matched_name, y=body_mass)) + geom_point() + scale_y_log10()

### LONGEVITY ###
# maximum recorded
long_amniote = amniote[,c('species','longevity_y')]
long_hbw = hbw[,c('species','longevity_maxrecorded')]
long_benoit = benoit[,c('species','Longevity_Species')]

# change column names to represent data source
colnames(long_amniote) [2] = "long_amniote"
colnames(long_hbw) [2] = "long_hbw"
colnames(long_benoit) [2] = "long_benoit"

## amniote data
x = merge(species, long_amniote, by = 'matched_name', all.x = T)
colnames(long_amniote)[1] = 'species'
colnames(long_amniote)[2] = 'long2'
y = merge(x, long_amniote, by = 'species', all.x = T)
y$long_amniote[is.na(y$long_amniote)] = y$long2[is.na(y$long_amniote)]
y = y[,c(1:4)]

## hbw data
z = merge(y, long_hbw, by = 'species', all.x = T)
long_hbw2 = hbw[,c(1,2,8)]
long_hbw2 = long_hbw2[complete.cases(long_hbw2[c(1,2)]),]
long_hbw2 = long_hbw2[,c(2,3)]
xx = merge(z, long_hbw2, by = 'matched_name', all.x = T)

xx$long_hbw[is.na(xx$long_hbw)] = xx$longevity_maxrecorded[is.na(xx$long_hbw)]
xx = xx[,c(1:5)]

## benoit data
zz = merge(xx, long_benoit, by = 'Tips', all.x = T)

# visualise data
ggplot(zz, aes(x=matched_name, y=value)) + geom_point(aes(colour = variable)) 
ggplot(zz, aes(x=matched_name, y=log(value))) + geom_point(aes(colour = variable)) 

# take the maximum value for each species
zz[,"max_longevity"] = apply(zz[,4:6], 1, max, na.rm=T)
zz[zz == (-Inf)] = NA

longevity = zz[,c(1,7)]

### AGE AT FIRST BREEDING ###

# only hbw data has age at first breeding (in months)
# the others have age (in days)
amniote_breeding = amniote[,c('species','female_maturity_d')]

hbw_breeding = hbw[,c('species','firstbreeding_minmonths', 'firstbreeding_maxmonths', 'firstbreeding_meanmonths')]
# change measurement into days
hbw_breeding$first_breeding_min_hbw = (hbw_breeding$firstbreeding_minmonths/12)*365
hbw_breeding$first_breeding_max_hbw = (hbw_breeding$firstbreeding_maxmonths/12)*365
hbw_breeding$first_breeding_mean_hbw = (hbw_breeding$firstbreeding_meanmonths/12)*365
# calculate mean 
hbw_breeding = hbw_breeding[,c(1,5:7)]
hbw_breeding$computed_mean = (hbw_breeding$first_breeding_min_hbw + hbw_breeding$first_breeding_max_hbw) / 2

benoit_breeding = benoit[,c('matched_name','SexMatF_Species')]

# change column names to represent data source
colnames(amniote_breeding) [2] = "breeding_amniote"
colnames(benoit_breeding) [2] = "breeding_benoit"

## amniote
x = merge(species, amniote_breeding, by = 'matched_name', all.x = T)
colnames(amniote_breeding)[1] = 'species'
colnames(amniote_breeding)[2] = 'gt2'
y = merge(x, amniote_breeding, by = 'species', all.x = T)
y$breeding_amniote[is.na(y$breeding_amniote)] = y$gt2[is.na(y$breeding_amniote)]
y = y[,c(1:4)]

## benoit
z = merge(y, benoit_breeding, by = 'Tips', all.x = T)
z$breeding_amniote[is.na(z$breeding_amniote)] = z$breeding_benoit[is.na(z$breeding_amniote)]
z = z[,c(1:4)]

## hbw
xx = merge(z, hbw_breeding, by = 'species', all.x =T)

# if mean data for hbw not available, use computed mean
xx$mean_hbw = xx$first_breeding_mean_hbw
xx$mean_hbw[is.na(xx$mean_hbw)] = xx$computed_mean[is.na(xx$mean_hbw)]
xx = xx[,c(1:4,9)]

xx$breeding_amniote[is.na(xx$breeding_amniote)] = xx$mean_hbw[is.na(xx$breeding_amniote)]
breeding_age = xx[,c(2,4)]
colnames(breeding_age)[2] = 'age_first_breeding'

# visually check for outliers
breeding_age1 = breeding_age[,c(1,2,3,8)]
data1 = melt(breeding_age1, id = c("matched_name"))
ggplot(data1, aes(x=matched_name, y=value)) + geom_point(aes(colour = variable)) 
ggplot(data1, aes(x=matched_name, y=log(value))) + geom_point(aes(colour = variable)) 

### CLUTCH SIZE ###
# amniote data, then hbw (not avail Benoit)

amniote_clutch_size = amniote[,c('species', 'litter_or_clutch_size_n')]
hbw_clutch_size = hbw[,c('species', 'clutchsize_min', 'clutchsize_max', 'clutchsize_mean.usual')]

# compute mean of min and max for hbw
hbw_clutch_size$computed_mean = (hbw_clutch_size$clutchsize_min + hbw_clutch_size$clutchsize_max) / 2
hbw_clutch_size$hbw_clutch_size = hbw_clutch_size$clutchsize_mean.usual
hbw_clutch_size$hbw_clutch_size[is.na(hbw_clutch_size$hbw_clutch_size)] = hbw_clutch_size$computed_mean[is.na(hbw_clutch_size$hbw_clutch_size)]
hbw_clutch_size = hbw_clutch_size[,c(1,6)]

# amniote
x = merge(species, amniote_clutch_size, by = 'matched_name', all.x = T)
colnames(amniote_clutch_size)[1] = 'species'
colnames(amniote_clutch_size)[2] = 'cs'
y = merge(x, amniote_clutch_size, by = 'species', all.x = T)
y$litter_or_clutch_size_n[is.na(y$litter_or_clutch_size_n)] = y$cs[is.na(y$litter_or_clutch_size_n)]
y = y[,c(1:4)]

## hbw
z = merge(y, hbw_clutch_size, by = 'species', all.x = T)
z$litter_or_clutch_size_n[is.na(z$litter_or_clutch_size_n)] = z$hbw_clutch_size[is.na(z$litter_or_clutch_size_n)]

clutch_size = z[,c(3,4)]

# visually check for outliers
clutch_size1 = clutch_size[,c(1,2,7)]
data1 = melt(clutch_size1, id = c("matched_name"))
ggplot(data1, aes(x=matched_name, y=value)) + geom_point(aes(colour = variable)) 
ggplot(data1, aes(x=matched_name, y=value)) + geom_point(aes(colour = variable)) + scale_y_log10()

### CLUTCHES PER YEAR ###
# amniote data, then hbw
clutches_year_am = amniote[,c('species','litters_or_clutches_per_y')]
clutches_year_hbw = hbw[,c('species','X.maxsucessful._clutches_per_year')]

colnames(clutches_year_am) [2] = "amniote_clutches"
colnames(clutches_year_hbw) [2] = "hbw_clutches"

# amniote
x = merge(species, clutches_year_am, by = 'matched_name', all.x = T)
colnames(clutches_year_am)[1] = 'species'
colnames(clutches_year_am)[2] = 'cy'
y = merge(x, clutches_year_am, by = 'species', all.x = T)
y$amniote_clutches[is.na(y$amniote_clutches)] = y$cy[is.na(y$amniote_clutches)]
y = y[,c(1:4)]

# hbw
z = merge(y, clutches_year_hbw, by = 'species', all.x = T)
z$amniote_clutches[is.na(z$amniote_clutches)] = z$hbw_clutches[is.na(z$amniote_clutches)]
z = z[,c(1:4)]

cs_hbw2 = hbw[,c(1,2,12)]
cs_hbw2 = cs_hbw2[complete.cases(cs_hbw2[c(1,2)]),]
cs_hbw2 = cs_hbw2[,c(2,3)]

xx = merge(z, cs_hbw2, by = 'matched_name', all.x = T)
xx$amniote_clutches[is.na(xx$amniote_clutches)] = xx$X.maxsuccesful._clutches_per_year[is.na(xx$amniote_clutches)]
xx = xx[,c(3,4)]

colnames(xx)[2] = 'clutches_per_year'

clutches_year = xx

# visually check for outliers
clutches_year1 = clutches_year[,c(1,2,3)]
data1 = melt(clutches_year1, id = c("matched_name"))
ggplot(data1, aes(x=matched_name, y=value)) + geom_point(aes(colour = variable)) 
ggplot(data1, aes(x=matched_name, y=value)) + geom_point(aes(colour = variable)) + scale_y_log10()


## merge together all LHTs by "accepted_name"
X = list(mass, longevity, breeding_age, clutch_size, clutches_year)
bird_LHT = Reduce(function(x,y) merge(x,y,by="Tips",all=T),X)

bird_LHT$fecundity = bird_LHT$litter_or_clutch_size_n * bird_LHT$clutches_per_year

write.csv(bird_LHT, file = "traits.csv")
