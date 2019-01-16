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
amniote_mass = amniote[,c(2,6)]
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
ben_mass = benoit[,c(1,5)]
# change mass column to specify the data
colnames(ben_mass) [2] = "body_mass_ben"

mass2 = merge(mass, ben_mass, by = 'Tips', all.x = T)

# if am not availiable use ben
mass2$body_mass = mass2$body_mass_am
mass2$body_mass[is.na(mass2$body_mass)] = mass2$body_mass_ben[is.na(mass2$body_mass)]

# keep only matched name and mass columns
mass = mass2[,c(1:3,6)]

# plot the data and visually check for outliers
ggplot(mass2, aes(x=matched_name, y=body_mass)) + geom_point() 
ggplot(mass2, aes(x=matched_name, y=body_mass)) + geom_point() + scale_y_log10()

# These species need checking
# Casuarius bennetti. our mass 30kgs looks more likely to be 18kgs

### LONGEVITY ###
# maximum recorded
long_amniote = amniote[,c(2,7)]
long_hbw = hbw[,c(1,8)]
long_benoit = benoit[,c(1,6)]

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

ggplot(zz, aes(x=matched_name, y=value)) + geom_point(aes(colour = variable)) 
ggplot(zz, aes(x=matched_name, y=log(value))) + geom_point(aes(colour = variable)) 

# species to check
# Apteryx owenii
# Grus rubicunda
# Grus virgo

# take the maximum value for each species
zz[,"max_longevity"] = apply(zz[,4:6], 1, max, na.rm=T)
zz[zz == (-Inf)] = NA

longevity = zz[,c(1,7)]

### AGE AT FIRST BREEDING ###
# amniote, then benoit, then hbw

# only hbw data has age at first breeding (in months)
# the others have age of female maturity maturity
amniote_breeding = amniote[,c(2,3)]
hbw_breeding = hbw[,c(1,16:18)]
# change measurement into days
hbw_breeding$first_breeding_min_hbw = (hbw_breeding$firstbreeding_minmonths/12)*365
hbw_breeding$first_breeding_max_hbw = (hbw_breeding$firstbreeding_maxmonths/12)*365
hbw_breeding$first_breeding_mean_hbw = (hbw_breeding$firstbreeding_meanmonths/12)*365
# calculate mean 
hbw_breeding = hbw_breeding[,c(1,5:7)]
hbw_breeding$computed_mean = (hbw_breeding$first_breeding_min_hbw + hbw_breeding$first_breeding_max_hbw) / 2

benoit_breeding = benoit[,c(1,7)]

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

# check these species
# Phoebastria nigripes
# Phoebastria immutabilis
# Procellaria cinerea
# Cathartes aura
# Diomedea chrysostoma
# Coturnix japonica
# Taeniopygia guttata

### CLUTCH SIZE ###
# amniote data, then hbw (not avail Benoit)

amniote_clutch_size = amniote[,c(2,4)]
hbw_clutch_size = hbw[,c(1,9:11)]

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

# species to check
# Perdix dauurica
# Rhea americana

### CLUTCHES PER YEAR ###
# amniote data, then hbw
clutches_year_am = amniote[,c(2,5)]
clutches_year_hbw = hbw[,c(1,12)]

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

# species to check
# Vidua chalybeata


## merge together all LHTs by "accepted_name"
X = list(mass, longevity, breeding_age, clutch_size, clutches_year)
bird_LHT = Reduce(function(x,y) merge(x,y,by="Tips",all=T),X)

bird_LHT$fecundity = bird_LHT$litter_or_clutch_size_n * bird_LHT$clutches_per_year

write.csv(bird_LHT, file = "traits.csv")

check = bird_LHT[complete.cases(bird_LHT),] ## just to check how many species have complete data













#####################################
bird_LHT2 = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/traits.csv", stringsAsFactors = F)

library(ggplot2)
ggplot(bird_LHT, aes(x=body_mass, y=pop_total_mean)) + geom_point(aes(colour = factor(flighted.flightless))) + scale_x_log10() + scale_y_log10() + geom_smooth(method=lm, se=FALSE)
ggplot(bird_LHT, aes(x=body_mass, y=distribution)) + geom_point() + scale_x_log10() + scale_y_log10() + geom_smooth(method=lm, se=FALSE)
ggplot(bird_LHT, aes(x=body_mass, y=egg_mass)) + geom_point(aes(colour = factor(flighted.flightless))) + scale_x_log10() + scale_y_log10()
ggplot(bird_LHT, aes(x=body_mass, y=max_longevity)) + geom_point(aes(colour = factor(flighted.flightless))) + scale_x_log10() + geom_smooth(method=lm, se=FALSE)
ggplot(bird_LHT, aes(x=body_mass, y=breeding_age)) + geom_point() + scale_x_log10() + geom_smooth(method="lm", formula= (y ~ exp(x)), se=FALSE, linetype = 1)
ggplot(bird_LHT, aes(x=body_mass, y=clutch_size)) + geom_point() + scale_x_log10() + geom_smooth(method=lm, se=FALSE)
ggplot(bird_LHT, aes(x=body_mass, y=clutches_year)) + geom_point(aes(colour = factor(flighted.flightless))) + scale_x_log10() + geom_smooth(method=lm, se=FALSE)

ggplot(bird_LHT, aes(x=egg_mass, y=clutches_year)) + geom_point(aes(colour = factor(flighted.flightless))) + scale_x_log10() + geom_smooth(method=lm, se=FALSE)

ggplot(bird_LHT, aes(x=clutch_size, y=clutches_year)) + geom_point() + geom_smooth(method=lm, se=FALSE) + scale_y_log10()

ggplot(bird_LHT, aes(x=egg_mass, y=clutch_size)) + geom_point(aes(colour = factor(flighted.flightless))) + scale_x_log10() + geom_smooth(method=lm, se=FALSE)
ggplot(bird_LHT, aes(x=egg_mass, y=breeding_age)) + geom_point() + scale_x_log10() + geom_smooth(method="lm", formula= (y ~ exp(x)), se=FALSE, linetype = 1)
ggplot(bird_LHT, aes(x=clutch_size, y=breeding_age)) + geom_point() + scale_x_log10() + geom_smooth(method="lm", formula= (y ~ exp(-x)), se=FALSE, linetype = 1)

ggplot(bird_LHT, aes(x=clutch_size, y=max_longevity)) + geom_point() + scale_x_log10() + geom_smooth(method=lm, se=FALSE)
ggplot(bird_LHT, aes(x=egg_mass, y=max_longevity)) + geom_point(aes(colour = factor(flighted.flightless))) + scale_x_log10() + geom_smooth(method=lm, se=FALSE)
ggplot(bird_LHT, aes(x=breeding_age, y=max_longevity)) + geom_point(aes(colour = factor(flighted.flightless))) + scale_x_log10() + geom_smooth(method=lm, se=FALSE)


ggplot(bird_LHT, aes(x=body_mass, y=accepted_name)) + geom_point(aes(colour = factor(flighted.flightless))) + scale_x_log10() 
ggplot(bird_LHT, aes(x=max_longevity, y=accepted_name)) + geom_point(aes(colour = factor(flighted.flightless)))  
ggplot(bird_LHT, aes(x=egg_mass, y=accepted_name)) + geom_point(aes(colour = factor(flighted.flightless))) + scale_x_log10() 

## tells you how correlated population size is with distribution
ggplot(bird_LHT, aes(x=distribution, y=pop_total_mean)) + geom_point(aes(colour = factor(flighted.flightless))) + scale_x_log10() + scale_y_log10() + geom_smooth(method=lm, se=FALSE)
summary(lm(distribution ~ pop_total_mean, data=bird_LHT))
#### should repeat this correcting for phylogenetic relatedness

# linear model of egg and body mass (forced through the origin)
summary(lm(egg_mass ~ 0 + body_mass, data=bird_LHT))

summary(lm(body_mass ~ pop_total_mean, data=bird_LHT))




###############
#OTHER TRAITS (NOT USING)

### FLIGHTED/NOT ###
flighted = hbw[,c(1,3)]

### SURELY THERE IS A BETTER WAY TO DO THIS??? ###
# add in missing values
flighted <- within(flighted, flighted.flightless[matched_name == "Amazilia versicolor"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Anhinga rufa"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Ardea alba"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Buteo burmanicus"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Butorides striata"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Chroicocephalus ridibundus"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Agelasticus cyanopus"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Agelasticus xanthophthalmus"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Hesperiphona vespertina"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Crossoptilon harmani"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Thalassarche melanophris"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Garrulax cineraceus"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Gavia pacifica"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Grus paradisea"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Grus virgo"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Chlorodrepanis flava"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Akialoa stejnegeri"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Chlorodrepanis virens"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Himantopus mexicanus"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Leiothrix argentauris"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Lophura swinhoii"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Tarsiger cyanurus"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Morus serrator"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Oedistoma iliolophum"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Pardaliparus venustulus"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Phaethon rubricauda"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Phoebastria immutabilis"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Phoebastria nigripes"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Porphyrio hochstetteri"] <- 0)
flighted <- within(flighted, flighted.flightless[matched_name == "Primolius couloni"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Pygoscelis antarcticus"] <- 0)
flighted <- within(flighted, flighted.flightless[matched_name == "Saxicola rubicola"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Stercorarius maccormicki"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Thalasseus acuflavidus"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Sylvia crassirostris"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Tetrastes sewerzowi"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Urocissa erythrorhyncha"] <- 1)
flighted <- within(flighted, flighted.flightless[matched_name == "Drepanis coccinea"] <- 1)



### POPULATION ###
# birdlife only
# no data availiable in either of the other datasets
population = birdlife[,c(1,3:5)]

population$pop_computed_mean = (population$pop_min + population$pop_max) / 2

population$pop_total_mean = population$pop_mean
population$pop_total_mean[is.na(population$pop_total_mean)] = population$pop_computed_mean[is.na(population$pop_total_mean)]
population = population[,c(1:3,6)]

# plot the data and visually check for outliers
ggplot(population, aes(x=matched_name, y=pop_total_mean)) + geom_point() 
ggplot(population, aes(x=matched_name, y=log(pop_total_mean))) + geom_point() 

# These species need checking
# Alauda arvensis. our population 6E08 looks more likely to be smaller
# Corvus frugilegus. our population 1.3E08 more likely smaller
# Loxia curvirostra
# Ara glaucogularis

### DISTRIBUTION ###
# birdlife only
# no data availiable in either of the other datasets
distribution = birdlife[,c(1,2)]

# plot the data and visually check for outliers
ggplot(distribution, aes(x=matched_name, y=distribution)) + geom_point() 
ggplot(distribution, aes(x=matched_name, y=log(distribution))) + geom_point() 

# These species need checking
# Ardea cinerea
# Phaethon rubricauda
# Procellaria cinerea
# Tyto alba

# Pseudonestor xanthophrys
# Phalacrocorax chalconotus
# Loxops caeruleirostris

### EGG MASS ###
# amniote first, then hbw (no data in Benoit)
amniote_egg = amniote[,c(1,12)]

# add in the HBW data (min, max and mean)
hbw_egg = hbw[,c(1,13:15)]

egg = merge(amniote_egg, hbw_egg, by="matched_name", all=T)
colnames(egg) [2] = "egg_mass_amniote"

# create new column to calculate the mean of min and max columns
egg$computed_mean = (egg$eggs_mass_ming + egg$eggs_mass_maxg) / 2

# if mean data is not available from hbw, use computed mean data
egg$eggs_mass_meang[is.na(egg$eggs_mass_meang)] = egg$computed_mean[is.na(egg$eggs_mass_meang)]

# if amniote data is not available, use hbw mean data
egg$egg_mass = egg$egg_mass_amniote
egg$egg_mass[is.na(egg$egg_mass)] = egg$eggs_mass_meang[is.na(egg$egg_mass)]

# plot data to check for outliers
egg1 = egg[,c(1,2,5)]
data1 = melt(egg1, id = c("matched_name"))
ggplot(data1, aes(x=matched_name, y=value)) + geom_point(aes(colour = variable)) 
ggplot(data1, aes(x=matched_name, y=log(value))) + geom_point(aes(colour = variable)) 

# values to check
# Casuarius casuarius
# Rhea americana

# Aegithalos caudatus
# Aegithalos glaucogularis
# Archilochus colubris

egg = egg[,c(1,7)]
