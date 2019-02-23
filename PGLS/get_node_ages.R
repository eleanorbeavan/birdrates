# A function for getting the age node of a tree
library(phytools)

t = read.tree('~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/calibration_2_new.tre')
t = read.tree('~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/calibration_sets_4.tre')
t = read.tree('~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/no_calibration.tre')
data = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/all_traits_genera/genera_mass_rates.csv", stringsAsFactors = F)
data = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/all_traits_genera/genera_gentime_rates.csv", stringsAsFactors = F)
data = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/all_traits_genera/genera_fecundity_rates.csv", stringsAsFactors = F)
data = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/all_traits_genera/genera_longevity_rates.csv", stringsAsFactors = F)

# get all pairs out of data frame
sister.pairs = list()
for (i in 1:nrow(data)) {
  x = as.character(data[i,1:2])
  sister.pairs = c(sister.pairs, list(x))
}

# extract node ages
node.ages = setNames(data.frame(matrix(ncol = 3, nrow = 1)), c('spp1', 'spp2', 'node.age'))
for (i in 1:nrow(data)) {
  x = sister.pairs[[i]]
  y = as.character(findMRCA(t, tips = x))
  age = branching.times(t)[[y]]
  node.ages = rbind(node.ages, c(x, age))
}

x = node.ages[,c(1,3)]
y = merge(data, x, by = 'spp1', all.x = T)

write.csv(y , file = '~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/ages_true_gnetime_cal2.csv')

#### ORDER ANALYSIS
## for order data
library(picante)
library(dplyr)

t = read.tree('~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/calibration_sets_2.tre')
t = read.tree('~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/calibration_sets_4.tre')
t = read.tree('~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/no_calibration.tre')

# get node ages for terminal branch for each species in the tree
phy.age = node.age(t)
BL.position = cbind(phy.age$edge,phy.age$age, t$edge.length)
dist.tip = max(phy.age$age)-BL.position[,3]
BL.positions = cbind(BL.position,dist.tip)
ages = BL.positions[,5]+BL.positions[,4]
BL.positions = cbind(BL.positions,ages)
node.ages = as.data.frame(BL.positions)
names(node.ages) = c("parental.node","daughter.node","dist.root","BL","dist.tip","mrca.age")
## node.ages is a data frame listing as variables the identity of parental and daughter nodes, the distance from the root and from the present of each node, the branch lenght and the age of the most recent common ancestor
species.ages = node.ages[node.ages[,2]<length(t$tip)+1,]

names = as.data.frame(cbind(t$tip.label, seq(1, 475)))
colnames(names) = c('Tips', 'daughter.node')

species.ages = merge(species.ages, names, by = 'daughter.node')
species.ages = species.ages[,c(6,7)]

data = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/order/rates_traits.csv", stringsAsFactors = F)

data.ages = merge(data, species.ages, by = 'Tips', all.x = T)

write.csv(data.ages, "~/Dropbox/PhD/bird_rates/processed_data/CSV_files/order/ages_cal2_order.csv")
