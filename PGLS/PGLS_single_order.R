# PGLS order analysis
library(ape)
library(phytools)
library(geiger)
library(caper)
library(ggplot2)
library(dplyr)
library(car)
library(cowplot)

# read in data
calibration.2 = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/order/ages_cal2_order.csv", stringsAsFactors = F)
calibration.4 = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/order/ages_cal4_order.csv", stringsAsFactors = F)
no.calibration = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/order/ages_nocal_order.csv", stringsAsFactors = F)

t2 = read.tree("~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/calibration_sets_2.tre")
t4 = read.tree("~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/calibration_sets_4.tre")
t = read.tree("~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/no_calibration.tre")

# divide dS and dN by node age
abs.ds.2 = calibration.2$ds/calibration.2$mrca.age
abs.dn.2 = calibration.2$dn/calibration.2$mrca.age
dnds = calibration.2$dn.ds

abs.ds.4 = calibration.4$ds/calibration.4$mrca.age
abs.dn.4 = calibration.4$dn/calibration.4$mrca.age

abs.ds.no = no.calibration$ds/no.calibration$mrca.age
abs.dn.no = no.calibration$dn/no.calibration$mrca.age

pgls.data.sorted = calibration.2[,c(1,6,10:12,15)]
pgls.data.sorted = cbind(pgls.data.sorted, abs.ds.2, abs.ds.4, abs.ds.no, abs.dn.2, abs.dn.4, abs.dn.no, dnds)

# need to subset the tree to include only the tips in the data
tips.to.drop = setdiff(t2$tip.label, calibration.2$Tips)
tree2 = drop.tip(t2, tips.to.drop)

# rows in dataframe need to be in the same order as tips in the tree
pgls.data.sorted = pgls.data.sorted[match(tree2$tip.label, calibration.2$Tips),]

# how many species for each trait?
x = calibration.2 %>% filter(!is.na(body_mass))
x = calibration.2 %>% filter(!is.na(max_longevity))
x = calibration.2 %>% filter(!is.na(age_first_breeding))
x = calibration.2 %>% filter(!is.na(fecundity))

# dS pgls model calibration set 2
pgls1 = comparative.data(tree2, pgls.data.sorted, Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds.2) ~ log(body_mass), pgls1, lambda = 'ML', kappa = 0.899, delta = 0.998)
summary(model1)

model1 = pgls(log(abs.ds.2) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 0.939, delta = 0.897)
summary(model1)

model1 = pgls(log(abs.ds.2) ~ log(age_first_breeding), pgls1, lambda = 'ML', kappa = 0.834, delta = 0.898)
summary(model1)

model1 = pgls(log(abs.ds.2) ~ log(fecundity), pgls1, lambda = 'ML', kappa = 0.850, delta = 0.832)
summary(model1)

# dS pgls model calibration set 4
pgls1 = comparative.data(tree2, pgls.data.sorted, Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds.4) ~ log(body_mass), pgls1, lambda = 'ML', kappa = 0.773, delta = 0.814)
summary(model1)

model1 = pgls(log(abs.ds.4) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 0.744, delta = 0.680)
summary(model1)

model1 = pgls(log(abs.ds.4) ~ log(age_first_breeding), pgls1, lambda = 'ML', kappa = 0.663, delta = 0.678)
summary(model1)

model1 = pgls(log(abs.ds.4) ~ log(fecundity), pgls1, lambda = 'ML', kappa = 0.670, delta = 0.621)
summary(model1)

# dS pgls model no calibration 
pgls1 = comparative.data(tree2, pgls.data.sorted, Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds.no) ~ log(body_mass), pgls1, lambda = 'ML', kappa = 0.903, delta = 0.982)
summary(model1)

model1 = pgls(log(abs.ds.no) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 0.958, delta = 0.891)
summary(model1)

model1 = pgls(log(abs.ds.no) ~ log(age_first_breeding), pgls1, lambda = 'ML', kappa = 0.850, delta = 0.889)
summary(model1)

model1 = pgls(log(abs.ds.no) ~ log(fecundity), pgls1, lambda = 'ML', kappa = 0.871, delta = 0.820)
summary(model1)

# omega pgls model 
pgls1 = comparative.data(tree2, pgls.data.sorted, Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(dnds) ~ log(body_mass), pgls1, lambda = 'ML', kappa = 0.799, delta = 1.025)
summary(model1)

model1 = pgls(log(dnds) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 0.746, delta = 1.004)
summary(model1)

model1 = pgls(log(dnds) ~ log(age_first_breeding), pgls1, lambda = 'ML', kappa = 0.776, delta = 1.044)
summary(model1)

model1 = pgls(log(dnds) ~ log(fecundity), pgls1, lambda = 'ML', kappa = 0.736, delta = 1.022)
summary(model1)

# dn pgls model calibration set 2
pgls.data.sorted = pgls.data.sorted %>% filter(abs.dn.2 > 0)
pgls1 = comparative.data(tree2, pgls.data.sorted, Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.dn.2) ~ log(body_mass), pgls1, lambda = 'ML', kappa = 0.817, delta = 0.828)
summary(model1)

model1 = pgls(log(abs.dn.2) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 0.885, delta = 0.782)
summary(model1)

model1 = pgls(log(abs.dn.2) ~ log(age_first_breeding), pgls1, lambda = 'ML', kappa = 0.798, delta = 0.819)
summary(model1)

model1 = pgls(log(abs.dn.2) ~ log(fecundity), pgls1, lambda = 'ML', kappa = 0.810, delta = 0.723)
summary(model1)

# dn pgls model calibration set 4
pgls1 = comparative.data(tree2, pgls.data.sorted, Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.dn.4) ~ log(body_mass), pgls1, lambda = 'ML', kappa = 0.695, delta = 0.639)
summary(model1)

model1 = pgls(log(abs.dn.4) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 0.701, delta = 0.562)
summary(model1)

model1 = pgls(log(abs.dn.4) ~ log(age_first_breeding), pgls1, lambda = 'ML', kappa = 0.625, delta = 0.599)
summary(model1)

model1 = pgls(log(abs.dn.4) ~ log(fecundity), pgls1, lambda = 'ML', kappa = 0.639, delta = 0.536)
summary(model1)

# dn pgls model no calibration 
pgls1 = comparative.data(tree2, pgls.data.sorted, Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.dn.no) ~ log(body_mass), pgls1, lambda = 'ML', kappa = 0.804, delta = 0.798)
summary(model1)

model1 = pgls(log(abs.dn.no) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 0.882, delta = 0.749)
summary(model1)

model1 = pgls(log(abs.dn.no) ~ log(age_first_breeding), pgls1, lambda = 'ML', kappa = 0.796, delta = 0.793)
summary(model1)

model1 = pgls(log(abs.dn.no) ~ log(fecundity), pgls1, lambda = 'ML', kappa = 0.810, delta = 0.692)
summary(model1)

