# PGLS analysis of bird rate data
library(ape)
library(phytools)
library(geiger)
library(caper)
library(ggplot2)

# load in data
calibration.2 = read.csv("~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/ages_mass_set_2.csv", stringsAsFactors = F)
calibration.4 = read.csv("~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/ages_mass_set_4.csv", stringsAsFactors = F)
no.calibration = read.csv("~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/ages_mass_no_cal.csv", stringsAsFactors = F)

t2 = read.tree("~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/calibration_sets_2.tre")
t4 = read.tree("~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/calibration_sets_4.tre")
t = read.tree("~/Dropbox/PhD/bird_rates/processed_data/calibration_analysis/no_calibration.tre")

# divide dS and dN by node age for each species in each sister pair
abs.ds.2.sp1 = calibration.2$sp1_ds/calibration.2$node.age
abs.ds.2.sp2 = calibration.2$sp2_ds/calibration.2$node.age

abs.dn.2.sp1 = calibration.2$sp1_dn/calibration.2$node.age
abs.dn.2.sp2 = calibration.2$sp2_dn/calibration.2$node.age

abs.ds.4.sp1 = calibration.4$sp1_ds/calibration.4$node.age
abs.ds.4.sp2 = calibration.4$sp2_ds/calibration.4$node.age

abs.dn.4.sp1 = calibration.4$sp1_dn/calibration.4$node.age
abs.dn.4.sp2 = calibration.4$sp2_dn/calibration.4$node.age

abs.ds.sp1 = no.calibration$sp1_ds/no.calibration$node.age
abs.ds.sp2 = no.calibration$sp2_ds/no.calibration$node.age

abs.dn.sp1 = no.calibration$sp1_dn/no.calibration$node.age
abs.dn.sp2 = no.calibration$sp2_dn/no.calibration$node.age

# dN/dS doesnt need to be calibrated
dn.ds.sp1 = calibration.2$sp1_w
dn.ds.sp2 = calibration.2$sp2_w

# combine them into one list
abs.ds.2 = c(abs.ds.2.sp1, abs.ds.2.sp2)
abs.dn.2 = c(abs.dn.2.sp1, abs.dn.2.sp2)

abs.ds.4 = c(abs.ds.4.sp1, abs.ds.4.sp2)
abs.dn.4 = c(abs.dn.4.sp1, abs.dn.4.sp2)

abs.ds = c(abs.ds.sp1, abs.ds.sp2)
abs.dn = c(abs.dn.sp1, abs.dn.sp2)
dn.ds = c(dn.ds.sp1, dn.ds.sp2)

# make a list of species names
species = c(calibration.2$spp1, calibration.2$spp2)

# make a list of traits
mass = c(calibration.2$mass_sp1, calibration.2$mass_sp2)

# combine them into one dataframe
pgls.data = data.frame(species, mass, abs.ds.2, abs.dn.2, abs.ds.4, abs.dn.4, abs.ds, abs.dn, dn.ds)

# need to subset the tree to include only the tips in the data
tips.to.drop = setdiff(t2$tip.label, species)
tree2 = drop.tip(t2, tips.to.drop)

# rows in dataframe need to be in the same order as tips in the tree
pgls.data.sorted = pgls.data[match(tree2$tip.label, pgls.data$species),]

# dS vs mass pgls model calibration set 2
pgls1 = comparative.data(tree2, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds.2) ~ log(mass), pgls1, lambda = 0.581, kappa = 'ML', delta = 0.370)
summary(model1)

# visualise the likelihood surfaces of lambda, kappa and delta
mod.l = pgls.profile(model1, 'lambda')
mod.k = pgls.profile(model1, 'kappa')
mod.d = pgls.profile(model1, 'delta')
par(mfrow=c(1,3))
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model1, which=c(1:6))

# visualise the data
ggplot(pgls.data, aes(log(mass), log(abs.ds.2))) + geom_point() + scale_x_continuous(limits = c(0, 15)) + scale_y_continuous(limits = c(-7.5, 0)) + geom_smooth(method = "lm", se = F) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

# dN vs mass pgls model calibration set 2
pgls1b = comparative.data(tree2, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model1b = pgls(log(abs.dn.2) ~ log(mass), pgls1b, lambda = 0.678, delta = 'ML', kappa = 0.432)
summary(model1b)

# visualise the likelihood surfaces of lambda, kappa and delta
mod.l = pgls.profile(model1b, 'lambda')
mod.k = pgls.profile(model1b, 'kappa')
mod.d = pgls.profile(model1b, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model1, which=c(1:6))

# visualise the data
ggplot(pgls.data, aes(log(mass), log(abs.dn.2))) + geom_point() + scale_x_continuous(limits = c(0, 15)) + scale_y_continuous(limits = c(-10, -5)) + geom_smooth(method = "lm", se = F) + xlab("log mass") + ylab("log dN (calibration set 2)") + theme_minimal()

#############
# need to subset the tree to include only the tips in the data
tips.to.drop = setdiff(t4$tip.label, species)
tree4 = drop.tip(t4, tips.to.drop)

# rows in dataframe need to be in the same order as tips in the tree
pgls.data.sorted = pgls.data[match(tree4$tip.label, pgls.data$species),]

# dS vs mass pgls model calibration set 4
pgls2 = comparative.data(tree4, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model2 = pgls(log(abs.ds.4) ~ log(mass), pgls2, lambda = 0.701, kappa = 'ML', delta = 0.495)
summary(model2)

# visualise the likelihood surfaces of lambda, kappa and delta
mod.l = pgls.profile(model2, 'lambda')
mod.k = pgls.profile(model2, 'kappa')
mod.d = pgls.profile(model2, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model2, which=c(1:6))

# visualise the data
ggplot(pgls.data, aes(log(mass), log(abs.ds.4))) + geom_point() + xlab("log mass") + ylab("log dS (calibration set 4)") + theme_minimal()

# dN vs mass pgls model calibration set 4
pgls2b = comparative.data(tree4, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model2b = pgls(log(abs.dn.4) ~ log(mass), pgls2, lambda = 0.737, kappa = 0.449, delta = 'ML')
summary(model2b)

# visualise the likelihood surfaces of lambda, kappa and delta
mod.l = pgls.profile(model2b, 'lambda')
mod.k = pgls.profile(model2b, 'kappa')
mod.d = pgls.profile(model2b, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model2, which=c(1:6))

# visualise the data
ggplot(pgls.data, aes(log(mass), log(abs.dn.4))) + geom_point() + scale_x_continuous(limits = c(0, 15)) + scale_y_continuous(limits = c(-10, -5)) + xlab("log mass") + ylab("log dN (calibration set 4)") + theme_minimal()

########
# need to subset the tree to include only the tips in the data
tips.to.drop = setdiff(t$tip.label, species)
tree = drop.tip(t, tips.to.drop)

# rows in dataframe need to be in the same order as tips in the tree
pgls.data.sorted = pgls.data[match(tree$tip.label, pgls.data$species),]

# dS vs mass pgls model calibration set 4
pgls3 = comparative.data(tree, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model3 = pgls(log(abs.ds) ~ log(mass), pgls3, lambda = 0.675, kappa = 'ML', delta = 0.412)
summary(model3)

# visualise the likelihood surfaces of lambda, kappa and delta
mod.l = pgls.profile(model3, 'lambda')
mod.k = pgls.profile(model3, 'kappa')
mod.d = pgls.profile(model3, 'delta')
par(mfrow=c(1,3))
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model3, which=c(1:6))

# visualise the data
ggplot(pgls.data, aes(log(mass), log(abs.ds))) + geom_point() + xlab("log mass") + ylab("log dS (no calibration)") + theme_minimal()

# dN vs mass pgls model calibration set 4
pgls3b = comparative.data(tree, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model3b = pgls(log(abs.dn) ~ log(mass), pgls3b, lambda = 0.828, kappa = 'ML', delta = 0.584)
summary(model3b)

# visualise the likelihood surfaces of lambda, kappa and delta
mod.l = pgls.profile(model3b, 'lambda')
mod.k = pgls.profile(model3b, 'kappa')
mod.d = pgls.profile(model3b, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model3b, which=c(1:6))

# visualise the data
ggplot(pgls.data, aes(log(mass), log(abs.dn))) + geom_point() + xlab("log mass") + ylab("log dN (no calibration)") + theme_minimal()

## DN/DS
# dN/dS vs mass pgls model 
pgls4 = comparative.data(tree2, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model4 = pgls(log(dn.ds) ~ log(mass), pgls4, lambda = 0.456, kappa = 0.697, delta = 'ML')
summary(model4)

# visualise the likelihood surfaces of lambda, kappa and delta
mod.l = pgls.profile(model4, 'lambda')
mod.k = pgls.profile(model4, 'kappa')
mod.d = pgls.profile(model4, 'delta')
par(mfrow=c(1,3))
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model4, which=c(1:6))

# visualise the data
ggplot(pgls.data, aes(log(mass), log(dn.ds))) + geom_point() + scale_x_continuous(limits = c(0, 15)) + scale_y_continuous(limits = c(-7.5, 0)) + geom_smooth(method = "lm", se = F) + xlab("log mass") + ylab("log dN/dS (calibration set 2)") + theme_minimal()

## if it's the branch lengths transformations that are making the trend significant
# I should be able to set them to 1 and get the PIC results?
model4 = pgls(log(abs.ds.2) ~ log(mass), pgls1, lambda = 1, kappa = 1, delta = 1)
summary(model4)

par(mfrow=c(1,3))
mod.l = pgls.profile(model4, 'lambda')
mod.k = pgls.profile(model4, 'kappa')
mod.d = pgls.profile(model4, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

model4b = pgls(log(abs.dn.2) ~ log(mass), pgls1b, lambda = 1, kappa = 1, delta = 1)
summary(model4b)

mod.l = pgls.profile(model4b, 'lambda')
mod.k = pgls.profile(model4b, 'kappa')
mod.d = pgls.profile(model4b, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

model5 = pgls(log(abs.ds.4) ~ log(mass), pgls2, lambda = 1, kappa = 1, delta = 1)
summary(model5)

mod.l = pgls.profile(model5, 'lambda')
mod.k = pgls.profile(model5, 'kappa')
mod.d = pgls.profile(model5, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

model5b = pgls(log(abs.dn.4) ~ log(mass), pgls2b, lambda = 1, kappa = 1, delta = 1)
summary(model5b)

mod.l = pgls.profile(model5b, 'lambda')
mod.k = pgls.profile(model5b, 'kappa')
mod.d = pgls.profile(model5b, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

model6 = pgls(log(abs.ds) ~ log(mass), pgls3, lambda = 1, kappa = 1, delta = 1)
summary(model6)

mod.l = pgls.profile(model6, 'lambda')
mod.k = pgls.profile(model6, 'kappa')
mod.d = pgls.profile(model6, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

model6b = pgls(log(abs.dn) ~ log(mass), pgls3b, lambda = 1, kappa = 1, delta = 1)
summary(model6b)

mod.l = pgls.profile(model6b, 'lambda')
mod.k = pgls.profile(model6b, 'kappa')
mod.d = pgls.profile(model6b, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)


## Order analysis 
pgls.data = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/order_mass.csv", stringsAsFactors = F)
pgls.data$Order = as.factor(as.character(pgls.data$Order))

# rows in dataframe need to be in the same order as tips in the tree
pgls.data.sorted = pgls.data[match(tree2$tip.label, pgls.data$species),]

# dS vs mass
pgls4 = comparative.data(tree2, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model7 = pgls(log(abs.ds.2) ~ log(mass) + Order, pgls4, lambda = 'ML', kappa = 0.0001, delta = 0.865)
summary(model7)

# visualise the likelihood surfaces of lambda, kappa and delta
par(mfrow=c(1,3))
mod.l = pgls.profile(model7, 'lambda')
mod.k = pgls.profile(model7, 'kappa')
mod.d = pgls.profile(model7, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model7, which=c(1:6))

# dN vs mass
model7b = pgls(log(abs.dn.2) ~ log(mass) + Order, pgls4, lambda = 0.0001, kappa = 'ML', delta = 1.275)
summary(model7b)

# visualise the likelihood surfaces of lambda, kappa and delta
par(mfrow=c(1,3))
mod.l = pgls.profile(model7b, 'lambda')
mod.k = pgls.profile(model7b, 'kappa')
mod.d = pgls.profile(model7b, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model7, which=c(1:6))

# SET 4 
# dS vs mass
pgls5 = comparative.data(tree4, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model8 = pgls(log(abs.ds.4) ~ log(mass) + Order, pgls5, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model8)

# visualise the likelihood surfaces of lambda, kappa and delta
par(mfrow=c(1,3))
mod.l = pgls.profile(model8, 'lambda')
mod.k = pgls.profile(model8, 'kappa')
mod.d = pgls.profile(model8, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model8, which=c(1:6))

# dN vs mass
model8b = pgls(log(abs.dn.4) ~ log(mass) + Order, pgls5, lambda = 0.0001, kappa = 0.733, delta = 'ML')
summary(model8b)

# visualise the likelihood surfaces of lambda, kappa and delta
par(mfrow=c(1,3))
mod.l = pgls.profile(model8b, 'lambda')
mod.k = pgls.profile(model8b, 'kappa')
mod.d = pgls.profile(model8b, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model8b, which=c(1:6))

## NO CAL
# dS vs mass
pgls6 = comparative.data(tree, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model9 = pgls(log(abs.ds) ~ log(mass) + Order, pgls6, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model9)

# visualise the likelihood surfaces of lambda, kappa and delta
par(mfrow=c(1,3))
mod.l = pgls.profile(model9, 'lambda')
mod.k = pgls.profile(model9, 'kappa')
mod.d = pgls.profile(model9, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model9, which=c(1:6))

# dN vs mass
model9b = pgls(log(abs.dn) ~ log(mass) + Order, pgls6, lambda = 0.0001, kappa = 0.464, delta = 'ML')
summary(model9b)

# visualise the likelihood surfaces of lambda, kappa and delta
par(mfrow=c(1,3))
mod.l = pgls.profile(model9b, 'lambda')
mod.k = pgls.profile(model9b, 'kappa')
mod.d = pgls.profile(model9b, 'delta')
plot(mod.l); plot(mod.k); plot(mod.d)

# visualise the model
par(mfrow=c(2,2))
plot(model9b, which=c(1:6))

# dS set 2
## to check which models better fit the data:
AIC(model1, model7)
# dN set 2
# the lower the better so model 7 - with order - is better
AIC(model1b, model7b)
# model 7b is lower - with order for dN

# dS set 4
AIC(model2, model8)
# model 8 is  better
# dN set 4
AIC(model2b, model8b)
# model 8b is better

# dS no cal
AIC(model3, model9)
# model 9 is better
# dN no cal
AIC(model3b, model9b)
# model 9b is better

# The order analysis is consistently better, suggesting order explains something that the other traits don't

# use melt and ggplot to plot all 3 calibration sets on one plot
library(reshape2)
library(viridis)

# dS
ds.data = pgls.data[,c(1,2,3,5,7)]
plot.ds.data = melt(ds.data, id.vars = c('species', 'mass'), measure.vars = c('abs.ds.2', 'abs.ds.4', 'abs.ds'))

p3 = ggplot(plot.ds.data, aes(log(mass), log(value), colour = variable)) + geom_point(size = 1) + xlab("log mass") + ylab("log dS") + stat_smooth(method = 'lm', se = F, size = 0.5) + theme_minimal() + theme(legend.position = "none")
p3 

dn.data = pgls.data[,c(1,2,4,6,8)]
plot.dn.data = melt(dn.data, id.vars = c('species', 'mass'), measure.vars = c('abs.dn.2', 'abs.dn.4', 'abs.dn'))

p4 = ggplot(plot.dn.data, aes(log(mass), log(value), colour = variable)) + geom_point(size = 1) + xlab("log mass") + ylab("log dN") + stat_smooth(method = 'lm', se = F, size = 0.5) + theme_minimal() + theme(legend.position = "none")
