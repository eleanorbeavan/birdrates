# PGLS analysis of bird rate data
library(ape)
library(phytools)
library(geiger)
library(caper)
library(ggplot2)

# load in data
calibration.2 = read.csv("name/of/output/from/get_node_ages.R", stringsAsFactors = F)
calibration.4 = read.csv("name/of/output/from/get_node_ages.R", stringsAsFactors = F)
no.calibration = read.csv("name/of/output/from/get_node_ages.R", stringsAsFactors = F)

t2 = read.tree("path/to/calibration_sets_2.tre")
t4 = read.tree("path/to/calibration_sets_4.tre")
t = read.tree("path/to/no_calibration.tre")

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
trait = c(calibration.2$trait_sp1, calibration.2$trait_sp2) ## change this to whichever trait you are looking at

# combine them into one dataframe
pgls.data = data.frame(species, trait, abs.ds.2, abs.dn.2, abs.ds.4, abs.dn.4, abs.ds, abs.dn, dn.ds)

# need to subset the tree to include only the tips in the data
tips.to.drop = setdiff(t2$tip.label, species)
tree2 = drop.tip(t2, tips.to.drop)

# rows in dataframe need to be in the same order as tips in the tree
pgls.data.sorted = pgls.data[match(tree2$tip.label, pgls.data$species),]

# dS vs trait pgls model calibration set 2
pgls1 = comparative.data(tree2, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds.2) ~ log(trait), pgls1, lambda = , kappa = , delta = ) ## fill in values of k, l and d
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
ggplot(pgls.data, aes(log(trait), log(abs.ds.2))) + geom_point() + scale_x_continuous(limits = c(0, 15)) + scale_y_continuous(limits = c(-7.5, 0)) + geom_smooth(method = "lm", se = F) + xlab("log trait") + ylab("log dS (calibration set 2)") + theme_minimal()

# dN vs triat pgls model calibration set 2
pgls1b = comparative.data(tree2, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model1b = pgls(log(abs.dn.2) ~ log(trait), pgls1b, lambda = , delta = , kappa = ) ## fill in values of k, l and d
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
ggplot(pgls.data, aes(log(trait), log(abs.dn.2))) + geom_point() + scale_x_continuous(limits = c(0, 15)) + scale_y_continuous(limits = c(-10, -5)) + geom_smooth(method = "lm", se = F) + xlab("log trait") + ylab("log dN (calibration set 2)") + theme_minimal()

#############
# need to subset the tree to include only the tips in the data
tips.to.drop = setdiff(t4$tip.label, species)
tree4 = drop.tip(t4, tips.to.drop)

# rows in dataframe need to be in the same order as tips in the tree
pgls.data.sorted = pgls.data[match(tree4$tip.label, pgls.data$species),]

# dS vs trait pgls model calibration set 4
pgls2 = comparative.data(tree4, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model2 = pgls(log(abs.ds.4) ~ log(trait), pgls2, lambda = , kappa = , delta = )
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
ggplot(pgls.data, aes(log(trait), log(abs.ds.4))) + geom_point() + xlab("log trait") + ylab("log dS (calibration set 4)") + theme_minimal()

# dN vs trait pgls model calibration set 4
pgls2b = comparative.data(tree4, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model2b = pgls(log(abs.dn.4) ~ log(trait), pgls2, lambda = , kappa = , delta = )
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
ggplot(pgls.data, aes(log(trait), log(abs.dn.4))) + geom_point() + scale_x_continuous(limits = c(0, 15)) + scale_y_continuous(limits = c(-10, -5)) + xlab("log trait") + ylab("log dN (calibration set 4)") + theme_minimal()

########
# need to subset the tree to include only the tips in the data
tips.to.drop = setdiff(t$tip.label, species)
tree = drop.tip(t, tips.to.drop)

# rows in dataframe need to be in the same order as tips in the tree
pgls.data.sorted = pgls.data[match(tree$tip.label, pgls.data$species),]

# dS vs trait pgls model calibration set 4
pgls3 = comparative.data(tree, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model3 = pgls(log(abs.ds) ~ log(trait), pgls3, lambda = , kappa = , delta = )
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
ggplot(pgls.data, aes(log(trait), log(abs.ds))) + geom_point() + xlab("log trait") + ylab("log dS (no calibration)") + theme_minimal()

# dN vs trait pgls model calibration set 4
pgls3b = comparative.data(tree, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model3b = pgls(log(abs.dn) ~ log(trait), pgls3b, lambda = , kappa = , delta = )
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
ggplot(pgls.data, aes(log(trait), log(abs.dn))) + geom_point() + xlab("log trait") + ylab("log dN (no calibration)") + theme_minimal()

## DN/DS
# dN/dS vs trait pgls model 
pgls4 = comparative.data(tree2, pgls.data.sorted, species, vcv= T, vcv.dim = 3)
model4 = pgls(log(dn.ds) ~ log(trait), pgls4, lambda = , kappa = , delta = )
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
ggplot(pgls.data, aes(log(trait), log(dn.ds))) + geom_point() + scale_x_continuous(limits = c(0, 15)) + scale_y_continuous(limits = c(-7.5, 0)) + geom_smooth(method = "lm", se = F) + xlab("log trait") + ylab("log dN/dS") + theme_minimal()
