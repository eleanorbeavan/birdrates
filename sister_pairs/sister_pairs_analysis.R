# a script to perform sister pairs analysis on bird data
library(ggplot2)

#data = read.csv("./cherries_mass_rates_consensus_pairs.csv", stringsAsFactors = F) # cherries analysis
mass = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/all_traits/genera_mass_rates.csv", stringsAsFactors = F) # mass genera analysis
fec = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/all_traits/genera_fecundity_rates.csv", stringsAsFactors = F) # fecundity genera analysis
long = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/all_traits/genera_longevity_rates.csv", stringsAsFactors = F) # longevity genera analysis
gen = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/all_traits/genera_gentime_rates.csv", stringsAsFactors = F) # generation time genera analysis
true_gen = read.csv("~/Dropbox/PhD/bird_rates/processed_data/CSV_files/all_traits/genera_true_gentime_rates.csv", stringsAsFactors = F) 
### NEW ###
order = read.csv("~/Dropbox/PhD/bird_rates/processed_data/order_analysis/rates_traits.csv", stringsAsFactors = F)

# transform data according to Freckleton test
# mass
ds.a.mass = log(mass$sp1_ds)
ds.b.mass = log(mass$sp2_ds)

dn.a.mass = log(mass$sp1_dn)
dn.b.mass = log(mass$sp2_dn)

w.a.mass = log(mass$sp1_w)
w.b.mass = log(mass$sp2_w)

mass.a = log(mass$mass_sp1)
mass.b = log(mass$mass_sp2)

# longevity
ds.a.long = log(long$sp1_ds)
ds.b.long = log(long$sp2_ds)

dn.a.long = log(long$sp1_dn)
dn.b.long = log(long$sp2_dn)

w.a.long = log(long$sp1_w)
w.b.long = log(long$sp2_w)

long.a = log(long$longevity_sp1)
long.b = log(long$longevity_sp2)

# age at first breeding
ds.a.gen = log(gen$sp1_ds)
ds.b.gen = log(gen$sp2_ds)

dn.a.gen = log(gen$sp1_dn)
dn.b.gen = log(gen$sp2_dn)

w.a.gen = log(gen$sp1_w)
w.b.gen = log(gen$sp2_w)

gen.a = log(gen$gentime_sp1)
gen.b = log(gen$gentime_sp2)

# fecundity
ds.a.fec = log(fec$sp1_ds)
ds.b.fec = log(fec$sp2_ds)

dn.a.fec = log(fec$sp1_dn)
dn.b.fec = log(fec$sp2_dn)

w.a.fec = log(fec$sp1_w)
w.b.fec = log(fec$sp2_w)

fec.a = log(fec$fecundity_sp1)
fec.b = log(fec$fecundity_sp2)

# true generation time
ds.a.gen.t = log(true_gen$sp1_ds)
ds.b.gen.t = log(true_gen$sp2_ds)

dn.a.gen.t = log(true_gen$sp1_dn)
dn.b.gen.t = log(true_gen$sp2_dn)

w.a.gen.t = log(true_gen$sp1_w)
w.b.gen.t = log(true_gen$sp2_w)

gen.t.a = log(true_gen$true_gentime_sp1)
gen.t.b = log(true_gen$true_gentime_sp2)

# standardise contrasts according to Garland test
# mass
ds.garl.correction.mass = sqrt(mass$sp1_ds + mass$sp2_ds)
dn.garl.correction.mass = sqrt(mass$sp1_dn + mass$sp2_dn)
w.garl.correction.mass = sqrt(mass$sp1_w + mass$sp2_w)

ds.mass = (ds.a.mass - ds.b.mass)/(ds.garl.correction.mass)
dn.mass = (dn.a.mass - dn.b.mass)/(dn.garl.correction.mass)
w.mass = (w.a.mass - w.b.mass)/(w.garl.correction.mass)

mass.corrected.ds = (mass.a - mass.b)/(ds.garl.correction.mass)
mass.corrected.dn = (mass.a - mass.b)/(dn.garl.correction.mass)
mass.corrected.w = (mass.a - mass.b)/(w.garl.correction.mass)

# longevity
ds.garl.correction.long = sqrt(long$sp1_ds + long$sp2_ds)
dn.garl.correction.long = sqrt(long$sp1_dn + long$sp2_dn)
w.garl.correction.long = sqrt(log(long$sp1_w*1000) + log(long$sp2_w*1000)) # see data_tests_longevity.R

ds.long = (ds.a.long - ds.b.long)/(ds.garl.correction.long)
dn.long = (dn.a.long - dn.b.long)/(dn.garl.correction.long)
w.long = (w.a.long - w.b.long)/(w.garl.correction.long)

longevity.corrected.ds = (long.a - long.b)/(ds.garl.correction.long)
longevity.corrected.dn = (long.a - long.b)/(dn.garl.correction.long)
longevity.corrected.w = (long.a - long.b)/(w.garl.correction.long)

# age at first breeding
ds.garl.correction.gen = sqrt(gen$sp1_ds + gen$sp2_ds)
dn.garl.correction.gen = sqrt(gen$sp1_dn + gen$sp2_dn)
w.garl.correction.gen = sqrt(gen$sp1_w + gen$sp2_w)

ds.gen = (ds.a.gen - ds.b.gen)/(ds.garl.correction.gen)
dn.gen = (dn.a.gen - dn.b.gen)/(dn.garl.correction.gen)
w.gen = (w.a.gen - w.b.gen)/(w.garl.correction.gen)

gen.time.corrected.ds = (gen.a - gen.b)/(ds.garl.correction.gen) 
gen.time.corrected.dn = (gen.a - gen.b)/(dn.garl.correction.gen)
gen.time.corrected.w = (gen.a - gen.b)/(w.garl.correction.gen)

# fecundity
ds.garl.correction.fec = sqrt(fec$sp1_ds + fec$sp2_ds)
dn.garl.correction.fec = sqrt(fec$sp1_dn + fec$sp2_dn)
w.garl.correction.fec = sqrt(fec$sp1_w + fec$sp2_w)

ds.fec = (ds.a.fec - ds.b.fec)/(ds.garl.correction.fec)
dn.fec = (dn.a.fec - dn.b.fec)/(dn.garl.correction.fec)
w.fec = (w.a.fec - w.b.fec)/(w.garl.correction.fec)

fecundity.corrected.ds = (fec.a - fec.b)/(ds.garl.correction.fec)
fecundity.corrected.dn = (fec.a - fec.b)/(dn.garl.correction.fec)
fecundity.corrected.w = (fec.a - fec.b)/(w.garl.correction.fec)

# true gentime
ds.garl.correction.gen.t = sqrt(true_gen$sp1_ds + true_gen$sp2_ds)
dn.garl.correction.gen.t = sqrt(true_gen$sp1_dn + true_gen$sp2_dn)
w.garl.correction.gen.t = sqrt(true_gen$sp1_w + true_gen$sp2_w)

ds.gen.t = (ds.a.gen.t - ds.b.gen.t)/(ds.garl.correction.gen.t)
dn.gen.t = (dn.a.gen.t - dn.b.gen.t)/(dn.garl.correction.gen.t)
w.gen.t = (w.a.gen.t - w.b.gen.t)/(w.garl.correction.gen.t)

gen.t.corrected.ds = (gen.t.a - gen.t.b)/(ds.garl.correction.gen.t)
gen.t.corrected.dn = (gen.t.a - gen.t.b)/(dn.garl.correction.gen.t)
gen.t.corrected.w = (gen.t.a - gen.t.b)/(w.garl.correction.gen.t)

# remove contrasts accroding to Welch and Waxman test
to.remove = 0.51 # only a problem for mass dS analysis (comes from data_tests_mass.R)

sum.ds.bl <- sqrt(mass$sp1_ds + mass$sp2_ds)
ds.points.to.exclude <- which((sum.ds.bl)<=to.remove)
for(i in ds.points.to.exclude){ds.mass[i] <- NA}
for(i in ds.points.to.exclude){mass.corrected.ds[i] <- NA}

# sign tests
# mass
ds.neg.mass <- length(which(ds.mass*mass.corrected.ds<0))
ds.pos.mass <- length(which(ds.mass*mass.corrected.ds>0))
ds.signtest.mass <- binom.test(c(ds.neg.mass, ds.pos.mass), p=0.5, alternative="two.sided")

dn.neg.mass <- length(which(dn.mass*mass.corrected.dn<0))
dn.pos.mass <- length(which(dn.mass*mass.corrected.dn>0))
dn.signtest.mass <- binom.test(c(dn.neg.mass, dn.pos.mass), p=0.5, alternative="two.sided")

w.neg.mass <- length(which(w.mass*mass.corrected.w<0))
w.pos.mass <- length(which(w.mass*mass.corrected.w>0))
w.signtest.mass <- binom.test(c(w.neg.mass, w.pos.mass), p=0.5, alternative="two.sided")

# gentime
ds.neg.gen <- length(which(ds.gen*gen.time.corrected.ds<0))
ds.pos.gen <- length(which(ds.gen*gen.time.corrected.ds>0))
ds.signtest.gen <- binom.test(c(ds.neg.gen, ds.pos.gen), p=0.5, alternative="two.sided")

dn.neg.gen <- length(which(dn.gen*gen.time.corrected.dn<0))
dn.pos.gen <- length(which(dn.gen*gen.time.corrected.dn>0))
dn.signtest.gen <- binom.test(c(dn.neg.gen, dn.pos.gen), p=0.5, alternative="two.sided")

w.neg.gen <- length(which(w.gen*gen.time.corrected.w<0))
w.pos.gen <- length(which(w.gen*gen.time.corrected.w>0))
w.signtest.gen <- binom.test(c(w.neg.gen, w.pos.gen), p=0.5, alternative="two.sided")

# longevity
ds.neg.long <- length(which(ds.long*longevity.corrected.ds<0))
ds.pos.long <- length(which(ds.long*longevity.corrected.ds>0))
ds.signtest.long <- binom.test(c(ds.neg.long, ds.pos.long), p=0.5, alternative="two.sided")

dn.neg.long <- length(which(dn.long*longevity.corrected.dn<0))
dn.pos.long <- length(which(dn.long*longevity.corrected.dn>0))
dn.signtest.long <- binom.test(c(dn.neg.long, dn.pos.long), p=0.5, alternative="two.sided")

w.neg.long <- length(which(w.long*longevity.corrected.w<0))
w.pos.long <- length(which(w.long*longevity.corrected.w>0))
w.signtest.long <- binom.test(c(w.neg.long, w.pos.long), p=0.5, alternative="two.sided")

# fecundity
ds.neg.fec <- length(which(ds.fec*fecundity.corrected.ds<0))
ds.pos.fec <- length(which(ds.fec*fecundity.corrected.ds>0))
ds.signtest.fec <- binom.test(c(ds.neg.fec, ds.pos.fec), p=0.5, alternative="two.sided")

dn.neg.fec <- length(which(dn.fec*fecundity.corrected.dn<0))
dn.pos.fec <- length(which(dn.fec*fecundity.corrected.dn>0))
dn.signtest.fec <- binom.test(c(dn.neg.fec, dn.pos.fec), p=0.5, alternative="two.sided")

w.neg.fec <- length(which(w.fec*fecundity.corrected.w<0))
w.pos.fec <- length(which(w.fec*fecundity.corrected.w>0))
w.signtest.fec <- binom.test(c(w.neg.fec, w.pos.fec), p=0.5, alternative="two.sided")

# true gentime
ds.neg.gen.t <- length(which(ds.gen.t*gen.t.corrected.ds<0))
ds.pos.gen.t <- length(which(ds.gen.t*gen.t.corrected.ds>0))
ds.signtest.gen.t <- binom.test(c(ds.neg.gen.t, ds.pos.gen.t), p=0.5, alternative="two.sided")

dn.neg.gen.t <- length(which(dn.gen.t*gen.t.corrected.dn<0))
dn.pos.gen.t <- length(which(dn.gen.t*gen.t.corrected.dn>0))
dn.signtest.gen.t <- binom.test(c(dn.neg.gen.t, dn.pos.gen.t), p=0.5, alternative="two.sided")

w.neg.gen.t <- length(which(w.gen.t*gen.t.corrected.w<0))
w.pos.gen.t <- length(which(w.gen.t*gen.t.corrected.w>0))
w.signtest.gen.t <- binom.test(c(w.neg.gen.t, w.pos.gen.t), p=0.5, alternative="two.sided")

# linear model
model.1a <- lm(ds.mass ~ mass.corrected.ds - 1)
model.1b <- lm(dn.mass ~ mass.corrected.dn - 1)
model.1c <- lm(w.mass ~ mass.corrected.w - 1)

model.2a <- lm(ds.long ~ longevity.corrected.ds - 1)
model.2b <- lm(dn.long ~ longevity.corrected.dn - 1)
model.2c <- lm(w.long ~ longevity.corrected.w - 1)

model.3a <- lm(ds.fec ~ fecundity.corrected.ds - 1)
model.3b <- lm(dn.fec ~ fecundity.corrected.dn - 1)
model.3c <- lm(w.fec ~ fecundity.corrected.w - 1)

model.4a <- lm(ds.gen ~ gen.time.corrected.ds - 1)
model.4b <- lm(dn.gen ~ gen.time.corrected.dn - 1)
model.4c <- lm(w.gen ~ gen.time.corrected.w - 1)

model.5a <- lm(ds.gen.t ~ gen.t.corrected.ds - 1)
model.5b <- lm(dn.gen.t ~ gen.t.corrected.dn - 1)
model.5c <- lm(w.gen.t ~ gen.t.corrected.w - 1)

# visualise the data
ggplot(mass, aes(mass.corrected.ds, ds.mass)) + geom_point(, colour = 'deepskyblue3') + scale_x_continuous(limits = c(-2, 3)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F, colour = 'gray50') + xlab("difference in ln(mass)") + ylab("difference in ln(dS)") + theme_minimal()
ggplot(mass, aes(mass.corrected.dn, dn.mass)) + geom_point() + scale_x_continuous(limits = c(-15, 20)) + scale_y_continuous(limits = c(-13, 13)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(mass)") + ylab("difference in ln(dN)") + theme_minimal()
ggplot(mass, aes(mass.corrected.w, w.mass)) + geom_point() + scale_x_continuous(limits = c(-10, 15)) + scale_y_continuous(limits = c(-5, 5)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(mass)") + ylab("difference in ln(W)") + theme_minimal()

ggplot(long, aes(longevity.corrected.ds, ds.long)) + geom_point() + scale_x_continuous(limits = c(-1.5, 1.5)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(longevity)") + ylab("difference in ln(dS)") + theme_minimal()
ggplot(long, aes(longevity.corrected.dn, dn.long)) + geom_point() + scale_x_continuous(limits = c(-10, 10)) + scale_y_continuous(limits = c(-10, 10)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(longevity)") + ylab("difference in ln(dN)") + theme_minimal()
ggplot(long, aes(longevity.corrected.w, w.long)) + geom_point() + scale_x_continuous(limits = c(-0.5, 0.5)) + scale_y_continuous(limits = c(-0.5, 0.5)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(longevity)") + ylab("difference in ln(W)") + theme_minimal()

ggplot(fec, aes(fecundity.corrected.ds, ds.fec)) + geom_point() + scale_x_continuous(limits = c(-2, 3)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(fecundity)") + ylab("difference in ln(dS)") + theme_minimal()
ggplot(fec, aes(fecundity.corrected.dn, dn.fec)) + geom_point() + scale_x_continuous(limits = c(-10, 20)) + scale_y_continuous(limits = c(-10, 10)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(fecundity)") + ylab("difference in ln(dN)") + theme_minimal()
ggplot(fec, aes(fecundity.corrected.w, w.fec)) + geom_point() + scale_x_continuous(limits = c(-10, 15)) + scale_y_continuous(limits = c(-5, 5)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(fecundity)") + ylab("difference in ln(W)") + theme_minimal()

ggplot(gen, aes(gen.time.corrected.ds, ds.gen)) + geom_point() + scale_x_continuous(limits = c(-2, 2)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(breeding age)") + ylab("difference in ln(dS)") + theme_minimal()
ggplot(gen, aes(gen.time.corrected.dn, dn.gen)) + geom_point() + scale_x_continuous(limits = c(-10, 15)) + scale_y_continuous(limits = c(-10, 10)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(breeding age)") + ylab("difference in ln(dN)") + theme_minimal()
ggplot(gen, aes(gen.time.corrected.w, w.gen)) + geom_point() + scale_x_continuous(limits = c(-10, 13)) + scale_y_continuous(limits = c(-5, 5)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(breeding age)") + ylab("difference in ln(W)") + theme_minimal()

ggplot(true_gen, aes(gen.t.corrected.ds, ds.gen.t)) + geom_point() + scale_x_continuous(limits = c(-2, 2)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(generation time)") + ylab("difference in ln(dS)") + theme_minimal()
ggplot(true_gen, aes(gen.t.corrected.dn, dn.gen.t)) + geom_point() + scale_x_continuous(limits = c(-10, 15)) + scale_y_continuous(limits = c(-10, 10)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(generation time)") + ylab("difference in ln(dN)") + theme_minimal()
ggplot(true_gen, aes(gen.t.corrected.w, w.gen.t)) + geom_point() + scale_x_continuous(limits = c(-10, 13)) + scale_y_continuous(limits = c(-5, 5)) + stat_smooth(method=lm, formula = y~x-1, se = F) + xlab("difference in ln(generation time)") + ylab("difference in ln(W)") + theme_minimal()

# for talk 
p1 = ggplot(mass, aes(mass.corrected.ds, ds.mass)) + geom_point(colour = 'mediumpurple4', size = 0.75) + scale_x_continuous(limits = c(-2, 3)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F, colour = 'gray50', size = 0.5) + xlab("difference in ln(mass)") + ylab("difference in ln(rate)") + theme_minimal()
p2 = ggplot(long, aes(longevity.corrected.ds, ds.long)) + geom_point(colour = 'orchid3', size = 0.75) + scale_x_continuous(limits = c(-1.5, 1.5)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F, colour = 'gray50', size = 0.5) + xlab("difference in ln(longevity)") + ylab("") + theme_minimal()
p3 = ggplot(fec, aes(fecundity.corrected.ds, ds.fec)) + geom_point(colour = 'indianred3', size = 0.75) + scale_x_continuous(limits = c(-2, 3)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F, colour = 'gray50', size = 0.5) + xlab("difference in ln(fecundity)") + ylab("difference in ln(rate)") + theme_minimal()
p4 = ggplot(gen, aes(gen.time.corrected.ds, ds.gen)) + geom_point(colour = 'plum4', size = 0.75) + scale_x_continuous(limits = c(-2, 2)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F, colour = 'gray50', size = 0.5) + xlab("difference in ln(breeding age)") + ylab("") + theme_minimal()

library(cowplot)
plot_grid(p1, p2, p3 ,p4)

# diagnostics
par(mfrow=c(2,3))
plot(model.1a, which=c(1:6))
plot(model.1b, which=c(1:6))
plot(model.1c, which=c(1:6))

plot(model.2a, which=c(1:6))
plot(model.2b, which=c(1:6))
plot(model.2c, which=c(1:6))

plot(model.3a, which=c(1:6))
plot(model.3b, which=c(1:6))
plot(model.3c, which=c(1:6))

plot(model.4a, which=c(1:6))
plot(model.4b, which=c(1:6))
plot(model.4c, which=c(1:6))

### MASS ####
################
# MODEL.1A
# cooks distance
cooksd = cooks.distance(model.1a)
# leverage
leverage = influence(model.1a)$hat

# remove data points with cooks distance > 0.05 and leverage > 0.05
influential = c(as.numeric(names(cooksd)[(cooksd > 0.05)]), as.numeric(names(leverage)[(leverage > 0.05)]))
influential = unique(influential)
influential = as.numeric(unlist(influential))

# remove influential data points 
mass.corrected.ds[influential] = NA
ds.mass[influential] = NA

# second linear model
model.5a <- lm(ds.mass ~ mass.corrected.ds - 1)
summary(model.5a)
par(mfrow=c(2,3))
plot(model.5a, which=c(1:6))

# MODEL.1B
# cooks distance
cooksd = cooks.distance(model.1b)
# leverage
leverage = influence(model.1b)$hat

# remove data points with cooks distance > 0.05 and leverage > 0.05
influential = c(as.numeric(names(cooksd)[(cooksd > 0.05)]), as.numeric(names(leverage)[(leverage > 0.05)]))
influential = unique(influential)
influential = as.numeric(unlist(influential))

# remove influential data points 
mass.corrected.dn[influential] = NA
dn.mass[influential] = NA

# second linear model
model.5b <- lm(dn.mass ~ mass.corrected.dn - 1)
summary(model.5b)
par(mfrow=c(2,3))
plot(model.5b, which=c(1:6))

# MODEL.1C
# cooks distance
cooksd = cooks.distance(model.1c)
# leverage
leverage = influence(model.1c)$hat

# remove data points with cooks distance > 0.05 and leverage > 0.05
influential = c(as.numeric(names(cooksd)[(cooksd > 0.05)]), as.numeric(names(leverage)[(leverage > 0.05)]))
influential = unique(influential)
influential = as.numeric(unlist(influential))

# remove influential data points 
mass.corrected.w[influential] = NA
w.mass[influential] = NA

# second linear model
model.5c <- lm(w.mass ~ mass.corrected.w - 1)
summary(model.5c)
par(mfrow=c(2,3))
plot(model.5c, which=c(1:6))

### LONGEVITY ####
################
# cooks distance
cooksd = cooks.distance(model.2a)
# leverage
leverage = influence(model.2a)$hat

# remove data points with cooks distance > 0.05 and leverage > 0.05
influential = c(as.numeric(names(cooksd)[(cooksd > 0.1)]), as.numeric(names(leverage)[(leverage > 0.1)]))
influential = unique(influential)
influential = as.numeric(unlist(influential))

# remove influential data points 
longevity.corrected.ds[influential] = NA
ds.long[influential] = NA

# second linear model
model.6a <- lm(ds.long ~ longevity.corrected.ds - 1)
summary(model.6a)
par(mfrow=c(2,3))
plot(model.6a, which=c(1:6))

# REPEAT FOR MODEL2 B and C

### FECUNDITY ####
################
# cooks distance
cooksd = cooks.distance(model.7a)
# leverage
leverage = influence(model.7a)$hat

# remove data points with cooks distance > 0.05 and leverage > 0.05
influential = c(as.numeric(names(cooksd)[(cooksd > 0.1)]), as.numeric(names(leverage)[(leverage > 0.1)]))
influential = unique(influential)
influential = as.numeric(unlist(influential))

# remove influential data points 
ds.fec[influential] = NA
fecundity.corrected.ds[influential] = NA

# second linear model
model.7a <- lm(ds ~ fecundity - 1)
summary(model.7a)
par(mfrow=c(2,3))
plot(model.7a, which=c(1:6))

# plots FIX THESE
ggplot(data, aes(mass, ds)) + geom_point() + scale_x_continuous(limits = c(-2, 2)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F, colour = "black") + xlab("difference in ln(mass)") + ylab("difference in ln(dS)") + theme_minimal()
ggplot(data, aes(longevity, ds)) + geom_point() + scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F, colour = "black") + xlab("difference in ln(longevity)") + ylab("difference in ln(dS)") + theme_minimal()
ggplot(data, aes(fecundity, ds)) + geom_point() + scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(-1, 1)) + stat_smooth(method=lm, formula = y~x-1, se = F, colour = "black") + xlab("difference in ln(fecundity)") + ylab("difference in ln(dS)") + theme_minimal()


###########################
# sign test
ds.neg <- length(which(ds.mass*mass<0))
ds.pos <- length(which(ds.mass*mass>0))
ds.signtest <- binom.test(c(ds.neg, ds.pos), p=0.5, alternative="two.sided")

