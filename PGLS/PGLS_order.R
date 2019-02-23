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
calibration.2 = read.csv("name/of/output/from/get_node_ages.R", stringsAsFactors = F)
calibration.4 = read.csv("name/of/output/from/get_node_ages.R", stringsAsFactors = F)
no.calibration = read.csv("name/of/output/from/get_node_ages.R", stringsAsFactors = F)

t2 = read.tree("/path/to/calibration_sets_2.tre")
t4 = read.tree("/path/to/calibration_sets_4.tre")
t = read.tree("/path/to/no_calibration.tre")

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

# dS pgls model calibration set 2
pgls1 = comparative.data(tree2, pgls.data.sorted, Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds.2) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity), pgls1, lambda = 'ML', kappa = 0.996, delta = 0.998)
summary(model1)

### include order as a factor and compare AIC scores
model2 = pgls(log(abs.ds.2) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity)+order, pgls1, lambda = 'ML', kappa = 1.005, delta = 1.500)
summary(model2)

AIC(model1, model2)

# dS pgls model calibration set 4
model1 = pgls(log(abs.ds.4) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity), pgls1, lambda = 'ML', kappa = 0.728, delta = 0.806)
summary(model1)

### include order as a factor and compare AIC scores
model2 = pgls(log(abs.ds.4) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity)+order, pgls1, lambda = 'ML', kappa = 1.055, delta = 1.220)
summary(model2)

AIC(model1, model2)

# dS pgls model no calibration set 
model1 = pgls(log(abs.ds.no) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity), pgls1, lambda = 'ML', kappa = 1.005, delta = 1.003)
summary(model1)

### include order as a factor and compare AIC scores
model2 = pgls(log(abs.ds.no) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity)+order, pgls1, lambda = 'ML', kappa = 0.993, delta = 1.500)
summary(model2)

AIC(model1, model2)


## dnds
# omega pgls model no calibration set 
model1 = pgls(log(dnds) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity), pgls1, lambda = 'ML', kappa = 0.805, delta = 1.102)
summary(model1)

### include order as a factor and compare AIC scores
model2 = pgls(log(dnds) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity)+order, pgls1, lambda = 0.0001, kappa = 'ML', delta = 1.584)
summary(model2)

AIC(model1, model2)


# dN pgls model calibration set 2
## it doesnt work if dn is 0 because you can't do log 0
pgls.data.sorted = pgls.data.sorted %>% filter(abs.dn.2 > 0)
pgls1 = comparative.data(tree2, pgls.data.sorted, Tips, vcv= T, vcv.dim = 3)
# dN pgls model calibration set 2
model1 = pgls(log(abs.dn.2) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity), pgls1, lambda = 0.447, kappa = 'ML', delta = 0.902)
summary(model1)

### include order as a factor and compare AIC scores
model2 = pgls(log(abs.dn.2) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity)+order, pgls1, lambda = 'ML', kappa = 0.898, delta = 1.252)
summary(model2)

AIC(model1, model2)

### dN pgls model calibration set 4
model1 = pgls(log(abs.dn.4) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity), pgls1, lambda = 'ML', kappa = 0.528, delta = 0.656)
summary(model1)

### include order as a factor and compare AIC scores
model2 = pgls(log(abs.dn.4) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity)+order, pgls1, lambda = 'ML', kappa = 0.898, delta = 1.252)
summary(model2)

AIC(model1, model2)

## dN pgls model no calibration set 
model1 = pgls(log(abs.dn.no) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity), pgls1, lambda = 0.0001, kappa = 0.878, delta = 1.291)
summary(model1)

### include order as a factor and compare AIC scores
model2 = pgls(log(abs.dn.no) ~ log(body_mass)+log(max_longevity)+log(age_first_breeding)+log(fecundity)+order, pgls1, lambda = 'ML', kappa = 0.875, delta = 1.240)
summary(model2)

AIC(model1, model2)



###########################
###### WITHIN ORDERS ######
###########################
### MASS
pgls.data.mass = pgls.data.sorted %>% filter(body_mass != 'NA')
number_order = pgls.data.mass %>% group_by(order) %>% summarise(n = n())
keep = number_order %>% filter(n >= 5)
pgls.data.mass = merge(pgls.data.mass, keep, by = 'order', all.x = T)
pgls.data.mass = pgls.data.mass %>% filter(n >= 5)
pgls.data.mass = pgls.data.mass %>% group_by(order)
pgls.data.mass = pgls.data.mass[,c(1:4,8:10)]

# loop to get separate dataframes for all orders > 5
orders = list()
for (i in 1:length(unique(pgls.data.mass$order))) {
  x = attr(pgls.data.mass, 'labels')[i,]
  x = pgls.data.mass %>% filter(order == x)
  orders = c(orders, list(as.data.frame(x)))
}

#Accipitriformes
pgls1 = comparative.data(tree2, orders[[1]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 0.802, kappa = 0.248, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 1, kappa = 0.271, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 1.176, delta = 'ML')
summary(model1)

ggplot(orders[[1]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Anseriformes
pgls1 = comparative.data(tree2, orders[[2]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 0.001, kappa = 'ML', delta = 2.191)
summary(model1)

orders[[2]] = orders[[2]] %>% filter(abs.dn > 0)
orders[[2]]$log.dn = log(orders[[2]]$abs.dn)
pgls1 = comparative.data(tree2, orders[[2]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log.dn ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 0.622, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[2]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Apodiformes
pgls1 = comparative.data(tree2, orders[[3]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 3, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 'ML', kappa = 1.200, delta = 0.990)
summary(model1)

ggplot(orders[[3]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Bucerotiformes
pgls1 = comparative.data(tree2, orders[[4]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 3, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 3, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 1, kappa = 'ML', delta = 1.878)
summary(model1)

ggplot(orders[[4]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Charadriiformes
pgls1 = comparative.data(tree2, orders[[5]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 1, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 0.127, delta = 1.228, kappa = 'ML')
summary(model1)

ggplot(orders[[5]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Columbiformes
pgls1 = comparative.data(tree2, orders[[6]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 1, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 1, kappa = 3, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 1, kappa = 3, delta = 'ML')
summary(model1)

ggplot(orders[[6]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Falconiformes
pgls1 = comparative.data(tree2, orders[[7]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 2.303, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 3, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 0.404, delta = 'ML')
summary(model1)

ggplot(orders[[7]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Galliformes
pgls1 = comparative.data(tree2, orders[[8]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 'ML', kappa = 0.0001, delta = 1.108)
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 0.001, delta = 'ML')
summary(model1)

ggplot(orders[[8]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Gruiformes
pgls1 = comparative.data(tree2, orders[[9]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 1, kappa = 0.791, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 'ML', kappa = 1.889, delta = 0.675)
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 0.709, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[9]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Passeriformes
pgls1 = comparative.data(tree2, orders[[10]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 0.546, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[10]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Pelecaniformes
pgls1 = comparative.data(tree2, orders[[11]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 1, kappa = 0.076, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 1, kappa = 0.351, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 0.522, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[11]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Procellariiformes
pgls1 = comparative.data(tree2, orders[[12]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 1, kappa = 1.658, delta = ,'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 'ML', kappa = 1.176, delta = 1.800)
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 0.00001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[12]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Psittaciformes
orders[[13]] = orders[[13]] %>% filter(abs.dn > 0)
pgls1 = comparative.data(tree2, orders[[13]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 1, kappa = 0.132, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 0.864, kappa = 1.921, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 1, kappa = 1.366, delta = 'ML')
summary(model1)

ggplot(orders[[13]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

#Sphenisciformes
pgls1 = comparative.data(tree2, orders[[14]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(body_mass), pgls1, lambda = 1, kappa = 1, delta = 0.079)
summary(model1)

model1 = pgls(log(abs.dn) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 1.616, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(body_mass), pgls1, lambda = 0.0001, kappa = 1.322, delta = 'ML')
summary(model1)

ggplot(orders[[14]], aes(log(body_mass), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log mass") + ylab("log dS (calibration set 2)") + theme_minimal()

### LONGEVITY
pgls.data.lon = pgls.data.sorted %>% filter(max_longevity != 'NA')
number_order = pgls.data.lon %>% group_by(order) %>% summarise(n = n())
keep = number_order %>% filter(n >= 5)
pgls.data.lon = merge(pgls.data.lon, keep, by = 'order', all.x = T)
pgls.data.lon = pgls.data.lon %>% filter(n >= 5)
pgls.data.lon = pgls.data.lon %>% group_by(order)
pgls.data.lon = pgls.data.lon[,c(1:3,5,8:10)]

# loop to get separate dataframes for all orders > 5
orders = list()
for (i in 1:length(unique(pgls.data.lon$order))) {
  x = attr(pgls.data.lon, 'labels')[i,]
  x = pgls.data.lon %>% filter(order == x)
  orders = c(orders, list(as.data.frame(x)))
}

#Accipitriformes
pgls1 = comparative.data(tree2, orders[[1]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(max_longevity), pgls1, lambda = 0.820, kappa = 0.159, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 0.322, delta = 2.036)
summary(model1)

model1 = pgls(log(dn.ds) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 1, delta ='ML')
summary(model1)

ggplot(orders[[1]], aes(log(max_longevity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log longevity") + ylab("log dS (calibration set 2)") + theme_minimal()

#Anseriformes
orders[[2]] = orders[[2]] %>% filter(abs.dn > 0)
pgls1 = comparative.data(tree2, orders[[2]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 'ML', delta = 1.5)
summary(model1)

model1 = pgls(log(abs.dn) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[2]], aes(log(max_longevity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log longevity") + ylab("log dS (calibration set 2)") + theme_minimal()

#Charadriiformes
pgls1 = comparative.data(tree2, orders[[3]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(max_longevity), pgls1, lambda = 1, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 0.5, delta = 1.500)
summary(model1)

ggplot(orders[[3]], aes(log(max_longevity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log longevity") + ylab("log dS (calibration set 2)") + theme_minimal()

#Columbiformes
pgls1 = comparative.data(tree2, orders[[4]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 0.0001, delta= 1.208)
summary(model1)

model1 = pgls(log(abs.dn) ~ log(max_longevity), pgls1, lambda = 1, kappa = 'ML', delta = 1.134)
summary(model1)

model1 = pgls(log(dn.ds) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 1.801, delta = 'ML')
summary(model1)

ggplot(orders[[4]], aes(log(max_longevity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log longevity") + ylab("log dS (calibration set 2)") + theme_minimal()

#Falconiformes
pgls1 = comparative.data(tree2, orders[[5]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(max_longevity), pgls1, lambda = 0.001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[5]], aes(log(max_longevity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log longevity") + ylab("log dS (calibration set 2)") + theme_minimal()

#Galliformes
pgls1 = comparative.data(tree2, orders[[6]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(max_longevity), pgls1, lambda = 0.568, kappa = 1.910, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(max_longevity), pgls1, lambda = 0.556, kappa = 1.889, delta = 'ML')
summary(model1)

ggplot(orders[[6]], aes(log(max_longevity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log longevity") + ylab("log dS (calibration set 2)") + theme_minimal()

#Gruiformes
pgls1 = comparative.data(tree2, orders[[7]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(max_longevity), pgls1, lambda = 1, kappa = 'ML', delta = 0.512)
summary(model1)

model1 = pgls(log(abs.dn) ~ log(max_longevity), pgls1, lambda = 0.720, kappa = 1.889, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(max_longevity), pgls1, lambda = 0.737, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[7]], aes(log(max_longevity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log longevity") + ylab("log dS (calibration set 2)") + theme_minimal()

# Passeriformes
pgls1 = comparative.data(tree2, orders[[8]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 0.0001, delta = 2.505)
summary(model1)

model1 = pgls(log(abs.dn) ~ log(max_longevity), pgls1, lambda = 0.398, kappa = 1.374, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(max_longevity), pgls1, lambda = 0.239, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[8]], aes(log(max_longevity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log longevity") + ylab("log dS (calibration set 2)") + theme_minimal()

#Pelecaniformes
pgls1 = comparative.data(tree2, orders[[9]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(max_longevity), pgls1, lambda = 1, kappa = 'ML', delta = 0.931)
summary(model1)

model1 = pgls(log(abs.dn) ~ log(max_longevity), pgls1, lambda = 0.643, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[9]], aes(log(max_longevity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log longevity") + ylab("log dS (calibration set 2)") + theme_minimal()

#Procellariiformes
pgls1 = comparative.data(tree2, orders[[10]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 0.321, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(max_longevity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[10]], aes(log(max_longevity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log longevity") + ylab("log dS (calibration set 2)") + theme_minimal()

#Psittaciformes
orders[[11]] = orders[[11]] %>% filter(abs.dn > 0)
pgls1 = comparative.data(tree2, orders[[11]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 0.624, delta = 1.794)
summary(model1)

model1 = pgls(log(abs.dn) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 2.790, delta = 0.956)
summary(model1)

model1 = pgls(log(dn.ds) ~ log(max_longevity), pgls1, lambda = 'ML', kappa = 1.587, delta = 3)
summary(model1)

ggplot(orders[[11]], aes(log(max_longevity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log longevity") + ylab("log dS (calibration set 2)") + theme_minimal()


### GEN TIME
pgls.data.gen = pgls.data.sorted %>% filter(age_first_breeding != 'NA')
number_order = pgls.data.gen %>% group_by(order) %>% summarise(n = n())
keep = number_order %>% filter(n >= 5)
pgls.data.gen = merge(pgls.data.gen, keep, by = 'order', all.x = T)
pgls.data.gen = pgls.data.gen %>% filter(n >= 5)
pgls.data.gen = pgls.data.gen %>% group_by(order)
pgls.data.gen = pgls.data.gen[,c(1:3,6,8:10)]

# loop to get separate dataframes for all orders > 5
orders = list()
for (i in 1:length(unique(pgls.data.gen$order))) {
  x = attr(pgls.data.gen, 'labels')[i,]
  x = pgls.data.gen %>% filter(order == x)
  orders = c(orders, list(as.data.frame(x)))
}

# Accipitriformes
pgls1 = comparative.data(tree2, orders[[1]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(age_first_breeding), pgls1, lambda = 0.833, kappa = 0.030, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(age_first_breeding), pgls1, lambda = 0.994, kappa = 0.259, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 1.681, delta = 'ML')
summary(model1)

ggplot(orders[[1]], aes(log(age_first_breeding), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log age first breeding") + ylab("log dS (calibration set 2)") + theme_minimal()

# Anseriformes
pgls1 = comparative.data(tree2, orders[[2]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 0.689, delta = 'ML')
summary(model1)

orders[[2]] = orders[[2]] %>% filter(abs.dn > 0)
pgls1 = comparative.data(tree2, orders[[2]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.dn) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 0.797, delta ='ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[2]], aes(log(age_first_breeding), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log age first breeding") + ylab("log dS (calibration set 2)") + theme_minimal()

#Charadriiformes
pgls1 = comparative.data(tree2, orders[[3]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(age_first_breeding), pgls1, lambda = 1, kappa = 'ML', delta = 2.281)
summary(model1)

model1 = pgls(log(abs.dn) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[3]], aes(log(age_first_breeding), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log age first breeding") + ylab("log dS (calibration set 2)") + theme_minimal()

#Falconiformes
pgls1 = comparative.data(tree2, orders[[4]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 1.118, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(age_first_breeding), pgls1, lambda = 0.001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[4]], aes(log(age_first_breeding), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log age first breeding") + ylab("log dS (calibration set 2)") + theme_minimal()

#Galliformes
pgls1 = comparative.data(tree2, orders[[5]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(age_first_breeding), pgls1, lambda = 0.921, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(age_first_breeding), pgls1, lambda = 0.712, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[5]], aes(log(age_first_breeding), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log age first breeding") + ylab("log dS (calibration set 2)") + theme_minimal()

#Gruiformes
pgls1 = comparative.data(tree2, orders[[6]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(age_first_breeding), pgls1, lambda = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(age_first_breeding), pgls1, lambda = 0.560, kappa = 2.557, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(age_first_breeding), pgls1, lambda = 0.613, kappa = 0.869, delta = 'ML')
summary(model1)

ggplot(orders[[6]], aes(log(age_first_breeding), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log age first breeding") + ylab("log dS (calibration set 2)") + theme_minimal()

# Passeriformes
pgls1 = comparative.data(tree2, orders[[7]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(age_first_breeding), pgls1, lambda = 0.468, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(age_first_breeding), pgls1, lambda = 0.308, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(age_first_breeding), pgls1, lambda = 0.001, kappa = 0.459, delta = 'ML')
summary(model1)

ggplot(orders[[7]], aes(log(age_first_breeding), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log age first breeding") + ylab("log dS (calibration set 2)") + theme_minimal()

#Pelecaniformes
pgls1 = comparative.data(tree2, orders[[8]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(age_first_breeding), pgls1, lambda = 1, kappa = 3, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(age_first_breeding), pgls1, lambda = 1, kappa = 1.271, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 2.079, delta = 'ML')
summary(model1)

ggplot(orders[[8]], aes(log(age_first_breeding), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log age first breeding") + ylab("log dS (calibration set 2)") + theme_minimal()

#Procellariiformes
pgls1 = comparative.data(tree2, orders[[9]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(age_first_breeding), pgls1, lambda = 1, kappa = 1.449, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 1.424, delta = 'ML')
summary(model1)

ggplot(orders[[9]], aes(log(age_first_breeding), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log age first breeding") + ylab("log dS (calibration set 2)") + theme_minimal()

#Sphenisciformes
pgls1 = comparative.data(tree2, orders[[10]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 1.897, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(age_first_breeding), pgls1, lambda = 0.0001, kappa = 2.495, delta = 'ML')
summary(model1)

ggplot(orders[[10]], aes(log(age_first_breeding), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log age first breeding") + ylab("log dS (calibration set 2)") + theme_minimal()

### FECUNDITY
pgls.data.fec = pgls.data.sorted %>% filter(fecundity != 'NA')
number_order = pgls.data.fec %>% group_by(order) %>% summarise(n = n())
keep = number_order %>% filter(n >= 5)
pgls.data.fec = merge(pgls.data.fec, keep, by = 'order', all.x = T)
pgls.data.fec = pgls.data.fec %>% filter(n >= 5)
pgls.data.fec = pgls.data.fec %>% group_by(order)
pgls.data.fec = pgls.data.fec[,c(1:3,7,8:10)]

# loop to get separate dataframes for all orders > 5
orders = list()
for (i in 1:length(unique(pgls.data.fec$order))) {
  x = attr(pgls.data.fec, 'labels')[i,]
  x = pgls.data.fec %>% filter(order == x)
  orders = c(orders, list(as.data.frame(x)))
}

# Accipitriformes
pgls1 = comparative.data(tree2, orders[[1]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(fecundity), pgls1, lambda = 1, kappa = 0.023, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(fecundity), pgls1, lambda = 0.971, kappa = 0.369, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[1]], aes(log(fecundity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log fecundity") + ylab("log dS (calibration set 2)") + theme_minimal()

# Anseriformes
pgls1 = comparative.data(tree2, orders[[2]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

orders[[2]] = orders[[2]] %>% filter(abs.dn > 0)
pgls1 = comparative.data(tree2, orders[[2]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.dn) ~ log(fecundity), pgls1, lambda = 0.351, kappa = 'ML', delta = 1.219)
summary(model1)

model1 = pgls(log(dn.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[2]], aes(log(fecundity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log fecundity") + ylab("log dS (calibration set 2)") + theme_minimal()

# Charadriiformes
pgls1 = comparative.data(tree2, orders[[3]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(fecundity), pgls1, lambda = 0.994, kappa = 1.406, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[3]], aes(log(fecundity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log fecundity") + ylab("log dS (calibration set 2)") + theme_minimal()

# Falconiformes
pgls1 = comparative.data(tree2, orders[[4]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 1.466, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[4]], aes(log(fecundity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log fecundity") + ylab("log dS (calibration set 2)") + theme_minimal()

# Galliformes
pgls1 = comparative.data(tree2, orders[[5]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 0.637, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[5]], aes(log(fecundity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log fecundity") + ylab("log dS (calibration set 2)") + theme_minimal()

# Gruiformes
pgls1 = comparative.data(tree2, orders[[6]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(fecundity), pgls1, lambda = 'ML', kappa = 0.0001, delta = 1.186)
summary(model1)

model1 = pgls(log(abs.dn) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 'ML', delta = 1.470)
summary(model1)

model1 = pgls(log(dn.ds) ~ log(fecundity), pgls1, lambda = 0.681, kappa = 0.108, delta = 'ML')
summary(model1)

ggplot(orders[[6]], aes(log(fecundity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log fecundity") + ylab("log dS (calibration set 2)") + theme_minimal()

# Passeriformes
pgls1 = comparative.data(tree2, orders[[7]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(fecundity), pgls1, lambda = 0.771, kappa = 0.0001, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(fecundity), pgls1, lambda = 0.446, kappa = 0.404, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(fecundity), pgls1, lambda = 0.267, kappa = 0.0001, delta = 'ML')
summary(model1)

ggplot(orders[[7]], aes(log(fecundity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log fecundity") + ylab("log dS (calibration set 2)") + theme_minimal()

# Pelecaniformes
pgls1 = comparative.data(tree2, orders[[8]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(fecundity), pgls1, lambda = 1, kappa = 1.550, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(fecundity), pgls1, lambda = 0.810, kappa = 1.279, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 'ML', delta = 1.378)
summary(model1)

ggplot(orders[[8]], aes(log(fecundity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log fecundity") + ylab("log dS (calibration set 2)") + theme_minimal()

# Procellariiformes
pgls1 = comparative.data(tree2, orders[[9]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(fecundity), pgls1, lambda = 1, kappa = 1.966, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(fecundity), pgls1, lambda = 0.238, kappa = 1, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 1, delta = 'ML')
summary(model1)

ggplot(orders[[9]], aes(log(fecundity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log fecundity") + ylab("log dS (calibration set 2)") + theme_minimal()

# Sphenisciformes
pgls1 = comparative.data(tree2, orders[[10]], Tips, vcv= T, vcv.dim = 3)
model1 = pgls(log(abs.ds) ~ log(fecundity), pgls1, lambda = 1, kappa = 3, delta = 'ML')
summary(model1)

model1 = pgls(log(abs.dn) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 1, delta = 'ML')
summary(model1)

model1 = pgls(log(dn.ds) ~ log(fecundity), pgls1, lambda = 0.0001, kappa = 1, delta = 'ML')
summary(model1)

ggplot(orders[[10]], aes(log(fecundity), log(abs.ds))) + geom_point() + geom_smooth(method = "lm", se = F, aes(linetype = order)) + xlab("log fecundity") + ylab("log dS (calibration set 2)") + theme_minimal()
