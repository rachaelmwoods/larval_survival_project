library(lme4)
library(boot)
source("R/functions.R")

# Load datasets
dat <- read.csv("data/data_table_factors - data.csv", as.is=TRUE)
water <- read.csv("data/water_samples.csv", as.is=TRUE)

# Data preparation done in separate file
source("R/data_prep.R")

## FERTILISATION - MODEL - GLM - FULL model - sediment, ammonium, phosphorous, copper, tributyltin, zinc, cadmium, salinity, salinity sq, nitrate

data <- dat[dat$life.stage == "fertilisation",]
data$rep <- 1:nrow(data)

# Seperate factors to be used
#data2 <- with(dat, data.frame(life.stage, sediment_mg_per_l, ammonium_microM, phosphorous_microM, copper_ug_per_l, tributyltin_ug_per_l, zinc_ug_per_l, cadmium_ug_per_l, salinity_psu, nitrate_microM, acidification_pH, tempertaure_degrees_celcius, tempertaure_degrees_celcius_sq))
# Subset main dataset
#data <- dat[apply(!is.na(data2[, -1]), 1, sum) > 0 & data2$life.stage == "fertilisation",]
# Add random observation variable for overdispersion

# Simple model without random effects to get feel for patterns



mod_fert_full <- glm(cbind(success, failure) ~ sediment_mg_per_l + ammonium_microM + phosphorous_microM + copper_ug_per_l + tributyltin_ug_per_l + zinc_ug_per_l + cadmium_ug_per_l + salinity_psu + salinity_psu_sq + nitrate_microM + acidification_pH + acidification_pH_sq + tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq, family=binomial, data)


# This doesn't work because overparametrised... 
mod_fert_full <- glmer(cbind(success, failure) ~ sediment_mg_per_l + ammonium_microM + phosphorous_microM + copper_ug_per_l + tributyltin_ug_per_l + zinc_ug_per_l + cadmium_ug_per_l + salinity_psu + salinity_psu_sq + nitrate_microM + acidification_pH + acidification_pH_sq + tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq + (1 | experiment) + (1 | rep), family=binomial, data)

# TRY with fewer variables -- remove ones we know not having effect
data2 <- with(dat2, data.frame(life.stage, sediment_mg_per_l, phosphorous_microM, copper_ug_per_l, zinc_ug_per_l, salinity_psu, nitrate_microM, tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq))
data = dat[apply(!is.na(data2[, -1]), 1, sum) > 0 & data2$life.stage == "fertilisation",]
data$rep <- 1:dim(data)[1]

# Works now, complains about rescaling some variables, but okay
mod_fert_full <- glmer(cbind(success, failure) ~ sediment_mg_per_l + phosphorous_microM + copper_ug_per_l + zinc_ug_per_l + salinity_psu + salinity_psu_sq + nitrate_microM + tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq +(1 | experiment) + (1 | rep), family=binomial, data)
# Below should be less than one (or else overdispersion)
sum(residuals(mod_fert_full, type="pearson")^2)/df.residual(mod_fert_full)

drop1(mod_fert_full, test="Chisq")

# Drop some based on above for final model
data2 <- with(dat2, data.frame(life.stage, sediment_mg_per_l, phosphorous_microM, copper_ug_per_l, salinity_psu))
data = dat[apply(!is.na(data2[, -1]), 1, sum) > 0 & data2$life.stage == "fertilisation",]
data$rep <- 1:dim(data)[1]

mod_fert_final <- glmer(cbind(success, failure) ~ sediment_mg_per_l + phosphorous_microM + copper_ug_per_l + salinity_psu + salinity_psu_sq + (1 | experiment) + (1 | rep), family=binomial, data, control=glmerControl(optimizer="bobyqa"))
# Note "rep" is needed to remove overdispersion; "experiment" not necessarily needed, but improves AIC and makes sense to keep.
sum(residuals(mod_fert_final, type="pearson")^2)/df.residual(mod_fert_final)

drop1(mod_fert_final, test="Chisq")
summary(mod_fert_final)


# PLOTS FOR FERTILISATION
# COPPER
par(mfrow=c(2,2))

ss <- seq(1, max(data$copper_ug_per_l), 1)
newdat <- expand.grid(sediment_mg_per_l=0, phosphorous_microM=0, copper_ug_per_l = ss, salinity_psu = 35, salinity_psu_sq = 35^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_final),newdat)

newdat$success <- mm %*% fixef(mod_fert_final)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_final),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_final)$experiment[1] + VarCorr(mod_fert_final)$rep[1]  ## must be adapted for more complex models

plot(data$copper_ug_per_l, data$mean_value_prop, xlab="Copper (?g/L)", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)

# SALINITY
ss <- seq(1, max(data$salinity_psu), 0.1)
newdat <- expand.grid(sediment_mg_per_l=0, phosphorous_microM=0, copper_ug_per_l = 0, salinity_psu = ss, salinity_psu_sq = 35^2, success=0, failure=0)
newdat$salinity_psu_sq = newdat$salinity_psu^2
mm <- model.matrix(terms(mod_fert_final),newdat)

newdat$success <- mm %*% fixef(mod_fert_final)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_final),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_final)$experiment[1] + VarCorr(mod_fert_final)$rep[1]  ## must be adapted for more complex models

plot(data$salinity_psu, data$mean_value_prop, xlab="Salinity (psu)", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)

# PHOSPHOROUS
ss <- seq(1, max(data$phosphorous_microM), 1)
newdat <- expand.grid(sediment_mg_per_l=0, phosphorous_microM=ss, copper_ug_per_l = 0, salinity_psu = 35, salinity_psu_sq = 35^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_final),newdat)

newdat$success <- mm %*% fixef(mod_fert_final)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_final),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_final)$experiment[1] + VarCorr(mod_fert_final)$rep[1]  ## must be adapted for more complex models

plot(data$phosphorous_microM, data$mean_value_prop, xlab="Phosphorous (?M)", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)

# SEDIMENT
ss <- seq(1, max(data$sediment_mg_per_l), 1)
newdat <- expand.grid(sediment_mg_per_l=ss, phosphorous_microM=0, copper_ug_per_l = 0, salinity_psu = 35, salinity_psu_sq = 35^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_final),newdat)

newdat$success <- mm %*% fixef(mod_fert_final)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_final),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_final)$experiment[1] + VarCorr(mod_fert_final)$rep[1]  ## must be adapted for more complex models

plot(data$sediment_mg_per_l, data$mean_value_prop, xlab="Sediment (mg/L)", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)



#SURVIVAL - MODEL - ammonium, copper, mercury, lead, salinity, sediment, acidification, temperature
# Seperate factors to be used
data2 <- with(dat2, data.frame(life.stage, ammonium_microM, copper_ug_per_l, mercury_ug._per_l, lead_ug_per_l, salinity_psu, acidification_pH, tempertaure_degrees_celcius))
                                                            
# Subset main dataset
data = dat[apply(!is.na(data2[, -1]), 1, sum) > 0 & data2$life.stage == "survivorship",]
# Add random observation variable for overdispersion
data$rep <- 1:dim(data)[1]

# This doesn't work because overparametrised... 
mod_surv_full <- glmer(cbind(success, failure) ~ ammonium_microM + copper_ug_per_l + mercury_ug._per_l + lead_ug_per_l + salinity_psu + salinity_psu_sq + acidification_pH + acidification_pH_sq + tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq + (1 | experiment) + (1 | rep), family=binomial, data)

# TRY with fewer variables -- remove ones we know not having effect
data2 <- with(dat2, data.frame(life.stage, copper_ug_per_l, lead_ug_per_l, salinity_psu, tempertaure_degrees_celcius))
data = dat[apply(!is.na(data2[, -1]), 1, sum) > 0 & data2$life.stage == "survivorship",]
data$rep <- 1:dim(data)[1]

# Works now, complains about rescaling some variables, but okay
mod_surv_full <- glmer(cbind(success, failure) ~ copper_ug_per_l + lead_ug_per_l + salinity_psu + salinity_psu_sq + tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq + (1 | experiment) + (1 | rep), family=binomial, data, control=glmerControl(optimizer="bobyqa"))
# Below should be less than one (or else overdispersion)
sum(residuals(mod_surv_full, type="pearson")^2)/df.residual(mod_surv_full)

drop1(mod_surv_full, test="Chisq")

# Drop some based on above for final model
data2 <- with(dat2, data.frame(life.stage, copper_ug_per_l, lead_ug_per_l, salinity_psu, tempertaure_degrees_celcius))
data = dat[apply(!is.na(data2[, -1]), 1, sum) > 0 & data2$life.stage == "survivorship",]
data$rep <- 1:dim(data)[1]


mod_surv_final <- glmer(cbind(success, failure) ~ copper_ug_per_l + lead_ug_per_l + salinity_psu + salinity_psu_sq + tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq + (1 | experiment) + (1 | rep), family=binomial, data, control=glmerControl(optimizer="bobyqa"))
# Note "rep" is needed to remove overdispersion; "experiment" not necessarily needed, but improves AIC and makes sense to keep.
sum(residuals(mod_surv_final, type="pearson")^2)/df.residual(mod_surv_final)

drop1(mod_surv_final, test="Chisq")
summary(mod_surv_final)


#PLOTS
par(mfrow=c(2,2))

# COPPER
ss <- seq(1, max(data$copper_ug_per_l), 1)
newdat <- expand.grid(copper_ug_per_l=ss, lead_ug_per_l=0, salinity_psu=35, salinity_psu_sq=35^2, tempertaure_degrees_celcius=29, tempertaure_degrees_celcius_sq=29^2, success=0, failure=0)
mm <- model.matrix(terms(mod_surv_final),newdat)

newdat$success <- mm %*% fixef(mod_surv_final)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_final),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_final)$experiment[1] + VarCorr(mod_surv_final)$rep[1]  ## must be adapted for more complex models

plot(data$copper_ug_per_l, data$mean_value_prop, xlab="Copper (?g/L)", ylab="Proportion of Larval Survivorship", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)

# LEAD
ss <- seq(1, max(data$lead_ug_per_l), 100)
newdat <- expand.grid(copper_ug_per_l=0, lead_ug_per_l=ss, salinity_psu=35, salinity_psu_sq=35^2, tempertaure_degrees_celcius=29, tempertaure_degrees_celcius_sq=29^2, success=0, failure=0)
mm <- model.matrix(terms(mod_surv_final),newdat)

newdat$success <- mm %*% fixef(mod_surv_final)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_final),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_final)$experiment[1] + VarCorr(mod_surv_final)$rep[1]  ## must be adapted for more complex models

plot(data$lead_ug_per_l, data$mean_value_prop, xlab="Lead (?g/L)", ylab="Proportion of Larval Survivorship", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)


# TEMPERATURE
ss <- seq(1, max(data$tempertaure_degrees_celcius), 0.1)
newdat <- expand.grid(copper_ug_per_l=0, lead_ug_per_l=0, salinity_psu=35, salinity_psu_sq=35^2, tempertaure_degrees_celcius=ss, success=0, failure=0)
newdat$tempertaure_degrees_celcius_sq = newdat$tempertaure_degrees_celcius^2
mm <- model.matrix(terms(mod_surv_final),newdat)

newdat$success <- mm %*% fixef(mod_surv_final)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_final),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_final)$experiment[1] + VarCorr(mod_surv_final)$rep[1]  ## must be adapted for more complex models

plot(data$tempertaure_degrees_celcius, data$mean_value_prop, xlab="Temperature (?C)", ylab="Proportion of Larval Survivorship", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)


# SALINITY
ss <- seq(1, max(data$salinity_psu), 0.1)
newdat <- expand.grid(copper_ug_per_l=0, lead_ug_per_l=0, salinity_psu=ss, tempertaure_degrees_celcius=29, tempertaure_degrees_celcius_sq=29^2, success=0, failure=0)
newdat$salinity_psu_sq = newdat$salinity_psu^2
mm <- model.matrix(terms(mod_surv_final),newdat)

newdat$success <- mm %*% fixef(mod_surv_final)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_final),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_final)$experiment[1] + VarCorr(mod_surv_final)$rep[1]  ## must be adapted for more complex models

plot(data$salinity_psu, data$mean_value_prop, xlab="Salinity (psu)", ylab="Proportion of Larval Survivorship", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)


# COMBINED MODEL

copper_store <- c()

for (cc in seq(0, 200, 1)) {

  
# p.fert <- predict(mod_fert_final, list(sediment_mg_per_l=0, phosphorous_microM=0, copper_ug_per_l=cc, salinity_psu=35, salinity_psu_sq = 35^2))

  newdat <- expand.grid(sediment_mg_per_l=0, phosphorous_microM=0, copper_ug_per_l = cc, salinity_psu = 35, salinity_psu_sq = 35^2, success=0, failure=0)
  mm <- model.matrix(terms(mod_fert_final), newdat)
  success_fert <- mm %*% fixef(mod_fert_final)
  pvar1_fert <- 2*sqrt(diag(mm %*% tcrossprod(vcov(mod_fert_final), mm)))

  # p.larv <- predict(mod_surv_final, list(copper_ug_per_l, lead_ug_per_l, tempertaure_degrees_celcius, tempertaure_degrees_celcius_sq = tempertaure_degrees_celcius^2, salinity_psu, salinity_psu_sq = salinity_psu^2), type="link", se.fit=TRUE)

  newdat <- expand.grid(copper_ug_per_l=cc, lead_ug_per_l=0, salinity_psu=35, salinity_psu_sq=35^2, tempertaure_degrees_celcius=29, tempertaure_degrees_celcius_sq=29^2, success=0, failure=0)
  mm <- model.matrix(terms(mod_surv_final),newdat)
  success_surv <- mm %*% fixef(mod_surv_final)
  pvar1_surv <- 2*sqrt(diag(mm %*% tcrossprod(vcov(mod_surv_final),mm)))

  # this part takes the mean expectation (and se) and generates 10000 random normally distributed 
  # samples for each of the two models and multiplies them together.
  # It then sorts these samples
  vars <- sort(inv.logit(rnorm(10000, success_fert, pvar1_fert)) * inv.logit(rnorm(10000, success_surv, pvar1_surv)))

  # vars[5000] # This will be the mean (median)
  # vars[250]  # This the lower 95% confidence interval
  # vars[9750] # The upper interval -- do you know why?

  copper_store <- rbind(copper_store, c(cc, vars[5000], vars[250], vars[9750]))

}

par(mfrow=c(1,1))

plot(copper_store[,1], copper_store[,2], xlab="Copper (?g/L)", ylab="Proportion of Successful Larvae", ylim=c(0, 1), type="l")
# "CI based on fixed-effects uncertainty ONLY"
# lines(copper_store[,1], copper_store[,3], lty=2)
# lines(copper_store[,1], copper_store[,4], lty=2)

polygon(c(copper_store[,1], rev(copper_store[,1])), c(copper_store[,3], rev(copper_store[,4])), col=rgb(0,0,0,0.2), border=NA)


# y25<-inv.logit(newdat$success-2*sqrt(pvar1))
# y97<-inv.logit(newdat$success+2*sqrt(pvar1))

# polygon(c(ss, ss[length(ss)], ss[length(ss):1], ss[1]), c(y25, y97[length(y97)], y97[length(y97):1], y25[1]), col=make.transparent('grey80', 0.5), border=NA)
# label(.05, 1.2, leg, font=2, xpd=NA)




# The mean above should be approximately the same as that caluclated directly form your models

inv.logit(p.fert$fit) * inv.logit(p.larv$fit)

#### SALINITY

salinity_store <- c()

for (ss in seq(18.4, 36.8, 0.1)) {

  # p.fert <- predict(mod_fert_final, list(sediment_mg_per_l=0, phosphorous_microM=0, copper_ug_per_l=cc, salinity_psu=35, salinity_psu_sq = 35^2))

  newdat <- expand.grid(sediment_mg_per_l=0, phosphorous_microM=0, copper_ug_per_l = 0, salinity_psu = ss, salinity_psu_sq = ss^2, success=0, failure=0)
  mm <- model.matrix(terms(mod_fert_final), newdat)
  success_fert <- mm %*% fixef(mod_fert_final)
  pvar1_fert <- 2*sqrt(diag(mm %*% tcrossprod(vcov(mod_fert_final), mm)))

  # p.larv <- predict(mod_surv_final, list(copper_ug_per_l, lead_ug_per_l, tempertaure_degrees_celcius, tempertaure_degrees_celcius_sq = tempertaure_degrees_celcius^2, salinity_psu, salinity_psu_sq = salinity_psu^2), type="link", se.fit=TRUE)

  newdat <- expand.grid(copper_ug_per_l=0, lead_ug_per_l=0, salinity_psu=ss, salinity_psu_sq=ss^2, tempertaure_degrees_celcius=29, tempertaure_degrees_celcius_sq=29^2, success=0, failure=0)
  mm <- model.matrix(terms(mod_surv_final),newdat)
  success_surv <- mm %*% fixef(mod_surv_final)
  pvar1_surv <- 2*sqrt(diag(mm %*% tcrossprod(vcov(mod_surv_final),mm)))

  # this part takes the mean expectation (and se) and generates 10000 random normally distributed 
  # samples for each of the two models and multiplies them together.
  # It then sorts these samples
  vars <- sort(inv.logit(rnorm(10000, success_fert, pvar1_fert)) * inv.logit(rnorm(10000, success_surv, pvar1_surv)))

  # vars[5000] # This will be the mean (median)
  # vars[250]  # This the lower 95% confidence interval
  # vars[9750] # The upper interval -- do you know why?

  salinity_store <- rbind(salinity_store, c(ss, vars[5000], vars[250], vars[9750]))

}

par(mfrow=c(1,1))

plot(salinity_store[,1], salinity_store[,2], xlab="Salinity (psu)", ylab="Proportion of Successful Larvae", ylim=c(0, 1), type="l")
# "CI based on fixed-effects uncertainty ONLY"
# lines(salinity_store[,1], salinity_store[,3], lty=2)
# lines(salinity_store[,1], salinity_store[,4], lty=2)
polygon(c(salinity_store[,1], rev(salinity_store[,1])), c(salinity_store[,3], rev(salinity_store[,4])), col=rgb(0,0,0,0.2), border=NA)


# y25<-inv.logit(newdat$success-2*sqrt(pvar1))
# y97<-inv.logit(newdat$success+2*sqrt(pvar1))

# polygon(c(ss, ss[length(ss)], ss[length(ss):1], ss[1]), c(y25, y97[length(y97)], y97[length(y97):1], y25[1]), col=make.transparent('grey80', 0.5), border=NA)
# label(.05, 1.2, leg, font=2, xpd=NA)



#LOCATIONS/WATER SAMPLES

#fertilisation model
water$salinity_g.l <- c(30, 31.5, 27.5)
water_fert <- data.frame(sediment_mg_per_l=water$suspended_solids_mg.l, 
                         phosphorous_microM=water$phosphorus_mg.l, 
                         copper_ug_per_l=water$copper_ug.l, 
                         salinity_psu=water$salinity_g.l, 
                         salinity_psu_sq =water$salinity_g.l^2, 
                         success=0, failure=0)

mm_fert <- model.matrix(terms(mod_fert_final), water_fert)

water_fert$success <- mm_fert %*% fixef(mod_fert_final)
pvar1 <- diag(mm_fert %*% tcrossprod(vcov(mod_fert_final), mm_fert))
tvar1 <- pvar1 + VarCorr(mod_fert_final)$experiment[1] + VarCorr(mod_fert_final)$rep[1]  
#bar plot
bp <- barplot(t(inv.logit(water_fert$success)), xlab="Location", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1), names.arg=water$sample)
#error bars
arrows(bp, inv.logit(water_fert$success-2*sqrt(pvar1)), bp, inv.logit(water_fert$success+2*sqrt(pvar1)), code=3, angle=90)


#survival model
 
