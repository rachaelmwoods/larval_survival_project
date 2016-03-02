library(lme4)
library(boot)
library(hier.part)
source("R/functions.R")

# Load datasets
dat <- read.csv("data/data_table_factors - data.csv", as.is=TRUE)
water <- read.csv("data/water_samples.csv", as.is=TRUE)

# Data preparation done in separate file
source("R/data_prep.R")


#######################
##FERTILISATION MODEL##
#######################
dat_fert <- dat[dat$life.stage == "fertilisation",]
dat_fert$rep <- 1:nrow(dat_fert)

par(mfrow=c(2,5))

#sediment
mod_fert <- glm(cbind(success, failure) ~ sediment_mg_per_l, family=binomial, data=dat_fert)
plot(mean_value_prop ~ sediment_mg_per_l, col = "dodgerblue", xlab = "Sediment (mg/L)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$sediment_mg_per_l)
lines(ss, predict(mod_fert, list(sediment_mg_per_l = ss), type="response"), lwd = 1)

#ammonium
mod_fert <- glm(cbind(success, failure) ~ ammonium_microM, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ ammonium_microM, col = "dodgerblue", xlab = "Ammonium (µM)", ylab = "",  pch=16, data=dat_fert)
ss <- sort(dat_fert$ammonium_microM)
lines(ss, predict(mod_fert, list(ammonium_microM = ss), type="response"))

#phosphate
mod_fert <- glm(cbind(success, failure) ~ phosphorous_microM, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ phosphorous_microM, col = "dodgerblue", xlab = "Phosphate (µM)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$phosphorous_microM)
lines(ss, predict(mod_fert, list(phosphorous_microM = ss), type="response"))

#copper
mod_fert <- glm(cbind(success, failure) ~ copper_ug_per_l, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ copper_ug_per_l, col = "dodgerblue", xlab = "Copper (µg/L)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$copper_ug_per_l)
lines(ss, predict(mod_fert, list(copper_ug_per_l = ss), type="response"))

#tributyltin_ug_per_l - NOT ENOUGH DATA

#zinc - REMOVE
mod_fert <- glm(cbind(success, failure) ~ zinc_ug_per_l, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ zinc_ug_per_l, col = "dodgerblue", xlab = "Zinc (µg/L)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$zinc_ug_per_l)
lines(ss, predict(mod_fert, list(zinc_ug_per_l = ss), type="response"))

#cadmium = VERY WEAK RELATOPNSHIP
mod_fert <- glm(cbind(success, failure) ~ cadmium_ug_per_l, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ cadmium_ug_per_l, col = "dodgerblue", xlab = "Cadmium (µg/L)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$cadmium_ug_per_l)
lines(ss, predict(mod_fert, list(cadmium_ug_per_l = ss), type="response"))
    
#nitrate_microM - REMOVE
mod_fert <- glm(cbind(success, failure) ~ nitrate_microM, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ nitrate_microM, col = "dodgerblue", xlab = "Nitrate (µM)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$nitrate_microM)
lines(ss, predict(mod_fert, list(nitrate_microM = ss), type="response"))

#salinity
mod_fert <- glm(cbind(success, failure) ~ salinity_psu + salinity_psu_sq, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ salinity_psu, col = "dodgerblue", xlab = "Salinity (psu)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$salinity_psu)
lines(ss, predict(mod_fert, list(salinity_psu = ss, salinity_psu_sq = ss^2), type="response"))

#acidification - VERY WEAK RESPONSE
mod_fert <- glm(cbind(success, failure) ~ acidification_pH + acidification_pH_sq, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ acidification_pH, col = "dodgerblue", xlab = "Acidification (pH)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$acidification_pH)
lines(ss, predict(mod_fert, list(acidification_pH = ss, acidification_pH_sq = ss^2), type="response"))
drop1(mod_fert, test="Chisq")

#tempertaure - REMOVE not enough data
mod_fert <- glm(cbind(success, failure) ~ tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ tempertaure_degrees_celcius, col = "dodgerblue", xlab = "Temperature (C)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$tempertaure_degrees_celcius)
lines(ss, predict(mod_fert, list(tempertaure_degrees_celcius = ss, tempertaure_degrees_celcius_sq = ss^2), type="response"))

title(ylab = "Some Values", outer = TRUE, line = 3)


## FULL MODEL
dat_fert2 <- dat_fert[c("success", "failure", "sediment_mg_per_l", "ammonium_microM", "phosphorous_microM", "copper_ug_per_l", "salinity_psu", "salinity_psu_sq", "experiment", "mean_value_prop", "spawn.brood")]
dat_fert2 <- dat_fert2[!(is.na(dat_fert2$sediment_mg_per_l) & is.na(dat_fert2$ammonium_microM) & is.na(dat_fert2$phosphorous_microM) & is.na(dat_fert2$copper_ug_per_l) & is.na(dat_fert2$salinity_psu)),]

dat_fert2$sediment_mg_per_l[is.na(dat_fert2$sediment_mg_per_l)] <- log10(1)
dat_fert2$ammonium_microM[is.na(dat_fert2$ammonium_microM)] <- log10(0.01391)
dat_fert2$phosphorous_microM[is.na(dat_fert2$phosphorous_microM)] <- log10(0.446)
dat_fert2$copper_ug_per_l[is.na(dat_fert2$copper_ug_per_l)] <- log10(0.9)

dat_fert2$salinity_psu[is.na(dat_fert2$salinity_psu)] <- 34
#square salinity because it has a quadratic response - because it decreases on both sides from 35psu
dat_fert2$salinity_psu_sq <- dat_fert2$salinity_psu^2

dat_fert2$rep <- 1:nrow(dat_fert2)

##MODEL##
mod_fert_full <- glm(cbind(success, failure) ~ sediment_mg_per_l + ammonium_microM + phosphorous_microM + copper_ug_per_l + salinity_psu + salinity_psu_sq, family=binomial, data=dat_fert2)
summary(mod_fert_full)
drop1(mod_fert_full, test="Chisq")

mod_fert_full <- glmer(cbind(success, failure) ~ sediment_mg_per_l + ammonium_microM + phosphorous_microM + copper_ug_per_l + salinity_psu + salinity_psu_sq + (1 | experiment) + (1 | rep), family=binomial, data=dat_fert2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e4)))
summary(mod_fert_full)

# Below should be less than one (or else overdispersion)
sum(residuals(mod_fert_full, type="pearson")^2)/df.residual(mod_fert_full)

drop1(mod_fert_full, test="Chisq")

##PLOTS FOR FERTILISATION##
par(mfrow=c(3,2))

# COPPER
ss <- seq(min(dat_fert2$copper_ug_per_l), max(dat_fert2$copper_ug_per_l), 0.05)
newdat <- expand.grid(sediment_mg_per_l=log10(1), ammonium_microM=log10(0.01391), phosphorous_microM=log10(0.446), copper_ug_per_l = ss, salinity_psu = 34, salinity_psu_sq = 34^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_full), newdat)

newdat$success <- mm %*% fixef(mod_fert_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_full),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1]  ## must be adapted for more complex models

plot(dat_fert2$copper_ug_per_l, dat_fert2$mean_value_prop, xlab="Copper", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)

# SEDIMENT
ss <- seq(min(dat_fert2$sediment_mg_per_l), max(dat_fert2$sediment_mg_per_l), 0.01)
newdat <- expand.grid(sediment_mg_per_l=ss, ammonium_microM=log10(0.01391), phosphorous_microM=log10(0.446), copper_ug_per_l = log10(0.9), salinity_psu = 34, salinity_psu_sq = 34^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_full),newdat)

newdat$success <- mm %*% fixef(mod_fert_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_full),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1]  ## must be adapted for more complex models

plot(dat_fert2$sediment_mg_per_l, dat_fert2$mean_value_prop, xlab="Sediment", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)

# AMMONIUM
ss <- seq(min(dat_fert2$ammonium_microM), max(dat_fert2$ammonium_microM), 0.05)
newdat <- expand.grid(sediment_mg_per_l=log10(1), ammonium_microM=ss, phosphorous_microM=log10(0.446), copper_ug_per_l = log10(0.9), salinity_psu = 34, salinity_psu_sq = 34^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_full),newdat)

newdat$success <- mm %*% fixef(mod_fert_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_full),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1]  ## must be adapted for more complex models

plot(dat_fert2$ammonium_microM, dat_fert2$mean_value_prop, xlab="Ammonium", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)

# PHOSPHOROUS
ss <- seq(min(dat_fert2$phosphorous_microM), max(dat_fert2$phosphorous_microM), 0.05)
newdat <- expand.grid(sediment_mg_per_l=log10(1), ammonium_microM=log10(0.01391), phosphorous_microM=ss, copper_ug_per_l = log10(0.9), salinity_psu = 34, salinity_psu_sq = 34^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_full),newdat)

newdat$success <- mm %*% fixef(mod_fert_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_full),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1]  ## must be adapted for more complex models

plot(dat_fert2$phosphorous_microM, dat_fert2$mean_value_prop, xlab="Phosphorous", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)

# SALINITY
ss <- seq(min(dat_fert2$salinity_psu), max(dat_fert2$salinity_psu), 0.1)
newdat <- expand.grid(sediment_mg_per_l=log10(1), ammonium_microM=log10(0.01391), phosphorous_microM=log10(0.446), copper_ug_per_l = log10(0.9), salinity_psu = ss, salinity_psu_sq = 0, success=0, failure=0)
newdat$salinity_psu_sq = newdat$salinity_psu^2
mm <- model.matrix(terms(mod_fert_full),newdat)

newdat$success <- mm %*% fixef(mod_fert_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_full),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1]  ## must be adapted for more complex models

plot(dat_fert2$salinity_psu, dat_fert2$mean_value_prop, xlab="Salinity", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
# "CI based on fixed-effects uncertainty ONLY"
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)
# "CI based on fixed-effects uncertainty + random-effects variance"
# lines(ss, inv.logit(newdat$success-2*sqrt(tvar1)), col="red", lty=2)
# lines(ss, inv.logit(newdat$success+2*sqrt(tvar1)), col="red", lty=2)


##################
##SURVIVAL MODEL##
##################
dat_surv <- dat[dat$life.stage == "survivorship",]
dat_surv$rep <- 1:nrow(dat_surv)

par(mfrow=c(1,1))

#ammonium - few points
mod_surv <- glm(cbind(success, failure) ~ ammonium_microM, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ ammonium_microM,  col = "dodgerblue", xlab = "Ammonium (µM)", ylab = "",  pch=16, las=0, data=dat_surv)
ss <- sort(dat_surv$ammonium_microM)
lines(ss, predict(mod_surv, list(ammonium_microM = ss), type="response"))

#copper
mod_surv <- glm(cbind(success, failure) ~ copper_ug_per_l, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ copper_ug_per_l,  col = "dodgerblue", xlab = "Copper (µg/L)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$copper_ug_per_l)
lines(ss, predict(mod_surv, list(copper_ug_per_l = ss), type="response"))

# mercury- REMOVE
mod_surv <- glm(cbind(success, failure) ~ mercury_ug._per_l, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ mercury_ug._per_l,  col = "dodgerblue", xlab = "Mercury (µg/L)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$mercury_ug._per_l)
lines(ss, predict(mod_surv, list(mercury_ug._per_l = ss), type="response"))

# lead
mod_surv <- glm(cbind(success, failure) ~ lead_ug_per_l, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ lead_ug_per_l,  col = "dodgerblue", xlab = "Lead (µg/L)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$lead_ug_per_l)
lines(ss, predict(mod_surv, list(lead_ug_per_l = ss), type="response"))

#salinity - REMOVE
mod_surv <- glm(cbind(success, failure) ~ salinity_psu + salinity_psu_sq, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ salinity_psu,  col = "dodgerblue", xlab = "Salinity (psu)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$salinity_psu)
lines(ss, predict(mod_surv, list(salinity_psu = ss, salinity_psu_sq = ss^2), type="response"))

drop1(mod_surv, test="Chisq")

mod_surv2 <- glm(cbind(success, failure) ~ salinity_psu, family=binomial, data=dat_surv)
plot(mean_value_prop ~ salinity_psu, col = "dodgerblue", xlab = "Salinity (psu)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$salinity_psu)
lines(ss, predict(mod_surv2, list(salinity_psu = ss), type="response"))

drop1(mod_surv2, test="Chisq")

# acidification - VERY WEAK
mod_surv <- glm(cbind(success, failure) ~ acidification_pH + acidification_pH_sq, family=binomial, data=dat_surv)
summary(mod_surv)
drop1(mod_surv, test="Chisq")
plot(mean_value_prop ~ acidification_pH,  col = "dodgerblue", xlab = "Acidification (pH)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$acidification_pH)
lines(ss, predict(mod_surv, list(acidification_pH = ss, acidification_pH_sq = ss^2), type="response"))

# Drop squared part
mod_surv2 <- glm(cbind(success, failure) ~ acidification_pH, family=binomial, data=dat_surv)
summary(mod_surv2)
plot(mean_value_prop ~ acidification_pH,  col = "dodgerblue", xlab = "Temperature (C)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$acidification_pH)
lines(ss, predict(mod_surv2, list(acidification_pH = ss), type="response"))

drop1(mod_surv2, test="Chisq")


# temperature
mod_surv <- glm(cbind(success, failure) ~ tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ tempertaure_degrees_celcius,  col = "dodgerblue", xlab = "Temperature (C)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$tempertaure_degrees_celcius)
lines(ss, predict(mod_surv, list(tempertaure_degrees_celcius = ss, tempertaure_degrees_celcius_sq = ss^2), type="response"))



##FULL MODEL##, "acidification_pH"   
dat_surv2 <- dat_surv[c("success", "failure", "copper_ug_per_l", "lead_ug_per_l", "salinity_psu", "tempertaure_degrees_celcius", "acidification_pH", "experiment", "mean_value_prop", "spawn.brood")]
dat_surv2 <- dat_surv2[!(is.na(dat_surv2$copper_ug_per_l) & is.na(dat_surv2$lead_ug_per_l) & is.na(dat_surv2$salinity_psu) & is.na(dat_surv2$acidification_pH) & is.na(dat_surv2$tempertaure_degrees_celcius)),]


dat_surv2$copper_ug_per_l[is.na(dat_surv2$copper_ug_per_l)] <- log10(0.9)
dat_surv2$lead_ug_per_l[is.na(dat_surv2$lead_ug_per_l)] <- log10(0.03)
dat_surv2$salinity_psu[is.na(dat_surv2$salinity_psu)] <- 34

dat_surv2$acidification_pH[is.na(dat_surv2$acidification_pH)] <- 8.1

dat_surv2$tempertaure_degrees_celcius[is.na(dat_surv2$tempertaure_degrees_celcius)] <- 28
#square salinity because it has a quadratic response - because it decreases on both sides from 35psu
dat_surv2$tempertaure_degrees_celcius_sq <- dat_surv2$tempertaure_degrees_celcius^2


dat_surv2$rep <- 1:nrow(dat_surv2)

# # Re-scaling data
# dat_surv3 <- dat_surv2
# dat_surv3[c("copper_ug_per_l", "lead_ug_per_l", "acidification_pH", "tempertaure_degrees_kelvin", "tempertaure_degrees_kelvin_sq")] <- scale(dat_surv3[c("copper_ug_per_l", "lead_ug_per_l", "acidification_pH", "tempertaure_degrees_kelvin", "tempertaure_degrees_kelvin_sq")])

dat_surv2$experiment <- factor(dat_surv2$experiment)
dat_surv2$rep <- factor(dat_surv2$rep)

## MODEL##
mod_surv_full <- glm(cbind(success, failure) ~ copper_ug_per_l + lead_ug_per_l + salinity_psu + acidification_pH + tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq, family=binomial, data=dat_surv2)
summary(mod_surv_full)
drop1(mod_surv_full, test="Chisq")



mod_surv_full <- glmer(cbind(success, failure) ~ copper_ug_per_l + lead_ug_per_l + salinity_psu + acidification_pH + tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq + (1 | experiment) + (1 | rep), family=binomial, data=dat_surv2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e4)))

# mod_surv_full <- glmer(cbind(success, failure) ~  tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq + (1 | experiment) + (1 | rep), family=binomial, data=dat_surv2)


summary(mod_surv_full)

# Below should be less than one (or else overdispersion)
sum(residuals(mod_surv_full, type="pearson")^2)/df.residual(mod_surv_full)

drop1(mod_surv_full, test="Chisq")

mod_surv_full <- glmer(cbind(success, failure) ~ copper_ug_per_l + lead_ug_per_l + salinity_psu + (1 | experiment) + (1 | rep), family=binomial, data=dat_surv2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e4)))

summary(mod_surv_full)
drop1(mod_surv_full, test="Chisq")


##PLOTS FOR SURVIVAL##
par(mfrow=c(2,3))
    

# COPPER
ss <- seq(min(dat_surv2$copper_ug_per_l), max(dat_surv2$copper_ug_per_l), 0.05)
newdat <- expand.grid(copper_ug_per_l = ss, lead_ug_per_l=log10(0.03), salinity_psu = 34, success=0, failure=0)
mm <- model.matrix(terms(mod_surv_full), newdat)

newdat$success <- mm %*% fixef(mod_surv_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_full),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1]  

plot(dat_surv2$copper_ug_per_l, dat_surv2$mean_value_prop, col = "dodgerblue", xlab="Copper", ylab="", pch=16, las=1, ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1)
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)


# LEAD
ss <- seq(min(dat_surv2$lead_ug_per_l), max(dat_surv2$lead_ug_per_l), 0.05)
newdat <- expand.grid(copper_ug_per_l=log10(0.9), lead_ug_per_l = ss, salinity_psu = 34, success=0, failure=0)
mm <- model.matrix(terms(mod_surv_full), newdat)

newdat$success <- mm %*% fixef(mod_surv_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_full),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1]  

plot(dat_surv2$copper_ug_per_l, dat_surv2$mean_value_prop, xlab="Lead", ylab="", col = "dodgerblue", pch=16, las=1,ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)

#SALINITY
ss <- seq(min(dat_surv2$salinity_psu), max(dat_surv2$salinity_psu), 0.2)
newdat <- expand.grid(copper_ug_per_l=log10(0.9), lead_ug_per_l = log10(0.03), salinity_psu=ss, success=0, failure=0)
mm <- model.matrix(terms(mod_surv_full),newdat)

newdat$success <- mm %*% fixef(mod_surv_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_full),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1] 

plot(dat_fert2$salinity_psu, dat_fert2$mean_value_prop, xlab="Temperature", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)


###################
#VARIANCE ANALYSIS#
###################

par(mfrow=c(1,2))

dat_fert = dat[apply(!is.na(dat[, -1]), 1, sum) > 0 & dat$life.stage == "fertilisation",]
dat_surv = dat[apply(!is.na(dat[, -1]), 1, sum) > 0 & dat$life.stage == "survivorship",]

factors <- data_fert[, c(17,12,14,15,21)]
colnames(factors) <- c("Copper", "Sediment","Ammonium", "Phosphate", "Salinity") 
hier.part(data_fert$mean_value_prop, factors, family = "binomial", gof = "logLik", barplot = TRUE)


factors <- data_surv[, c(17,27,21)]
colnames(factors) <- c("Copper", "Lead", "Salinity" ) 
hier.part(data_surv$mean_value_prop, factors, family = "binomial", gof = "logLik", barplot = TRUE)




###########################
##LOCATIONS/WATER SAMPLES##
############################

##Fertilisation Model##
water_fert <- data.frame(sediment_mg_per_l=log10(water$suspended_solids_mg.l), 
                         copper_ug_per_l=log10(water$copper_ug.l),
                         ammonium_microM=log10(water$ammonia_mg.l),
                         phosphorous_microM=log10(water$phosphorus_mg.l), 
                         salinity_psu=water$salinity_g.l, 
                         salinity_psu_sq =water$salinity_g.l^2, 
                         success=0, failure=0)

mm_fert <- model.matrix(terms(mod_fert_full), water_fert)

water_fert$success <- mm_fert %*% fixef(mod_fert_full)
pvar1 <- diag(mm_fert %*% tcrossprod(vcov(mod_fert_full), mm_fert))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1]  
#bar plot
bp <- barplot(t(inv.logit(water_fert$success)), xlab="Location", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1), names.arg=water$sample)
#error bars
arrows(bp, inv.logit(water_fert$success-2*sqrt(pvar1)), bp, inv.logit(water_fert$success+2*sqrt(pvar1)), code=3, angle=90)


##Survival Model##
water_surv <- data.frame(copper_ug_per_l=log10(water$copper_ug.l), 
                         lead_ug_per_l=log10(water$lead_ug.l),
                         salinity_psu=water$salinity_g.l,
                         success=0, failure=0)

mm_surv <- model.matrix(terms(mod_surv_full), water_surv)

water_surv$success <- mm_surv %*% fixef(mod_surv_full)
pvar1 <- diag(mm_surv %*% tcrossprod(vcov(mod_surv_full), mm_surv))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1]  
#bar plot
bp <- barplot(t(inv.logit(water_surv$success)), xlab="Location", ylab="Proportion of Larvae Survived", ylim=c(0, 1), names.arg=water$sample)
#error bars
arrows(bp, inv.logit(water_surv$success-2*sqrt(pvar1)), bp, inv.logit(water_surv$success+2*sqrt(pvar1)), code=3, angle=90)



##################
##COMBINED MODEL##
##################

combined_store <- c()

for (cc in 1:3) {

  temp_fert <- water_fert[cc,]
  mm_fert <- model.matrix(terms(mod_fert_full), temp_fert)
  temp_fert$success <- mm_fert %*% fixef(mod_fert_full)
  pvar1_fert <- diag(mm_fert %*% tcrossprod(vcov(mod_fert_full), mm_fert))


  temp_surv <- water_surv[cc,]
  mm_surv <- model.matrix(terms(mod_surv_full), temp_surv)
  temp_surv$success <- mm_surv %*% fixef(mod_surv_full)
  pvar1_surv <- diag(mm_surv %*% tcrossprod(vcov(mod_surv_full), mm_surv))


vars <- sort(inv.logit(rnorm(10000, temp_fert$success, pvar1_fert)) * inv.logit(rnorm(10000, temp_surv$success, pvar1_surv)))

combined_store <- rbind(combined_store, c(cc, vars[5000], vars[250], vars[9750]))

}

bp <- barplot(combined_store[,2], xlab="Location", ylab="Proportion of Succesful Larvae", ylim=c(0, 1), names.arg=water$sample)
#error bars
arrows(bp, combined_store[,3], bp, combined_store[,4], code=3, angle=90)


inv.logit(p.fert$fit) * inv.logit(p.larv$fit)

