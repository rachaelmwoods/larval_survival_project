library(lme4)
library(boot)
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
summary(mod_fert)
plot(mean_value_prop ~ sediment_mg_per_l, data=dat_fert)
ss <- sort(dat_fert$sediment_mg_per_l)
lines(ss, predict(mod_fert, list(sediment_mg_per_l = ss), type="response"))

#ammonium
mod_fert <- glm(cbind(success, failure) ~ ammonium_microM, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ ammonium_microM, data=dat_fert)
ss <- sort(dat_fert$ammonium_microM)
lines(ss, predict(mod_fert, list(ammonium_microM = ss), type="response"))

#phosphorous
mod_fert <- glm(cbind(success, failure) ~ phosphorous_microM, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ phosphorous_microM, data=dat_fert)
ss <- sort(dat_fert$phosphorous_microM)
lines(ss, predict(mod_fert, list(phosphorous_microM = ss), type="response"))

#copper
mod_fert <- glm(cbind(success, failure) ~ copper_ug_per_l, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ copper_ug_per_l, data=dat_fert)
ss <- sort(dat_fert$copper_ug_per_l)
lines(ss, predict(mod_fert, list(copper_ug_per_l = ss), type="response"))

#tributyltin_ug_per_l - REMOVE not enough data!!

#zinc - REMOVE
mod_fert <- glm(cbind(success, failure) ~ zinc_ug_per_l, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ zinc_ug_per_l, data=dat_fert)
ss <- sort(dat_fert$zinc_ug_per_l)
lines(ss, predict(mod_fert, list(zinc_ug_per_l = ss), type="response"))

#cadmium = VERY WEAK RELATOPNSHIP
mod_fert <- glm(cbind(success, failure) ~ cadmium_ug_per_l, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ cadmium_ug_per_l, data=dat_fert)
ss <- sort(dat_fert$cadmium_ug_per_l)
lines(ss, predict(mod_fert, list(cadmium_ug_per_l = ss), type="response"))

#nitrate_microM - REMOVE
mod_fert <- glm(cbind(success, failure) ~ nitrate_microM, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ nitrate_microM, data=dat_fert)
ss <- sort(dat_fert$nitrate_microM)
lines(ss, predict(mod_fert, list(nitrate_microM = ss), type="response"))

#salinity
mod_fert <- glm(cbind(success, failure) ~ salinity_psu + salinity_psu_sq, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ salinity_psu, data=dat_fert)
ss <- sort(dat_fert$salinity_psu)
lines(ss, predict(mod_fert, list(salinity_psu = ss, salinity_psu_sq = ss^2), type="response"))

#acidification - VERY WEAK RESPONSE
mod_fert <- glm(cbind(success, failure) ~ acidification_pH + acidification_pH_sq, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ acidification_pH, data=dat_fert)
ss <- sort(dat_fert$acidification_pH)
lines(ss, predict(mod_fert, list(acidification_pH = ss, acidification_pH_sq = ss^2), type="response"))
drop1(mod_fert, test="Chisq")

#tempertaure - REMOVE not enough data
mod_fert <- glm(cbind(success, failure) ~ tempertaure_degrees_kelvin, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ tempertaure_degrees_kelvin, data=dat_fert)


## FULL MODEL
dat_fert2 <- dat_fert[c("success", "failure", "sediment_mg_per_l", "ammonium_microM", "phosphorous_microM", "copper_ug_per_l", "salinity_psu", "salinity_psu_sq", "experiment", "mean_value_prop")]
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

mod_fert_full <- glmer(cbind(success, failure) ~ sediment_mg_per_l + ammonium_microM + phosphorous_microM + copper_ug_per_l + salinity_psu + salinity_psu_sq + (1 | experiment) + (1 | rep), family=binomial, data=dat_fert2)

summary(mod_fert_full)

# Below should be less than one (or else overdispersion)
sum(residuals(mod_fert_full, type="pearson")^2)/df.residual(mod_fert_full)

drop1(mod_fert_full, test="Chisq")

##PLOTS FOR FERTILISATION##
par(mfrow=c(5,2))

# COPPER
ss <- seq(min(dat_fert2$copper_ug_per_l), max(dat_fert2$copper_ug_per_l), 0.05)
newdat <- expand.grid(sediment_mg_per_l=log10(1), ammonium_microM=log10(0.01391), phosphorous_microM=log10(0.446), copper_ug_per_l = ss, salinity_psu = 34, salinity_psu_sq = 34^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_full), newdat)

newdat$success <- mm %*% fixef(mod_fert_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_full),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1]  ## must be adapted for more complex models

plot(dat_fert2$copper_ug_per_l, dat_fert2$mean_value_prop, xlab="Copper (?g/L)", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
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

plot(dat_fert2$sediment_mg_per_l, dat_fert2$mean_value_prop, xlab="Copper (?g/L)", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
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

plot(dat_fert2$ammonium_microM, dat_fert2$mean_value_prop, xlab="Copper (?g/L)", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
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

plot(dat_fert2$phosphorous_microM, dat_fert2$mean_value_prop, xlab="Copper (?g/L)", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
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

plot(dat_fert2$salinity_psu, dat_fert2$mean_value_prop, xlab="Copper (?g/L)", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
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

#ammonium - few points
mod_surv <- glm(cbind(success, failure) ~ ammonium_microM, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ ammonium_microM, data=dat_surv)
ss <- sort(dat_surv$ammonium_microM)
lines(ss, predict(mod_surv, list(ammonium_microM = ss), type="response"))

#copper
mod_surv <- glm(cbind(success, failure) ~ copper_ug_per_l, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ copper_ug_per_l, data=dat_surv)
ss <- sort(dat_surv$copper_ug_per_l)
lines(ss, predict(mod_surv, list(copper_ug_per_l = ss), type="response"))

# mercury- REMOVE
mod_surv <- glm(cbind(success, failure) ~ mercury_ug._per_l, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ mercury_ug._per_l, data=dat_surv)
ss <- sort(dat_surv$mercury_ug._per_l)
lines(ss, predict(mod_surv, list(mercury_ug._per_l = ss), type="response"))

# lead
mod_surv <- glm(cbind(success, failure) ~ lead_ug_per_l, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ lead_ug_per_l, data=dat_surv)
ss <- sort(dat_surv$lead_ug_per_l)
lines(ss, predict(mod_surv, list(lead_ug_per_l = ss), type="response"))

#salinity - REMOVE
mod_surv <- glm(cbind(success, failure) ~ salinity_psu + salinity_psu_sq, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ salinity_psu, data=dat_surv)
ss <- sort(dat_surv$salinity_psu)
lines(ss, predict(mod_surv, list(salinity_psu = ss, salinity_psu_sq = ss^2), type="response"))

# acidification - VERY WEAK
mod_surv <- glm(cbind(success, failure) ~ acidification_pH + acidification_pH_sq, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ acidification_pH, data=dat_surv)
ss <- sort(dat_surv$acidification_pH)
lines(ss, predict(mod_surv, list(acidification_pH = ss, acidification_pH_sq = ss^2), type="response"))

# temperature
mod_surv <- glm(cbind(success, failure) ~ tempertaure_degrees_kelvin + tempertaure_degrees_kelvin_sq, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ tempertaure_degrees_kelvin, data=dat_surv)
ss <- sort(dat_surv$tempertaure_degrees_kelvin)
lines(ss, predict(mod_surv, list(tempertaure_degrees_kelvin = ss, tempertaure_degrees_kelvin_sq = ss^2), type="response"))


##FULL MODEL##
dat_surv2 <- dat_surv[c("success", "failure", "ammonium_microM", "copper_ug_per_l", "lead_ug_per_l", "salinity_psu", "salinity_psu_sq", "acidification_pH", "acidification_pH_sq", "tempertaure_degrees_kelvin", "tempertaure_degrees_kelvin_sq", "experiment", "mean_value_prop")]
dat_surv2 <- dat_surv2[!(is.na(dat_surv2$ammonium_microM) & is.na(dat_surv2$copper_ug_per_l) & is.na(dat_surv2$lead_ug_per_l) & is.na(dat_surv2$salinity_psu) & is.na(dat_surv2$acidification_pH) & is.na(dat_surv2$tempertaure_degrees_kelvin)),]


dat_surv2$ammonium_microM[is.na(dat_surv2$ammonium_microM)] <- log10(0.01391)
dat_surv2$copper_ug_per_l[is.na(dat_surv2$copper_ug_per_l)] <- log10(0.9)
dat_surv2$lead_ug_per_l[is.na(dat_surv2$lead_ug_per_l)] <- log10(0.03)

dat_surv2$salinity_psu[is.na(dat_surv2$salinity_psu)] <- 34
#square salinity because it has a quadratic response - because it decreases on both sides from 35psu
dat_surv2$salinity_psu_sq <- dat_surv2$salinity_psu^2

dat_surv2$acidification_pH[is.na(dat_surv2$acidification_pH)] <- 8.1
#square salinity because it has a quadratic response - because it decreases on both sides from 35psu
dat_surv2$acidification_pH <- dat_surv2$acidification_pH^2

dat_surv2$tempertaure_degrees_kelvin[is.na(dat_surv2$tempertaure_degrees_kelvin)] <- 301
#square salinity because it has a quadratic response - because it decreases on both sides from 35psu
dat_surv2$salinity_psu_sq <- dat_surv2$tempertaure_degrees_kelvin^2


dat_surv2$rep <- 1:nrow(dat_surv2)

## MODEL##
mod_surv_full <- glm(cbind(success, failure) ~ ammonium_microM + copper_ug_per_l + lead_ug_per_l + salinity_psu + salinity_psu_sq + acidification_pH + acidification_pH_sq + tempertaure_degrees_kelvin + tempertaure_degrees_kelvin_sq, family=binomial, data=dat_surv2)
summary(mod_surv_full)
drop1(mod_surv_full, test="Chisq")


mod_surv_full <- glmer(cbind(success, failure) ~ ammonium_microM + copper_ug_per_l + lead_ug_per_l + salinity_psu + salinity_psu_sq + acidification_pH + acidification_pH_sq + tempertaure_degrees_kelvin + tempertaure_degrees_kelvin_sq + (1 | experiment) + (1 | rep), family=binomial, data=dat_surv2)

summary(mod_surv_full)

# Below should be less than one (or else overdispersion)
sum(residuals(mod_surv_full, type="pearson")^2)/df.residual(mod_surv_full)

drop1(mod_surv_full, test="Chisq")


##PLOTS FOR SURVIVAL##
par(mfrow=c(2,2))

# COPPER
ss <- seq(min(dat_surv2$copper_ug_per_l), max(dat_surv2$copper_ug_per_l), 0.05)
newdat <- expand.grid(copper_ug_per_l = ss, lead_ug_per_l=log10(0.03), salinity_psu = 34, salinity_psu_sq = 34^2, tempertaure_degrees_kelvin = 301, tempertaure_degrees_kelvin_sq = 301^2 success=0, failure=0)
mm <- model.matrix(terms(mod_surv_full), newdat)

newdat$success <- mm %*% fixef(mod_surv_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_full),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1]  

plot(dat_surv2$copper_ug_per_l, dat_surv2$mean_value_prop, xlab="Copper (?g/L)", ylab="Proportion of Larvae Survived", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)

# LEAD
ss <- seq(min(dat_surv2$lead_ug_per_l), max(dat_surv2$lead_ug_per_l), 0.05)
newdat <- expand.grid(lead_ug_per_l = ss, copper_ug_per_l=log10(0.9), salinity_psu = 34, salinity_psu_sq = 34^2, tempertaure_degrees_kelvin = 301, tempertaure_degrees_kelvin_sq = 301^2 success=0, failure=0)
mm <- model.matrix(terms(mod_surv_full), newdat)

newdat$success <- mm %*% fixef(mod_surv_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_full),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1]  

plot(dat_surv2$copper_ug_per_l, dat_surv2$mean_value_prop, xlab="Lead", ylab="Proportion of Larvae Survived", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)

#TEMPERATURE
ss <- seq(min(dat_surv2$tempertaure_degrees_kelvin), max(dat_surv2$tempertaure_degrees_kelvin), 0.1)
newdat <- expand.grid(lead_ug_per_l = 0.03, copper_ug_per_l=log10(0.9), salinity_psu = 34, salinity_psu_sq = 34^2, tempertaure_degrees_kelvin = ss, tempertaure_degrees_kelvin_sq = 0, success=0, failure=0)
newdat$salinity_psu_sq = newdat$tempertaure_degrees_kelvin^2
mm <- model.matrix(terms(mod_surv_full),newdat)

newdat$success <- mm %*% fixef(mod_surv_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_full),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1] 

plot(dat_fert2$tempertaure_degrees_kelvin, dat_fert2$mean_value_prop, xlab="Sailinity", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)

#SALINITY
ss <- seq(min(dat_surv2$salinity_psu), max(dat_surv2$salinity_psu), 0.1)
newdat <- expand.grid(lead_ug_per_l = 0.03, copper_ug_per_l=log10(0.9), salinity_psu = ss, salinity_psu_sq = 0, tempertaure_degrees_kelvin = 301, tempertaure_degrees_kelvin_sq = 301^2, success=0, failure=0)
newdat$salinity_psu_sq = newdat$salinity_psu^2
mm <- model.matrix(terms(mod_surv_full),newdat)

newdat$success <- mm %*% fixef(mod_surv_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_full),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1] 

plot(dat_fert2$salinity_psu, dat_fert2$mean_value_prop, xlab="Sailinity", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1))
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
lines(ss, inv.logit(newdat$success-2*sqrt(pvar1)), lty=2)
lines(ss, inv.logit(newdat$success+2*sqrt(pvar1)), lty=2)

###########################
##LOCATIONS/WATER SAMPLES##
############################

##Fertilisation Model##
water_fert <- data.frame(sediment_mg_per_l=water$suspended_solids_mg.l, 
                         copper_ug_per_l=water$copper_ug.l,
                         ammonium_microM=water$ammonia_mg/l,
                         phosphorous_microM=water$phosphorus_mg.l, 
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

water_surv <- data.frame(copper_ug_per_l=water$copper_ug.l, 
                         lead_ug_per_l=water$lead_ug/l,
                         tempertaure_degrees_kelvin=water$temperature_degrees_celcius,
                         tempertaure_degrees_kelvin_sq=water$temperature_degrees_celcius^2
                         salinity_psu=water$salinity_g.l, 
                         salinity_psu_sq =water$salinity_g.l^2, 
                         success=0, failure=0)

mm_surv <- model.matrix(terms(mod_surv_full), water_surv)

water_surv$success <- mm_fert %*% fixef(mod_surv_full)
pvar1 <- diag(mm_fert %*% tcrossprod(vcov(mod_surv_full), mm_surv))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1]  
#bar plot
bp <- barplot(t(inv.logit(water_surv$success)), xlab="Location", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1), names.arg=water$sample)
#error bars
arrows(bp, inv.logit(water_surv$success-2*sqrt(pvar1)), bp, inv.logit(water_surv$success+2*sqrt(pvar1)), code=3, angle=90)






##################
##COMBINED MODEL##
##################

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




 
