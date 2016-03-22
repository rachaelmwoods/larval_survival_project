library(lme4)
library(boot)
library(hier.part)
source("R/functions.R")

# Load datasets
dat <- read.csv("data/data_table_factors - data.csv", as.is=TRUE)

# change mean value fertilisation from percentage to proportion (for plotting)
dat$mean_value_prop <- dat$mean_value/100

#########################
## FERTILISATION MODEL ##
#########################
dat_fert <- dat[dat$life.stage == "fertilisation",]
dat_fert$rep <- 1:nrow(dat_fert)

pdf("figures/figure_S1.pdf", 4.5, 10)

par(mfrow=c(5,2), oma=c(0, 2, 0, 0), mar=c(5, 4, 1, 1))

#sediment
mod_fert <- glm(cbind(success, failure) ~ sediment_mg_per_l, family=binomial, data=dat_fert)
plot(mean_value_prop ~ sediment_mg_per_l, col = "dodgerblue", xlab = "Sediment (mg/L)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$sediment_mg_per_l)
pred_fert <- predict(mod_fert, list(sediment_mg_per_l = ss), type="response", se.fit = TRUE)
lines(ss, pred_fert$fit)

#ammonium
mod_fert <- glm(cbind(success, failure) ~ ammonium_microM, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ ammonium_microM, col = "dodgerblue", xlab = "Ammonium (µM)", ylab = "",  pch=16, data=dat_fert)
ss <- sort(dat_fert$ammonium_microM)
pred_fert <- predict(mod_fert, list(ammonium_microM = ss), type="response", se.fit = TRUE)
lines(ss, pred_fert$fit)

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
pred_fert <- predict(mod_fert, list(cadmium_ug_per_l = ss), type="response", se.fit = TRUE)
lines(ss, pred_fert$fit)
  
#nitrate_microM - REMOVE
mod_fert <- glm(cbind(success, failure) ~ nitrate_microM, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ nitrate_microM, col = "dodgerblue", xlab = "Nitrate (µM)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$nitrate_microM)
lines(ss, predict(mod_fert, list(nitrate_microM = ss), type="response"))

#salinity
dat_fert$salinity_psu_sq <- dat_fert$salinity_psu^2

mod_fert <- glm(cbind(success, failure) ~ salinity_psu + salinity_psu_sq, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ salinity_psu, col = "dodgerblue", xlab = "Salinity (psu)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$salinity_psu)
lines(ss, predict(mod_fert, list(salinity_psu = ss, salinity_psu_sq = ss^2), type="response"))

#acidification - VERY WEAK RESPONSE
dat_fert$acidification_pH_sq <- dat_fert$acidification_pH^2
mod_fert <- glm(cbind(success, failure) ~ acidification_pH + acidification_pH_sq, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ acidification_pH, col = "dodgerblue", xlab = "Acidification (pH)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$acidification_pH)
lines(ss, predict(mod_fert, list(acidification_pH = ss, acidification_pH_sq = ss^2), type="response"))
drop1(mod_fert, test="Chisq")

#tempertaure - REMOVE not enough data
dat_fert$tempertaure_degrees_celcius_sq <- dat_fert$tempertaure_degrees_celcius^2
mod_fert <- glm(cbind(success, failure) ~ tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq, family=binomial, data=dat_fert)
summary(mod_fert)
plot(mean_value_prop ~ tempertaure_degrees_celcius, col = "dodgerblue", xlab = "Temperature (deg C)", ylab = "",  pch=16, las=1, data=dat_fert)
ss <- sort(dat_fert$tempertaure_degrees_celcius)
lines(ss, predict(mod_fert, list(tempertaure_degrees_celcius = ss, tempertaure_degrees_celcius_sq = ss^2), type="response"))

mtext("Proportion fertilised", 2, line=0, outer=TRUE)

dev.off()



##################
##SURVIVAL MODEL##
##################
dat_surv <- dat[dat$life.stage == "survivorship",]
dat_surv$rep <- 1:nrow(dat_surv)

pdf("figures/figure_S2.pdf", 4.5, 10)

par(mfrow=c(4,2), oma=c(0, 2, 0, 0), mar=c(5, 4, 1, 1))

#ammonium - few points
mod_surv <- glm(cbind(success, failure) ~ ammonium_microM, family=binomial, data=dat_surv)
plot(mean_value_prop ~ ammonium_microM,  col = "dodgerblue", xlab = "Ammonium (µM)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$ammonium_microM)
pred_surv <- predict(mod_surv, list(ammonium_microM = ss), type="response", se.fit = TRUE)
lines(ss, pred_surv$fit)

#copper
mod_surv <- glm(cbind(success, failure) ~ copper_ug_per_l, family=binomial, data=dat_surv)
plot(mean_value_prop ~ copper_ug_per_l,  col = "dodgerblue", xlab = "Copper (µg/L)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$copper_ug_per_l)
pred_surv <- predict(mod_surv, list(copper_ug_per_l = ss), type="response", se.fit = TRUE)
lines(ss, pred_surv$fit)

# mercury- REMOVE
mod_surv <- glm(cbind(success, failure) ~ mercury_ug._per_l, family=binomial, data=dat_surv)
plot(mean_value_prop ~ mercury_ug._per_l,  col = "dodgerblue", xlab = "Mercury (µg/L)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$mercury_ug._per_l)
pred_surv <- predict(mod_surv, list(mercury_ug._per_l = ss), type="response", se.fit = TRUE)
lines(ss, pred_surv$fit)

# lead
mod_surv <- glm(cbind(success, failure) ~ lead_ug_per_l, family=binomial, data=dat_surv)
plot(mean_value_prop ~ lead_ug_per_l,  col = "dodgerblue", xlab = "Lead (µg/L)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$lead_ug_per_l)
pred_surv <- predict(mod_surv, list(lead_ug_per_l = ss), type="response", se.fit = TRUE)
lines(ss, pred_surv$fit)

#salinity - REMOVE - 
dat_surv$salinity_psu_sq <- dat_surv$salinity_psu^2
mod_surv <- glm(cbind(success, failure) ~ salinity_psu + salinity_psu_sq, family=binomial, data=dat_surv)
plot(mean_value_prop ~ salinity_psu, col = "dodgerblue", xlab = "Salinity (psu)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$salinity_psu)
pred_surv <- predict(mod_surv, list(salinity_psu = ss, salinity_psu_sq = ss^2), type="response", se.fit = TRUE)
lines(ss, pred_surv$fit)
drop1(mod_surv, test="Chisq")

mod_surv <- glm(cbind(success, failure) ~ salinity_psu, family=binomial, data=dat_surv)


# acidification - VERY WEAK
dat_surv$acidification_pH_sq <- dat_surv$acidification_pH^2
mod_surv <- glm(cbind(success, failure) ~ acidification_pH + acidification_pH_sq, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ acidification_pH,  col = "dodgerblue", xlab = "Acidification (pH)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$acidification_pH)
pred_surv <- predict(mod_surv, list(acidification_pH = ss, acidification_pH_sq= ss^2), type="response", se.fit = TRUE)
lines(ss, pred_surv$fit)
drop1(mod_surv, test="Chisq")

# temperature
dat_surv$tempertaure_degrees_celcius_sq <- dat_surv$tempertaure_degrees_celcius^2
mod_surv <- glm(cbind(success, failure) ~ tempertaure_degrees_celcius + tempertaure_degrees_celcius_sq, family=binomial, data=dat_surv)
summary(mod_surv)
plot(mean_value_prop ~ tempertaure_degrees_celcius,  col = "dodgerblue", xlab = "Temperature (deg C)", ylab = "",  pch=16, las=1, data=dat_surv)
ss <- sort(dat_surv$tempertaure_degrees_celcius)
pred_surv <- predict(mod_surv, list(tempertaure_degrees_celcius = ss, tempertaure_degrees_celcius_sq = ss^2), type="response", se.fit = TRUE)
lines(ss, pred_surv$fit)


mtext("Proportion survived", 2, line=0, outer=TRUE)

dev.off()


############################
## FULL MODEL FERTILISATION
##############################

# Load datasets
dat <- read.csv("data/data_table_factors - data.csv", as.is=TRUE)

dat$success <- round(dat$success)
dat$failure <- round(dat$failure)

# change mean value fertilisation from percentage to proportion (for plotting)
dat$mean_value_prop <- dat$mean_value/100
dat_fert2 <- dat_fert <- dat[dat$life.stage == "fertilisation",]

dat_fert <- dat_fert[c("success", "failure", "sediment_mg_per_l", "ammonium_microM", "phosphorous_microM", "copper_ug_per_l", "salinity_psu", "experiment", "mean_value_prop", "spawn.brood")]

dat_fert <- dat_fert[!(is.na(dat_fert$sediment_mg_per_l) & is.na(dat_fert$ammonium_microM) & is.na(dat_fert$phosphorous_microM) & is.na(dat_fert$copper_ug_per_l) & is.na(dat_fert$salinity_psu)),]

dat_fert$sediment_mg_per_l[is.na(dat_fert$sediment_mg_per_l) | dat_fert$sediment_mg_per_l == 0] <- 0.001 * max(dat_fert2$sediment_mg_per_l, na.rm=TRUE)
dat_fert$ammonium_microM[is.na(dat_fert$ammonium_microM) | dat_fert$ammonium_microM == 0] <- 0.001 * max(dat_fert2$ammonium_microM, na.rm=TRUE)
dat_fert$phosphorous_microM[is.na(dat_fert$phosphorous_microM) | dat_fert$phosphorous_microM == 0] <- 0.001 * max(dat_fert2$phosphorous_microM, na.rm=TRUE)
dat_fert$copper_ug_per_l[is.na(dat_fert$copper_ug_per_l) | dat_fert$copper_ug_per_l == 0] <- 0.01 * max(dat_fert2$copper_ug_per_l, na.rm=TRUE)

dat_fert$salinity_psu[is.na(dat_fert$salinity_psu)] <- 34

dat_fert$rep <- 1:nrow(dat_fert)

# Rescale
dat_fert$rs_sediment_mg_per_l <- rescale(dat_fert$sediment_mg_per_l)
dat_fert$rs_ammonium_microM <- rescale(dat_fert$ammonium_microM)
dat_fert$rs_phosphorous_microM <- rescale(dat_fert$phosphorous_microM)
dat_fert$rs_copper_ug_per_l <- rescale(dat_fert$copper_ug_per_l)
dat_fert$rs_salinity_psu <- rescale(dat_fert$salinity_psu)
dat_fert$rs_salinity_psu_sq <- dat_fert$rs_salinity_psu^2

##MODEL##
mod_fert_full <- glm(cbind(success, failure) ~ rs_sediment_mg_per_l + rs_ammonium_microM + rs_phosphorous_microM + rs_copper_ug_per_l + rs_salinity_psu + rs_salinity_psu_sq, family=binomial, data=dat_fert)
summary(mod_fert_full)
drop1(mod_fert_full, test="Chisq")

dat_fert$experiment <- factor(dat_fert$experiment)
dat_fert$rep <- factor(dat_fert$rep)

mod_fert_full <- glmer(cbind(success, failure) ~ rs_sediment_mg_per_l + rs_ammonium_microM + rs_phosphorous_microM + rs_copper_ug_per_l + rs_salinity_psu + rs_salinity_psu_sq + (1 | experiment) + (1 | rep), family=binomial, data=dat_fert, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000000)))
summary(mod_fert_full)

sum(residuals(mod_fert_full, type="pearson")^2)/df.residual(mod_fert_full)
drop1(mod_fert_full, test="Chisq")

##PLOTS FOR FERTILISATION##

# TESTS, should all equal 0, accept salinity
dat_fert$rs_sediment_mg_per_l[dat_fert$sediment_mg_per_l==0.01 * max(dat_fert2$sediment_mg_per_l, na.rm=TRUE)][1]
dat_fert$rs_ammonium_microM[dat_fert$ammonium_microM==0.001 * max(dat_fert2$ammonium_microM, na.rm=TRUE)][1]
dat_fert$rs_phosphorous_microM[dat_fert$phosphorous_microM==0.001 * max(dat_fert2$phosphorous_microM, na.rm=TRUE)][1]
dat_fert$rs_copper_ug_per_l[dat_fert$copper_ug_per_l==0.01 * max(dat_fert2$copper_ug_per_l, na.rm=TRUE)][1]
dat_fert$rs_salinity_psu[dat_fert$salinity_psu==34][1]


pdf("figures/figure_1.pdf", 5.5, 8)

par(mfrow=c(3,2), oma=c(0,2,0,0), mar=c(4, 4, 2, 1))

# COPPER
# ss <- seq(min(dat_fert2$copper_ug_per_l), max(dat_fert2$copper_ug_per_l), 0.05)
ss <- seq(0, 1, 0.01)
newdat <- expand.grid(rs_sediment_mg_per_l=0, rs_ammonium_microM=0, rs_phosphorous_microM=0, rs_copper_ug_per_l = ss, rs_salinity_psu = 0.8478261, rs_salinity_psu_sq = 0.8478261^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_full), newdat)
newdat$success <- mm %*% fixef(mod_fert_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_full),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1]  

plot(dat_fert$rs_copper_ug_per_l, dat_fert$mean_value_prop, xlab="Copper (µg/L)", ylab="", ylim=c(0, 1), axes=FALSE)
axis(1)
axis(2, las=2)
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
polygon(c(ss, rev(ss)), c(inv.logit(newdat$success+2*sqrt(pvar1)), rev(inv.logit(newdat$success-2*sqrt(pvar1)))), col=rgb(0,0,0,0.2), border=NA)
mtext("A", side=3, line=0, adj=0, cex=1.2)

# SEDIMENT
# ss <- seq(min(dat_fert2$sediment_mg_per_l), max(dat_fert2$sediment_mg_per_l), 0.01)
newdat <- expand.grid(rs_sediment_mg_per_l=ss, rs_ammonium_microM=0, rs_phosphorous_microM=0, rs_copper_ug_per_l = 0, rs_salinity_psu = 0.8478261, rs_salinity_psu_sq = 0.8478261^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_full),newdat)
newdat$success <- mm %*% fixef(mod_fert_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_full),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1]  

plot(dat_fert$rs_sediment_mg_per_l, dat_fert$mean_value_prop, xlab="Sediment (mg/L)", ylab="", ylim=c(0, 1), axes=FALSE)
axis(1)
axis(2, las=2)
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
polygon(c(ss, rev(ss)), c(inv.logit(newdat$success+2*sqrt(pvar1)), rev(inv.logit(newdat$success-2*sqrt(pvar1)))), col=rgb(0,0,0,0.2), border=NA)
mtext("B", side=3, line=0, adj=0, cex=1.2)

# AMMONIUM
# ss <- seq(min(dat_fert$ammonium_microM), max(dat_fert$ammonium_microM), 0.05)
newdat <- expand.grid(rs_sediment_mg_per_l=0, rs_ammonium_microM=ss, rs_phosphorous_microM=0, rs_copper_ug_per_l = 0, rs_salinity_psu = 0.8478261, rs_salinity_psu_sq = 0.8478261^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_full),newdat)
newdat$success <- mm %*% fixef(mod_fert_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_full),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1]  

plot(dat_fert$rs_ammonium_microM, dat_fert$mean_value_prop, xlab="Ammonium (µM)", ylab="", ylim=c(0, 1), axes=FALSE)
axis(1)
axis(2, las=2)
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
polygon(c(ss, rev(ss)), c(inv.logit(newdat$success+2*sqrt(pvar1)), rev(inv.logit(newdat$success-2*sqrt(pvar1)))), col=rgb(0,0,0,0.2), border=NA)
mtext("C", side=3, line=0, adj=0, cex=1.2)

# PHOSPHOROUS
# ss <- seq(min(dat_fert$phosphorous_microM), max(dat_fert$phosphorous_microM), 0.05)
newdat <- expand.grid(rs_sediment_mg_per_l=0, rs_ammonium_microM=0, rs_phosphorous_microM=ss, rs_copper_ug_per_l = 0, rs_salinity_psu = 0.8478261, rs_salinity_psu_sq = 0.8478261^2, success=0, failure=0)
mm <- model.matrix(terms(mod_fert_full),newdat)
newdat$success <- mm %*% fixef(mod_fert_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_full),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1] 

plot(dat_fert$rs_phosphorous_microM, dat_fert$mean_value_prop, xlab="Phosphorous (µM)", ylab="", ylim=c(0, 1), axes=FALSE)
axis(1)
axis(2, las=2)
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
polygon(c(ss, rev(ss)), c(inv.logit(newdat$success+2*sqrt(pvar1)), rev(inv.logit(newdat$success-2*sqrt(pvar1)))), col=rgb(0,0,0,0.2), border=NA)
mtext("D", side=3, line=0, adj=0, cex=1.2)

# SALINITY
# ss <- seq(min(dat_fert$salinity_psu), max(dat_fert$salinity_psu), 0.1)
newdat <- expand.grid(rs_sediment_mg_per_l=0, rs_ammonium_microM=0, rs_phosphorous_microM=0, rs_copper_ug_per_l = 0, rs_salinity_psu = ss, success=0, failure=0)
newdat$rs_salinity_psu_sq = newdat$rs_salinity_psu^2
mm <- model.matrix(terms(mod_fert_full),newdat)
newdat$success <- mm %*% fixef(mod_fert_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_fert_full),mm))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1] 

plot(dat_fert$rs_salinity_psu, dat_fert$mean_value_prop, xlab="Salinity (psu)", ylab="", ylim=c(0, 1), axes=FALSE)
axis(1)
axis(2, las=2)
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
polygon(c(ss, rev(ss)), c(inv.logit(newdat$success+2*sqrt(pvar1)), rev(inv.logit(newdat$success-2*sqrt(pvar1)))), col=rgb(0,0,0,0.2), border=NA)
mtext("E", side=3, line=0, adj=0, cex=1.2)


mtext("Proportion fertilised", 2, line=0, outer=TRUE)

dev.off()

#######################
##FULL MODEL SURVIVAL##   
#######################

dat <- read.csv("data/data_table_factors - data.csv", as.is=TRUE)

dat$success <- round(dat$success)
dat$failure <- round(dat$failure)

# change mean value fertilisation from percentage to proportion (for plotting)
dat$mean_value_prop <- dat$mean_value/100
dat_surv2 <- dat_surv <- dat[dat$life.stage == "survivorship",]

dat_surv <- dat_surv[c("success", "failure", "copper_ug_per_l", "lead_ug_per_l", "salinity_psu", "tempertaure_degrees_celcius", "acidification_pH", "experiment", "mean_value_prop", "spawn.brood")]

dat_surv <- dat_surv[!(is.na(dat_surv$copper_ug_per_l) & is.na(dat_surv$lead_ug_per_l) & is.na(dat_surv$salinity_psu) & is.na(dat_surv$acidification_pH) & is.na(dat_surv$tempertaure_degrees_celcius)),]

dat_surv$copper_ug_per_l[is.na(dat_surv$copper_ug_per_l) | dat_surv$copper_ug_per_l == 0] <- 0.01 * max(dat_surv2$copper_ug_per_l, na.rm=TRUE)
dat_surv$lead_ug_per_l[is.na(dat_surv$lead_ug_per_l) | dat_surv$lead_ug_per_l == 0] <- 0.001 * max(dat_surv2$lead_ug_per_l, na.rm=TRUE)

dat_surv$salinity_psu[is.na(dat_surv$salinity_psu)] <- 34
dat_surv$acidification_pH[is.na(dat_surv$acidification_pH)] <- 8.1
dat_surv$tempertaure_degrees_celcius[is.na(dat_surv$tempertaure_degrees_celcius)] <- 28

dat_surv$rep <- 1:nrow(dat_surv)


# Rescale
dat_surv$rs_copper_ug_per_l <- rescale(dat_surv$copper_ug_per_l)
dat_surv$rs_lead_ug_per_l <- rescale(dat_surv$lead_ug_per_l)
dat_surv$rs_salinity_psu <- rescale(dat_surv$salinity_psu)
dat_surv$rs_salinity_psu_sq <- dat_surv$rs_salinity_psu^2
dat_surv$rs_acidification_pH <- rescale(dat_surv$acidification_pH)
# dat_surv$rs_acidification_pH_sq <- dat_surv$rs_acidification_pH^2
dat_surv$rs_tempertaure_degrees_celcius <- rescale(dat_surv$tempertaure_degrees_celcius)
dat_surv$rs_tempertaure_degrees_celcius_sq <- dat_surv$rs_tempertaure_degrees_celcius^2


dat_surv$experiment <- factor(dat_surv$experiment)
dat_surv$rep <- factor(dat_surv$rep)

## MODEL##
mod_surv_full <- glm(cbind(success, failure) ~ rs_copper_ug_per_l + rs_lead_ug_per_l + rs_salinity_psu + rs_salinity_psu_sq + rs_acidification_pH + rs_tempertaure_degrees_celcius + rs_tempertaure_degrees_celcius_sq, family=binomial, data=dat_surv)
summary(mod_surv_full)
drop1(mod_surv_full, test="Chisq")

mod_surv_full <- glmer(cbind(success, failure) ~ rs_copper_ug_per_l + rs_lead_ug_per_l + rs_salinity_psu + rs_salinity_psu_sq + rs_acidification_pH + rs_tempertaure_degrees_celcius + rs_tempertaure_degrees_celcius_sq + (1 | experiment) + (1 | rep), family=binomial, data=dat_surv, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e4)))

summary(mod_surv_full)

sum(residuals(mod_surv_full, type="pearson")^2)/df.residual(mod_surv_full)

drop1(mod_surv_full, test="Chisq")

mod_surv_full <- glmer(cbind(success, failure) ~ rs_copper_ug_per_l + rs_lead_ug_per_l + rs_salinity_psu + rs_salinity_psu_sq + (1 | experiment) + (1 | rep), family=binomial, data=dat_surv, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e4)))

summary(mod_surv_full)
drop1(mod_surv_full, test="Chisq")

# TESTS, should all equal 0, accept salinity
dat_surv$rs_copper_ug_per_l[dat_surv$copper_ug_per_l==0.01 * max(dat_surv2$copper_ug_per_l, na.rm=TRUE)][1]
dat_surv$rs_lead_ug_per_l[dat_surv$lead_ug_per_l==0.001 * max(dat_surv2$lead_ug_per_l, na.rm=TRUE)][1]
dat_surv$rs_salinity_psu[dat_surv$salinity_psu==34][1]

##PLOTS FOR SURVIVAL##


pdf("figures/figure_2.pdf", 5.5, 8)

par(mfrow=c(3,1), oma=c(0,2,0,0), mar=c(5, 4, 2, 7))

# COPPER

ss <- seq(0, 1, 0.01)
# ss <- seq(min(dat_surv2$copper_ug_per_l), max(dat_surv2$copper_ug_per_l), 0.05)
newdat <- expand.grid(rs_copper_ug_per_l = ss, rs_lead_ug_per_l=0, rs_salinity_psu = 0.8478261, success=0, failure=0)
mm <- model.matrix(terms(mod_surv_full), newdat)
newdat$success <- mm %*% fixef(mod_surv_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_full),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1]  

plot(dat_surv$rs_copper_ug_per_l, dat_surv$mean_value_prop, xlab="Copper (µg/L)", ylab="", ylim=c(0, 1),  axes=FALSE)
axis(1)
axis(2)
lines(ss, inv.logit(newdat$success), lwd=1, lty=1)
polygon(c(ss, rev(ss)), c(inv.logit(newdat$success+2*sqrt(pvar1)), rev(inv.logit(newdat$success-2*sqrt(pvar1)))), col=rgb(0,0,0,0.2), border=NA)
mtext("A", side=3, line=0, adj=0, cex=1.2)

# LEAD
# ss <- seq(min(dat_surv2$lead_ug_per_l), max(dat_surv2$lead_ug_per_l), 0.05)
newdat <- expand.grid(rs_copper_ug_per_l=0, rs_lead_ug_per_l = ss, rs_salinity_psu = 0.8478261, success=0, failure=0)
mm <- model.matrix(terms(mod_surv_full), newdat)
newdat$success <- mm %*% fixef(mod_surv_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_full),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1]  

plot(dat_surv$rs_copper_ug_per_l, dat_surv$mean_value_prop, xlab="Lead (µg/L)", ylab="", ylim=c(0, 1), axes=FALSE)
axis(1)
axis(2)
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) # inf
polygon(c(ss, rev(ss)), c(inv.logit(newdat$success+2*sqrt(pvar1)), rev(inv.logit(newdat$success-2*sqrt(pvar1)))), col=rgb(0,0,0,0.2), border=NA)
mtext("B", side=3, line=0, adj=0, cex=1.2)

#SALINITY
# ss <- seq(min(dat_surv2$salinity_psu), max(dat_surv2$salinity_psu), 0.2)
newdat <- expand.grid(rs_copper_ug_per_l=0, rs_lead_ug_per_l = 0, rs_salinity_psu=ss, success=0, failure=0)
newdat$rs_salinity_psu_sq <- newdat$rs_salinity_psu^2
mm <- model.matrix(terms(mod_surv_full),newdat)
newdat$success <- mm %*% fixef(mod_surv_full)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod_surv_full),mm))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1] 

plot(dat_surv$rs_salinity_psu, dat_surv$mean_value_prop, xlab="Salinity (psu)", ylab="", ylim=c(0, 1), axes=FALSE)
axis(1)
axis(2)
lines(ss, inv.logit(newdat$success), lwd=1, lty=1) 
polygon(c(ss, rev(ss)), c(inv.logit(newdat$success+2*sqrt(pvar1)), rev(inv.logit(newdat$success-2*sqrt(pvar1)))), col=rgb(0,0,0,0.2), border=NA)
mtext("C", side=3, line=0, adj=0, cex=1.2)

mtext("Proportion survived", 2, line=0, outer=TRUE)

dev.off()


###################
#VARIANCE ANALYSIS#
###################

pdf("figures/figure_3.pdf", 5.5, 8)

par(mfrow=c(2,1), oma=c(0,2,0,0), mar=c(4, 2, 2, 1))

#dat_fert = dat[apply(!is.na(dat[, -1]), 1, sum) > 0 & dat$life.stage == "fertilisation",]
#dat_surv = dat[apply(!is.na(dat[, -1]), 1, sum) > 0 & dat$life.stage == "survivorship",]

factors <- dat_fert[, c(17,12,14,15,21)]
colnames(factors) <- c("Copper", "Sediment","Ammonium", "Phosphate", "Salinity") 
hier.part(dat_fert$mean_value_prop, factors, family = "binomial", gof = "logLik", barplot = TRUE)

factors <- dat_surv[, c(17,27,21)]
colnames(factors) <- c("Copper", "Lead", "Salinity") 
hier.part(dat_surv$mean_value_prop, factors, family = "binomial", gof = "logLik", barplot = TRUE)

dev.off()

###########################
##LOCATIONS/WATER SAMPLES##
############################

water <- read.csv("data/water_samples.csv", as.is=TRUE)

# Water data preperation, deadling with non-numeric characters like "<"
water[water == "<1"] <- "0.5"
water[water == "<0.05"] <- "0.025"
water[water == "<0.1"] <- "0.05"
water[water == "<0.25"] <- "0.125"
water[water == "<0.5"] <- "0.25"
water[water == "<0.005"] <- "0.0025"
water[water == "<5"] <- "2.5"

water[,2:14] <- apply(water[,2:14], 2, as.numeric)
water <- water[1:3,]



##Fertilisation Model##
pdf("figures/figure_4.pdf", 5.5, 8)
par(mfrow=c(2,1), oma=c(0,2,0,0), mar=c(5, 5, 3, 1))

mod_sediment_mg_per_l <- lm(rs_sediment_mg_per_l ~ sediment_mg_per_l, dat_fert)
mod_copper_ug_per_l <- lm(rs_copper_ug_per_l ~ copper_ug_per_l, dat_fert)
mod_ammonium_microM <- lm(rs_ammonium_microM ~ ammonium_microM, dat_fert)
mod_phosphorous_microM <- lm(rs_phosphorous_microM ~ phosphorous_microM, dat_fert)
mod_salinity_psu <- lm(rs_salinity_psu ~ salinity_psu, dat_fert)


water_fert <- data.frame(
  rs_sediment_mg_per_l=predict(mod_sediment_mg_per_l, list(sediment_mg_per_l=water$suspended_solids_mg.l)), 
  rs_copper_ug_per_l=predict(mod_copper_ug_per_l, list(copper_ug_per_l=water$copper_ug.l)),
  rs_ammonium_microM=predict(mod_ammonium_microM, list(ammonium_microM=water$ammonia_mg.l)),
  rs_phosphorous_microM=predict(mod_phosphorous_microM, list(phosphorous_microM=water$phosphorus_mg.l)), 
  rs_salinity_psu=predict(mod_salinity_psu, list(salinity_psu=water$salinity_g.l)), 
  success=0, failure=0
)

water_fert$rs_salinity_psu_sq <- water_fert$rs_salinity_psu^2

mm_fert <- model.matrix(terms(mod_fert_full), water_fert)

water_fert$success <- mm_fert %*% fixef(mod_fert_full)
pvar1 <- diag(mm_fert %*% tcrossprod(vcov(mod_fert_full), mm_fert))
tvar1 <- pvar1 + VarCorr(mod_fert_full)$experiment[1] + VarCorr(mod_fert_full)$rep[1]  
#bar plot
bp <- barplot(t(inv.logit(water_fert$success)), xlab="Location", ylab="Proportion of Larvae Fertilised", ylim=c(0, 1), names.arg=water$sample)
#error bars
arrows(bp, inv.logit(water_fert$success-2*sqrt(pvar1)), bp, inv.logit(water_fert$success+2*sqrt(pvar1)), code=3, angle=90)
mtext("E", side=3, line=0, adj=0, cex=1.2)

##Survival Model##

mod_copper_ug_per_l <- lm(rs_copper_ug_per_l ~ copper_ug_per_l, dat_surv)
mod_lead_ug_per_l <- lm(rs_lead_ug_per_l ~ lead_ug_per_l, dat_surv)

water_surv <- data.frame(
  rs_copper_ug_per_l=predict(mod_copper_ug_per_l, list(copper_ug_per_l=water$copper_ug.l)),
  rs_lead_ug_per_l=predict(mod_lead_ug_per_l, list(lead_ug_per_l=water$lead_ug.l)),
  rs_salinity_psu=predict(mod_salinity_psu, list(salinity_psu=water$salinity_g.l)), 
  success=0, failure=0
)
water_surv$rs_salinity_psu_sq <- water_surv$rs_salinity_psu^2

mm_surv <- model.matrix(terms(mod_surv_full), water_surv)

water_surv$success <- mm_surv %*% fixef(mod_surv_full)
pvar1 <- diag(mm_surv %*% tcrossprod(vcov(mod_surv_full), mm_surv))
tvar1 <- pvar1 + VarCorr(mod_surv_full)$experiment[1] + VarCorr(mod_surv_full)$rep[1]  
#bar plot
bp <- barplot(t(inv.logit(water_surv$success)), xlab="Location", ylab="Proportion of Larvae Survived", ylim=c(0, 1), names.arg=water$sample)
#error bars
arrows(bp, inv.logit(water_surv$success-2*sqrt(pvar1)), bp, inv.logit(water_surv$success+2*sqrt(pvar1)), code=3, angle=90)

dev.off()

##################
##COMBINED MODEL##
##################

pdf("figures/figure_5.pdf", 5.5, 8)
par(mfrow=c(2,1), oma=c(0,2,0,0), mar=c(5, 5, 3, 1))

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

dev.off()




