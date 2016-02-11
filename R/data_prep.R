# Data preparation

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

# Look at data

dat$tempertaure_degrees_kelvin <- log10(dat$tempertaure_degrees_celcius + 273.15)

#remove spawn/brood data
dat$spawnbrood[dat$spawn.brood==""] <- NA

# #filling all NAs in factors (because not measured in expt) with real data
dat$sediment_mg_per_l[dat$sediment_mg_per_l == 0] <- 1
dat$ammonium_microM[dat$ammonium_microM == 0] <- 0.01391
dat$phosphorous_microM[dat$phosphorous_microM == 0] <- 0.446
dat$copper_ug_per_l[dat$copper_ug_per_l == 0] <- 0.9
dat$mercury_ug._per_l[dat$mercury_ug._per_l == 0] <- 0.15
dat$tributyltin_ug_per_l[dat$tributyltin_ug_per_l == 0] <- 0.01
dat$zinc_ug_per_l[dat$zinc_ug_per_l == 0] <- 5
dat$cadmium_ug_per_l[dat$cadmium_ug_per_l == 0] <- 0.11
dat$nitrate_microM[dat$nitrate_microM == 0] <- 0.254
dat$lead_ug_per_l[dat$lead_ug_per_l == 0] <- 0.03
dat$nickel_ug_per_l[dat$nickel_ug_per_l == 0] <- 6.6


# Log10 continuous variables based on histograms above. Salinity, Temperature and pH not log-transformed
dat$sediment_mg_per_l <- log10(dat$sediment_mg_per_l)
dat$ammonium_microM <- log10(dat$ammonium_microM)
dat$phosphorous_microM <- log10(dat$phosphorous_microM)
dat$copper_ug_per_l <- log10(dat$copper_ug_per_l)
dat$mercury_ug._per_l <- log10(dat$mercury_ug._per_l)
dat$tributyltin_ug_per_l <- log10(dat$tributyltin_ug_per_l)
dat$zinc_ug_per_l <- log10(dat$zinc_ug_per_l)
dat$cadmium_ug_per_l <- log10(dat$cadmium_ug_per_l)
dat$nitrate_microM <- log10(dat$nitrate_microM)
dat$lead_ug_per_l <- log10(dat$lead_ug_per_l)
dat$nickel_ug_per_l <- log10(dat$nickel_ug_per_l)

#where there is no salinity we're making it 35
#dat$salinity_psu[is.na(dat$salinity_psu)] <- 35
#square salinity because it has a quadratic response - because it decreases on both sides from 35psu
dat$salinity_psu_sq <- dat$salinity_psu^2

#where there is no acidification_pH we're making it 8.1
#dat$acidification_pH[is.na(dat$acidification_pH)] <- 8.1
dat$acidification_pH <- dat$acidification_pH
dat$acidification_pH_sq <- dat$acidification_pH^2

#where there is no tempertaure_degrees_celcius we're making it 28
#dat$tempertaure_degrees_kelvin[is.na(dat$tempertaure_degrees_kelvin)] <- 28
dat$tempertaure_degrees_celcius_sq <- dat$tempertaure_degrees_celcius^2
dat$tempertaure_degrees_kelvin_sq <- dat$tempertaure_degrees_kelvin^2

# change mean value fertilisation from percentage to proportion (for plotting)
dat$mean_value_prop <- dat$mean_value/100

dat$success <- round(dat$success)
dat$failure <- round(dat$failure)


