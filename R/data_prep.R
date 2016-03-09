# # Data preparation


# # Look at data
# dat$tempertaure_degrees_kelvin <- log10(dat$tempertaure_degrees_celcius + 273.15)

# #remove spawn/brood data
# dat$spawnbrood[dat$spawn.brood==""] <- NA

# #filling all NAs in factors (because not measured in expt) with real data
# dat$sediment_mg_per_l[is.na(dat$sediment_mg_per_l) | dat$sediment_mg_per_l == 0] <- 0.01 * max(dat$sediment_mg_per_l, na.rm=TRUE)
# dat$ammonium_microM[is.na(dat$ammonium_microM) | dat$ammonium_microM == 0] <- 0.01 * max(dat$ammonium_microM, na.rm=TRUE)
# dat$phosphorous_microM[is.na(dat$phosphorous_microM) | dat$phosphorous_microM == 0] <- 0.01 * max(dat$phosphorous_microM, na.rm=TRUE)
# dat$copper_ug_per_l[is.na(dat$copper_ug_per_l) | dat$copper_ug_per_l == 0] <- 0.01 * max(dat$copper_ug_per_l, na.rm=TRUE)
# dat$mercury_ug._per_l[is.na(dat$mercury_ug._per_l) | dat$mercury_ug._per_l == 0] <- 0.01 * max(dat$mercury_ug._per_l, na.rm=TRUE)
# dat$tributyltin_ug_per_l[is.na(dat$tributyltin_ug_per_l) | dat$tributyltin_ug_per_l == 0] <- 0.01 * max(dat$tributyltin_ug_per_l, na.rm=TRUE)
# dat$zinc_ug_per_l[is.na(dat$zinc_ug_per_l) | dat$zinc_ug_per_l == 0] <- 0.01 * max(dat$zinc_ug_per_l, na.rm=TRUE)
# dat$cadmium_ug_per_l[is.na(dat$cadmium_ug_per_l) | dat$cadmium_ug_per_l == 0] <- 0.01 * max(dat$cadmium_ug_per_l, na.rm=TRUE)
# dat$nitrate_microM[is.na(dat$nitrate_microM) | dat$nitrate_microM == 0] <- 0.01 * max(dat$nitrate_microM, na.rm=TRUE)
# dat$lead_ug_per_l[is.na(dat$lead_ug_per_l) | dat$lead_ug_per_l == 0] <- 0.001 * max(dat$lead_ug_per_l, na.rm=TRUE)
# dat$nickel_ug_per_l[is.na(dat$nickel_ug_per_l) | dat$nickel_ug_per_l == 0] <- 0.01 * max(dat$nickel_ug_per_l, na.rm=TRUE)

# # Log10 + 1 continuous variables based on histograms above. Salinity, Temperature and pH not log-transformed
# dat$sediment_mg_per_l <- log10(dat$sediment_mg_per_l)
# dat$ammonium_microM <- log10(dat$ammonium_microM)
# dat$phosphorous_microM <- log10(dat$phosphorous_microM)
# dat$copper_ug_per_l <- log10(dat$copper_ug_per_l)
# dat$mercury_ug._per_l <- log10(dat$mercury_ug._per_l)
# dat$tributyltin_ug_per_l <- log10(dat$tributyltin_ug_per_l)
# dat$zinc_ug_per_l <- log10(dat$zinc_ug_per_l)
# dat$cadmium_ug_per_l <- log10(dat$cadmium_ug_per_l)
# dat$nitrate_microM <- log10(dat$nitrate_microM)
# dat$lead_ug_per_l <- log10(dat$lead_ug_per_l)
# dat$nickel_ug_per_l <- log10(dat$nickel_ug_per_l)

# #where there is no salinity we're making it 34, from the raw data fit
# dat$salinity_psu[is.na(dat$salinity_psu)] <- 34
# #square salinity because it has a quadratic response - because it decreases on both sides from 35psu
# dat$salinity_psu_sq <- dat$salinity_psu^2

# dat$log_salinity_psu <- log10(dat$salinity_psu)
# dat$log_salinity_psu_sq <- log10(dat$salinity_psu_sq)


# #where there is no acidification_pH we're making it 8.1
# dat$acidification_pH[is.na(dat$acidification_pH)] <- 8.1
# dat$acidification_pH <- dat$acidification_pH
# dat$acidification_pH_sq <- dat$acidification_pH^2

# #where there is no tempertaure_degrees_celcius we're making it 28
# dat$tempertaure_degrees_celcius[is.na(dat$tempertaure_degrees_celcius)] <- 28
# dat$tempertaure_degrees_celcius_sq <- dat$tempertaure_degrees_celcius^2
# dat$tempertaure_degrees_kelvin_sq <- dat$tempertaure_degrees_kelvin^2


# dat$success <- round(dat$success)
# dat$failure <- round(dat$failure)


