###############################################################################
# Imputation of the mRS at Baseline for all treated patients:
###############################################################################
set.seed(123)
variables_to_imp <- c("gender", "age_Consultation", "Positive_familial_history",
                      "prf_hypertension_awareness", "prf_smokedtobacco", 
                      "aneurisk", "maxdiam", "type", "side",
                      "recommended", "perprotocol", "modality",
                      "mR_Baseline", "mR_1", "residence",
                      "pseudonym","mR_Discharge")

dat_imp <- dat_s[, c( variables_to_imp)]
dat_imp$modality <- as.factor(dat_imp$modality)
dat_imp$mR_Baseline <- as.numeric(as.character(dat_imp$mR_Baseline))
table(dat_imp$mR_Baseline, useNA = "always")

# mR_Baseline for all treated patients is unknown, set to NA
dat_imp$mR_Baseline[dat_imp$modality !="Observation"] <- NA
dat_imp$mR_Baseline <-  as.ordered(dat_imp$mR_Baseline) #model it as ordered factor
table(dat_imp$mR_Baseline, useNA = "always")

# Covariate-Shift model via Propensity Score:
ps_model <- glm(modality == "Observation" ~ gender + scale(age_Consultation) +
                  Positive_familial_history + prf_hypertension_awareness +
                  prf_smokedtobacco + aneurisk + scale(maxdiam) + type + side +
                  perprotocol + residence+mR_1,
                data = dat_imp, family = binomial)
# summary(ps_model)

dat_imp$pscore <- predict(ps_model, type = "response")
dat_imp$p_observed <- mean(dat_imp$modality == "Observation")
dat_imp$ipw_stab <- ifelse(dat_imp$modality == "Observation", 
                           dat_imp$p_observed / dat_imp$pscore, 
                           (1 - dat_imp$p_observed) / (1 - dat_imp$pscore))
dat_imp$ipw_stab <- pmin(dat_imp$ipw_stab, 5)  # Trimmen


# Set the Imputation Methods
meth <- make.method(dat_imp)
meth["mR_Baseline"] <- "pmm"  #predictive mean matching
meth["ipw_stab"] <- ""
meth["pseudonym"] <- ""  
meth["mR_Discharge"] <- ""  

predMat <- make.predictorMatrix(dat_imp)
predMat["mR_Baseline", ] <- 1
predMat["mR_Baseline", "mR_Baseline"] <- 0  # can't predict itself
predMat["mR_Baseline", "ipw_stab"] <- 0
predMat["pseudonym", ] <- 0  # pseudonym is not imputed
predMat[, "pseudonym"] <- 0  # pseudonym is not used as predictor as it is an identifier
predMat["mR_Discharge", ] <- 0  # mR_D is not impute
predMat[, "mR_Discharge"] <- 0  # mRS at Discharge is not used as predictor, 
# as it is too close to mRS at Baseline

# weighted imputation with MICE:
imp <- mice(dat_imp, m=10, method=meth, predictorMatrix=predMat, 
            weights=dat_imp$ipw_stab, seed=123, printFlag=TRUE)
names(imp$data)[names(imp$data) == "mR_Baseline"] <- "mR_BaselineMICE"


imp_long <- complete(imp, "long")  
# table(imp_long$pseudonym) #santiy check, should be equal to m, for all
# sum((table(imp_long$pseudonym))!=10) #sanity check, should be 0

## compute the mean resulting mRS at Baseline per patient across all imputations:
imp_res <- imp_long %>%
  select(.imp, pseudonym, modality, mR_BaselineMICE)

imp_mean <- imp_res %>%
  group_by(pseudonym, modality) %>%
  summarise(
    mR_Baseline_mean = round(mean(as.numeric(as.character(mR_BaselineMICE)), na.rm = TRUE)),
    .groups = "drop")

dat_s$mR_BaselinePredMICEmean <- imp_mean$mR_Baseline_mean[match(dat_s$pseudonym, imp_mean$pseudonym)]
dat_s$mR_BaselinePredMICEmean <- as.factor(dat_s$mR_BaselinePredMICEmean)

table(imp_mean$modality, imp_mean$mR_Baseline_mean, useNA = "always")

rm(meth, predMat, variables_to_imp, ps_model, imp_res)

###################################################################
# Generate a long-data set
###################################################################
## make a long data format:
datlong <- tidyr::pivot_longer(dat_s,
                               cols = c(mR_BaselinePredMICEmean,mR_Discharge, mR_1),
                               names_to = "timepoint",
                               values_to = "mRS") %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("mR_BaselinePredMICEmean", "mR_Discharge", "mR_1"),
                            labels=c("Baseline", "Discharge", "1-year FU"))) %>% 
  mutate(age_tp = case_when(
    timepoint == "Baseline" & modality == "Observation" ~ age_e - 5/52,
    timepoint == "Baseline" & modality != "Observation" ~ age_e - 3/12 - 5/52,
    timepoint == "Discharge" & modality == "Observation"~ age_e,
    timepoint == "Discharge" & modality != "Observation"~ age_e+6/365,
    timepoint == "1-year FU" & modality == "Observation"~ age_e +1,
    timepoint == "1-year FU" & modality != "Observation"~ age_e +1+6/365,
    TRUE ~ NA_real_
  )) 
