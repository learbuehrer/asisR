##############################################################################
# Analysis of the Management Effect:
# 1. Causal Regression
# 2. abn and MCMCabn
##############################################################################
# Causal Regression Models for the MICE data
##############################################################################
# 0. preparation of the data set:
##############################################################################
dat_i_list <- vector("list", length = imp$m)
dat_long_list <- vector("list", length = imp$m)

for (i in 1:imp$m) {
  dat_i <- complete(imp, i)
  
  ## Additional code chunk for sensitivity analysis:
  # start of the code chunk:
  ## sensitivity analysis
  # table(dat_i$mR_Baseline, dat_i$modality)
  # dat_i$mR_Baseline <- as.numeric(as.character(dat_i$mR_Baseline))
  # dat_i$mR_Baseline <- ifelse(dat_i$modality == "Observation", dat_i$mR_Baseline, dat_i$mR_Baseline-1)
  # dat_i$mR_Baseline <- ifelse(dat_i$mR_Baseline < 0, 0, dat_i$mR_Baseline)
  # dat_i$mR_Baseline <- ifelse(dat_i$mR_Baseline > 6, 6, dat_i$mR_Baseline)
  # table(dat_i$mR_Baseline, dat_i$modality)
  # dat_i$mR_Baseline <- factor(dat_i$mR_Baseline, ordered=T)
  
  fitmodality <- multinom(modality ~
                            recommended +
                            gender +
                            scale(age_Consultation) +
                            Positive_familial_history +
                            prf_hypertension_awareness +
                            prf_smokedtobacco +
                            aneurisk +
                            scale(maxdiam) +
                            type +
                            side,
                          data = dat_i)
  
  dat_i$propensity <- predict(fitmodality, type = "prob")
  dat_i$ipw_causal <- 1 / dat_i$propensity[cbind(1:nrow(dat_i), as.numeric(dat_i$modality))]
  dat_i$ipw_causal <- ifelse(dat_i$ipw_causal > 5, 5, dat_i$ipw_causal)
  
  modality_probs <- predict(fitmodality, newdata = dat_i, type = "probs")
  dat_i$modality_hat <- factor(apply(modality_probs, 1, function(x) which.max(x)),
                               labels = c("Observation", "Endovascular", "Microneurosurgery"))
  
  dat_i <- dat_i %>%
    mutate(modality_hat = factor(modality_hat, levels = c("Observation", "Endovascular", "Microneurosurgery")),
           mR_BaselineMICE = factor(mR_BaselineMICE, levels = 0:6, ordered = FALSE))
  
  dat_long_i <- tidyr::pivot_longer(dat_i,
                                    cols = c(mR_BaselineMICE, mR_Discharge, mR_1),
                                    names_to = "timepoint",
                                    values_to = "mRS") %>%
    mutate(timepoint = factor(timepoint,
                              levels = c("mR_BaselineMICE", "mR_Discharge", "mR_1"),
                              labels = c("Baseline", "Discharge", "1-year FU")),
           mRS = as.numeric(mRS) - 1,
           age_tp = case_when(
             timepoint == "Baseline" & modality == "Observation" ~ age_Consultation - 5/52,
             timepoint == "Baseline" & modality != "Observation" ~ age_Consultation - 5/52,
             timepoint == "Discharge" & modality == "Observation" ~ age_Consultation,
             timepoint == "Discharge" & modality != "Observation" ~ age_Consultation + 6/365+3/12,
             timepoint == "1-year FU" & modality == "Observation" ~ age_Consultation + 1,
             timepoint == "1-year FU" & modality != "Observation" ~ age_Consultation + 1 + 6/365+3/12,
             TRUE ~ NA_real_
           ),
           pseudonym = droplevels(pseudonym))
  
  # Speichern in Listen
  dat_i_list[[i]] <- dat_i
  dat_long_list[[i]] <- dat_long_i
}

###############################################################################
### Copmute the ESS:
###############################################################################
w <- dat_i$ipw_causal
ESS <- (sum(w)^2) / sum(w^2)
ESS/nrow(dat_i)

ESS_treated  <- (sum(dat_i$ipw_causal[!(dat_i$modality == "Observation")])^2) / sum(dat_i$ipw_causal[!(dat_i$modality == "Observation")]^2)
ESS_endo  <- (sum(dat_i$ipw_causal[(dat_i$modality == "Endovascular")])^2) / sum(dat_i$ipw_causal[(dat_i$modality == "Endovascular")]^2)
ESS_micro  <- (sum(dat_i$ipw_causal[(dat_i$modality == "Microneurosurgery")])^2) / sum(dat_i$ipw_causal[(dat_i$modality == "Microneurosurgery")]^2)
ESS_obs  <- (sum(dat_i$ipw_causal[dat_i$modality == "Observation"])^2) / sum(dat_i$ipw_causal[dat_i$modality == "Observation"]^2)

table(dat_i$modality)
ESS_treated/(table(dat_i$modality)["Microneurosurgery"] + table(dat_i$modality)["Endovascular"])
ESS_obs/table(dat_i$modality)["Observation"]

ESS_endo/table(dat_i$modality)["Endovascular"]
ESS_micro/table(dat_i$modality)["Microneurosurgery"]

##################################################################################
# Plot of the PS weights and IPW
##################################################################################
# check the distribution of the IPW weights and the propensity scores:
dat_i <- dat_i_list[[1]]
a <- ggplot(data=dat_i)+
  geom_density(aes(x=propensity[, "Endovascular"]), fill="#005A8C", alpha=0.5) +
  geom_density(aes(x=propensity[, "Microneurosurgery"]), fill="#007F5F", alpha=0.5) +
  geom_density(aes(x=propensity[, "Observation"]), fill="#B07C00", alpha=0.5) +
  theme_classic()+
  xlab("Propensity Score") +
  ylab("Density") +
  theme(text = element_text(size=12),
        axis.text = element_text(size=12))

dat_i$ipw_causal <- 1/dat_i$propensity[cbind(1:nrow(dat_i), as.numeric(dat_i$modality))]
dat_i$ipw_causal <- ifelse(dat_i$ipw_causal >5, 5, dat_i$ipw_causal)

b <- ggplot(data=dat_i)+
  geom_density(aes(x=ipw_causal, fill=modality), alpha=0.5) +
  theme_classic()+
  xlab("Inverse Probability Weights") +
  ylab("Density") +
  xlim(c(0,5))+
  theme(legend.position = "bottom")+
  scale_fill_manual(values=darkColorPalette[c(1,5,3)])+
  theme(text = element_text(size=12),
        axis.text = element_text(size=12),
        legend.text = element_text(size=12))


p <- ggarrange(a,b, ncol=1, nrow=2)
p

###############################################################################
# Love plot:
###############################################################################
# Step 1: specify the covariates used in the propensity score model
covariates <- c("gender","age_Consultation",
                "Positive_familial_history",
                "prf_hypertension_awareness",
                "prf_smokedtobacco",
                "aneurisk", "maxdiam",
                "type","side",
                "recommended")

# Step 2: calculate weighted balance using cobalt
bal_tab <- bal.tab(
  x = dat_i[, covariates],
  treat = dat_i$modality,
  weights = dat_i$ipw_causal,
  method = "weighting",
  estimand = "ATE",
  un= T, 
  continuous="std",
  binary="std"
)

# Step 3: generate the love plot
love.plot(
  bal_tab,
  stats = "mean.diffs",      # standardized mean differences
  
  abs = TRUE,
  var.order = "unadjusted",  # order covariates by unweighted imbalance
  threshold = 0.2,
  stars = "std",             # distinguishes unweighted vs weighted
  colors = c("grey40", "black"), # red = before, blue = after
  shapes = c(15, 17),
  title = "Covariate Balance (SMDs) Before and After IPW",
  line = T
)


###############################################################################
# Test Proportional Odds assumption for ordinal mRS
##############################################################################
# Ensure mRS is ordered
datlong$mRS <- factor(datlong$mRS, ordered = TRUE)

# Fit a proportional odds model with one predictor of interest
# Replace 'predictor' with your actual variable
model <- ordinal::clm(mRS ~ age_Consultation, data = datlong)

# Generate values for predictor across its range
x_vals <- seq(min(datlong$age_Consultation), max(datlong$age_Consultation), length.out = 100)

# Get cumulative probabilities P(Y >= j) for each threshold
thresholds <- model$alpha  # cumulative intercepts
beta <- coef(model)["age_Consultation"]

cum_probs <- sapply(1:length(thresholds), function(j) {
  plogis(-thresholds[j] + beta * x_vals)
})

# Convert to data frame for ggplot
plot_df <- data.frame(
  age_Consultation = rep(x_vals, times = length(thresholds)),
  Threshold = factor(rep(paste0("≥", 1:length(thresholds)), each = length(x_vals))),
  CumProb = as.vector(cum_probs)
)

# Plot
ggplot(plot_df, aes(x = age_Consultation, y = CumProb, color = Threshold)) +
  geom_line(size = 1.2) +
  labs(x = "Predictor", y = "P(Y ≥ j)", 
       title = "Cumulative Probability Curves (Proportional Odds Model)") +
  theme_minimal()


model_polr <- polr(mRS ~ gender+
                     age_Consultation+
                     Positive_familial_history+
                     prf_hypertension_awareness+
                     prf_smokedtobacco+
                     maxdiam+
                     type+
                     side+
                     recommended+
                     modality,
                   data = datlong, Hess=TRUE)

# Brant test
brant_test <- brant::brant(model_polr)
print(brant_test)


################################################################################
# 1. Unweighted model for comparison:
################################################################################
fitsStandard <- lapply(dat_long_list, function(dat_long_i) {
  glmmTMB(mRS ~
            gender +
            Positive_familial_history +
            prf_hypertension_awareness +
            prf_smokedtobacco +
            aneurisk +
            scale(maxdiam) +
            type +
            side +
            timepoint +
            modality * scale(age_tp) +
            (1 | pseudonym),
          data = dat_long_i,
          family = poisson(),
          dispformula = ~1)})


fitStandard <- pool(fitsStandard)
ForestPlot_MICE(fitStandard, title ="Standard Model, Poisson (MICE)" )

## test for underdispersion:
fitstandard10 <- glmmTMB(mRS ~
                           gender +
                           Positive_familial_history +
                           prf_hypertension_awareness +
                           prf_smokedtobacco +
                           aneurisk +
                           scale(maxdiam) +
                           type +
                           side +
                           timepoint +
                           modality * scale(age_tp) +
                           (1 | pseudonym),
                         data = dat_long_i,
                         family = poisson(),
                         dispformula = ~1)

summary(fitstandard10)
testDispersion(fitstandard10)  
# Here for the synthetic data we do not obtain any underdispersion


#COMPOIS-Model, to adjust for the underdispersion
fitsCOMpoisStandard <- lapply(dat_long_list, function(dat_long_i) {
  glmmTMB(mRS ~
            gender +
            Positive_familial_history +
            prf_hypertension_awareness +
            prf_smokedtobacco +
            aneurisk +
            scale(maxdiam) +
            type +
            side +
            modality * scale(age_tp) +
            timepoint +
            (1 | pseudonym),
          data = dat_long_i,
          family = compois(),
          dispformula = ~1)})

fitCOMpoisStandard <- pool(fitsCOMpoisStandard)
bic_summary(fitsCOMpoisStandard) # find the bic per model

# fit a model without interaction between age and modality:
fitsCOMpoisStandardwoI <- lapply(dat_long_list, function(dat_long_i) {
  glmmTMB(mRS ~
            gender +
            Positive_familial_history +
            prf_hypertension_awareness +
            prf_smokedtobacco +
            aneurisk +
            scale(maxdiam) +
            type +
            side +
            modality + scale(age_tp) +
            timepoint +
            (1 | pseudonym),
          data = dat_long_i,
          family = compois(),
          dispformula = ~1)})

fitCOMpoisStandardwoI <- pool(fitsCOMpoisStandardwoI)

# D1-test for the interaction:
anova_pool <- mice::D1(
  fit0 = fitsCOMpoisStandardwoI,
  fit1 = fitsCOMpoisStandard)

## save the results:
# saveRDS(fitsStandard, file = "./results/fits/fitsStandard.rds")
# saveRDS(fitStandard,  file = "./results/fits/fitStandard.rds")
# saveRDS(fitCOMpoisStandardwoI,  file = "./results/fits/fitCOMpoisStandardwoI.rds")
# saveRDS(fitsCOMpoisStandardwoI, file = "./results/fits/fitsCOMpoisStandardwoI.rds")
# saveRDS(fitsCOMpoisStandard, file = "./results/fits/fitsCOMpoisStandard.rds")
# saveRDS(fitCOMpoisStandard,  file = "./results/fits/fitCOMpoisStandard_pooled.rds")

# or simply load all fits, and then print the summaries and forestplots:
fitsStandard <- readRDS("./results/fits/fitsStandard.rds")
fitStandard <- readRDS("./results/fits/fitStandard.rds")
fitsCOMpoisStandardwoI <- readRDS("./results/fits/fitsCOMpoisStandardwoI.rds")
fitCOMpoisStandardwoI <- readRDS("./results/fits/fitCOMpoisStandardwoI.rds")
fitsCOMpoisStandard <- readRDS("./results/fits/fitsCOMpoisStandard.rds")
fitCOMpoisStandard <- readRDS("./results/fits/fitCOMpoisStandard_pooled.rds")


# Forestplots for the individual models:
ForestPlot_MICE(fitStandard, title ="Standard Model, Poisson (MICE)", col_def = "black")
ForestPlot_MICE(fitCOMpoisStandard, title = "Standard Model, COMPOIS (MICE)", col_def = "black")

# Forestplot comparing both models:
ForestPlot_MICE_multi(
  models = list(fitCOMpoisStandard, fitStandard),
  model_names = c( "COMPoisson", "Poisson"),
  title = "Standard, Model Comparison (MICE)")

## Forest Plots for the individual and pooled models per variable:
res <- forestplot_mice_param(
  fits_list = fitsCOMpoisStandard,
  param = "modalityEndovascular:scale(age_tp)", # change the variable here to see other effects
  title = "",
  exponentiate = TRUE,
  xlim = c(0,2.5))

res$plot

################################################################################
### 2. IPW model:
################################################################################
# Poisson Model with interaction between modality and age_tp
# Note that we do not use "Discharge" here anymore

fits_IPW <- lapply(dat_long_list, function(dat_long_i) {
  fitIPW <- glmmTMB(mRS ~
                      gender +
                      Positive_familial_history +
                      prf_hypertension_awareness +
                      prf_smokedtobacco +
                      aneurisk +
                      scale(maxdiam) +
                      type +
                      side +
                      modality * scale(age_tp) +
                      timepoint +
                      (1 | pseudonym),
                    data = dat_long_i[dat_long_i$timepoint != "Discharge",],
                    family = poisson,
                    weights = ipw_causal)
})

fitIPW <- pool(fits_IPW)

# check for underdispersion in one of the imputed data sets:
fitIPW10 <- glmmTMB(mRS ~
                      gender +
                      Positive_familial_history +
                      prf_hypertension_awareness +
                      prf_smokedtobacco +
                      aneurisk +
                      scale(maxdiam) +
                      type +
                      side +
                      modality * scale(age_tp) +
                      timepoint +
                      (1 | pseudonym),
                    data = dat_long_i[dat_long_i$timepoint != "Discharge",],
                    family = poisson,
                    weights = ipw_causal)

summary(fitIPW10)
testDispersion(fitIPW10)  #Yep, we have underdispersion

# COMPOIS Model, to adjust for the underdispersion
fits_IPW_COMPOIS <- lapply(dat_long_list, function(dat_long_i) {
  fitIPW <- glmmTMB(mRS ~
                      gender +
                      Positive_familial_history +
                      prf_hypertension_awareness +
                      prf_smokedtobacco +
                      aneurisk +
                      scale(maxdiam) +
                      type +
                      side +
                      timepoint +
                      modality * scale(age_tp) +
                      (1 | pseudonym),
                    data = dat_long_i[dat_long_i$timepoint != "Discharge",],
                    family = compois,
                    weights = ipw_causal)
})


fitIPW_COMPOIS <- pool(fits_IPW_COMPOIS)
bic_summary(fits_IPW_COMPOIS)

# fit a model without interaction between age and modality:
fits_IPW_COMPOISwoI <- lapply(dat_long_list, function(dat_long_i) {
  fitIPW <- glmmTMB(mRS ~
                      gender +
                      Positive_familial_history +
                      prf_hypertension_awareness +
                      prf_smokedtobacco +
                      aneurisk +
                      scale(maxdiam) +
                      type +
                      side +
                      timepoint +
                      modality + scale(age_tp) +
                      (1 | pseudonym),
                    data = dat_long_i[dat_long_i$timepoint != "Discharge",],
                    family = compois,
                    weights = ipw_causal)
})

fitIPW_COMPOISwoI <- pool(fits_IPW_COMPOISwoI)

# D1-test for the interaction:
anova_pool <- mice::D1(
  fit0 = fits_IPW_COMPOISwoI,
  fit1 = fits_IPW_COMPOIS)


## save the results:
# saveRDS(fits_IPW, file = "./results/fits/fitsIPW.rds")
# saveRDS(fitIPW, file = "./results/fits/fitIPW.rds")
# saveRDS(fits_IPW_COMPOISwoI, file = "./results/fits/fitsIPW_COMPOISwoI.rds")
# saveRDS(fitIPW_COMPOISwoI,  file = "./results/fits/fitIPW_COMPOISwoI.rds")
# saveRDS(fits_IPW_COMPOIS, file = "./results/fits/fitsIPW_COMPOIS.rds")
# saveRDS(fitIPW_COMPOIS,  file = "./results/fits/fitIPW_COMPOIS.rds")

# or simply load all fits, and then print the summaries and forestplots:
fits_IPW <- readRDS("./results/fits/fitsIPW.rds")
fitIPW <- readRDS("./results/fits/fitIPW.rds")
fits_IPW_COMPOISwoI <-  readRDS(file = "./results/fits/fitsIPW_COMPOISwoI.rds")
fitIPW_COMPOISwoI<- readRDS(file = "./results/fits/fitIPW_COMPOISwoI.rds")
fits_IPW_COMPOIS <- readRDS("./results/fits/fitsIPW_COMPOIS.rds")
fitIPW_COMPOIS <- readRDS("./results/fits/fitIPW_COMPOIS.rds")

## Forest Plots:
ForestPlot_MICE(fitIPW, title ="IPW Model, Poisson (MICE)" )
ForestPlot_MICE(fitIPW_COMPOIS, title ="IPW Model, COMPOIS (MICE)" )

ForestPlot_MICE_multi(
  models = list(fitIPW_COMPOIS, fitIPW),
  model_names = c( "COMPoisson", "Poisson"),
  title = "IPW, Model Comparison (MICE)")

## Forest plots for the individual and pooled models:
res <- forestplot_mice_param(
  fits_list = fits_IPW_COMPOIS,
  param = "modalityMicroneurosurgery:scale(age_tp)", # change the variable here to see other effects
  title = "",
  exponentiate = TRUE,
  xlim = c(0,2.5))

res$plot


################################################################################
# 3. All models in one plot
################################################################################
ForestPlot_MICE_multi(
  models = list(fitIPW_COMPOIS, fitCOMpoisStandard ),
  model_names = c("COMPois, IPW","COMPois, unweighted"),
  title = "",
  col_defs = c("grey40", "grey0"))

## Forest Plot for multiple models at once, e.g. for the sensitity analysis performed:
ForestPlot_MICE_multi_n(
  models = list(fitCOMpoisStandard, 
                fitIPW_COMPOIS,                 
                fitIPW_COMPOISwoI),
  model_names = c("0","IPW",
                  "IPW w/ interaction"),
  title = "",
  col_defs = c("#041E42","#F2A900", "#005EB8"))


################################################################################
# 4. Get the average treatment effects (ATE)
################################################################################
pool_contrasts(fitsCOMpoisStandard, "modality")
pool_contrasts(fits_IPW_COMPOIS, "modality")



##############################################################################
# 5. ABN
#############################################################################
## Poisson model: 
keep <- c("gender","age_Consultation",
          "Positive_familial_history",
          "prf_hypertension_awareness",
          "prf_smokedtobacco",
          "aneurisk", "maxdiam",
          "type","side",
          "recommended",
          "modality",
          "mR_Baseline",
          "mR_Discharge",
          "mR_1")

datPoisson <- dat_s %>%
  mutate(Positive_familial_history = 
           factor(Positive_familial_history, levels=c("Yes", "No"))) %>%
  mutate(side = factor(side, levels = c("Left", "Midline", "Right"))) %>%
  mutate(age_Consultation = scale(age_Consultation)) %>% 
  mutate(maxdiam = scale(maxdiam)) %>%
  mutate(mR_Baseline = as.numeric(mR_BaselinePredMICEmean)-1,
         mR_Discharge = as.numeric(mR_Discharge)-1,
         mR_1 = as.numeric(mR_1)-1) %>%
  select(all_of(keep)) %>%
  na.omit()

dists <- list(gender =  "binomial",
              age_Consultation = "gaussian",
              Positive_familial_history = "binomial",
              prf_hypertension_awareness = "binomial",
              prf_smokedtobacco ="multinomial",
              aneurisk = "multinomial",
              maxdiam = "gaussian",
              type= "binomial",
              side="multinomial",
              recommended = "multinomial",
              modality = "multinomial",
              mR_Baseline = "poisson",
              mR_Discharge = "poisson",
              mR_1 = "poisson")


allowmatASIS<- readxl::read_excel("data/prior_knowledge.xlsx", sheet="allowmat") %>%
  tibble::column_to_rownames('ChildParent') %>%
  as.matrix()

## create banmat: 1-allowmat
banmatASIS <- 1-allowmatASIS 
banmatASIS <- banmatASIS[c(keep),
                         c(keep)]

retainmatASIS <- readxl::read_excel("data/prior_knowledge.xlsx", sheet="retainmat") %>%
  tibble::column_to_rownames('ChildParent') %>%
  as.matrix()

retainmatASIS <- retainmatASIS[c(keep),
                               c(keep)]

## single ABN
set.seed("0815")
starttime <- Sys.time()
chacheMLE <- buildScoreCache(data.df = datPoisson,
                             data.dists = dists,
                             method = "mle",
                             dag.banned = banmatASIS,
                             dag.retained = retainmatASIS,
                             centre = TRUE,
                             max.parents = length(datPoisson)-1)

endtime <- Sys.time()
cat(paste("\n****************************\nEnd ABN., Time used [h]:",
          round(difftime(endtime, starttime, units = "hours"), 4)))

dagMP <- mostProbable(score.cache = chacheMLE, prior.choice = 1, score="bic")
par(mfrow=c(1,1))
plot(dagMP)

# CPDAG
par(mfrow=c(1,2))
dag <- bnlearn::empty.graph(names(datPoisson))
bnlearn::amat(dag) <- t(dagMP$dag)
plot(dag, main ="DAG")
cons.cpdag <- bnlearn::cpdag(dag)
cons.amat <- bnlearn::amat(cons.cpdag)
# # Remove all edges that are present in the banned matrix (also an adjacency matrix)
cons.amat[cons.amat == 1 & t(banmatASIS) == 1] <- 0
cpdag.cons.amat <- bnlearn::empty.graph(names(datPoisson))
bnlearn::amat(cpdag.cons.amat) <- cons.amat
plot(cpdag.cons.amat, main="CPDAG II")

################################################################################
# 6. Imputations and Consensus DAG with MCMCabn
################################################################################
# Access individual models and consensus
modelmaxparents <- findMaxParents(dat=datPoisson, dists=dists, banmat = banmatASIS,
                                  retainmat = retainmatASIS,
                                  SCORE = "bic",
                                  FILENAMEbase =  paste0("./results/management_effect/"),
                                  PRIOR = 1)

# save the maxparents:
saveRDS(modelmaxparents,
        file = paste0("./results/management_effect/modelmaxparents.rds"))

modelmaxparents <- readRDS(file = paste0("./results/management_effect/modelmaxparents.rds"))

# Run models & save
results <- run_MCMC_ABN_all_imputations(
  dat_i_list = dat_i_list[1:10],
  runname = "management_effect",
  dists = dists,
  banmatASIS = banmatASIS,
  retainmatASIS = retainmatASIS,
  modelmaxparents = modelmaxparents,
  MCMC.SEEDS=c("081514", "123456", "2025", "3141"),
  MCMC.SCHEME = c(10000,0,0),
  burnin.length = 1000,
  thinningsteps = 2,
  keep = keep,
  baseline = 0
)

## or
results<- load_MCMC_ABN_results("management_effect")

# Later or in a new session:
consensus_results <- compute_consensus_dag_from_saved("management_effect")

# Access final DAG:
consensus_results$consensus


# Run the analysis for each imputation
for (i in 1:10) {
  analyze_MCMC_ABN_results(
    modelMCMC = results$models[[paste0("imp", i)]],
    dag = results$dags[[paste0("imp", i)]],
    dat = dat_i_list[[i]],
    dists = dists,
    banmatASIS = banmatASIS,
    retainmatASIS = retainmatASIS,
    modelmaxparents = modelmaxparents,
    runname = paste0("management_effect"),
    impnumber = paste0("imp_", i),
    burnin.length = 2000,
    thinningsteps = 2
  )
}

# And finally for the global consensus
analyze_MCMC_ABN_results(
  modelMCMC = results$models[[1]],  # reuse any model structure
  dag = consensus_results$consensus,
  dat = dat_i_list[[1]],            # representative data set
  dists = dists,
  banmatASIS = banmatASIS,
  retainmatASIS = retainmatASIS,
  modelmaxparents = modelmaxparents,
  runname = "management_effect",
  impnumber = "consensus",
  burnin.length = 2000,
  thinningsteps = 2
)

###############################################################################
# 7. Compute the parameters for the individual nodes using the
#    imputed data sets and pooling the results with Rubin's rules
###############################################################################

## Effects in the individual models:
results$models$imp1$fabn.maxpar
results$models$imp2$fabn.maxpar

fit <- glm(mR_1 ~mR_Discharge+maxdiam, data=datPoisson, family=poisson)
biostatUZH::tableRegression(fit)

#fit for each imputed data set:
fitsmR1 <- lapply(dat_i_list, function(dat_i) {
  glm(I(as.numeric(mR_1)-1) ~
        I(as.numeric(mR_Discharge)-1)+
        maxdiam,
      data = dat_i,
      family = poisson)})
fitmR1 <- pool(fitsmR1)
summary(fitmR1)
ForestPlot_MICE(fitmR1)


fitsmRdischarge <- lapply(dat_i_list, function(dat_i) {
  glm(I(as.numeric(mR_Discharge)-1) ~
        I(as.numeric(mR_BaselineMICE)-1),
      data = dat_i,
      family = poisson)})
fitmRdischarge <- pool(fitsmRdischarge)
summary(fitmRdischarge)
ForestPlot_MICE(fitmRdischarge)

fitsmodality <- lapply(dat_i_list, function(dat_i) {
  multinom(modality ~
             maxdiam+aneurisk+type,
           data = dat_i)})
fitmodality <- pool(fitsmodality)
summary(fitmodality)


