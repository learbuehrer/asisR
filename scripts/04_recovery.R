##############################################################################
# Analysis of the Recovery
##############################################################################

############################################################################
# Regression models:
############################################################################
# Poisson model w/o any interaction
fitrecovery <- glmmTMB(as.numeric(mR_1) ~
                      gender +
                      Positive_familial_history +
                      prf_hypertension_awareness +
                      prf_smokedtobacco +
                      aneurisk +
                      scale(maxdiam) +
                      type +
                      side +
                      modality + scale(age_Consultation) +
                      as.numeric(mR_Discharge),
                    data = dat_s,
                    family = poisson(),
                    dispformula = ~1)

# Poisson model with interaction term for modality and age
fitrecoveryinteraction <- glmmTMB(as.numeric(mR_1) ~
                                 gender +
                                 Positive_familial_history +
                                 prf_hypertension_awareness +
                                 prf_smokedtobacco +
                                 aneurisk +
                                 scale(maxdiam) +
                                 type +
                                 side +
                                 modality * scale(age_Consultation) +
                                 as.numeric(mR_Discharge),
                               data = dat_s,
                               family = poisson(),
                               dispformula = ~1)

anova(fitrecovery, fitrecoveryinteraction)
# Conclusion: Interaction does not add information to the model.

# Check for an interaction with mRS at Discharge and applied management strategy
fitrecoveryinteractionmrsD <- glmmTMB(as.numeric(mR_1) ~
                                     gender +
                                     scale(age_Consultation) +
                                     Positive_familial_history +
                                     prf_hypertension_awareness +
                                     prf_smokedtobacco +
                                     aneurisk +
                                     scale(maxdiam) +
                                     type +
                                     side +
                                     modality * 
                                     as.numeric(mR_Discharge),
                                   data = dat_s,
                                   family = poisson(),
                                   dispformula = ~1)
anova(fitrecovery, fitrecoveryinteractionmrsD)
# Conclusion: Interaction does in this case add information to the model
# this is not the case for the original data.

###############################################################################
# Check dispersion in the Poisson model:
###############################################################################
# Calculate Pearson residuals
pearson_resid <- residuals(fitrecovery, type = "pearson")
sum_sq_pearson <- sum(pearson_resid^2)
rdf <- df.residual(fitrecovery)
dispersion <- sum_sq_pearson / rdf
cat("Estimated dispersion:", dispersion, "\n")

# Simulate residuals from your glmmTMB model
res <- simulateResiduals(fitrecovery)
plot(res)
testDispersion(res)

# Conclusion: we obtain underdispersion in the Poisson model.

###############################################################################
# Compois model:
###############################################################################
fitrecoverycompois <- glmmTMB(as.numeric(mR_1) ~
                             gender +
                             Positive_familial_history +
                             prf_hypertension_awareness +
                             prf_smokedtobacco +
                             aneurisk +
                             scale(maxdiam) +
                             type +
                             side +
                             modality +
                             scale(age_Consultation) +
                             as.numeric(mR_Discharge),
                           data = dat_s,
                           family = compois(),
                           dispformula = ~1)

summary(fitrecoverycompois)
SummaryTable(fitrecoverycompois, n=16, digits = 2)

ForestPlot_GLM(fitrecoverycompois, title = "",
               exponentiate = T, term1 = "Better \nCondition",
               term2 = "Worse \nCondition", col_def  = "black",
               measure = "IRR")


############################################################################
# ABN model:
############################################################################
## Poisson model: 
keep <- c("gender","age_Consultation",
          "Positive_familial_history",
          "prf_hypertension_awareness",
          "prf_smokedtobacco",
          "aneurisk", "maxdiam",
          "type","side",
          "recommended",
          "modality",
          "mR_Discharge",
          "mR_1")

datPoisson <- dat_s %>%
  mutate(Positive_familial_history = 
           factor(Positive_familial_history, levels=c("Yes", "No"))) %>%
  mutate(side = factor(side, levels = c("Left", "Midline", "Right"))) %>%
  mutate(age_Consultation = scale(age_Consultation)) %>% 
  mutate(maxdiam = scale(maxdiam)) %>%
  mutate(mR_Discharge = as.numeric(mR_Discharge)-1,
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

#############################################################################
### fit a single ABN:
#############################################################################
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

fitBN <- fitAbn(object= dagMP, method = "mle")

#############################################################################
### fit the mcmcABN for recovery:
#############################################################################
modelmaxparents <- findMaxParents(dat=datPoisson,
                                  dists=dists, 
                                  banmat = banmatASIS,
                                  retainmat = retainmatASIS,
                                  SCORE = "bic",
                                  FILENAMEbase =  paste0("./results/recovery/"),
                                  PRIOR = 1)


runname ="recovery" #only adapt path here and not later
MCMC.SEEDS=c("081514", "123456", "2025", "3141")
MCMC.SCHEME = c(10000,0,0)
burnin.length = 2000
thinningsteps = 2

ScorePlotMaxParent(modelmaxparents$net.scores, FILENAMEbase =  
                     paste0("./results/",runname, "/"))

DAGMaxParent(net.scores.dags = modelmaxparents$net.scores.dags, net.scores = modelmaxparents$net.scores,
             par1=3, par2=4,
             retainmat = retainmatASIS,
             FILENAMEbase = paste0("./results/",runname, "/"),
             filename = "DAGMaxParents")

## MCMC runs:
modelMCMC <- mcmcmABN(FILENAMEbase =  paste0("./results/", runname, "/"),
                      dat = datPoisson, dists= dists,
                      banmat = banmatASIS, 
                      retainmat = retainmatASIS,
                      FILENAME = runname,
                      max.par = modelmaxparents$max.par,
                      MCMC.SEEDS = MCMC.SEEDS,
                      MCMC.SCHEME = MCMC.SCHEME,
                      PROB.REV = 0.05, PROB.MBR = 0, MCMC.PRIOR=1,
                      SCORE="bic")

saveRDS(modelMCMC,
        file = paste0("./results/", runname, "/",
                      "modelMCMC_", runname, ".rds"))

modelMCMC <- readRDS(
  file = paste0("./results/", runname, "/",
                "modelMCMC_", runname, ".rds"))

## generate all plots:
## postpreprocessing - THINNING and BURNIN
mcmc.out.list.burn <- postBURNin(mcmc.out.list = modelMCMC$mcmc.out.list, burnin.length = burnin.length)
mcmc.out.list.thin <- postTHINN(mcmc.out.list = modelMCMC$mcmc.out.list, thinningsteps = thinningsteps)
mcmc.out.list.burn.thin <- postTHINN(mcmc.out.list = mcmc.out.list.burn, thinningsteps = thinningsteps)

# original list:
mc.out.score.list <- lapply(modelMCMC$mcmc.out.list, function(x) mcmc(x$scores))
list.mc.out.score <- mcmc.list(mc.out.score.list)


# thinned only
mc.out.thin.1 <- mcmc.out.list.thin[[1]]
mc.out.thin.2 <- mcmc.out.list.thin[[2]]
mc.out.thin.3 <- mcmc.out.list.thin[[3]]
mc.out.thin.4 <- mcmc.out.list.thin[[4]]

mc.out.thin.score.1 <- mcmc(mc.out.thin.1$scores)
mc.out.thin.score.2 <- mcmc(mc.out.thin.2$scores)
mc.out.thin.score.3 <- mcmc(mc.out.thin.3$scores)
mc.out.thin.score.4 <- mcmc(mc.out.thin.4$scores)

list.mc.out.thin.score <- mcmc.list(mc.out.thin.score.1, mc.out.thin.score.2,
                                    mc.out.thin.score.3, mc.out.thin.score.4)


# burned and thinned
mc.out.burn.thin.1 <- mcmc.out.list.burn.thin[[1]]
mc.out.burn.thin.2 <- mcmc.out.list.burn.thin[[2]]
mc.out.burn.thin.3 <- mcmc.out.list.burn.thin[[3]]
mc.out.burn.thin.4 <- mcmc.out.list.burn.thin[[4]]

mc.out.burn.thin.score.1 <- mcmc(mc.out.burn.thin.1$scores)
mc.out.burn.thin.score.2 <- mcmc(mc.out.burn.thin.2$scores)
mc.out.burn.thin.score.3 <- mcmc(mc.out.burn.thin.3$scores)
mc.out.burn.thin.score.4 <- mcmc(mc.out.burn.thin.4$scores)

list.mc.out.burn.thin.score <- mcmc.list(mc.out.burn.thin.score.1, mc.out.burn.thin.score.2,
                                         mc.out.burn.thin.score.3, mc.out.burn.thin.score.4)

mc.out.burn.thin.dag.1 <- mc.out.burn.thin.1$dags
mc.out.burn.thin.dag.2 <- mc.out.burn.thin.2$dags
mc.out.burn.thin.dag.3 <- mc.out.burn.thin.3$dags
mc.out.burn.thin.dag.4 <- mc.out.burn.thin.4$dags

list.mc.out.burn.thin.dag <- abind::abind(mc.out.burn.thin.dag.1, mc.out.burn.thin.dag.2,
                                          mc.out.burn.thin.dag.3, mc.out.burn.thin.dag.4)


## Quality Assessment:
MCMCQualityAssessment(mcmcoutput = list.mc.out.score,
                      mcmcoutputburnedthinned = list.mc.out.burn.thin.score,
                      FILENAMEbase = paste0("./results/", runname, "/"), 
                      filename = "GelmanPlot")


TracePlot(mcouthin1 = mc.out.burn.thin.1,
          mcoutthin2 = mc.out.burn.thin.2,
          mcoutthin3 = mc.out.burn.thin.3,
          mcoutthin4 = mc.out.burn.thin.4,
          MCMC.SEEDS=c("0815", "0510", "1107", "2611"),
          dat=datPoisson,dists= dists,
          METHOD="mle", modeloutput = modelMCMC,
          filename=paste0("TracePlot_thinned", thinningsteps, "burned", burnin.length),
          FILENAMEbase=paste0("./results/", runname, "/"))

## Consensus DAGs:
dag <- apply(list.mc.out.burn.thin.dag, 1:2, mean)


plotDAGthresholds(dag=dag,
                  dists = dists,
                  dat=datPoisson,
                  threshold = c(0.5, 0.6),
                  FILENAMEbase = paste0("./results/", runname, "/"),
                  filename = "DAG_different_thresholds")


# Arc strength significance threshold
arc.stren.sign.threshold <-arc.stren.threshold(dag,
                                               method = "l1")


ecdf <-ecdf_plot_ggplot(strengthdag=dag,
                        arc.stren.sign.threshold = arc.stren.sign.threshold,
                        filename = "ECDF_arc_sten_threshold2",
                        FILENAMEbase = paste0("./results/", runname, "/"))


plotDAGthresholds(threshold = arc.stren.sign.threshold,
                  dag=dag,
                  dists = dists,
                  dat=datPoisson,
                  FILENAMEbase = paste0("./results/", runname, "/"),
                  filename = "DAG_arcstreng_threshold")


ThresholdHeatmap(dag=dag,banmat = banmatASIS,
                 threshold = arc.stren.sign.threshold,
                 FILENAMEbase = paste0("./results/", runname, "/"),
                 filename = "ThresholdHeatmap")



AdjacencyMatrix(banmat = banmatASIS,
                dag=dag,
                threshold = arc.stren.sign.threshold,
                filename = "AdjacencyMatrix",
                FILENAMEbase = paste0("./results/", runname, "/"))

# Comparison Plot:
abnMCMCabnComparision(threshold = arc.stren.sign.threshold+0.01,
                      dag=dag,
                      dists = dists,
                      dat=datPoisson,
                      prior=1,
                      banmat=banmatASIS,
                      retainmat=retainmatASIS,
                      max.parents = modelmaxparents$max.par,
                      FILENAMEbase = paste0("./results/", runname, "/"),
                      filename = "DAG_Comparison")


# CPDAG for the MCMCDAG
dagMCMC <- dag
dagMCMC[dagMCMC > arc.stren.sign.threshold] <- 1
dagMCMC[dagMCMC <= arc.stren.sign.threshold] <- 0
colnames(dagMCMC) = rownames(dagMCMC) <- names(datPoisson)

# Open a PNG device
png(filename = paste0(paste0("results/", runname,"/CPDAG.png")), 
    width = 14*600, height = 9*600, res = 600)
par(mfrow=c(1,2))

dagBN <- bnlearn::empty.graph(names(datPoisson))
bnlearn::amat(dagBN) <- t(dagMCMC)
plot(dagBN, main ="DAG")
cons.cpdag <- bnlearn::cpdag(dagBN)

# Convert the CPDAG to an adjacency matrix
cons.amat <- bnlearn::amat(cons.cpdag)
# # Remove all edges that are present in the banned matrix (also an adjacency matrix)
cons.amat[cons.amat == 1 & t(banmatASIS) == 1] <- 0
cpdag.cons.amat <- bnlearn::empty.graph(names(datPoisson))
bnlearn::amat(cpdag.cons.amat) <- cons.amat
plot(cpdag.cons.amat, main="CPDAG II")
dev.off()

# SVG plot:
dag_threshold <- dag
dag_threshold[dag_threshold > arc.stren.sign.threshold+0.01] <- 1
dag_threshold[dag_threshold <= arc.stren.sign.threshold+0.01] <- 0
edgestren <- dag
edgestren[edgestren <= arc.stren.sign.threshold+0.01] <- 0
cons.dag.plt.edgestrength <- plotAbn(dag = dag_threshold,
                                     data.df = dat,
                                     data.dists = dists,
                                     digits = 2,
                                     edge.strength = edgestren,
                                     plot = T, 
                                     main=arc.stren.sign.threshold)

## Effects for the recovery model:
modelMCMC$fabn.maxpar
fit <- glm(mR_1 ~mR_Discharge+maxdiam+Positive_familial_history+prf_hypertension_awareness, data=datPoisson, family=poisson)
biostatUZH::tableRegression(fit)

