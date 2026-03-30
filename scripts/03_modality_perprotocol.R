#################################################################################
# Management and Per protocol Analysis
#################################################################################
keep = c("pseudonym","gender", "age_e", "age_Consultation",
         "Positive_familial_history",
         "prf_hypertension_awareness",
         "prf_smokedtobacco",
         "aneurisk", "type", "side", "maxdiam",
         "recommended", "modality", "perprotocol",
         "mR_Discharge", "mR_1")

# create a new dedicated data set for this analysis, as 
# a finer granularity of the recommendation is needed
# an the unlogarithmized IA size is wanted

dat_pp <- dat_s %>%
  mutate(maxdiam = exp(maxdiam)) %>%
  select(all_of(keep)) %>%
  na.omit()%>% 
  mutate(modality_treatment = ifelse(modality %in% c("Endovascular", "Microneurosurgery"),
                                     "Applied: Treatment", "Applied: Observation")) %>%
  mutate(recommendation_treatment = ifelse(recommended %in% c("Microneurosurgery", "Endovascular",
                                                              "Either Treatment"),
                                           "treatment", "observation")) %>% 
  mutate(treat_decision = ifelse(modality == "Observation" &
                                   recommendation_treatment == "treatment",
                                 "decision to observe",
                                 ifelse(modality %in% c("Endovascular", "Microneurosurgery")&
                                          recommendation_treatment =="observation",
                                        "decision to treat", 
                                        ifelse(modality == "Observation"&
                                                 recommendation_treatment =="observation",
                                               "per protocol observed", "per protocol treated")))) %>% 
  mutate(pp = ifelse(modality == "Observation" &
                       recommendation_treatment == "observation",
                     "per protocol",
                     ifelse(modality %in% c("Endovascular", "Microneurosurgery")&
                              recommendation_treatment =="treatment",
                            "per protocol", "not per protocol")))

table(dat_pp$treat_decision)

# Frequency table of recommended vs applied management per modality:
table(dat_pp$recommended, dat_pp$perprotocol, dat_pp$modality)

###############################################################################
# Binomial Regression Model for the Management Strategy:
###############################################################################
# Create binary variable for modality:
dat_pp$mod <- ifelse(dat_pp$modality=="Observation", 0, 1)

## Regression for the modality choosen:
fit_modality <- glm(mod~gender+
                      Positive_familial_history+
                      prf_hypertension_awareness+
                      prf_smokedtobacco+
                      aneurisk+
                      scale(maxdiam)+
                      type+
                      side+
                      scale(age_Consultation), data = dat_pp,
                    family="binomial")
summary(fit_modality)
xtable(SummaryTable(fit_modality, n=13, digits = 3))


ForestPlot_GLM(fit_modality, title = "",
               exponentiate = T, term1 = "Towards \nObservation",
               term2 = "Towards \nTreatment", col_def  = "black")


###############################################################################
# Descriptives - Per-Protocol:
##############################################################################
dat_pp %>% 
  mutate(recommended = factor(recommended, levels = c("Either Treatment","Microneurosurgery", 
                                                      "Endovascular", 
                                                      "Observation"),
                              labels = c("Either \nTreatment",
                                         "Microneurosurgery",
                                         "Endovascular",
                                         "Observation"))) %>%
  mutate(modality = factor(modality, levels = c("Microneurosurgery","Endovascular", "Observation"),
                           labels = c("Microneuro-\nsurgery", "Endovascular",
                                      "Observation"))) %>% 
  ggplot(aes(axis1 = recommended, axis2 = modality)) +
  geom_alluvium(aes(fill = modality), show.legend = F) +
  geom_stratum(fill = c(darkColorPalette[c(1,6,3,8,3,1,6)]), alpha=1) +
  scale_fill_manual(values = c(darkColorPalette[c(3,6,1)]), name="",
                    guide = guide_legend(reverse = TRUE)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), col="white") +
  scale_x_discrete(limits = c("Recommended \nManagement", "Applied \nManagement"))+
  theme_classic()+
  ylab("")+
  theme(legend.position = "bottom") +
  theme(text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


# Scatter plot of IA size vs Age at Consultation, colored by Per-Protocol
p <- dat_pp %>%
  mutate(perprotocol = factor(perprotocol, levels=c("Yes", "No"))) %>%
  ggplot(aes(x = age_Consultation, y = maxdiam)) +
  geom_point(aes(color = perprotocol, shape = modality_treatment)) +
  geom_smooth(aes(color = perprotocol, fill = perprotocol)) +
  scale_color_manual(values = c("darkgreen", "darkred"), name = "Per Protocol") +
  scale_fill_manual(values = c("darkgreen", "darkred"), name = "Per Protocol") +
  # scale_size_manual(values = c(1,1.5,2,2.5,3,3.5,4),name = "mRS at \nDischarge") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    plot.margin = margin(t = 10, r = 10, b = 50, l = 10) # extra bottom space
  )+
  xlab("Age at Consultation [Years]")+
  ylab("IA Size [mm]")+
  guides(color = guide_legend(order = 1, nrow = 2, byrow = TRUE),
         fill = guide_legend(order = 1, nrow = 2, byrow = TRUE),
         size = guide_legend(order = 2, nrow = 2, byrow = TRUE, reverse = FALSE))

# Add marginal histograms
ggMarginal(
  p,
  type = "histogram",
  margins = "both",
  size = 5,      # size of the marginal plot
  groupColour = F,
  groupFill = F
)

# plot stratified by applied management strategy: 
p <- dat_pp %>%
  ggplot(aes(x = age_Consultation, y = maxdiam)) +
  geom_point(aes(x = age_Consultation, y = maxdiam,
                 color = pp, pch=recommended), size=2) +
  geom_smooth(aes(color = pp, fill = pp))+
  scale_color_manual(
    values = c(
      "per protocol" = "darkgreen",
      "not per protocol" = "darkred"
    ),
    name = "Per Protocol"
  ) +
  scale_fill_manual(
    values = c(
      "per protocol" = "darkgreen",
      "not per protocol" = "darkred"
    ),
    name = "Per Protocol"
  ) +
  scale_shape_manual(values = c(4,1, 15,17,18),
                     name = "Recommended \nManagement",
                     labels=c("Observation",
                              "Endovascular",
                              "Microneurosurgery", "Either Treatment"))+
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    plot.margin = margin(t = 10, r = 10, b = 50, l = 10)
  ) +
  xlab("Age at Consultation [Years]") +
  ylab("IA Size [mm]") +
  guides(
    color = guide_legend(order = 1, nrow = 2, byrow = TRUE),
    fill = guide_legend(order = 1, nrow = 2, byrow = TRUE),
    size = guide_legend(order = 2, nrow = 2, byrow = TRUE, reverse = FALSE),
    shape = guide_legend(order = 2, nrow = 2, byrow = TRUE)
  )+
  theme(legend.title = element_text(size=12),
        strip.text =  element_text(size=16),
        legend.text = element_text(size=12))+
  facet_wrap(~modality_treatment)

p


###############################################################################
## Regression:
###############################################################################
fit_pp <- glm(perprotocol~gender+
                Positive_familial_history+
                prf_hypertension_awareness+
                prf_smokedtobacco+
                aneurisk+
                scale(maxdiam)+
                type+
                side+
                scale(age_Consultation),
              data = dat_pp,
              family="binomial")
summary(fit_pp)
xtable(SummaryTable(fit_pp, n=13, digits = 3))

ForestPlot_GLM(fit_pp, title = "",
               exponentiate = T, term1 = "Not Per-Protocol",
               term2 = "Per-Protocol", col_def  = "black")

