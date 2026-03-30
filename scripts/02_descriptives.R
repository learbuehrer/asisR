###############################################################################
# Descriptives for the imputation of mR_Baseline
###############################################################################
# check the propensity scores and IPW distributions of the Covariate Shift:
a <- ggplot(data=dat_imp)+
  geom_density(aes(x=pscore, fill=modality), alpha=0.5, show.legend = T) +
  theme_classic()+
  xlab("Propensity Score") +
  ylab("Density") +
  scale_fill_manual(values=darkColorPalette[c(1,5,3)])


b <- ggplot(data=dat_imp)+
  geom_density(aes(x=ipw_stab, fill=modality), alpha=0.5) +
  theme_classic()+
  xlab("Stabilized Inverse Probability Weights (IPW)") +
  ylab("Density") +
  xlim(c(0,5))+
  theme(legend.position = "bottom") +
  scale_fill_manual(values=darkColorPalette[c(1,5,3)])


p <- ggarrange(a,b, ncol=1, nrow=2, common.legend = T, legend="bottom")
p

## Check the final imputed mRS at Baseline:
densityplot(imp, ~ mR_BaselineMICE, xlab="mRS at Baseline",
            ylab="Density")
stripplot(imp, mR_BaselineMICE ~ .imp)

# redo the density plot:
complete_cases <- imp$data %>%
  filter(!is.na(mR_BaselineMICE))

ggplot()+
  # densities for each imputation
  geom_histogram(data = complete_cases,
                 aes(x = as.numeric(mR_BaselineMICE)-1, y = after_stat(density)),
                 color = "gray", fill = "gray", bins = 5) +
  geom_density(data = imp_long,
               aes(x = as.numeric(mR_BaselineMICE)-1, group = .imp, color = .imp),
               alpha = 0.3, show.legend = F, bw = 0.25) +
  labs(x = "mRS at Baseline",
       y = "Density",
       color = "Imputation") +
  theme_classic()+
  theme(element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12))

## show the mean imputed mRS at Baseline per patient across all imputations:
ggplot(imp_mean, aes(x = modality, fill = factor(mR_Baseline_mean))) +
  geom_bar() +
  theme_classic()+
  scale_fill_manual(values = colorPalette[c(4,5,3,9,10,8,11)],
                    name = "mRS at Baseline") +
  theme(legend.position = "bottom")+
  ylab("Patients [#]")+
  xlab("")+
  theme(element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12))

##############################################################################
# Descriptive Plots
##############################################################################
# simpler descriptive plots, according to PB:
a <- ggplot(dat_s)+
  geom_bar(aes(x=gender, fill=modality))+
  theme_classic()+
  scale_fill_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")+
  ylab("Patients [#]")+
  xlab("Gender")

b <- ggplot(dat_s)+
  geom_histogram(aes(x=age_Consultation, fill=modality), bins=20, position="identity", alpha=0.5)+
  geom_density(aes(x=age_Consultation, y=..count..*78/20, fill=modality, color=modality), alpha=0.2)+
  theme_classic()+
  scale_fill_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")+
  scale_color_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")+
  xlab("Age at Consultation [Years]")+
  ylab("Patients [#]")

c <- ggplot(dat_s)+
  geom_bar(aes(x=Positive_familial_history, fill=modality))+
  theme_classic()+
  scale_fill_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")+
  xlab("Positive Familial History")+
  ylab("Patients [#]")

d <- ggplot(dat_s)+
  geom_bar(aes(x=prf_hypertension_awareness, fill=modality))+
  theme_classic()+
  scale_fill_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")+
  xlab("Hypertension Awareness")+
  ylab("Patients [#]")

e <- ggplot(dat_s)+
  geom_bar(aes(x=prf_smokedtobacco, fill=modality))+
  theme_classic()+
  scale_fill_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")+
  xlab("Smoking Status")+
  ylab("Patients [#]")

f <- ggplot(dat_s)+
  geom_bar(aes(x=aneurisk, fill=modality))+
  theme_classic()+
  scale_fill_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")+
  xlab("IA Location")+
  ylab("Patients [#]")

g <- ggplot(dat_s)+
  geom_histogram(aes(x=exp(maxdiam), fill=modality), bins=20, position="identity", alpha=0.5)+
  geom_density(aes(x=exp(maxdiam), y=..count..*19.4/20, fill=modality, color=modality), alpha=0.2)+
  theme_classic()+
  scale_fill_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")+
  scale_color_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")+
  xlab("IA Size [mm]")+
  ylab("Patients [#]")

h <- ggplot(dat_s)+
  geom_bar(aes(x=type, fill=modality))+
  theme_classic()+
  scale_fill_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")+
  xlab("IA Type")+
  ylab("Patients [#]")

i <- ggplot(dat_s)+
  geom_bar(aes(x=side, fill=modality))+
  theme_classic()+
  scale_fill_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")+
  xlab("IA Side")+
  ylab("Patients [#]")

ggarrange(a,b,c,d,e,f,g,h,i,
          nrow=3, ncol=3,
          labels = c("A","B","C","D","E","F","G","H","I"),
          common.legend = T, legend = "bottom")

##############################################################################
# Distribution of the numeric variables:
###############################################################################
a <-dat_s %>%
  ggplot(aes(x = age_Consultation, fill = modality)) +
  geom_histogram(aes(y = ..count..), bins = 20, position = "stack", alpha = 0.7) +
  stat_function(
    fun = function(x) {
      binwidth <- (max(dat_s$age_Consultation, na.rm = TRUE) - min(dat_s$age_Consultation, na.rm = TRUE)) / 20
      length(dat_s$age_Consultation[!is.na(dat_s$age_Consultation)]) * binwidth * 
        dnorm(x, mean = mean(dat_s$age_Consultation, na.rm = TRUE), sd = sd(dat_s$age_Consultation, na.rm = TRUE))
    },
    color = "black", size = 1.2
  ) +
  theme_classic() +
  ylab("Patients [#]") +
  xlab("Age at Consultation [Years]") +
  scale_fill_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")

b <-dat_s %>%
  ggplot(aes(x = maxdiam, fill = modality)) +
  geom_histogram(aes(y = ..count..), bins = 20, position = "stack", alpha = 0.7) +
  stat_function(
    fun = function(x) {
      binwidth <- (max(dat_s$maxdiam, na.rm = TRUE) - min(dat_s$maxdiam, na.rm = TRUE)) / 20
      length(dat_s$maxdiam[!is.na(dat_s$maxdiam)]) * binwidth * 
        dnorm(x, mean = mean(dat_s$maxdiam, na.rm = TRUE), sd = sd(dat_s$maxdiam, na.rm = TRUE))
    },
    color = "black", size = 1.2
  ) +
  theme_classic() +
  ylab("Patients [#]") +
  xlab("IA size [log mm]") +
  scale_fill_manual(values = darkColorPalette[c(1, 6,3)], name = "Modality")


ggarrange(a,b,
          nrow=1, ncol=2, common.legend = T, legend = "bottom")


##############################################################################
# Distribution of the IA location:
###############################################################################
# aneulocation
dat_s %>%
  filter(pseudonym %in% dat_s$pseudonym) %>%
  ggplot(aes(x = aneurisk, fill = side)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  theme_classic() +
  xlab("IA Location") +
  ylab("Number of Patients") +
  scale_fill_manual(values = c("grey20", "grey80", "grey50"),
                    name = "Side") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))


aneurisk_colors <- c(
  Low = "#007F5F",
  Medium = "#B07C00",
  High = "#A04000")

side_alpha <- c(
  Left = 1,
  Midline = 0.5,
  Right = 0.2)


ggplot(dat_s, aes(x = aneurisk)) +
  ggpattern::geom_bar_pattern(
    aes(fill = aneurisk, alpha = side),
    position = position_stack(reverse = TRUE),
    stat = "count",
    color = "black",
    pattern = "none") +
  theme_classic() +
  xlab("IA Location") +
  ylab("Number of Patients") +
  scale_fill_manual(values = aneurisk_colors, name = "IA Location") +
  scale_alpha_manual(values = side_alpha, name = "Side") +
  guides(
    fill = guide_legend(nrow = 1),
    alpha = guide_legend(nrow = 1),
    pattern = guide_legend(
      nrow = 1,
      override.aes = list(
        pattern = c("none", "circle", "none"),
        fill = c("black", "black", "black"),
        alpha = c(1, 0.3, 0.2),
        pattern_fill = "black"))) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y = element_text(size=14),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    axis.title = element_text(size=14))

##############################################################################
## ASIS-cohort - allulium plot - mR_BaselineMICEmean, Discharge and 1-year FU:
##############################################################################
plot_dat <- dat_s %>%
  select(mR_BaselinePredMICEmean, mR_Discharge, mR_1, modality) %>%
  group_by(mR_BaselinePredMICEmean, mR_Discharge, mR_1,modality) %>%
  count() %>%
  na.omit() %>%
  mutate(mR_BaselinePredMICEmean = factor(mR_BaselinePredMICEmean, levels=(c("0","1","2","3","4")))) %>%
  mutate(mR_Discharge = factor(mR_Discharge, levels=(levels(mR_Discharge)))) %>%
  mutate(mR_1 = factor(mR_1, levels=(levels(mR_1)))) %>%
  group_by(modality) %>%
  mutate(ngroup = sum(n)) %>%
  mutate(prop=n/ngroup*100)
fn = sum(plot_dat$n)
nObservat = sum(plot_dat$n[plot_dat$modality=="Observation"])
nEndo = sum(plot_dat$n[plot_dat$modality=="Endovascular"])
nMicro = sum(plot_dat$n[f$modality=="Microneurosurgery"])

p<- plot_dat %>%
  ggplot(aes(y=prop,
             axis1= mR_BaselinePredMICEmean,
             axis2=mR_Discharge,
             axis3=mR_1))+
  geom_alluvium(aes(fill=mR_1), width = 1/12)+
  theme_classic()+
  scale_y_reverse(breaks = seq(100, 0, by = -20), labels = paste0(seq(0, 100, by = 20)))+
  scale_fill_manual(values=colorPalette[(c(4,5,3,9,10,8,11))], name ="mRS")+
  scale_x_continuous(breaks = c(1,2,3),
                     labels=c("1"="Baseline",
                              "2" = "Discharge",
                              "3"= "1-year FU"))+
  guides(fill=guide_legend(nrow=1, reverse = F, byrow=T))+
  theme(legend.position = "bottom")+
  facet_grid(~modality)+
  labs(caption = paste0("n = ", fn))+
  geom_stratum(width = 1/12, alpha=1,
               fill = c(colorPalette[c(4,5,3,9,10)],
                        colorPalette[c(4,5,3,9,10)],
                        colorPalette[c(4,5,3,9,11)],
                        # EndoVasc
                        colorPalette[c(4,5,3)], 
                        colorPalette[c(4,5,3,9,11)],
                        colorPalette[c(4,5,3,9,8,11)],
                        # MicroNeuroSurg
                        colorPalette[c(4,5,3)], 
                        colorPalette[c(4,5,3,9)],
                        colorPalette[c(4,5,3)])) +
  ylab("Patients [%]")+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.y = element_text(size=14),
        plot.caption = element_text(size=14),
        strip.text = element_text(size=21))
p


##############################################################################
## Check of outliers:
##############################################################################
ggplot(dat_s)+
  geom_boxplot(aes(x=modality, y=exp(maxdiam), fill=modality))+
  theme_classic()+
  scale_fill_manual(values=colorPalette[c(1,3,9)],
                    name="Modality")+
  ylab("IA Size [mm]")+
  xlab("Modality")+
  theme(axis.text.x = element_text(angle=45, hjust=1))


ggplot(dat_s)+
  geom_boxplot(aes(x=aneurisk, y=(maxdiam), fill=modality))+
  theme_classic()+
  scale_fill_manual(values=colorPalette[c(1,3,9)],
                    name="Modality")+
  ylab("Log IA Size [log mm]")+
  xlab("Modality")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "bottom")


rm(a,b,c,d,e,f,g,h,i,plot_dat,p,
   complete_cases, nEndo, nMicro, nObservat, fn)

