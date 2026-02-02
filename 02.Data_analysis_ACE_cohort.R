library(openxlsx)
library(data.table)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(gridExtra)
library(tidyr)
library(ggpubr)
library(VennDiagram)
library(ggplot2)
library(gridExtra)
library(forestplot)
library(glmnet)
library(survival)
library(survival)
library("survminer")
library(RColorBrewer)
library(meta)

## Load input data
load("omics_clinical_data/Data_for_analysis_V1.RData")
load("omics_clinical_data/Metadata_cohort.RData")

## Demographics full cohort

harpone.ids <- read.xlsx("omics_clinical_data/HARPONE_REDUX_integrated_data_20240416.xlsx")$Codigo_ACE_LCR
csf.for.demog <- csf[which(csf$Codigo_ACE_LCR%in%harpone.ids),]

csf.for.demog$Baseline_Diagnostic<-csf.for.demog$Diagnostic_Sindromic_PL %>% 
  recode_factor(., "1"="Control_QSM", "2"="MCI", "3"="Dementia", "4"="Control_QSM") %>%
  replace_na("Control_QSM") ## Creo que estos NAs son todos FACEHBI. Pero recomprobar en su momento.
csf.for.demog$e4_carrier <- ifelse(is.na(csf.for.demog$APOE_merge), NA, grepl("e4", csf.for.demog$APOE_merge))

csf.for.demog$AT_Status <- paste0(ifelse(csf.for.demog$A_1pos_0neg, "A+", "A-"),
                                  ifelse(csf.for.demog$T_1pos_0neg|csf.for.demog$N_1pos_0neg, "T+", "T-")) %>% ifelse(grepl("NA",.),NA,.)

by.group.demog<-csf.for.demog %>% group_by_("Baseline_Diagnostic") %>% 
  summarise(Sample_Size = n(),
            A_T_Minus = sum(AT_Status == "A-T-", na.rm = TRUE),
            A_T_Plus = sum(AT_Status == "A-T+", na.rm = TRUE),
            A_Plus_T_Minus = sum(AT_Status == "A+T-", na.rm = TRUE),
            A_Plus_T_Plus = sum(AT_Status == "A+T+", na.rm = TRUE),
            Mean_Age = mean(Age_LP, na.rm = TRUE),
            SD_Age = sd(Age_LP, na.rm = TRUE),
            Percent_Female = mean(Sex_1m_2f == 2, na.rm = TRUE) * 100,
            Percent_APOE_Carriers = mean(e4_carrier, na.rm = TRUE) * 100,
            Mean_MMSE = mean(as.numeric(MMSE_PL), na.rm=T),
            SD_MMSE = sd(as.numeric(MMSE_PL), na.rm=T),
            Mean_Education=mean(as.numeric(Years_School), na.rm=T),
            SD_Education=sd(as.numeric(Years_School), na.rm=T))

all.demog<-csf.for.demog %>% 
  summarise(Baseline_Diagnostic = "All",
            A_T_Minus = sum(AT_Status == "A-T-", na.rm = TRUE),
            A_T_Plus = sum(AT_Status == "A-T+", na.rm = TRUE),
            A_Plus_T_Minus = sum(AT_Status == "A+T-", na.rm = TRUE),
            A_Plus_T_Plus = sum(AT_Status == "A+T+", na.rm = TRUE),
            Sample_Size = n(),
            Mean_Age = mean(Age_LP, na.rm = TRUE),
            SD_Age = sd(Age_LP, na.rm = TRUE),
            Percent_Female = mean(Sex_1m_2f == 2, na.rm = TRUE) * 100,
            Percent_APOE_Carriers = mean(e4_carrier, na.rm = TRUE) * 100,
            Mean_MMSE = mean(as.numeric(MMSE_PL), na.rm=T),
            SD_MMSE = sd(as.numeric(MMSE_PL), na.rm=T),
            Mean_Education=mean(as.numeric(Years_School), na.rm=T),
            SD_Education=sd(as.numeric(Years_School), na.rm=T)) 

demog<-rbind(by.group.demog,all.demog)
demog[,-1] <- sapply(demog[,-1], round, digits=1)


## Supp. table with OPCML tertiles (Table 3)

csf.for.demog <- data.for.analysis.v1
csf.for.demog$Baseline_Diagnostic<-csf.for.demog$Diagnostic_Sindromic_PL %>% 
  recode_factor(., "1"="Control_QSM", "2"="MCI", "3"="Dementia", "4"="Control_QSM") %>%
  replace_na("Control_QSM") ## Creo que estos NAs son todos FACEHBI. Pero recomprobar en su momento.
csf.for.demog$e4_carrier <- ifelse(is.na(csf.for.demog$APOE_merge), NA, grepl("e4", csf.for.demog$APOE_merge))

csf.for.demog$AT_Status <- paste0(ifelse(csf.for.demog$A_1pos_0neg, "A+", "A-"),
                                  ifelse(csf.for.demog$T_1pos_0neg|csf.for.demog$N_1pos_0neg, "T+", "T-")) %>% ifelse(grepl("NA",.),NA,.)


csf.for.demog <- csf.for.demog %>% mutate(Tertiles_OBCAM=cut(OBCAM_seq.15622.13,
                                                              breaks = quantile(OBCAM_seq.15622.13, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                                                              include.lowest = TRUE, 
                                                              labels = c("Low", "Medium", "High"))) %>% filter(!is.na(Tertiles_OBCAM))

supptable <- lapply(c("Control_QSM", "MCI", "Dementia"), function(group){
  
  test.data <- csf.for.demog[csf.for.demog$Baseline_Diagnostic == group, ] 
  
  # summarized table
  demog.x <- test.data %>%
    group_by(Tertiles_OBCAM) %>%
    summarise(
      Sample_Size = n(),
      A_T_Minus = sum(AT_Status == "A-T-", na.rm = TRUE),
      A_T_Plus = sum(AT_Status == "A-T+", na.rm = TRUE),
      A_Plus_T_Minus = sum(AT_Status == "A+T-", na.rm = TRUE),
      A_Plus_T_Plus = sum(AT_Status == "A+T+", na.rm = TRUE),
      Mean_Abeta = median(Abeta_42_LCR, na.rm = TRUE),    ## es median en realidad
      SD_Abeta = sd(Abeta_42_LCR, na.rm = TRUE),
      Mean_ptau = median(P_tau_LCR, na.rm = TRUE),        ## es median en realidad
      SD_ptau = sd(P_tau_LCR, na.rm = TRUE),
      Mean_Age = mean(Age_absolute, na.rm = TRUE),
      SD_Age = sd(Age_absolute, na.rm = TRUE),
      Percent_Female = mean(Sex_1m_2f == 2, na.rm = TRUE) * 100,
      Percent_APOE_Carriers = mean(e4_carrier, na.rm = TRUE) * 100,
      Mean_MMSE = mean(as.numeric(MMSE_PL), na.rm=TRUE),
      SD_MMSE = sd(as.numeric(MMSE_PL), na.rm=TRUE),
      Mean_CDR = mean(as.numeric(CDR_PL), na.rm=TRUE),
      SD_CDR = sd(as.numeric(CDR_PL), na.rm=TRUE),
      Mean_Education=mean(as.numeric(Years_School), na.rm=TRUE),
      SD_Education=sd(as.numeric(Years_School), na.rm=TRUE)
    )
  
  # statistical tests:
  p_abeta <- kruskal.test(Abeta_42_LCR ~ Tertiles_OBCAM, data=test.data)$p.value
  p_ptau <- kruskal.test(P_tau_LCR ~ Tertiles_OBCAM, data=test.data)$p.value
  p_age <- kruskal.test(Age_absolute ~ Tertiles_OBCAM, data=test.data)$p.value
  p_sex <- chisq.test(table(test.data$Tertiles_OBCAM, test.data$Sex_1m_2f))$p.value
  p_apoe <- chisq.test(table(test.data$Tertiles_OBCAM, test.data$e4_carrier))$p.value
  p_mmse <- kruskal.test(as.numeric(MMSE_PL) ~ Tertiles_OBCAM, data=test.data)$p.value
  p_edu <- kruskal.test(as.numeric(Years_School) ~ Tertiles_OBCAM, data=test.data)$p.value
  
  # for AT_Status distribution
  at_table <- table(test.data$Tertiles_OBCAM, test.data$AT_Status)
  p_atstatus <- suppressWarnings(chisq.test(at_table)$p.value) %>% format(., digits = 3, scientific = T)
  
  # Chisq in each AT level
  p_atgroups <- lapply(1:ncol(at_table), function(i){suppressWarnings(chisq.test(at_table[,i])$p.value)}) %>% 
    unlist() %>% format(., digits = 3, scientific = T)
  
  p_values <- c(
    p_atstatus,
    p_atgroups,  
    format(p_abeta, digits = 3, scientific = T),
    format(p_ptau, digits = 3, scientific = T),
    round(p_age, 3),
    round(p_sex, 3),
    round(p_apoe, 3),
    round(p_mmse, 3),
    round(p_edu, 3)
  )
  
  demog.x <- t(demog.x[,-1]) %>% round(., 1)
  colnames(demog.x) <- c("Low", "Medium", "High")
  
  meansd <- lapply(c("Abeta", "ptau", "Age", "MMSE", "Education"), function(x){
    paste0(demog.x[paste0("Mean_", x), ], " ;", demog.x[paste0("SD_", x), ])
  }) %>% do.call("rbind", .)
  row.names(meansd) <- c("Abeta", "ptau", "Age", "MMSE", "Education")
  
  demog.x <- rbind(demog.x, meansd)
  
  # final ordering
  demog.x <- demog.x[c("Sample_Size", "A_T_Minus", "A_T_Plus", "A_Plus_T_Minus", "A_Plus_T_Plus", "Abeta", "ptau",
                       "Age", "Percent_Female", "Percent_APOE_Carriers", "MMSE", "Education"), ]
  
  # add p-values as a column
  demog.x <- cbind(demog.x, p_value = p_values)
  
  demog.x
}) %>% do.call("cbind", .) %>% as.data.frame() 

write.xlsx(supptable, "Supp_Table_OCBAMtertiles.xlsx")

## Barplot tertile distribution by tertile and AT group (Fig. 5b)

plotdata <- csf.for.demog %>%
  filter(!is.na(Tertiles_OBCAM)) %>%
  group_by(Baseline_Diagnostic, AT_Status, Tertiles_OBCAM) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Baseline_Diagnostic, AT_Status) %>%
  mutate(Percent = N / sum(N) * 100)

# plot
ggplot(plotdata, aes(x = AT_Status, y = Percent, fill = Tertiles_OBCAM)) +
  geom_bar(stat = "identity", position = position_stack()) +
  facet_wrap(~Baseline_Diagnostic) +
  labs(
    x = "AT Status",
    y = "Percent of Tertiles",
    fill = "GAGE2A Tertile",
    title = "Distribution of OPCML Tertiles by AT Status and Baseline Diagnostic Group"
  ) +
  theme_minimal(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave("Barplot_OBCAM_tertiles_by_AT_and_Strata.png", width=2200, height=1000, units="px")



## Plots proteomics vs lipidomcis intensity (Fig. 1e):
ggplot(data.for.analysis.v1, aes(Lipometrix_intensity, Somalogic_intensity)) + geom_point(alpha=0.4) + geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)
lm_model <- lm(Somalogic_intensity ~ Lipometrix_intensity, data = data.for.analysis.v1)


## Boxplot with Wilcoxon test by syndromic status:

my_comparisons<-list(c(1,2), c(1,3), c(2,3))

ggboxplot(data.for.analysis.v1 %>% 
            mutate(Syndromic_Status=recode_factor(data.for.analysis.v1$Syndromic_Status, "1"="CU", "2"="MCI", "3"="Dementia") %>%
                     replace_na("CU")), 
          x = "Syndromic_Status", y = "Somalogic_intensity",
          color = "Syndromic_Status", palette = "jco",
          add = "jitter", add.params = list(alpha = 0.3)) + theme(axis.text.x=element_blank()) + stat_compare_means(comparisons = my_comparisons, step.increase=0.20) + theme(legend.position="none") + scale_y_continuous(limits=c(-4,6))
# ggsave("boxplot_syndromic_Wilcox_Somascan_int.png", width=900, height=850, units="px") 

ggboxplot(data.for.analysis.v1 %>% 
            mutate(Syndromic_Status=recode_factor(data.for.analysis.v1$Syndromic_Status, "1"="CU", "2"="MCI", "3"="Dementia") %>%
                     replace_na("CU")), 
          x = "Syndromic_Status", y = "Lipometrix_intensity",
          color = "Syndromic_Status", palette = "jco",
          add = "jitter", add.params = list(alpha = 0.3)) + theme(axis.text.x=element_blank()) + stat_compare_means(comparisons = my_comparisons, step.increase=0.20) + theme(legend.position="none") + scale_y_continuous(limits=c(-4,6))
# ggsave("boxplot_syndromic_Wilcox_Lipometrix_int.png", width=900, height=850, units="px")


ggboxplot(data.for.analysis.v1 %>% 
            mutate(Syndromic_Status=recode_factor(data.for.analysis.v1$Syndromic_Status, "1"="CU", "2"="MCI", "3"="Dementia") %>%
                     replace_na("CU")), 
          x = "Syndromic_Status", y = "Somalogic_PC.1",
          color = "Syndromic_Status", palette = "jco",
          add = "jitter", add.params = list(alpha = 0.3)) + theme(axis.text.x=element_blank()) + stat_compare_means(comparisons = my_comparisons, step.increase=0.20) + theme(legend.position="none") + scale_y_continuous(limits=c(-4,6))
# ggsave("boxplot_syndromic_Wilcox_SomascanPC1.png", width=900, height=850, units="px")

ggboxplot(data.for.analysis.v1 %>% 
            mutate(Syndromic_Status=recode_factor(data.for.analysis.v1$Syndromic_Status, "1"="CU", "2"="MCI", "3"="Dementia") %>%
                     replace_na("CU")), 
          x = "Syndromic_Status", y = "Lipometrix_PC.1",
          color = "Syndromic_Status", palette = "jco",
          add = "jitter", add.params = list(alpha = 0.3)) + theme(axis.text.x=element_blank()) + stat_compare_means(comparisons = my_comparisons, step.increase=0.20) + theme(legend.position="none") + scale_y_continuous(limits=c(-4,6))
# ggsave("boxplot_syndromic_Wilcox_LipometrixPC1.png", width=900, height=850, units="px")


ggboxplot(data.for.analysis.v1 %>% 
            mutate(Syndromic_Status=recode_factor(data.for.analysis.v1$Syndromic_Status, "1"="CU", "2"="MCI", "3"="Dementia") %>%
                     replace_na("CU")), 
          x = "Syndromic_Status", y = "Somalogic_PC.2",
          color = "Syndromic_Status", palette = "jco",
          add = "jitter", add.params = list(alpha = 0.3)) + theme(axis.text.x=element_blank()) + stat_compare_means(comparisons = my_comparisons, step.increase=0.20) + theme(legend.position="none") + scale_y_continuous(limits=c(-4,6))
# ggsave("boxplot_syndromic_Wilcox_SomascanPC2.png", width=900, height=850, units="px")

ggboxplot(data.for.analysis.v1 %>% 
            mutate(Syndromic_Status=recode_factor(data.for.analysis.v1$Syndromic_Status, "1"="CU", "2"="MCI", "3"="Dementia") %>%
                     replace_na("CU")), 
          x = "Syndromic_Status", y = "Lipometrix_PC.2",
          color = "Syndromic_Status", palette = "jco",
          add = "jitter", add.params = list(alpha = 0.3)) + theme(axis.text.x=element_blank()) + stat_compare_means(comparisons = my_comparisons, step.increase=0.20) + theme(legend.position="none") + scale_y_continuous(limits=c(-4,6))
# ggsave("boxplot_syndromic_Wilcox_LipometrixPC2.png", width=900, height=850, units="px")


### PC variance barplot (Fig. 2a)

lip.variance <- read.table("data/Variance_PCs_lipidomics.txt")[,1]
prot.variance <- read.table("data/Variance_PCs_proteomics.txt")[,1]

n.pcs=4
pve.pcs <- data.frame(Dataset=c(rep("Lipometrix", n.pcs),rep("SOMAscan", n.pcs)),
                      PC=rep(paste0("PC", sprintf("%02d",1:n.pcs)), 2),
                      Pve=c(lip.variance[1:n.pcs]/sum(lip.variance), prot.variance[1:n.pcs]/sum(prot.variance)))

ggplot(pve.pcs, aes(PC, 100*Pve, fill=Dataset)) + geom_bar(position="dodge", stat="identity") + 
  ylab("Variance explained (%)") + xlab("Principal component")


### BOXPLOTS tertiles vs AD biomarkers (Fig. 5c)

### PC1
## OPCML tertiles

ggplot(data.for.analysis.v1 %>% 
         mutate(Syndromic_Status = recode_factor(data.for.analysis.v1$Syndromic_Status, 
                                                 "1"="CU", "2"="MCI", "3"="Dementia") %>%
                  replace_na("CU")) %>% 
         mutate(Tertiles_OPCML = cut(OBCAM_seq.15622.13, 
                                     breaks = quantile(OBCAM_seq.15622.13, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                                     include.lowest = TRUE, 
                                     labels = c("Low", "Medium", "High"))) %>% 
         filter(!is.na(Tertiles_OPCML)), 
       aes(x = Syndromic_Status, y = P_tau_LCR, color = Tertiles_OPCML)) + 
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75), 
              size = 1.0, alpha = 0.3) +
  geom_hline(yintercept = 54, linetype = "dashed", color = "red", size = 0.4)
# ggsave("pTau_distrib_OPCML_tertiles.png", width=1600, height=1200, units="px")


ggplot(data.for.analysis.v1 %>% 
         mutate(Syndromic_Status = recode_factor(data.for.analysis.v1$Syndromic_Status, 
                                                 "1"="CU", "2"="MCI", "3"="Dementia") %>%
                  replace_na("CU")) %>% 
         mutate(Tertiles_OPCML = cut(OBCAM_seq.15622.13, 
                                     breaks = quantile(OBCAM_seq.15622.13, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                                     include.lowest = TRUE, 
                                     labels = c("Low", "Medium", "High"))) %>% 
         filter(!is.na(Tertiles_OPCML)), 
       aes(x = Syndromic_Status, y = Abeta_42_LCR, color = Tertiles_OPCML)) + 
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75), 
              size = 1.0, alpha = 0.3) +
  geom_hline(yintercept = 796, linetype = "dashed", color = "red", size = 0.4)
# ggsave("Abeta42_distrib_OPCML_tertiles.png", width=1600, height=1200, units="px")


### Boxplots clinical history (Supp. Fig. 16) #### 

plot_boxplots_with_jitter <- function(data, group_var, y_vars) {
  
  plot_list <- lapply(y_vars, function(y) {
    ggplot(data %>% subset(!is.na(data[,group_var])), aes(x = as.factor(.data[[group_var]]), y = .data[[y]], color = as.factor(.data[[group_var]]))) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.4, shape = 16) +
      stat_compare_means(method = "t.test", label = "p.format") +
      ggtitle(y) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      )
  })
  grid.arrange(grobs = plot_list, ncol = 2)
}

y_vars <- c("Lipometrix_PC.1", "Lipometrix_PC.2", "Somalogic_PC.1", "Somalogic_PC.2")
x_vars <- c("Artrosis", "AVC", "Cardiopatia", "Dislipidemia", "Diabetes",
            "Hipertension", "EPOC", "Insuficiencia_renal", "Depresion", "Esquizofrenia", "Epilepsia", "Alcoholismo", "Tabaquismo")

## Make and save plots:

lapply(x_vars, function(x){ 
  a <- plot_boxplots_with_jitter(data.for.analysis.v1, x, y_vars)
  ggsave(plot = a,paste0("Boxplot_antecedentes.", x, ".png"), width=1800, height=1800, units="px")
})


## Check effect of stroke and hypertension PRSs (Supp. Fig. 17)

y_vars=c("PRS_Keaton_BloodPressure_DBP_2024", 
         "PRS_Keaton_BloodPressure_PP_2024", "PRS_Keaton_BloodPressure_SBP_2024")

x_vars <- c("Artrosis", "AVC", "Cardiopatia", "Dislipidemia", "Diabetes",
            "Hipertension", "EPOC", "Insuficiencia_renal", "Depresion", "Esquizofrenia", "Epilepsia", "Alcoholismo", "Tabaquismo")

lapply(x_vars, function(x){ 
  plot_boxplots_with_jitter(data.for.analysis.v1, x, y_vars)
})


## Not so with stroke PRSs and cerebrovascular accident
y_vars=c("PRS_Mishra_Stroke_2022_AIS", "PRS_Mishra_Stroke_2022_AS", "PRS_Mishra_Stroke_2022_CES", 
         "PRS_Mishra_Stroke_2022_LAS", "PRS_Mishra_Stroke_2022_SVS")

lapply(x_vars, function(x){ 
  plot_boxplots_with_jitter(data.for.analysis.v1, x, y_vars)
})


plot_and_save_prs_boxplots <- function(data, prs_vars, y_var = "Somalogic_PC.2", output_dir = ".") {
  for (prs in prs_vars) {
    # Create tertiles
    data_with_tertiles <- data %>%
      mutate(Tertiles_PRS = cut(.data[[prs]], 
                                breaks = quantile(.data[[prs]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                                include.lowest = TRUE, 
                                labels = c("Low", "Medium", "High"))) %>%
      filter(!is.na(Tertiles_PRS))
    
    # Generate the plot with t-tests
    p <- ggplot(data_with_tertiles, aes(x = Tertiles_PRS, y = .data[[y_var]], color = Tertiles_PRS)) +
      geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75), 
                  size = 1.0, alpha = 0.3) +
      stat_compare_means(method = "t.test", 
                         comparisons = list(c("Low", "Medium"), c("Medium", "High"), c("Low", "High")),
                         label = "p.format") +
      ggtitle(prs) + theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
    
    # Save the plot
    filename <- file.path(output_dir, paste0(prs, "_boxplot.png"))
    ggsave(filename, plot = p, width = 4, height = 4, dpi = 300)
  }
}

plot_and_save_prs_boxplots(data.for.analysis.v1, grep("^PRS", names(data.for.analysis.v1), value=T))

# No parece que el PRS de hipertension se asocie mucho al PC2
ggplot(data.for.analysis.v1 %>% 
         mutate(Tertiles_PRS = cut(PRS_Keaton_BloodPressure_DBP_2024, 
                                   breaks = quantile(PRS_Keaton_BloodPressure_DBP_2024, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                                   include.lowest = TRUE, 
                                   labels = c("Low", "Medium", "High"))) %>% 
         filter(!is.na(Tertiles_PRS)), 
       aes(x = Tertiles_PRS, y = Somalogic_PC.2, color = Tertiles_PRS)) + 
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75), 
              size = 1.0, alpha = 0.3) 


## Test hypertension vs PRSs
summary(glm(Hipertension ~ PRS_Keaton_BloodPressure_DBP_2024 + Sex_1m_2f + Age_LP, data=data.for.analysis.v1, family=binomial))
summary(glm(Hipertension ~ PRS_Keaton_BloodPressure_PP_2024 + Sex_1m_2f + Age_LP, data=data.for.analysis.v1, family=binomial))
summary(glm(Hipertension ~ PRS_Keaton_BloodPressure_SBP_2024 + Sex_1m_2f + Age_LP, data=data.for.analysis.v1, family=binomial))

## El PRS no se asocia al PC2
summary(lm(Somalogic_PC.2 ~ PRS_Keaton_BloodPressure_SBP_2024 + Sex_1m_2f + Age_LP, data=data.for.analysis.v1))
summary(lm(Somalogic_PC.2 ~ PRS_Keaton_BloodPressure_DBP_2024 + Sex_1m_2f + Age_LP, data=data.for.analysis.v1))
summary(lm(Somalogic_PC.2 ~ PRS_Keaton_BloodPressure_PP_2024 + Sex_1m_2f + Age_LP, data=data.for.analysis.v1))

## El PRS no se asocia al PC2, pero la hipertension si
summary(lm(Somalogic_PC.2 ~ Hipertension + PRS_Keaton_BloodPressure_SBP_2024 + Sex_1m_2f + Age_LP, data=data.for.analysis.v1))
summary(lm(Somalogic_PC.2 ~ Hipertension + PRS_Keaton_BloodPressure_DBP_2024 + Sex_1m_2f + Age_LP, data=data.for.analysis.v1))
summary(lm(Somalogic_PC.2 ~ Hipertension + PRS_Keaton_BloodPressure_PP_2024 + Sex_1m_2f + Age_LP, data=data.for.analysis.v1))

## Si meto la total protein CSF se lleva casi todo, pero deja algo de hipertension
summary(lm(Somalogic_PC.2 ~ Hipertension + PRS_Keaton_BloodPressure_SBP_2024 + Sex_1m_2f + Age_LP + Total_protein_CSF, data=data.for.analysis.v1))
summary(lm(Somalogic_PC.2 ~ Hipertension + PRS_Keaton_BloodPressure_DBP_2024 + Sex_1m_2f + Age_LP + Total_protein_CSF,, data=data.for.analysis.v1))
summary(lm(Somalogic_PC.2 ~ Hipertension + PRS_Keaton_BloodPressure_PP_2024 + Sex_1m_2f + Age_LP + Total_protein_CSF,, data=data.for.analysis.v1))

## RBC count and PCs
data.for.analysis.v1$RBC_count_CSF <- gsub("[^0-9.]", "", data.for.analysis.v1$RECUENTO.CELULAR.Hematíes..mm3.Contaje.celular.N) %>% as.numeric() 

data.for.analysis.v1$RBC_count_CSF_factor <- cut(data.for.analysis.v1$RBC_count_CSF, breaks = c(-Inf, 10.1, 100, Inf), 
                                                 labels = c("<10", "10-100", ">100"), right = FALSE) 

ggplot(data.for.analysis.v1 %>% subset(!is.na(RBC_count_CSF_factor)), aes(x=RBC_count_CSF_factor, y=Somalogic_PC.2, color=RBC_count_CSF_factor)) + 
  geom_boxplot() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75)) +
  stat_compare_means(method = "t.test", label = "p.format", comparisons = list(c(1,2), c(1,3), c(2,3)))

### ASSOCIATIONS AND FOREST PLOTS (Fig. 2c)

multiassociation <- function(input.data,analytes,independent.vars, make.volcano=F, volcano.var=NULL, tag="Volcano plot", legend.pos="bottomleft", output.type="full"){
  
  # Matriz de datos
  input.data <- input.data
  
  # Especificar analitos a mano
  analytes <- analytes
  
  # Vector de variables independientes con covariables 
  independent.vars <- independent.vars
  
  # Adaptar variables categoricas
  input.data$Sex_1m_2f <- as.factor(input.data$Sex_1m_2f)
  
  report.assoc <- data.frame()
  
  for(i in 1:length(analytes)){
    
    temp.model <- input.data[,c(analytes[i],independent.vars)]
    temp.model <- temp.model[complete.cases(temp.model),]
    mymodel <- as.formula(paste(analytes[i], paste(independent.vars, collapse=" + "), sep=" ~ "))
    temp.assoc <- lm(mymodel, data=temp.model)
    indep.assocs <- complete.cases(confint(temp.assoc))
    n.assoc <- rep(nrow(temp.model), length(indep.assocs))
    marker <- rep(analytes[i], length(indep.assocs))
    indep.var <- c("Intercept", independent.vars)[indep.assocs] # me sobra 1 aqui no se por que 
    r1 <- cbind(marker,indep.var,summary(temp.assoc)$coefficients,confint(temp.assoc, level = 0.95)[complete.cases(confint(temp.assoc)),], n.assoc)
    row.names(r1) <- paste(analytes[i], row.names((r1)), sep=".")
    
    report.assoc <- rbind(report.assoc, r1)
  }
  
  
  report.assoc <- report.assoc[grep("Intercept", row.names(report.assoc), invert=T),]
  names(report.assoc)[names(report.assoc) == "Pr(>|t|)"] <- "pval"
  names(report.assoc)[names(report.assoc) == "n.assoc"] <- "N"
  
  # Se me guarda ahora como factor, convierto estas cols a numericas:
  report.assoc[,c("Estimate", "Std. Error", "t value", "pval", "2.5 %", "97.5 %")] <- sapply(report.assoc[,c("Estimate", "Std. Error", "t value", "pval", "2.5 %", "97.5 %")], function(x) {as.numeric(as.character(x))})
  
  if(make.volcano){
    
    report.volcano <- report.assoc[which(report.assoc$indep.var==volcano.var),]
    
    report.volcano$Padj_FDR <- p.adjust(report.volcano$pval, method="fdr", n=length(analytes))
    report.volcano$logPnegativo <- -log10(report.volcano$pval)
    
    plot(report.volcano$Estimate, report.volcano$logPnegativo, xlab="Estimate", ylab="-log10(pval)", main=tag)
    points(report.volcano$Estimate[which(report.volcano$pval<0.05)], report.volcano$logPnegativo[which(report.volcano$pval<0.05)], col="blue")
    points(report.volcano$Estimate[which(report.volcano$Padj<0.05)], report.volcano$logPnegativo[which(report.volcano$Padj<0.05)], col="red")
    abline(v=0, lty=2)
    report.volcano$chisq <- qchisq(1-report.volcano$pval,1)
    
    legend(legend.pos, c("P>=0.05", "P<0.05", "Padj<0.05"), col=c("black", "blue", "red"), pch=c(1,1,1),cex = 0.8)
    
    upregulated <- sum(report.volcano$Padj_FDR<0.05&report.volcano$Estimate>0)
    downregulated <- sum(report.volcano$Padj_FDR<0.05&report.volcano$Estimate<0)
    lambda.gc <- median(qchisq(1-report.volcano$pval,1))/qchisq(0.5,1) # guardar lambda
    print(paste("Lambda Inflation Factor =", lambda.gc))
    print(paste(upregulated, "Upregulated,", downregulated, "Downregulated"))
    btest.result <- binom.test(x = c(upregulated,downregulated), alternative = "two.sided", conf.level = 0.95)
    print(paste("Binomial test p =",format(btest.result$p.value, scientific = T, digits=2)))
    gaston::qqplot.pvalues(report.volcano$pval, col.abline = "red",CB = F)
    
    report.volcano <- report.volcano[order(report.volcano$pval),]
    
  }else{
    report.volcano <- report.assoc[which(report.assoc$indep.var==volcano.var),]
    report.volcano$Padj_FDR <- p.adjust(report.volcano$pval, method="fdr", n=length(analytes))
    report.volcano$logPnegativo <- -log10(report.volcano$pval)
    report.volcano <- report.volcano[order(report.volcano$pval),]
  }
  
  if(output.type=="full"){return(report.assoc)}
  if(output.type=="summary"){return(report.volcano)}
}

forest.pcs <- function(input.data,main.variables, xlim=c(-2,2), tag.plot, fig.width=1500, fig.height=2000){
  
  input.data <- input.data
  main.variables <- main.variables
  xlim <- xlim
  tag.plot <- tag.plot
  
  for.forest.v1 <- input.data[input.data$marker%in%main.variables,]
  tags <- unique(for.forest.v1$marker)
  for.forest.v1 <- for.forest.v1[!for.forest.v1$indep.var=="Tecnica_Empleada_LCR",]
  
  for.forest.v1$Summary <- NA
  for.forest.v1$Significance <- NA
  for.forest.v1$Significance[for.forest.v1$pval < 0.05] <- "*"
  for.forest.v1$Significance[for.forest.v1$pval < 0.01] <- "**"
  for.forest.v1$Significance[for.forest.v1$pval < 0.001] <- "***"
  
  for.forest.v1 <- for.forest.v1[,c("marker", "indep.var", "Estimate", "Significance", "Summary", "Estimate", "2.5 %", "97.5 %", "pval")]
  
  for.forest.v1$Estimate <- paste(round(for.forest.v1$Estimate, digits=2), "; p=", format(for.forest.v1$pval, scientific = T, digits = 2), sep="")
  for.forest.v1$pval <- NULL
  
  names(for.forest.v1) <- c("marker", "indep.var", "Beta", "Significance", "Summary", "Estimate", "CI2.5", "CI97.5")
  for.forest.v1 <- rbind(c("PC", "Indep_variable", "Beta; p", "Significance", TRUE, NA, NA, NA),for.forest.v1)
  
  ## Forest plots:
  for.forest.2 <- for.forest.v1[1,]
  for(i in tags){
    
    temp.edit <- for.forest.v1[for.forest.v1$marker==i,]
    temp.edit$marker <- NA
    
    for.forest.2 <- rbind(for.forest.2,c(i, NA, NA, NA,NA, NA, NA, NA, NA))
    for.forest.2<-rbind(for.forest.2,temp.edit)
  }
  
  for.forest.2$Summary <- as.logical(for.forest.2$Summary)
  for.forest.2[,c("Estimate", "CI2.5", "CI97.5")] <- sapply(for.forest.2[,c("Estimate", "CI2.5", "CI97.5")], as.numeric)
  
  png(paste("Forest_", tag.plot, ".png", sep=""), fig.width, fig.height)
  # trellis.device(device="windows", height = 25, width = 40, color=TRUE)
  plot <- for.forest.2 %>% forestplot(labeltext=c(marker,indep.var,Significance,Beta), is.summary=Summary, 
                                      mean=Estimate, lower=CI2.5, upper=CI97.5, clip = xlim,
                                      # hrzl_lines = list("8.5" = gpar(lwd=460, lineend="butt",  col="#99999922"),
                                      #                   "32.5" = gpar(lwd=460, lineend="butt",  col="#99999922"),
                                      #                   "2" = gpar(lty = 1)),
                                      vertices=T, col = fpColors(box = "black",
                                                                 line = "black", 
                                                                 summary = "royalblue"),graph.pos=3,
                                      txt_gp = fpTxtGp(cex=2,ticks=gpar(cex=2)),xlog=F)
  print(plot)
  dev.off()
  
}

## Run models and forest
covariates.all <- c("Age_LP", "Sex_1m_2f", "pTau", "Abeta42", "Tecnica_Empleada_LCR", "Total_protein_CSF", "QAlb", "BMI", "Sample_longevity_Years",  "e4_alleles", "e2_alleles", "PRS_Bellenguez_etal_2020_dosage")
missingness <- sapply(data.for.analysis.v1[,covariates.all], function(x){sum(is.na(x))})

## Normal models

lipidomics.results.pc.1.20 <- multiassociation(data.for.analysis.v1,c("Lipometrix_intensity", paste("Lipometrix_PC", 1:20, sep=".")),covariates.all)
proteomics.results.pc.1.20 <- multiassociation(data.for.analysis.v1,c("Somalogic_intensity", paste("Somalogic_PC", 1:20, sep=".")),covariates.all)

rename.vars <- data.frame(VarName=c("Age_LP", "Sex_1m_2f", "pTau", "Abeta42", "Tecnica_Empleada_LCR", "Total_protein_CSF", "QAlb", "BMI", "Sample_longevity_Years",  "e4_alleles", "e2_alleles", "PRS_Bellenguez_etal_2020_dosage"),
                          Rename=c("Age", "Sex", "CSF p-tau181", "CSF A??42", "AD Biomarker technique", "CSF turbidimetry", "QAlb index", "BMI",  "Storage time", "APOE ??4", "APOE ??2", "AD PRS (Bellenguez et al. 2022)"))
for(i in 1:nrow(rename.vars)){
  lipidomics.results.pc.1.20$indep.var <- gsub(rename.vars$VarName[i], rename.vars$Rename[i],lipidomics.results.pc.1.20$indep.var)
  proteomics.results.pc.1.20$indep.var <- gsub(rename.vars$VarName[i], rename.vars$Rename[i],proteomics.results.pc.1.20$indep.var)
}

### Include Intensity and PC1-3
forest.pcs(lipidomics.results.pc.1.20,c("Lipometrix_intensity", paste("Lipometrix_PC", 1:3, sep=".")), xlim = c(-2,2),tag.plot = "Lipometrix_Int_PC1-3")
forest.pcs(proteomics.results.pc.1.20,c("Somalogic_intensity", paste("Somalogic_PC", 1:3, sep=".")), xlim = c(-2,2),tag.plot = "Somascan_Int_PC1-3")


### LASSO simulations (intensity vs PCs, Supp Fig. 2) 

y <- data.for.analysis.v1$Lipometrix_intensity[!is.na(data.for.analysis.v1$Avail_Lipometrix)]
x <- as.matrix(data.for.analysis.v1[!is.na(data.for.analysis.v1$Avail_Lipometrix), paste("Lipometrix_PC", 1:20, sep=".")])

cv_model <- cv.glmnet(x, y, alpha = 1)
best_lambda <- cv_model$lambda.min
best_lambda
plot(cv_model) 

best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
coef(best_model)
barplot(coef(best_model)[-1]/sum(coef(best_model)[-1]), names.arg = paste("PC", 1:20, sep=""), las=2, main="Lipometrix PCs, LASSO coefficients")
boxplot(coef(best_model)[-1]/sum(coef(best_model)[-1]), main="Lipometrix, Distribution of LASSO coefficients")


y <- data.for.analysis.v1$Somalogic_intensity[!is.na(data.for.analysis.v1$Avail_Somascan)]
x <- as.matrix(data.for.analysis.v1[!is.na(data.for.analysis.v1$Avail_Somascan), paste("Somalogic_PC", 1:20, sep=".")])


cv_model <- cv.glmnet(x, y, alpha = 1)
best_lambda <- cv_model$lambda.min
best_lambda
plot(cv_model) 

best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
coef(best_model)
barplot(coef(best_model)[-1]/sum(coef(best_model)[-1]), names.arg = paste("PC", 1:20, sep=""), las=2, main="SOMAscan PCs, LASSO coefficients")
boxplot(coef(best_model)[-1]/sum(coef(best_model)[-1]),  main="SOMAscan, Distribution of LASSO coefficients")

## Standardized data pheatmaps (Fig.1 )

check.raw.data.intensities <- function(input.data,select.variables,sort.variables=NULL,sort.analytes=NULL,cluster.rows=F,cluster.cols=F,show_rownames=F,show_colnames=F){
  
  input.data <- input.data
  select.variables <- select.variables
  sort.variables=sort.variables
  cluster.rows=cluster.rows
  cluster.cols=cluster.cols
  
  row.names(input.data) <- input.data$Codigo_ACE_LCR
  if(length(sort.variables)>0){input.data<-input.data[order(input.data[,sort.variables]),]}
  
  if(length(sort.analytes)>0){input.data<-input.data[,sort.analytes]}else{input.data <- input.data[,select.variables]}
  
  
  
  pheatmap(input.data, cluster_rows = cluster.rows, cluster_cols = cluster.cols, show_rownames=show_rownames, show_colnames=show_colnames)
}

## Sort by MSV and z-score

### Lipometrix 
intraplatform.report <- read.xlsx("data/IntraplatformQC_report.xlsx")
lipid.corresp <- read.table("data/lipid_corresp.txt", h=T)
lipid.species <- lipid.corresp[which(lipid.corresp$V1%in%intraplatform.report$Species[which(intraplatform.report$Quality=="Good"|intraplatform.report$Quality=="Caution")]),2]
selected.species.lipometrix <- lipid.species[lipid.species%in%names(data.for.analysis.v1)]

indep.var <- "Lipometrix_intensity"
covars <- NULL
get.tvals <- multiassociation(input.data = data.for.analysis.v1,
                              analytes = selected.species.lipometrix,
                              independent.vars=c(indep.var,covars),
                              make.volcano=T,
                              volcano.var=indep.var,
                              tag=indep.var,
                              legend.pos="bottomleft",
                              output.type="summary")

check.raw.data.intensities(data.for.analysis.v1[!is.na(data.for.analysis.v1$Avail_Lipometrix),],selected.species.lipometrix, sort.variables="Lipometrix_intensity", sort.analytes = get.tvals$marker[order(get.tvals$`t value`)])

## Somascan
selected.proteins.somascan <- grep("seq",names(data.for.analysis.v1), value=T)

indep.var <- "Somalogic_intensity"
covars <- NULL
get.tvals <- multiassociation(input.data = data.for.analysis.v1,
                              analytes = selected.proteins.somascan,
                              independent.vars=c(indep.var,covars),
                              make.volcano=T,
                              volcano.var=indep.var,
                              tag=indep.var,
                              legend.pos="bottomleft",
                              output.type="summary")


check.raw.data.intensities(data.for.analysis.v1[!is.na(data.for.analysis.v1$Avail_Somascan),],selected.proteins.somascan, sort.variables="Somalogic_intensity", sort.analytes = get.tvals$marker[order(get.tvals$`t value`)])


### Clustered

check.correlations <- function(input.data, transpose=T, cluster.rows=F, cluster.cols=F, show.names=F, grep.cols=NULL,display.numbers=F){
  
  input.data <- input.data
  transpose <- transpose
  cluster.rows <- cluster.rows
  cluster.cols <- cluster.cols
  
  input.data <- input.data[,-1]
  
  if(transpose){input.data <- t(input.data)}
  
  for(i in 1:ncol(input.data)){input.data[,i] <- scale(as.numeric(input.data[,i]))} 
  correlation.matrix <- cor(input.data, use = "pairwise.complete.obs", method = "pearson")
  correlation.matrix[is.na(correlation.matrix)] <- 0
  
  if(length(grep.cols)>0){correlation.matrix<-correlation.matrix[grep(grep.cols[1], colnames(correlation.matrix)),
                                                                 grep(grep.cols[2], rownames(correlation.matrix))]}
  
  plot.heatmap <- pheatmap(correlation.matrix, cluster_rows = cluster.rows, cluster_cols = cluster.cols, 
                           show_rownames=show.names, show_colnames=show.names, display_numbers = display.numbers,breaks=seq(-1, 1, length.out=101))
  print(plot.heatmap)
  
}

# Lipometrix
check.correlations(data.for.analysis.v1[!is.na(data.for.analysis.v1$Avail_Lipometrix),c("Codigo_ACE_LCR", selected.species.lipometrix)], 
                   cluster.rows = T, cluster.cols = T)

# Somascan
check.correlations(data.for.analysis.v1[!is.na(data.for.analysis.v1$Avail_Somascan),c("Codigo_ACE_LCR", selected.proteins.somascan)], 
                   cluster.rows = T, cluster.cols = T)



#### Correlations between MSV (intensity) and PCs (Fig. 2b)

check.correlations <- function(input.data, transpose=T, cluster.rows=F, cluster.cols=F, show.names=F, grep.cols=NULL,display.numbers=F, number_size=8){
  
  input.data <- input.data
  transpose <- transpose
  cluster.rows <- cluster.rows
  cluster.cols <- cluster.cols
  
  input.data <- input.data[,-1]
  
  if(transpose){input.data <- t(input.data)}
  
  input.data <- sapply(input.data, scale)
  correlation.matrix <- cor(input.data, use = "pairwise.complete.obs", method = "pearson")
  correlation.matrix[is.na(correlation.matrix)] <- 0
  
  if(length(grep.cols)>0){correlation.matrix<-correlation.matrix[grep(grep.cols[1], colnames(correlation.matrix)),grep(grep.cols[2], rownames(correlation.matrix))]}
  
  plot.heatmap <- pheatmap(correlation.matrix, cluster_rows = cluster.rows, cluster_cols = cluster.cols, show_rownames=show.names, show_colnames=show.names, display_numbers = display.numbers,breaks=seq(-1, 1, length.out=101), fontsize_number = number_size)
  print(plot.heatmap)
  
}

check.correlations(data.for.analysis.v1[,c("Codigo_ACE_LCR","Lipometrix_intensity","Somalogic_intensity",
                                           paste("Lipometrix_PC", 1:4, sep="."), paste("Somalogic_PC", 1:4, sep="."))], 
                   transpose = F,show.names=T, 
                   grep.cols=c("Lipometrix_PC|Lipometrix_intensity|Somalogic_intensity", "Somalogic_PC|Lipometrix_intensity|Somalogic_intensity"), 
                   display.numbers = T)


### MCI to AD-dementia conversion #### 

conversion.models <- function(input.data, independent.vars, analytes, conversion, print.summary=F, return.all=F){
  
  input.data <- input.data
  conversion <- conversion
  analytes <- analytes
  independent.vars <- independent.vars
  
  if(conversion=="Conv_Dementia"){input.data[,conversion] <- input.data$MCI_dementia_conversion_1y_0n}
  if(conversion=="Conv_AD"){
    input.data[,conversion] <- input.data$MCI_dementia_conversion_1y_0n
    input.data[,conversion][!is.na(input.data$MCI_conversion_1AD_0DementiaNoAD)] <- input.data$MCI_conversion_1AD_0DementiaNoAD[!is.na(input.data$MCI_conversion_1AD_0DementiaNoAD)]
  }
  
  report.assoc <- data.frame()
  for(i in 1:(length(analytes))){
    
    temp.model <- input.data[,c(analytes[i],conversion,"Years_FollowUp", independent.vars)]
    temp.model <- temp.model[complete.cases(temp.model),]
    mymodel <- as.formula(paste(paste("Surv(Years_FollowUp,", conversion,")", sep=""), paste(c(analytes[i],independent.vars), collapse=" + "), sep=" ~ "))
    temp.assoc <- coxph(mymodel, data=temp.model)
    if(print.summary){print(summary(temp.assoc))}
    n.assoc <- rep(nrow(temp.model), length(independent.vars)+1)
    marker <- gsub("^L_", "", rep(analytes[i], length(independent.vars)+1))
    indep.var <- c(analytes[i],independent.vars)
    r1 <- cbind(marker,indep.var,summary(temp.assoc)$coefficients,confint(temp.assoc, level = 0.95), n.assoc)
    row.names(r1) <- paste(analytes[i], row.names((r1)), sep=".")
    
    report.assoc <- rbind(report.assoc, r1)
  }
  
  names(report.assoc)[names(report.assoc) == "Pr(>|z|)"] <- "pval"
  names(report.assoc)[names(report.assoc) == "n.assoc"] <- "N"
  names(report.assoc)[names(report.assoc) == "coef"] <- "logHR"
  names(report.assoc)[names(report.assoc) == "exp(coef)"] <- "HR"
  names(report.assoc)[names(report.assoc) == "se(coef)"] <- "SE"
  
  report.assoc[, c("logHR", "HR", "SE", "z", "pval", "2.5 %", "97.5 %", "N")] <- sapply(report.assoc[, c("logHR", "HR", "SE", "z", "pval", "2.5 %", "97.5 %", "N")], function(x) {as.numeric(as.character(x))})
  
  main.variables <- report.assoc[report.assoc$indep.var%in%analytes,]
  main.variables$P_bonf <- p.adjust(main.variables$pval, method="bonferroni", n=length(analytes))
  
  if(return.all){return(report.assoc)}else{return(main.variables)}
  return(main.variables)
}

## Run Adjusted Cox models (Supp. Table 5)
convAD.lipid.PCs  <- conversion.models(input.data = data.for.analysis.v1, 
                                       independent.vars=c("Sex_1m_2f", "Age_LP"), 
                                       analytes=c("Lipometrix_intensity", paste("Lipometrix_PC", 1:20, sep=".")), 
                                       conversion="Conv_AD")
convAD.prot.PCs  <- conversion.models(input.data = data.for.analysis.v1, 
                                      independent.vars=c("Sex_1m_2f", "Age_LP"), 
                                      analytes=c("Somalogic_intensity", paste("Somalogic_PC", 1:20, sep=".")), 
                                      conversion="Conv_AD")

convAD.lipid.PCs[,8:9] <- exp(convAD.lipid.PCs[,8:9])
convAD.prot.PCs[,8:9] <- exp(convAD.prot.PCs[,8:9])

# write.xlsx(convAD.lipid.PCs, "convAD_lipidomics.xlsx")
# write.xlsx(convAD.prot.PCs, "convAD_proteomics.xlsx")

## Kaplan Meier plots (Fig. 2e)
kmplot <- function(input.data, conversion, variable.conv){

  input.data <- input.data
  conversion <- conversion
  variable.conv <- variable.conv
  
  print(variable.conv)
  print(conversion)
  
  if(conversion=="Conv_Dementia"){input.data[,conversion] <- input.data$MCI_dementia_conversion_1y_0n}
  if(conversion=="Conv_AD"){
    input.data$conversion <- input.data$MCI_dementia_conversion_1y_0n
    input.data$conversion[!is.na(input.data$MCI_conversion_1AD_0DementiaNoAD)] <- input.data$MCI_conversion_1AD_0DementiaNoAD[!is.na(input.data$MCI_conversion_1AD_0DementiaNoAD)]
  }
  
  input.data <- input.data[!is.na(input.data$conversion)&!is.na(input.data[,variable.conv]),]
  
  print(nrow(input.data))
  
  ## Tertiles
  # Create tertiles
  input.data$Tertile <- cut(input.data[,variable.conv], 
                            breaks = quantile(input.data[,variable.conv], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                            include.lowest = TRUE, 
                            labels = c("Low", "Medium", "High")) 
  # levels(input.data$Tertile) <- c("High", "Medium", "Low")
  
  # View(input.data[,c("Years_FollowUp",conversion,variable.conv)])
  
  fit <- survfit(Surv(Years_FollowUp, conversion) ~ Tertile, data = input.data)
  ggsurvplot(fit, data = input.data, conf.int = T, ) + ggtitle(paste(conversion, variable.conv)) + xlab("Time (years)")
}

# Lipidomics
kmplot(input.data=data.for.analysis.v1, conversion="Conv_AD", variable.conv="Lipometrix_intensity")
ggsave("KMPlot_ConvAD_Lipometrix_intensity.png", width=1100, height=1250, units="px")

kmplot(input.data=data.for.analysis.v1, conversion="Conv_AD", variable.conv="Lipometrix_PC.1")
ggsave("KMPlot_ConvAD_Lipometrix_PC.1.png", width=1100, height=1250, units="px")

kmplot(input.data=data.for.analysis.v1, conversion="Conv_AD", variable.conv="Lipometrix_PC.2")
ggsave("KMPlot_ConvAD_Lipometrix_PC.2.png", width=1100, height=1250, units="px")

# Proteomics
kmplot(input.data=data.for.analysis.v1, conversion="Conv_AD", variable.conv="Somalogic_intensity")
ggsave("KMPlot_ConvAD_Somalogic_intensity.png", width=1100, height=1250, units="px")

kmplot(input.data=data.for.analysis.v1, conversion="Conv_AD", variable.conv="Somalogic_PC.1")
ggsave("KMPlot_ConvAD_Somalogic_PC.1.png", width=1100, height=1250, units="px")

kmplot(input.data=data.for.analysis.v1, conversion="Conv_AD", variable.conv="Somalogic_PC.2")
ggsave("KMPlot_ConvAD_Somalogic_PC.2.png", width=1100, height=1250, units="px")


####### VENTRICULAR VOLUME INCREASING ALLELES ####### 

### GWAs association pheatmap (Fig 4a.)

sort.cols <- function(input.data){
  input.data<-input.data
  names(input.data)[grep("PC\\.\\d$", names(input.data))] <- gsub("PC\\.", "PC\\.0", names(input.data)[grep("PC\\.\\d$", names(input.data))])
  input.data <- input.data[,c(grep("Somalogic|Lipometrix", names(input.data), invert = T), order(names(input.data))[grep("Somalogic|Lipometrix", names(input.data))])]
  names(input.data) <- gsub("PC\\.0", "PC\\.", names(input.data))
  return(input.data)
}

summary.results <- read.table("data/Summary_Results_Ventricular_Volume_Vojinovic_2018.txt",h=T)
row.names(summary.results) <- summary.results$ID

external.sumstats <- read.xlsx("data/Ventricular_Volume_Vojinovic_2018.xlsx")
external.sumstats$RA <- external.sumstats$A1
external.sumstats$RA[which(external.sumstats$Zscore<0)] <- external.sumstats$A2[which(external.sumstats$Zscore<0)]

# Flip summary results to align Risk Alleles

for(i in 1:nrow(external.sumstats)){
  
  tested.snp <- external.sumstats$ID[i]
  match.snp <- which(summary.results$ID==tested.snp)
  
  if(external.sumstats$RA[i]!=summary.results$REF[match.snp]&external.sumstats$RA[i]!=summary.results$ALT[match.snp]){
    print(paste("Alleles dont match for", tested.snp, ", Exiting"))
    break
  }
  
  if(external.sumstats$RA[i]==summary.results$REF[match.snp]){
    print(i)
    summary.results[match.snp,grep("^BETA\\.|^T_STAT\\.", names(summary.results))] <- (-1)*summary.results[match.snp,grep("^BETA\\.|^T_STAT\\.", names(summary.results))]
  }
}

tested.phenos <- c("Abeta40", "Abeta42", "pTau", "Total_protein_CSF", "Lipometrix_intensity", paste("Lipometrix_PC", 1:2, sep="."), 
                   "Somalogic_intensity", paste("Somalogic_PC", 1:2, sep="."))

effects.matrix <- as.matrix(summary.results[,paste("T_STAT", tested.phenos, sep=".")])
colnames(effects.matrix) <- tested.phenos
p.matrix <- as.matrix(summary.results[,paste("P", tested.phenos, sep=".")])

effects.matrix<-effects.matrix[external.sumstats$ID,]
p.matrix<-p.matrix[external.sumstats$ID,]

one.ast <- which(p.matrix<0.05)
two.ast <- which(p.matrix<0.01)
three.ast <- which(p.matrix<0.001)
p.matrix[one.ast] <- "*"
p.matrix[two.ast] <- "**"
p.matrix[three.ast] <- "***"
p.matrix[-c(one.ast,two.ast, three.ast)] <- ""

color.range<-max(abs(effects.matrix),na.rm = T)
breaksList = seq(-color.range, color.range, by = 0.1)

pheatmap(effects.matrix,cluster_rows = F, cluster_cols = F, show_rownames=T, show_colnames=T,display_numbers = p.matrix,fontsize_number=15,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),breaks = breaksList)


#### VVIA proteomic and lipidomics signatures:

## Flip to align with increased ventricular volume.
flip.snps <- c("chr16.87191495.G.A_G", "chr7.2718304.C.T_C", "chr12.106083027.T.C_T")
data.for.analysis.v1[,flip.snps] <- sapply(data.for.analysis.v1[,flip.snps], function(x){x<-2-x})

test.snps <- c("chr3.190902357.G.A_G",
               "chr16.87191495.G.A_G",
               "chr7.2718304.C.T_C",
               "chr12.106083027.T.C_T",
               "chr22.37714443.C.T_C",
               "chr10.21589215.T.A_T",
               "chr11.111207001.A.G_A")


#### Lipidomics
selected.analytes <- selected.species.lipometrix
covars <- c("Age_LP", "Sex_1m_2f","BMI", "Sample_longevity_Years")

lipidomics.snp.assoc.results <- lapply(test.snps, function(indep.var){
  multiassociation(input.data = data.for.analysis.v1,
                   analytes = selected.analytes,
                   independent.vars=c(indep.var,covars),
                   make.volcano=T,
                   volcano.var=indep.var,
                   tag=indep.var,
                   legend.pos="bottomleft",
                   output.type="summary")
})

tags <- c("GMNC", "C16orf95", "AMZ1", "NUAK1", "TRIOBP", "MLLT10", "ARHGAP20")
keep.cols <- c("marker", "Estimate", "t value", "Padj_FDR")

for.merge <- list()
for(i in 1:length(lipidomics.snp.assoc.results)){
  temp.x <- lipidomics.snp.assoc.results[[i]][,keep.cols]
  temp.x$DE <-""
  temp.x$DE [which(temp.x$Padj_FDR<0.05&temp.x$Estimate<0)] <- "Downregulated"
  temp.x$DE [which(temp.x$Padj_FDR<0.05&temp.x$Estimate>0)] <- "Upregulated"
  names(temp.x)[-1] <- paste(names(temp.x)[-1], tags[i], sep=".")
  for.merge[[i]]<-temp.x
}

merged.results.lipidomics <- Reduce(function(x, y) merge(x,y,by="marker"), for.merge)

merged.results.lipidomics$Overall_regulation <- colSums(t(merged.results.lipidomics[,grep("^DE", names(merged.results.lipidomics))]=="Upregulated"))-colSums(t(merged.results.lipidomics[,grep("^DE", names(merged.results.lipidomics))]=="Downregulated"))
merged.results.lipidomics$Overall_eff_direction <- colSums(t(merged.results.lipidomics[,grep("^Estimate", names(merged.results.lipidomics))]>0))-colSums(t(merged.results.lipidomics[,grep("^Estimate", names(merged.results.lipidomics))]<0))


## Upregulated lipids Venn (Supp. Fig. 6)

protein.list <- list(merged.results.lipidomics$marker[which(merged.results.lipidomics$DE.GMNC=="Upregulated")],
                     merged.results.lipidomics$marker[which(merged.results.lipidomics$DE.C16orf95=="Upregulated")],
                     merged.results.lipidomics$marker[which(merged.results.lipidomics$DE.AMZ1=="Upregulated")])


myCol <- brewer.pal(length(protein.list), "Pastel2")
plt<-venn.diagram(protein.list, category.names=c("GMNC", "C16orf95", "AMZ1/GNA12"), filename=NULL, fill = myCol,print.mode	=c("raw", "percent"))
grid.newpage()
grid::grid.draw(plt)


## Downregulated lipids Venn (Supp. Fig. 6)

protein.list <- list(merged.results.lipidomics$marker[which(merged.results.lipidomics$DE.GMNC=="Downregulated")],
                     merged.results.lipidomics$marker[which(merged.results.lipidomics$DE.C16orf95=="Downregulated")],
                     merged.results.lipidomics$marker[which(merged.results.lipidomics$DE.AMZ1=="Downregulated")])

myCol <- brewer.pal(length(protein.list), "Pastel2")
plt<-venn.diagram(protein.list, category.names=c("GMNC", "C16orf95", "AMZ1/GNA12"), filename=NULL, fill = myCol, print.mode	=c("raw", "percent"))
grid.newpage()
grid::grid.draw(plt)


##### Proteomics 

#### Lipidomics
selected.analytes <- selected.proteins.somascan
covars <- c("Age_LP", "Sex_1m_2f","BMI", "Sample_longevity_Years")
list.good.proteins <- read.xlsx("data/SOMAscan_GoodProteins_QC.xlsx")

proteomics.snp.assoc.results <- lapply(test.snps, function(indep.var){
  multiassociation(input.data = data.for.analysis.v1,
                   analytes = selected.analytes,
                   independent.vars=c(indep.var,covars),
                   make.volcano=T,
                   volcano.var=indep.var,
                   tag=indep.var,
                   legend.pos="bottomleft",
                   output.type="summary")
})



tags <- c("GMNC", "C16orf95", "AMZ1", "NUAK1", "TRIOBP", "MLLT10", "ARHGAP20")
keep.cols <- c("marker", "Estimate", "Std. Error", "t value", "pval", "Padj_FDR")

for.merge <- list()
for(i in 1:length(proteomics.snp.assoc.results)){
  temp.x <- proteomics.snp.assoc.results[[i]][,keep.cols]
  temp.x$DE <-""
  temp.x$DE [which(temp.x$Padj_FDR<0.05&temp.x$Estimate<0)] <- "Downregulated"
  temp.x$DE [which(temp.x$Padj_FDR<0.05&temp.x$Estimate>0)] <- "Upregulated"
  names(temp.x)[-1] <- paste(names(temp.x)[-1], tags[i], sep=".")
  for.merge[[i]]<-temp.x
}

merged.results.proteomics <- Reduce(function(x, y) merge(x,y,by="marker"), for.merge)
names(merged.results.proteomics) <- gsub("Std. Error", "SE", names(merged.results.proteomics)) %>% gsub(" ", "_", .)

merged.results.proteomics$Overall_regulation <- colSums(t(merged.results.proteomics[,grep("^DE", names(merged.results.proteomics))]=="Upregulated"))-colSums(t(merged.results.proteomics[,grep("^DE", names(merged.results.proteomics))]=="Downregulated"))
merged.results.proteomics$Overall_eff_direction <- colSums(t(merged.results.proteomics[,grep("^Estimate", names(merged.results.proteomics))]>0))-colSums(t(merged.results.proteomics[,grep("^Estimate", names(merged.results.proteomics))]<0))

annotated.merged.results.proteomics <- merge(list.good.proteins[,1:10],merged.results.proteomics, by.x="TargetID_SS", by.y="marker")

## Export
# write.xlsx(annotated.merged.results.proteomics, "Annotated_Somascan_results_VentVolSNPs.xlsx")

## Venn diagrams
protein.list <- list(merged.results.proteomics$marker[which(merged.results.proteomics$DE.GMNC=="Upregulated")],
                     merged.results.proteomics$marker[which(merged.results.proteomics$DE.C16orf95=="Upregulated")],
                     merged.results.proteomics$marker[which(merged.results.proteomics$DE.AMZ1=="Upregulated")])

myCol <- brewer.pal(length(protein.list), "Pastel2")
plt<-venn.diagram(protein.list, category.names=c("GMNC", "C16orf95", "AMZ1/GNA12"), filename=NULL, fill = myCol,print.mode	=c("raw", "percent"))
grid.newpage()
grid::grid.draw(plt)

protein.list <- list(merged.results.proteomics$marker[which(merged.results.proteomics$DE.GMNC=="Downregulated")],
                     merged.results.proteomics$marker[which(merged.results.proteomics$DE.C16orf95=="Downregulated")],
                     merged.results.proteomics$marker[which(merged.results.proteomics$DE.AMZ1=="Downregulated")])

myCol <- brewer.pal(length(protein.list), "Pastel2")
plt<-venn.diagram(protein.list, category.names=c("GMNC", "C16orf95", "AMZ1/GNA12"), filename=NULL, fill = myCol, print.mode	=c("raw", "percent"))
grid.newpage()
grid::grid.draw(plt)

### Write tables for enrichment analysis
somascan.prepare.for.GSEA <- function(results.table.proteomics){
  results.table.proteomics <- results.table.proteomics
  
  results.table.proteomics <- merge(list.good.proteins[,1:10],results.table.proteomics, by.x="TargetID_SS", by.y="marker")
  results.table.proteomics <- results.table.proteomics[order(results.table.proteomics$logPnegativo, decreasing = T),]
  results.table.proteomics <- results.table.proteomics[!duplicated(results.table.proteomics$EntrezGeneSymbol),] # Eliminar duplicados, dejando el mas significativo
  results.table.proteomics$signedlogP <- results.table.proteomics$logPnegativo*sign(results.table.proteomics$Estimate)
  results.table.proteomics <- results.table.proteomics[order(results.table.proteomics$signedlogP, decreasing=T),]
  return(results.table.proteomics)
}

# proteomics.GMNC.normal.forGSEA <- somascan.prepare.for.GSEA(proteomics.GMNC.normal)
# write.table(proteomics.GMNC.normal.forGSEA[,c("EntrezGeneSymbol", "signedlogP")], "GMNC_proteomics_forGSEA.rnk", row.names=F, col.names=F, quote=F, sep="\t")
# write.table(proteomics.GMNC.normal.forGSEA$EntrezGeneSymbol[which(proteomics.GMNC.normal.forGSEA$Padj_FDR<0.05&proteomics.GMNC.normal.forGSEA$Estimate>0)],"Upregulated_GMNC_forNTA.txt", row.names=F, col.names=F, quote=F)
# write.table(proteomics.GMNC.normal.forGSEA$EntrezGeneSymbol[which(proteomics.GMNC.normal.forGSEA$Padj_FDR<0.05&proteomics.GMNC.normal.forGSEA$Estimate<0)],"Downregulated_GMNC_forNTA.txt", row.names=F, col.names=F, quote=F)
# 
# proteomics.C16orf95.normal.forGSEA <- somascan.prepare.for.GSEA(proteomics.C16orf95.normal)
# write.table(proteomics.C16orf95.normal.forGSEA[,c("EntrezGeneSymbol", "signedlogP")], "C16orf95_proteomics_forGSEA.rnk", row.names=F, col.names=F, quote=F, sep="\t")
# write.table(proteomics.C16orf95.normal.forGSEA$EntrezGeneSymbol[which(proteomics.C16orf95.normal.forGSEA$Padj_FDR<0.05&proteomics.C16orf95.normal.forGSEA$Estimate>0)],"Upregulated_C16orf95_forNTA.txt", row.names=F, col.names=F, quote=F)
# write.table(proteomics.C16orf95.normal.forGSEA$EntrezGeneSymbol[which(proteomics.C16orf95.normal.forGSEA$Padj_FDR<0.05&proteomics.C16orf95.normal.forGSEA$Estimate<0)],"Downregulated_C16orf95_forNTA.txt", row.names=F, col.names=F, quote=F)
# 
# proteomics.AMZ1.normal.forGSEA <- somascan.prepare.for.GSEA(proteomics.AMZ1.normal)
# write.table(proteomics.AMZ1.normal.forGSEA[,c("EntrezGeneSymbol", "signedlogP")], "AMZ1_proteomics_forGSEA.rnk", row.names=F, col.names=F, quote=F, sep="\t")
# write.table(proteomics.AMZ1.normal.forGSEA$EntrezGeneSymbol[which(proteomics.AMZ1.normal.forGSEA$Padj_FDR<0.05&proteomics.AMZ1.normal.forGSEA$Estimate>0)],"Upregulated_AMZ1_forNTA.txt", row.names=F, col.names=F, quote=F)
# write.table(proteomics.AMZ1.normal.forGSEA$EntrezGeneSymbol[which(proteomics.AMZ1.normal.forGSEA$Padj_FDR<0.05&proteomics.AMZ1.normal.forGSEA$Estimate<0)],"Downregulated_AMZ1_forNTA.txt", row.names=F, col.names=F, quote=F)


### Save vulcanos (Supp. Fig. 6-7)
tested.variants <- data.frame(Variant=c("chr3:190902357:G:A", "chr16:87191495:G:A", "chr7:2718304:C:T", 
                                        "chr12:106083027:T:C", "chr22:37714443:C:T", "chr10:21589215:T:A", 
                                        "chr11:111207001:A:G"),
                              Gene=c("GMNC", "C16orf95", "AMZ1/GNA12", "NUAK1", "TRIOBP", "MLLT10", "ARHGAP20/C11orf53"))
### Lipidomics

volcano.results<- lapply(1:length(lipidomics.snp.assoc.results), function(i){
  
  temp.table <- lipidomics.snp.assoc.results[[i]]
  temp.table$Significance <- ifelse(
    temp.table$Padj_FDR  < 0.05, "FDR<0.05",
    ifelse(temp.table$pval < 0.05, "p<0.05", "Non-significant")
  )
  temp.table$Significance <- factor(
    temp.table$Significance,
    levels = c("Non-significant", "p<0.05", "FDR<0.05")
  )
  
  lambda.gc <- median(qchisq(1-temp.table$pval,1))/qchisq(0.5,1) # guardar lambda
  
  ggplot(temp.table, aes(Estimate, logPnegativo , color=Significance)) + geom_point() +
    scale_color_manual(values = c("Non-significant" = "black",
                                  "p<0.05" = "blue",
                                  "FDR<0.05" = "red")) + 
    ggtitle(tested.variants$Gene[i])  +theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
    annotate("label", x = Inf, y = Inf,
             label = paste0("?? = ", round(lambda.gc, 2)),
             hjust = 1.1, vjust = 2,
             size = 5,                   # smaller text size
             fill = alpha("white", 0.7), # semi-transparent background
             label.size = NA) + 
    labs(x = "Beta", y = "-log10(P)")

})

export.plots <- do.call(grid.arrange, c(volcano.results, ncol = 3))

ggsave("VolcanoPlot_Lipidomics_VVIAs_ACE.png", plot=export.plots, width=3000, height=2500, units="px", )

### Proteomics

volcano.results<- lapply(1:length(proteomics.snp.assoc.results), function(i){
  
  temp.table <- proteomics.snp.assoc.results[[i]]
  temp.table$Significance <- ifelse(
    temp.table$Padj_FDR  < 0.05, "FDR<0.05",
    ifelse(temp.table$pval < 0.05, "p<0.05", "Non-significant")
  )
  temp.table$Significance <- factor(
    temp.table$Significance,
    levels = c("Non-significant", "p<0.05", "FDR<0.05")
  )
  
  lambda.gc <- median(qchisq(1-temp.table$pval,1))/qchisq(0.5,1) # guardar lambda
  
  ggplot(temp.table, aes(Estimate, logPnegativo , color=Significance)) + geom_point() +
    scale_color_manual(values = c("Non-significant" = "black",
                                  "p<0.05" = "blue",
                                  "FDR<0.05" = "red")) + 
    ggtitle(tested.variants$Gene[i])  +theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
    annotate("label", x = Inf, y = Inf,
             label = paste0("?? = ", round(lambda.gc, 2)),
             hjust = 1.1, vjust = 2,
             size = 5,                   # smaller text size
             fill = alpha("white", 0.7), # semi-transparent background
             label.size = NA) + 
    labs(x = "Beta", y = "-log10(P)")
  
})

export.plots <- do.call(grid.arrange, c(volcano.results, ncol = 3))
ggsave("VolcanoPlot_Proteomics_VVIAs_ACE.png", plot=export.plots, width=3000, height=2500, units="px", )


#### Integrate proteogenomic signature

### Run raw models
run.raw.models <- function(pheno, analytes){
  
  
  all.results <- lapply(analytes, function(x){
    model.res <- summary(lm(data.for.analysis.v1[,pheno] ~ data.for.analysis.v1[,x]))
    c(model.res$coefficients[2,], model.res$r.squared, model.res$adj.r.squared)
  })
  
  temp.df <- cbind(analytes, do.call("rbind", all.results)) %>% as.data.frame()
  names(temp.df) <- c("marker", "Beta", "SE", "zscore", "pval", "rsq", "adj.rsq")
  temp.df[,-1] <- sapply(temp.df[,-1], as.numeric)
  names(temp.df)[-1] <- paste(names(temp.df)[-1], pheno, sep=".")
  return(temp.df)
  
}


## Metaanalyze Lipidomics and proteomics PCs
do.metaanalysis <- function(predictors,phenos.formeta){
  
  
  table.formeta <- lapply(phenos.formeta, function(x){
    temp <- predictors[,c("marker", paste(c("Beta", "SE"),x,sep="."))]
    names(temp) <- c("marker","Beta", "SE")
    temp
  }) %>% do.call("rbind", .) %>% as.data.frame()
  
  meta.results <- lapply(selected.analytes, function(x){
    
    m1 <- metagen(Beta, SE, data=table.formeta[table.formeta$marker==x,-1], comb.fixed = T, comb.random = T, prediction=TRUE, sm="SMD")
    meta.summary <- summary(m1)
    result.statistics <- c(x,
                           meta.summary$TE.fixed,
                           meta.summary$seTE.fixed,
                           meta.summary$pval.fixed,
                           meta.summary$TE.random,
                           meta.summary$seTE.random,
                           meta.summary$pval.random,
                           meta.summary$I2)
    
  }) %>% do.call("rbind", .) %>% as.data.frame()
  
  names(meta.results) <- c("marker", "Beta_Fixed", "SE_Fixed", "pval_Fixed","Beta_Random", "SE_Random", "pval_Random", "I2")
  meta.results[,-1] <- sapply(meta.results[,-1], as.numeric)
  return(meta.results)
  
}


## Run raw models
test.phenotypes <- c("Somalogic_PC.1", "Somalogic_PC.2", "Lipometrix_PC.1", "Lipometrix_PC.2")

raw.models.results <- lapply(c("Somalogic_PC.1", "Somalogic_PC.2", "Lipometrix_PC.1", "Lipometrix_PC.2"), run.raw.models, selected.analytes)
proteomics.predictors <- Reduce(function(x, y) merge(x,y,by="marker"), raw.models.results)

## Metaanalyze
meta.pc1<- do.metaanalysis(proteomics.predictors,c("Somalogic_PC.1", "Lipometrix_PC.1"))
meta.pc2<- do.metaanalysis(proteomics.predictors,c("Somalogic_PC.2", "Lipometrix_PC.2"))

## Merge all
names(meta.pc1)[-1] <- paste(names(meta.pc1)[-1], "META_PC1",  sep=".")
names(meta.pc2)[-1] <- paste(names(meta.pc2)[-1], "META_PC2",  sep=".")

proteomics.pc.predictors <- merge(proteomics.predictors,meta.pc1,by="marker") %>% 
  merge(., meta.pc2, by="marker") %>% 
  merge(list.good.proteins[,1:10],.,by.x="TargetID_SS", by.y="marker")
names(proteomics.pc.predictors)[1] <- "marker"

## Merge VentVol SNPs results

vent.vol.snps.results <- annotated.merged.results.proteomics

proteomics.pc.predictors <- merge(proteomics.pc.predictors,vent.vol.snps.results[,-c(2:10)], by.x="marker", by.y="TargetID_SS")

### Run Cox models - MCI to AD

convAD.proteomics  <- conversion.models(input.data = data.for.analysis.v1, 
                                        independent.vars=c("Sex_1m_2f", "Age_LP"), 
                                        analytes=selected.analytes, 
                                        conversion="Conv_AD")


proteomics.pc.predictors.v1 <- merge(proteomics.pc.predictors,convAD.proteomics[,c("marker", "HR", "pval")], by="marker")

## Run ANOVAs

table.for.anova <- data.for.analysis.v1 %>% 
  mutate(Syndromic_Status=recode_factor(data.for.analysis.v1$Syndromic_Status, "1"="CU", "2"="MCI", "3"="Dementia") %>%
           replace_na("CU")) %>% 
  filter(!is.na(Syndromic_Status)&!is.na(Avail_Somascan))


proteomics.pc.predictors.v1$ANOVA_Syndromic_Status <- lapply(proteomics.pc.predictors.v1$marker, function(x){
  summary(aov(table.for.anova[!is.na(table.for.anova[,x]),x] ~ table.for.anova$Syndromic_Status[!is.na(table.for.anova[,x])]))[[1]]$`Pr(>F)`[1]
}) %>% unlist()

## Assoc sex

results.sex <- lapply(selected.analytes, function(x){
  summary(lm(data.for.analysis.v1$Sex_1m_2f ~ data.for.analysis.v1[,x]))$coefficients[-1,c(1,3,4)]
}) %>% do.call("rbind", .) %>% 
  as.data.frame() %>% cbind(selected.analytes, .) %>% as.data.frame()
names(results.sex) <- c("marker", "Beta_sex", "z_sex", "pval_sex")
results.sex[,-1] <- sapply(results.sex[,-1], as.numeric)

## Assoc age

results.age <- lapply(selected.analytes, function(x){
  summary(lm(data.for.analysis.v1$Age_LP ~ data.for.analysis.v1[,x]))$coefficients[-1,c(1,3,4)]
}) %>% do.call("rbind", .) %>% 
  as.data.frame() %>% cbind(selected.analytes, .) %>% as.data.frame()

names(results.age) <- c("marker", "Beta_age", "z_age", "pval_age")
results.age[,-1] <- sapply(results.age[,-1], as.numeric)

proteomics.pc.predictors.v1 <- merge(proteomics.pc.predictors.v1, results.sex, by="marker") %>% 
  merge(., results.age, by="marker")

proteomics.pc.predictors.v1[, c("Pbonf", "ANOVA_bonf", "Pbonf_sex", "Pbonf_age")] <- 
  sapply(proteomics.pc.predictors.v1[, c("pval", "ANOVA_Syndromic_Status", "pval_sex", "pval_age")], function(x){ifelse(x*2395>1,1,x*2395)})

## Add AUCs 2 years 

### DETERMINED AUCS for ptau abeta42 normal
# Mean SD 2Years AUC normal biomarkers (10-fold, 100 repetitions): 82.26193, 4.992148
# Mean SD 5Years AUC normal biomarkers (10-fold, 100 repetitions): 88.9011, 7.13667

auc.2years <- read.table("data/AUC_summary_2Years_K10_Repeats100_SomaScan2395.txt", h=T)
auc.2years$marker <- gsub(".*,", "", auc.2years$Included_variables) %>% gsub("^Abeta42_adj_", "",.)
auc.2years <- auc.2years[c("marker", "mean_AUC", "sd_AUC")]
names(auc.2years) <- c("marker", "Mean_AUC_2Y", "SD_AUC_2Y")
auc.2years$Improves_AUC2Y <- auc.2years$Mean_AUC_2Y > 83.26 # 1 point increase

proteomics.pc.predictors.v2<-merge(proteomics.pc.predictors.v1,auc.2years,by="marker", all.x=T)

## Add lambdas after adjustment

lambda.results <- c("data/adjusted_abeta42_models_lambda.txt", "data/adjusted_ptau_models_lambda.txt")
inflation.results <- lapply(lambda.results, read.table, h=T) %>% Reduce(function(x, y) merge(x,y,by="marker"),.)
names(inflation.results)[-1] <- c("lambda_abeta42", "lambda_ptau")
proteomics.pc.predictors.v2<-merge(proteomics.pc.predictors.v2,inflation.results,by="marker", all.x=T)

## Combined rank
proteomics.pc.predictors.v2$comb_rank <- rank(proteomics.pc.predictors.v2$lambda_abeta42) + rank(proteomics.pc.predictors.v2$lambda_ptau) + rank(-proteomics.pc.predictors.v2$Mean_AUC_2Y)

# write.xlsx(proteomics.pc.predictors.v2, "Proteogenomic_signature.xlsx")

## Export with signed logP values for webgestalt
write.table(snp.proteome.assocs.predictors[,c("EntrezGeneSymbol", "META_PC1.Signed_random_logP")], "PC1_MetaRandomlogP_proteomics_forGSEA.rnk", row.names=F, col.names=F, quote=F, sep="\t")
write.table(snp.proteome.assocs.predictors[,c("EntrezGeneSymbol", "META_PC2.Signed_random_logP")], "PC2_MetaRandomlogP_proteomics_forGSEA.rnk", row.names=F, col.names=F, quote=F, sep="\t")


#### Plots

## Supp Fig 11.

## PC1
pc1.top.random.effect <- order(proteomics.pc.predictors.v2$Beta_Random.META_PC1, decreasing = TRUE)[1:50]
pc1.top.fixed.effect <- order(proteomics.pc.predictors.v2$Beta_Fixed.META_PC1, decreasing = TRUE)[1:50]
pc1.bot.random.effect <- order(proteomics.pc.predictors.v2$Beta_Random.META_PC1, decreasing = FALSE)[1:50]
pc1.bot.fixed.effect <- order(proteomics.pc.predictors.v2$Beta_Fixed.META_PC1, decreasing = FALSE)[1:50]

png("Random_effect_PC1.png", width=350, height=350)
plot(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.1, proteomics.pc.predictors.v2$Beta.Somalogic_PC.1, main="PC1, Random Effect", xlab="Beta Lipometrix PC1", ylab="Beta SomaScan PC1", xlim=c(-0.7,0.95), ylim=c(-0.7,0.95))
points(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.1[pc1.top.random.effect], proteomics.pc.predictors.v2$Beta.Somalogic_PC.1[pc1.top.random.effect], col="red")
points(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.1[pc1.bot.random.effect], proteomics.pc.predictors.v2$Beta.Somalogic_PC.1[pc1.bot.random.effect], col="blue")
dev.off()

png("Fixed_effect_PC1.png", width=350, height=350)
plot(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.1, proteomics.pc.predictors.v2$Beta.Somalogic_PC.1, main="PC1, Fixed Effect", xlab="Beta Lipometrix PC1", ylab="Beta SomaScan PC1", xlim=c(-0.7,0.95), ylim=c(-0.7,0.95))
points(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.1[pc1.top.fixed.effect], proteomics.pc.predictors.v2$Beta.Somalogic_PC.1[pc1.top.fixed.effect], col="red")
points(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.1[pc1.bot.fixed.effect], proteomics.pc.predictors.v2$Beta.Somalogic_PC.1[pc1.bot.fixed.effect], col="blue")
dev.off()

## PC2
pc2.top.random.effect <- order(proteomics.pc.predictors.v2$Beta_Random.META_PC2, decreasing = TRUE)[1:50]
pc2.top.fixed.effect <- order(proteomics.pc.predictors.v2$Beta_Fixed.META_PC2, decreasing = TRUE)[1:50]
pc2.bot.random.effect <- order(proteomics.pc.predictors.v2$Beta_Random.META_PC2, decreasing = FALSE)[1:50]
pc2.bot.fixed.effect <- order(proteomics.pc.predictors.v2$Beta_Fixed.META_PC2, decreasing = FALSE)[1:50]

png("Random_effect_PC2.png", width=350, height=350)
plot(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.2, proteomics.pc.predictors.v2$Beta.Somalogic_PC.2, main="PC2, Random Effect", xlab="Beta Lipometrix PC2", ylab="Beta SomaScan PC2", xlim=c(-0.7,0.95), ylim=c(-0.7,0.95))
points(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.2[pc2.top.random.effect], proteomics.pc.predictors.v2$Beta.Somalogic_PC.2[pc2.top.random.effect], col="red")
points(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.2[pc2.bot.random.effect], proteomics.pc.predictors.v2$Beta.Somalogic_PC.2[pc2.bot.random.effect], col="blue")
dev.off()

png("Fixed_effect_PC2.png", width=350, height=350)
plot(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.2, proteomics.pc.predictors.v2$Beta.Somalogic_PC.2, main="PC2, Fixed Effect", xlab="Beta Lipometrix PC2", ylab="Beta SomaScan PC2", xlim=c(-0.7,0.95), ylim=c(-0.7,0.95))
points(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.2[pc2.top.fixed.effect], proteomics.pc.predictors.v2$Beta.Somalogic_PC.2[pc2.top.fixed.effect], col="red")
points(proteomics.pc.predictors.v2$Beta.Lipometrix_PC.2[pc2.bot.fixed.effect], proteomics.pc.predictors.v2$Beta.Somalogic_PC.2[pc2.bot.fixed.effect], col="blue")
dev.off()

## N FDR upregulated PC1 random effect
sum(proteomics.pc.predictors.v2$Beta_Fixed.META_PC1>0&p.adjust(proteomics.pc.predictors.v2$pval_Random.META_PC1, method="fdr", n = nrow(proteomics.pc.predictors.v2))<0.05) # 1859
## N FDR downregulated PC1 random effect
sum(proteomics.pc.predictors.v2$Beta_Fixed.META_PC1<0&p.adjust(proteomics.pc.predictors.v2$pval_Random.META_PC1, method="fdr", n = nrow(proteomics.pc.predictors.v2))<0.05) # 261

## N FDR upregulated PC2 random effect
sum(proteomics.pc.predictors.v2$Beta_Fixed.META_PC2>0&p.adjust(proteomics.pc.predictors.v2$pval_Random.META_PC2, method="fdr", n = nrow(proteomics.pc.predictors.v2))<0.05) # 548
## N FDR downregulated PC2 random effect
sum(proteomics.pc.predictors.v2$Beta_Fixed.META_PC2<0&p.adjust(proteomics.pc.predictors.v2$pval_Random.META_PC2, method="fdr", n = nrow(proteomics.pc.predictors.v2))<0.05) # 199

padj.filter=5e-2

beta.re.pc1 <- proteomics.pc.predictors.v2$Beta_Random.META_PC1
beta.re.pc2 <- proteomics.pc.predictors.v2$Beta_Random.META_PC2
padj.re.pc1 <- p.adjust(proteomics.pc.predictors.v2$pval_Random.META_PC1, method="fdr", n = nrow(proteomics.pc.predictors.v2))
padj.re.pc2 <- p.adjust(proteomics.pc.predictors.v2$pval_Random.META_PC2, method="fdr", n = nrow(proteomics.pc.predictors.v2))

protein.list <- list(proteomics.pc.predictors.v2$marker[beta.re.pc1>0&padj.re.pc1<padj.filter],
                     proteomics.pc.predictors.v2$marker[beta.re.pc1<0&padj.re.pc1<padj.filter],
                     proteomics.pc.predictors.v2$marker[beta.re.pc2>0&padj.re.pc2<padj.filter],
                     proteomics.pc.predictors.v2$marker[beta.re.pc2<0&padj.re.pc2<padj.filter])

myCol <- brewer.pal(length(protein.list), "Pastel2")
plt<-venn.diagram(protein.list, category.names=c("PC1_up", "PC1_down", "PC2_up", "PC2_down"), filename=NULL, fill = myCol,print.mode	=c("raw", "percent"), cex=1)
grid.newpage()
grid::grid.draw(plt)
ggsave(plot=plt, "Venn_verysignificant_somascan.png", width=1700, height=1700, units="px")

lapply(protein.list, length)


### Export lists for ORA
## FDR 0.05
padj.filter=5e-02

protein.list <- list(proteomics.pc.predictors.v2$EntrezGeneSymbol[beta.re.pc1>0&padj.re.pc1<padj.filter],
                     proteomics.pc.predictors.v2$EntrezGeneSymbol[beta.re.pc1<0&padj.re.pc1<padj.filter],
                     proteomics.pc.predictors.v2$EntrezGeneSymbol[beta.re.pc2>0&padj.re.pc2<padj.filter],
                     proteomics.pc.predictors.v2$EntrezGeneSymbol[beta.re.pc2<0&padj.re.pc2<padj.filter]) %>% lapply(., unique)

filenames <- paste(padj.filter, c("Upreg_PC1","Downreg_PC1", "Upreg_PC2", "Downreg_PC2"), "txt", sep=".")

### Check GTEX (Supp. Fig 13)

write.table(paste(unique(proteomics.pc.predictors.v2$EntrezGeneSymbol[pc1.top.random.effect[1:25]]), collapse=", "), "ForGTEX_Top50Upreg_PC1.txt", quote=F, col.names=F, row.names=F)
write.table(paste(unique(proteomics.pc.predictors.v2$EntrezGeneSymbol[pc1.bot.random.effect[1:25]]), collapse=", "), "ForGTEX_Top50Downreg_PC1.txt", quote=F, col.names=F, row.names=F)
write.table(paste(unique(proteomics.pc.predictors.v2$EntrezGeneSymbol[pc2.top.random.effect[1:25]]), collapse=", "), "ForGTEX_Top50Upreg_PC2.txt", quote=F, col.names=F, row.names=F)
write.table(paste(unique(proteomics.pc.predictors.v2$EntrezGeneSymbol[pc2.bot.random.effect[1:25]]), collapse=", "), "ForGTEX_Top50Downreg_PC2.txt", quote=F, col.names=F, row.names=F)

### Fig. 5a 

proteomics.pc.predictors.v2 <- read.xlsx("data/Proteogenomic_signature.xlsx")
proteomics.pc.predictors.v2$lambda_ptau[proteomics.pc.predictors.v2$lambda_ptau==Inf] <- 151.14915

a <- ggplot(proteomics.pc.predictors.v2, 
            aes(x = Beta_Random.META_PC1^2, y = Mean_AUC_2Y)) +
  geom_point(color = "lightblue", alpha = 0.4, size = 2) +
  geom_smooth(color = "darkblue", fill = "blue", alpha = 0.25, size = 1.2) +
  geom_hline(yintercept = 82.4, linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = 83.0, linetype = "dashed", alpha = 0.6, color="red") +
  theme_classic2() +
  labs(x="PC1 Squared Beta", y="2-year AUC")

b <- ggplot(proteomics.pc.predictors.v2, 
            aes(x = Beta_Random.META_PC1^2, y = lambda_ptau)) +
  geom_point(color = "lightblue", alpha = 0.4, size = 2) +
  geom_smooth(color = "darkblue", fill = "blue", alpha = 0.25, size = 1.2) +
  geom_hline(yintercept = 151.14915, linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = 14.9, linetype = "dashed", alpha = 0.6, color="red") +
  theme_classic2()+
  labs(x="PC1 Squared Beta", y="p-tau ??")

c <- ggplot(proteomics.pc.predictors.v2, 
            aes(x = Beta_Random.META_PC1^2, y = lambda_abeta42)) +
  geom_point(color = "lightblue", alpha = 0.4, size = 2) +
  geom_smooth(color = "darkblue", fill = "blue", alpha = 0.25, size = 1.2) +
  geom_hline(yintercept = 62.3, linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = 10.4, linetype = "dashed", alpha = 0.6, color="red") +
  theme_classic2()+
  labs(x="PC1 Squared Beta", y="AB42 ??")

d <- ggplot(proteomics.pc.predictors.v2, 
            aes(x = Beta_Random.META_PC2^2, y = Mean_AUC_2Y)) +
  geom_point(color = "lightblue", alpha = 0.4, size = 2) +
  geom_smooth(color = "darkblue", fill = "blue", alpha = 0.25, size = 1.2) +
  geom_hline(yintercept = 82.4, linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = 83.0, linetype = "dashed", alpha = 0.6, color="red") +
  theme_classic2()+
  labs(x="PC2 Squared Beta", y="2-year AUC")

e <- ggplot(proteomics.pc.predictors.v2, 
            aes(x = Beta_Random.META_PC2^2, y = lambda_ptau)) +
  geom_point(color = "lightblue", alpha = 0.4, size = 2) +
  geom_smooth(color = "darkblue", fill = "blue", alpha = 0.25, size = 1.2) +
  geom_hline(yintercept = 151.14915, linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = 14.9, linetype = "dashed", alpha = 0.6, color="red") +
  theme_classic2()+
  labs(x="PC2 Squared Beta", y="p-tau ??")

f <- ggplot(proteomics.pc.predictors.v2, 
            aes(x = Beta_Random.META_PC2^2, y = lambda_abeta42)) +
  geom_point(color = "lightblue", alpha = 0.4, size = 2) +
  geom_smooth(color = "darkblue", fill = "blue", alpha = 0.25, size = 1.2) +
  geom_hline(yintercept = 62.3, linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = 10.4, linetype = "dashed", alpha = 0.6, color="red") +
  theme_classic2()+
  labs(x="PC2 Squared Beta", y="AB42 ??")

ggarrange(plotlist=list(a,d,b,e,c,f), nrow=3, ncol = 2)
ggsave("PC1_and_PC2_trend_with_AUC_and_lambdas.png", height=1500, width=1250, units="px")


#### Plot WebGestalt Results (Supp. Fig. 14)
## Variants
path.enrichments <- "C:/Users/pablo/OneDrive - FUNDACIO ACE/FACE/Escritorio/PROYECTOS/HARPONE/QTLs/Analisis/cosas/Webgestalt/01.feb2025/Ventvol"
gsea.variants <- lapply(dir(path.enrichments), function(x){
  read.table(paste(path.enrichments,x,grep("enrichment_results", dir(paste(path.enrichments,x,sep="/")), value=T),sep="/"), h=T, sep="\t") %>% mutate(EnrichmentDb=x)
})

## PC1
path.enrichments <- "C:/Users/pablo/OneDrive - FUNDACIO ACE/FACE/Escritorio/PROYECTOS/HARPONE/QTLs/Analisis/cosas/Webgestalt/03.apr2025/PC1"
gsea.pc1 <- lapply(dir(path.enrichments), function(x){
  read.table(paste(path.enrichments,x,grep("enrichment_results", dir(paste(path.enrichments,x,sep="/")), value=T),sep="/"), h=T, sep="\t") %>% mutate(EnrichmentDb=x)
})

## PC2
path.enrichments <- "C:/Users/pablo/OneDrive - FUNDACIO ACE/FACE/Escritorio/PROYECTOS/HARPONE/QTLs/Analisis/cosas/Webgestalt/03.apr2025/PC2"
gsea.pc2 <- lapply(dir(path.enrichments), function(x){
  read.table(paste(path.enrichments,x,grep("enrichment_results", dir(paste(path.enrichments,x,sep="/")), value=T),sep="/"), h=T, sep="\t") %>% mutate(EnrichmentDb=x)
})

gsea.variants.df <- do.call("rbind", gsea.variants)
gsea.pc1.df <- do.call("rbind", gsea.pc1)
gsea.pc2.df <- do.call("rbind", gsea.pc2)


### Plot enrichments (FDR significant results)

ggplot(gsea.variants.df %>% 
         mutate(Signed_logP=abs(log10(pValue))*enrichmentScore,
                Geneset=paste(ifelse(grepl("^GO", geneSet), "GO:", "KEGG:"), description)) %>%
         arrange(Signed_logP) %>% 
         mutate(Geneset = factor(Geneset, levels = Geneset))) + 
  geom_col(aes(Signed_logP, Geneset,fill = Signed_logP > 0)) + 
  scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "#0066ff")) + 
  theme(legend.position = "none")
ggsave("GSEA_VVIAs.png", width=2000, height=2500, units="px")


ggplot(gsea.pc1.df %>% 
         mutate(Geneset=paste(ifelse(grepl("^GO", geneSet), "GO:", "KEGG:"), description)) %>%
         arrange(normalizedEnrichmentScore) %>% 
         filter(FDR<0.05) %>% 
         mutate(Geneset = factor(Geneset, levels = Geneset))) + 
  geom_col(aes(normalizedEnrichmentScore, Geneset,fill = normalizedEnrichmentScore > 0)) + 
  scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "#0066ff")) + 
  theme(legend.position = "none")

ggsave("GSEA_PC1.png", width=2000, height=2500, units="px")

ggplot(gsea.pc2.df %>% 
         mutate(Geneset=paste(ifelse(grepl("^GO", geneSet), "GO:", "KEGG:"), description)) %>%
         arrange(normalizedEnrichmentScore) %>% 
         filter(FDR<0.05) %>% 
         mutate(Geneset = factor(Geneset, levels = Geneset))) + 
  geom_col(aes(normalizedEnrichmentScore, Geneset,fill = normalizedEnrichmentScore > 0)) + 
  scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "#0066ff")) + 
  theme(legend.position = "none")

ggsave("GSEA_PC2.png", width=2000, height=1000, units="px")

## ptau abeta ORA enrichments
path.enrichments <- "C:/Users/pablo/OneDrive - FUNDACIO ACE/FACE/Escritorio/PROYECTOS/HARPONE/QTLs/Analisis/cosas/Webgestalt/03.apr2025/ORA_ADbiom/enrichment_results_wg_result1746527741.txt"

oraptauabeta <- read.table(path.enrichments, h=T, sep="\t")
ggplot(oraptauabeta %>% 
         mutate(Geneset=paste(ifelse(grepl("^GO", geneSet), "GO:", "KEGG:"), description)) %>%
         arrange(enrichmentRatio) %>% 
         filter(FDR<0.05) %>% tail(., 30) %>% 
         mutate(Geneset = factor(Geneset, levels = Geneset))) + 
  geom_col(aes( enrichmentRatio, Geneset,fill =  enrichmentRatio > 0)) + 
  scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "#0066ff")) + 
  theme(legend.position = "none")


save.image("temp_image_script02.RData")
##### Generate Somacan associations with p-tau and AB42 (Figure 6 and Supp Table 11)

find.analyte.profile <- function(input.data,dep.var,analytes,independent.vars, tag="Volcano plot", legend.pos="bottomleft", full.report=F){
  
  # Matriz de datos
  input.data <- input.data
  
  # Especificar analitos a mano
  analytes <- analytes
  
  # Vector de variables independientes con covariables 
  independent.vars <- independent.vars
  
  # Adaptar variables categoricas
  input.data$Sex_1m_2f <- as.factor(input.data$Sex_1m_2f)
  
  report.assoc <- data.frame()
  
  for(i in 1:length(analytes)){
    
    if(sum(duplicated(c(dep.var, analytes[i],independent.vars)))>0){next} 
    
    if(i%%1000==0){print(i)}
    temp.model <- input.data[,c(dep.var, analytes[i],independent.vars)]
    temp.model <- temp.model[complete.cases(temp.model),]
    mymodel <- as.formula(paste(dep.var, paste(c(analytes[i],independent.vars), collapse=" + "), sep=" ~ "))
    temp.assoc <- lm(mymodel, data=temp.model)
    indep.assocs <- complete.cases(confint(temp.assoc))
    n.assoc <- rep(nrow(temp.model), length(indep.assocs))
    marker <- rep(analytes[i], length(indep.assocs))
    indep.var <- c("Intercept", analytes[i], independent.vars)[indep.assocs] 
    r1 <- cbind(marker,indep.var,summary(temp.assoc)$coefficients,confint(temp.assoc, level = 0.95)[complete.cases(confint(temp.assoc)),], n.assoc)
    row.names(r1) <- paste(analytes[i], row.names((r1)), sep=".")
    
    report.assoc <- rbind(report.assoc, r1)
  }
  
  
  report.assoc <- report.assoc[grep("Intercept", row.names(report.assoc), invert=T),]
  names(report.assoc)[names(report.assoc) == "Pr(>|t|)"] <- "pval"
  names(report.assoc)[names(report.assoc) == "n.assoc"] <- "N"
  
  # Se me guarda ahora como factor, convierto estas cols a numericas:
  report.assoc[,c("Estimate", "Std. Error", "t value", "pval", "2.5 %", "97.5 %")] <- sapply(report.assoc[,c("Estimate", "Std. Error", "t value", "pval", "2.5 %", "97.5 %")], function(x) {as.numeric(as.character(x))})
  
  if(full.report){return(report.assoc)}else{
    
    report.volcano <- report.assoc[which(report.assoc$indep.var%in%analytes),]
    
    # Filter out adjustment proteins
    adj.protein <- analytes[analytes%in%independent.vars]
    if(length(adj.protein)>0){report.volcano<-report.volcano[!report.volcano$indep.var==adj.protein,]}
    
    report.volcano$Padj_FDR <- p.adjust(report.volcano$pval, method="fdr", n=nrow(report.volcano))
    report.volcano$logPnegativo <- -log10(report.volcano$pval)
    
    plot(report.volcano$Estimate, report.volcano$logPnegativo, xlab="Estimate", ylab="-log10(pval)", main=tag)
    points(report.volcano$Estimate[which(report.volcano$pval<0.05)], report.volcano$logPnegativo[which(report.volcano$pval<0.05)], col="blue")
    points(report.volcano$Estimate[which(report.volcano$Padj<0.05)], report.volcano$logPnegativo[which(report.volcano$Padj<0.05)], col="red")
    abline(v=0, lty=2)
    report.volcano$chisq <- qchisq(1-report.volcano$pval,1)
    
    legend(legend.pos, c("P>=0.05", "P<0.05", "Padj<0.05"), col=c("black", "blue", "red"), pch=c(1,1,1),cex = 0.8)
    
    upregulated <- sum(report.volcano$Padj_FDR<0.05&report.volcano$Estimate>0)
    downregulated <- sum(report.volcano$Padj_FDR<0.05&report.volcano$Estimate<0)
    lambda.gc <- median(qchisq(1-report.volcano$pval,1))/qchisq(0.5,1) # guardar lambda
    print(paste("Lambda Inflation Factor =", lambda.gc))
    print(paste(upregulated, "Upregulated,", downregulated, "Downregulated"))
    btest.result <- binom.test(x = c(upregulated,downregulated), alternative = "two.sided", conf.level = 0.95)
    print(paste("Binomial test p =",format(btest.result$p.value, scientific = T, digits=2)))
    gaston::qqplot.pvalues(report.volcano$pval, col.abline = "red",CB = F)
    
    report.volcano <- report.volcano[order(report.volcano$pval),]
    names(report.volcano) <- c("marker", "indep.var", "Beta", "SE", "zscore", 
                               "pval", "CI2.5", "CI97.5", "N", "Padj_FDR", "logPnegativo", "chisq")
    
    return(report.volcano)
  }
}


make.volcano <- function(report.volcano, tag="Volcano plot", legend.pos="bottomleft", make.qqlot=T){
  
  plot(report.volcano$Beta, report.volcano$logPnegativo, xlab="Beta", ylab="-log10(pval)", main=tag)
  points(report.volcano$Beta[which(report.volcano$pval<0.05)], report.volcano$logPnegativo[which(report.volcano$pval<0.05)], col="blue")
  points(report.volcano$Beta[which(report.volcano$Padj<0.05)], report.volcano$logPnegativo[which(report.volcano$Padj<0.05)], col="red")
  abline(v=0, lty=2)
  report.volcano$chisq <- qchisq(1-report.volcano$pval,1)
  
  legend(legend.pos, c("P>=0.05", "P<0.05", "Padj<0.05"), col=c("black", "blue", "red"), pch=c(1,1,1),cex = 0.8)
  
  upregulated <- sum(report.volcano$Padj_FDR<0.05&report.volcano$Beta>0, na.rm=T)
  downregulated <- sum(report.volcano$Padj_FDR<0.05&report.volcano$Beta<0, na.rm=T)
  lambda.gc <- median(qchisq(1-report.volcano$pval,1))/qchisq(0.5,1) # guardar lambda
  print(paste("Lambda Inflation Factor =", lambda.gc))
  print(paste(upregulated, "Upregulated,", downregulated, "Downregulated"))
  btest.result <- binom.test(x = c(upregulated,downregulated), alternative = "two.sided", conf.level = 0.95)
  print(paste("Binomial test p =",format(btest.result$p.value, scientific = T, digits=2)))
  if(make.qqlot){gaston::qqplot.pvalues(report.volcano$pval, col.abline = "red",CB = F)}
  
  
}

somascan.prepare.for.ORA <- function(results.table.proteomics){
  results.table.proteomics <- results.table.proteomics
  
  results.table.proteomics <- merge(list.good.proteins[,1:10],results.table.proteomics, by.x="TargetID_SS", by.y="marker")
  results.table.proteomics <- results.table.proteomics[order(results.table.proteomics$logPnegativo, decreasing = T),]
  results.table.proteomics <- results.table.proteomics[!duplicated(results.table.proteomics$EntrezGeneSymbol),] # Eliminar duplicados, dejando el mas significativo
  results.table.proteomics <- results.table.proteomics[order(results.table.proteomics$logPnegativo, decreasing=T),]
  return(results.table.proteomics)
}


proteomics.dataset <- as.data.frame(fread("omics_clinical_data/HARPONE_Somascan_preprocessed_20240704.txt", h=T))
dim(proteomics.dataset) # 1330 7215 
row.names(proteomics.dataset) <- proteomics.dataset$Codigo_ACE_LCR
proteomics.dataset[,-1] <- scale(proteomics.dataset[,-1]) ## Scale somamers

# Merge, removing old values
data.for.analysis.v2 <- data.for.analysis.v1 %>% select(-matches("_seq\\.")) %>% left_join(., proteomics.dataset, by="Codigo_ACE_LCR")

selected.analytes <-grep("_seq\\.", names(data.for.analysis.v2), value=T)
covariates.all <- c("Age_LP", "Sex_1m_2f", "Tecnica_Empleada_LCR", "QAlb", "Sample_longevity_Years", "MMSE_PL")

run.models <- list(covariates.all,
                   c(covariates.all, "GAG2A_seq.18268.5"),
                   c(covariates.all, "GAG2A_seq.18268.5", "e4_alleles", "e2_alleles"),
                   c(covariates.all, "OBCAM_seq.15622.13"),
                   c(covariates.all, "OBCAM_seq.15622.13", "e4_alleles", "e2_alleles"),
                   c(covariates.all, "Nectin.like.protein.3_seq.16907.3"),
                   c(covariates.all, "Nectin.like.protein.3_seq.16907.3", "e4_alleles", "e2_alleles"),
                   c(covariates.all, "NPTN_seq.7194.36"),
                   c(covariates.all, "NPTN_seq.7194.36", "e4_alleles", "e2_alleles"))

names(run.models) <- c("model1", "model2_GAGE2A", "model3_GAGE2A", "model2_OPCML", "model3_OPCML",
                       "model2_CADM2", "model3_CADM2", "model2_NPTN", "model3_NPTN")

ptau.models <- lapply(run.models,function(covs){
  find.analyte.profile(data.for.analysis.v2,
                       dep.var="pTau",
                       analytes=selected.analytes,
                       independent.vars=covs)
})
names(ptau.models) <- paste("pTau", names(run.models), sep="_")

abeta.models <- lapply(run.models,function(covs){
  find.analyte.profile(data.for.analysis.v2,
                       dep.var="Abeta42",
                       analytes=selected.analytes,
                       independent.vars=covs)
})
names(abeta.models) <- paste("Abeta42", names(run.models), sep="_")

ad.biom.models <- c(abeta.models, ptau.models)
# save(ad.biom.models, file="AD_biomarker_models_20251001.RData")




