# 
# Crick/UCLH Legacy Cohort Neutralisation Data Analysis
# AZD1222 vaccine recipients
# Against SARS-CoV-2 Variants of Concern B.1.617.2, B.1.351, and B.1.1.7
#
# 24 June 2021
#
# David LV Bauer
# RNA Virus Replication Laboratory
# The Francis Crick Institute
#
# This work is licensed under a CC BY 4.0 license and is free to use with attribution
# 

### Remove all stored data from environment & reset plot ###
rm(list = ls())
dev.off()

### Load Legacy Study Data table ###
load("Crick_Legacy_2021-12-06_B-1-617-2_AZ_PUBLIC.Rda")
legacyAZ <- dtHashed
rm(dtHashed)

load("pastData/Crick_Legacy_2021-24-05_B-1-617-2_PUBLIC.Rda")
legacyPf <- dtHashed
rm(dtHashed)


### Load required packages ### 
# library(plotly)
library(tidyverse)
library(readxl)
library(khroma)
library(boot)
library(extrafont)
library(svglite)
library(ggpubr)
library(stats)
library(cowplot)
library(broom)
library(furniture)
library(rms)
# source("geom_split_violin.R")



### Set constants for various functions / plots
runBootstrapStats <- FALSE  # Set to TRUE in order to calculate bootstrap stats
runFontImport <- FALSE # Set to TRUE in order to import fonts into R environment
vaccineOrder <- c("AZD1222", "BNT162b2")
strainOrder <- c("Wildtype", "D614G", "B.1.1.7", "B.1.351", "B.1.617.2")
datasetForCorr <- c("D614G", "B.1.1.7", "B.1.351", "B.1.617.2")
referenceIC50 <- 2^9.802672




######################################################
#          Set up functions for later use            #
######################################################

if(runFontImport){
  font_import(path = "fonts/", prompt = FALSE)
}

boot_foldMedian <- function(d, i){
  d2 <- d[i,]
  firstCol <- d2 %>% pull(2)
  secondCol <- d2 %>% pull(3)
  foldMedian <- (median(firstCol, na.rm=TRUE) / median(secondCol, na.rm=TRUE))
  return(foldMedian)
}

CIbootstrap <- function(testTable, strainComparison, rpt=5000) {
  statcheck <- testTable %>% filter(strain %in% strainComparison) %>% pivot_wider(values_from = ic50, names_from = strain)
  results_boot <- boot(statcheck, boot_foldMedian, R=rpt)
  results_ci <- boot.ci(results_boot, type="basic")
  print(results_boot)
  print(results_ci)
  print("######################################################")
  print(strainComparison)
  print(c("Median1", "Median2", "Fold-decrease"))
  print(c(median(pull(statcheck, 2), na.rm=TRUE),
          median(pull(statcheck, 3), na.rm=TRUE),
          median(pull(statcheck, 2), na.rm=TRUE)/ median(pull(statcheck, 3), na.rm=TRUE) )
  )
}

boot_foldMedian_unpaired <- function(d, i){
  d2 <- d[i,]
  firstCol <- d2 %>% filter(COVID_vaccineName == "AZD1222") %>% pull("ic50")
  secondCol <- d2 %>% filter(COVID_vaccineName == "BNT162b2") %>% pull("ic50")
  foldMedian <- (median(secondCol , na.rm=TRUE) / median(firstCol, na.rm=TRUE))
  return(foldMedian)
}


CIbootstrap_unpaired <- function(testTable, strainComparison, rpt=1000, logTransform=TRUE) {
  if (logTransform == TRUE) {testTable$ic50 <- log2(testTable$ic50)}
  statcheck <- testTable %>% filter(strain %in% strainComparison)
  grp <- match(statcheck$COVID_vaccineName, unique(statcheck$COVID_vaccineName))
  results_boot <- boot(statcheck, boot_foldMedian_unpaired, R=rpt, strata=grp)
  results_ci <- boot.ci(results_boot, type="basic", conf = 0.95)
  print(results_boot)
  print(results_ci)
  if(logTransform == TRUE){
    print("Log-transformed.... antilog 2^...")
    print(2^(results_boot$t0))
    print(2^(results_ci$basic[c(4,5)]))
  }
  print("######################################################")
  print(strainComparison)
  # print(c("Median1", "Median2", "Fold-decrease"))
  # print(c(median(pull(statcheck, 2), na.rm=TRUE),
  #         median(pull(statcheck, 3), na.rm=TRUE),
  #         median(pull(statcheck, 2), na.rm=TRUE)/ median(pull(statcheck, 3), na.rm=TRUE) )
  # )
}

boot_foldMedian_unpaired_symp <- function(d, i){
  d2 <- d[i,]
  firstCol <- d2 %>% filter(COVID_symptoms == FALSE) %>% pull("ic50")
  secondCol <- d2 %>% filter(COVID_symptoms == TRUE) %>% pull("ic50")
  foldMedian <- (median(secondCol , na.rm=TRUE) / median(firstCol, na.rm=TRUE))
  return(foldMedian)
}

CIbootstrap_unpaired_symp <- function(testTable, strainComparison, rpt=1000, logTransform=TRUE) {
  if (logTransform == TRUE) {testTable$ic50 <- log2(testTable$ic50)}
  statcheck <- testTable %>% filter(strain %in% strainComparison)
  grp <- match(statcheck$COVID_symptoms, unique(statcheck$COVID_symptoms))
  results_boot <- boot(statcheck, boot_foldMedian_unpaired_symp, R=rpt, strata=grp)
  results_ci <- boot.ci(results_boot, type="basic", conf = 0.95)
  print(results_boot)
  print(results_ci)
  if(logTransform == TRUE){
    print("Log-transformed.... antilog 2^...")
    print(2^(results_boot$t0))
    print(2^(results_ci$basic[c(4,5)]))
  }
  print("######################################################")
  print(strainComparison)
  # print(c("Median1", "Median2", "Fold-decrease"))
  # print(c(median(pull(statcheck, 2), na.rm=TRUE),
  #         median(pull(statcheck, 3), na.rm=TRUE),
  #         median(pull(statcheck, 2), na.rm=TRUE)/ median(pull(statcheck, 3), na.rm=TRUE) )
  # )
}



pTtest12 <- function(relevantData, strainSel="Wuhan1", resultTable = 0){
  testTable <- relevantData %>% filter(strain ==strainSel) %>% ungroup() %>% 
    select(sample_barcode, bc_visit, ic50, ic50_log2) %>%
    mutate(bc_visit=replace(bc_visit, bc_visit > 1, 2))
  visit1 <- testTable %>% filter(bc_visit == 1) %>% pull(ic50_log2)
  visit2 <- testTable %>% filter(bc_visit == 2) %>% pull(ic50_log2)
  visitTable <- tibble(visit1, visit2) %>% drop_na()
  broomTable <- t.test(visitTable$visit1, visitTable$visit2, paired = TRUE) %>% tidy() %>% add_column(strain=strainSel, .before=1)
  if (is.numeric(resultTable)){
    return(broomTable)
  } else {
    return( bind_rows (resultTable, broomTable))
  }
}

pWtest12 <- function(relevantData, strainSel="Wuhan1", resultTable = 0){
  testTable <- relevantData %>% filter(strain ==strainSel) %>% ungroup() %>% 
    select(sample_barcode, bc_visit, ic50, ic50_log2) %>%
    mutate(bc_visit=replace(bc_visit, bc_visit > 1, 2))
  visit1 <- testTable %>% filter(bc_visit == 1) %>% pull(ic50_log2)
  visit2 <- testTable %>% filter(bc_visit == 2) %>% pull(ic50_log2)
  broomTable <- wilcox.test(visit1, visit2, paired = TRUE) %>% tidy() %>% add_column(strain=strainSel, .before=1)
  if (is.numeric(resultTable)){
    return(broomTable)
  } else {
    return( bind_rows (resultTable, broomTable))
  }
}

# From Paul McMurdie
# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}


######################################################
#    Define cohorts to analyse in this study         #
######################################################

### Merge data for further analysis
legacy <- bind_rows(legacyAZ, legacyPf)

### Drop bc_pad columns
legacy <- legacy %>% select(-bc_pad)

### Ensure if prior COVID symptoms was left blank that it is counted as "No"
legacy$COVID_symptoms <- legacy$COVID_symptoms %>% recode(.missing = 0)

### Reorder vaccine factors
legacy$COVID_vaccineName <- fct_relevel(legacy$COVID_vaccineName, vaccineOrder)

### Make readable factors
legacy$COVID_vaccStatus_pretty <- recode(legacy$COVID_vaccStatus, "1" = "One dose", "2" = "Two doses")
legacy$COVID_symptoms_pretty <- recode(legacy$COVID_symptoms, "0" = "Prior symptoms: NO", "1" = "Prior symptoms: YES")

### Pick only those that are the first visit following their jabs
dose1cohort <- legacy %>% filter(COVID_vaccStatus == 1,  sampleOrderInVaccStatus == 1)
dose2cohort <- legacy %>% filter(COVID_vaccStatus == 2,  sampleOrderInVaccStatus == 1)
doseBOTHcohort <- bind_rows(dose1cohort, dose2cohort)

### Subset those for AZ/Pfizer
dose1cohort_AZ <- dose1cohort %>% filter(COVID_vaccineName == "AZD1222")
dose1cohort_Pf <- dose1cohort %>% filter(COVID_vaccineName == "BNT162b2")
dose2cohort_AZ <- dose2cohort %>% filter(COVID_vaccineName == "AZD1222")
dose2cohort_Pf <- dose2cohort %>% filter(COVID_vaccineName == "BNT162b2")
doseBOTHcohort_AZ <- doseBOTHcohort %>% filter(COVID_vaccineName == "AZD1222")
doseBOTHcohort_Pf <- doseBOTHcohort %>% filter(COVID_vaccineName == "BNT162b2")

### Resample Pfizer cohort for Crick site, age, & long vaccine interval
dose1cohort_Pf_resam <- dose1cohort_Pf %>% filter(site == "Crick", age < 50)
dose2cohort_Pf_resam <- dose2cohort_Pf %>% filter(site == "Crick", age < 50, COVID_daysBetweenJabs >= 41)
dose1cohort_resam <- bind_rows(dose1cohort_AZ, dose1cohort_Pf_resam)
dose2cohort_resam <- bind_rows(dose2cohort_AZ, dose2cohort_Pf_resam)
doseBOTHcohort_resam <- bind_rows(dose1cohort_resam, dose2cohort_resam)

# Generate a merged set of original and resampled data for plotting comparisons (Figure 3)
doseBOTHcohort_Pf_resam_annot <- bind_rows(dose1cohort_Pf_resam, dose2cohort_Pf_resam) %>% 
                                  mutate(COVID_vaccineName = "BNT162b2 (matched)")
doseBOTHcohort_Pf_resam_annot$COVID_vaccineName <- fct_relevel(doseBOTHcohort_Pf_resam_annot$COVID_vaccineName,
                                                               c("AZD1222", "BNT162b2", "BNT162b2 (matched)"))

doseBOTHcohort_resam_ALL <- bind_rows(doseBOTHcohort, doseBOTHcohort_Pf_resam_annot) 
dose1cohort_resam_ALL <- doseBOTHcohort_resam_ALL %>% filter(COVID_vaccStatus == 1)
dose2cohort_resam_ALL <- doseBOTHcohort_resam_ALL %>% filter(COVID_vaccStatus == 2)

dose2cohort_Pf$WTtoDelta_ratio_ic50 <- dose2cohort_Pf$Wildtype_ic50 / dose2cohort_Pf$B.1.617.2_ic50
dose2cohort_Pf_resam$WTtoDelta_ratio_ic50 <- dose2cohort_Pf_resam$Wildtype_ic50 / dose2cohort_Pf_resam$B.1.617.2_ic50


dose2cohort_Pf_UCLH_short <- dose2cohort_Pf %>% filter(COVID_daysBetweenJabs == 21)
dose2cohort_Pf_UCLH_long <- dose2cohort_Pf %>% filter(COVID_daysBetweenJabs > 40, site == "UCLH")
dose2cohort_Pf_Crick_longOld <- dose2cohort_Pf %>% filter(COVID_daysBetweenJabs > 40, site == "Crick", age > 50)
dose2cohort_Pf_Crick_longYoung <- dose2cohort_Pf_resam


dose2cohort_Pf_UCLH_short %>% select(ends_with("_ic50")) %>% write_csv("dose2cohort_Pf_UCLH_short.csv")
dose2cohort_Pf_UCLH_long %>% select(ends_with("_ic50")) %>% write_csv("dose2cohort_Pf_UCLH_long.csv")
dose2cohort_Pf_Crick_longOld %>% select(ends_with("_ic50")) %>% write_csv("dose2cohort_Pf_Crick_longOld.csv")
dose2cohort_Pf_Crick_longYoung %>% select(ends_with("_ic50")) %>% write_csv("dose2cohort_Pf_Crick_longYoung.csv")


######################################################
#    Table S1. Generate table of participant data    #
######################################################

# Generate table with single row per participant
participantsAZ <- bind_rows(dose1cohort_AZ, dose2cohort_AZ) %>% group_by(bc_participant_id) %>% slice_head(n=1)
# Annotate if in each dosing cohort
participantsAZ$inDose1cohort <- participantsAZ$bc_participant_id %in% dose1cohort_AZ$bc_participant_id
participantsAZ$inDose2cohort <- participantsAZ$bc_participant_id %in% dose2cohort_AZ$bc_participant_id

# Tidy up data structure
participantsAZ$sex <- factor(participantsAZ$sex)
participantsAZ$site <- factor(participantsAZ$site)
# participantsAZ$COVID_sumJabs_niceLabels <- recode  c("One dose", "Two doses")[participantsAZ$inDose1cohort]
participantsAZ <- participantsAZ %>% ungroup()
# participantsAZ <- participantsAZ %>% group_by(COVID_sumJabs_niceLabels)

# Generate Supp Table 1 - a
participantsAZ_dose1 <- filter(participantsAZ, inDose1cohort == TRUE)
tabA <- table1((participantsAZ_dose1),
               site,
               age,
               sex,
               BMI,
               ethnicity_3,
               test = TRUE,
               # assume parametric distribution:
               param = TRUE, 
               na.rm = FALSE,
               output = "text2", export = "supplementary_table_Adose")

# Generate Supp Table 1 - b
participantsAZ_dose2 <- filter(participantsAZ, inDose2cohort == TRUE)
tabB <- table1((participantsAZ_dose2),
               site,
               age,
               sex,
               BMI,
               ethnicity_3,
               test = TRUE,
               # assume parametric distribution:
               param = TRUE, 
               na.rm = FALSE,
               output = "text2", export = "supplementary_table_Bdose")

# Generate Supp Table 1 - c

# Generate table with single row per participant
participantsPf2_resam <- dose2cohort_Pf_resam %>% group_by(bc_participant_id) %>% slice_head(n=1)

# Tidy up data structure
participantsPf2_resam$sex <- factor(participantsPf2_resam$sex)
participantsPf2_resam$site <- factor(participantsPf2_resam$site)
participantsPf2_resam <- participantsPf2_resam %>% ungroup()

tabC <- table1((participantsPf2_resam),
               site,
               age,
               sex,
               BMI,
               ethnicity_3,
               test = TRUE,
               # assume parametric distribution:
               param = TRUE, 
               na.rm = FALSE,
               output = "text2", export = "supplementary_table_Cdose")



# Narrative in manuscript...
    # Using a high-throughput live-virus SARS-CoV-2 neutralisation assay, we determined NAb titres (NAbTs) in
    # 106 participants (median age 34 years, [IQR 29-42]) 
    nrow(participantsAZ)
    fivenum(participantsAZ$age)
    # following either 1 dose (n = 50, median time after first dose = 41 days [IQR 30-51]) 
    nrow(dose1cohort_AZ)
    fivenum(dose1cohort_AZ$COVID_daysSinceJab1)
    # or 2 doses (n = 63, median time after second dose = 31 days [IQR 19.5-46],  
    nrow(dose2cohort_AZ)
    fivenum(dose2cohort_AZ$COVID_daysSinceJab2)
    # median interval between doses = 63 days [IQR 62-69.5])
    fivenum(dose2cohort_AZ$COVID_daysBetweenJabs)
    # of AZD1222 (Oxford-AstraZeneca) against five SARS-CoV-2 strains.

# Compare AZ participants' age to Pf participants
participantsPf <- bind_rows(dose1cohort_Pf, dose2cohort_Pf) %>% group_by(bc_participant_id) %>% slice_head(n=1)
t.test(participantsAZ$age, participantsPf$age)

# Compare AZ participants' age to resampled Pf participants
participantsPf_resam <- bind_rows(dose1cohort_Pf_resam, dose2cohort_Pf_resam) %>% group_by(bc_participant_id) %>% slice_head(n=1)
t.test(participantsAZ$age, participantsPf_resam$age)


### Random questions to consider...

  # Do 2-dose ChAdOx recipients' ic50 show any effect of prior symptoms?
  data_no <- dose2cohort_AZ %>% filter(COVID_symptoms == 0) %>% pull(Wildtype_ic50)
  data_yes <- dose2cohort_AZ %>% filter(COVID_symptoms == 1) %>% pull(Wildtype_ic50)
  median(data_no) - median(data_yes)
  wilcox.test(data_no, data_yes)
  
  # Do 2-dose Pfizer recipients ic50 show any effect of prior symptoms?
  data_no <- dose2cohort_Pf %>% filter(COVID_symptoms == 0) %>% pull(Wildtype_ic50)
  data_yes <- dose2cohort_Pf %>% filter(COVID_symptoms == 1) %>% pull(Wildtype_ic50)
  median(data_no) - median(data_yes)
  wilcox.test(data_no, data_yes)
  
  # Do 1-dose Pfizer recipients ic50 show any effect of prior symptoms? (NO!)
  data_no <- dose1cohort_Pf %>% filter(COVID_symptoms == 0) %>% pull(Wildtype_ic50)
  data_yes <- dose1cohort_Pf %>% filter(COVID_symptoms == 1) %>% pull(Wildtype_ic50)
  median(data_no) - median(data_yes)
  wilcox.test(data_no, data_yes)
  
  # Do resampled 1-dose Pfizer recipients ic50 show any effect of prior symptoms? (NO!)
  data_no <- dose1cohort_Pf_resam %>% filter(COVID_symptoms == 0) %>% pull(Wildtype_ic50)
  data_yes <- dose1cohort_Pf_resam %>% filter(COVID_symptoms == 1) %>% pull(Wildtype_ic50)
  median(data_no) - median(data_yes)
  wilcox.test(data_no, data_yes)
  
  # Do resampled 2-dose Pfizer recipients show any effect of age? (YES!)
  data_no <- dose2cohort_Pf_resam %>% filter(COVID_symptoms == 0) %>% pull(Wildtype_ic50)
  data_yes <- dose2cohort_Pf_resam %>% filter(COVID_symptoms == 1) %>% pull(Wildtype_ic50)
  median(data_no) - median(data_yes)
  wilcox.test(data_no, data_yes)
  
  # Do resampled 2-dose Pfizer recipients increase NAbTs against B.1.617.2? (YES!)
  wilcox.test(dose2cohort_Pf_resam$B.1.617.2_ic50, dose2cohort_Pf$B.1.617.2_ic50)
  


    
###############################################################################
#   Panel A. Vaccine responses per strain following 2nd dose  AZ vs Pfizer    #
###############################################################################
    
    relevantData <- dose2cohort %>%
      pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50")
    relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
    relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)
    
    relevantData %>% group_by(strain) %>% summarise(counts = n(), qtile = c(0.25, 0.5, 0.75),  value = quantile(ic50, c(0.25, 0.5, 0.75), na.rm=TRUE))
    
    
    outplot <- ggplot(relevantData, aes(x=COVID_vaccineName, y=ic50, color = strain, label = sample_barcode)) + 
      scale_colour_muted() +
      # scale_color_brewer(palette="Set1") + 
      geom_hline(yintercept = referenceIC50,linetype = 3 ) +
      geom_violin(trim=TRUE) + 
      # stat_summary(fun.data = give.n, geom = "text", fun.y = median,
      #              position = position_dodge(width = 0.75)) +
      # scale_y_log10(limits = c(0.1, 2^12))
      scale_y_continuous(trans='log2',
                         breaks=c(5, 10, 64, 256, 1024, 5120), 
                         labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
                         minor_breaks=c(32, 128, 512, 2048)) +
      ylab(bquote('Virus Neutralisation, '~IC[50]~ '')) +
      geom_jitter(shape=20, position=position_jitter(0.2), alpha=0.3) + 
      # facet_wrap(~ COVID_vaccStatus, labeller = label_both, nrow = 1) +
      stat_summary(fun=median, geom = "point", color="black",  shape=5, size=1, stroke=1) + 
      # stat_summary(fun=gm_mean, geom = "point", color="black",  shape=3, size=1, stroke=1) + 
      theme_bw(base_family = "Helvetica Neue Thin") +
      theme(legend.position="none")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + 
      theme(axis.text.y = element_text(size=12))  + 
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y=element_text(size=15),
        axis.title.x=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_blank()
        ) + facet_wrap(~strain, nrow=1) +
      ggtitle("Two Doses", subtitle="Oxford-AstraZeneca (AZD1222) vs. Pfizer-BioNTech (BNT162b2)")

    # outplot
    
    PanelA <- outplot
    
    ggsave("FIGURE-PanelA.jpg", outplot, width=40, height=15, units="cm")
    
    # ggplotly(outplot)
    
    ### Stats for Panel 1 ###
    grouped <- relevantData %>% group_by(strain) %>% group_by(COVID_vaccineName, .add = TRUE)
    summaryTable <- grouped %>% summarise(count=n(), median=median(ic50, na.rm = T), .groups = "keep")
    summaryTable
    # strain    COVID_vaccineName count median
    # <fct>     <fct>             <int>  <dbl>
    # 1 Wildtype  AZD1222              32  419. 
    # 2 Wildtype  BNT162b2            160  893. 
    # 3 D614G     AZD1222              32  113. 
    # 4 D614G     BNT162b2            160  393. 
    # 5 B.1.1.7   AZD1222              32  112. 
    # 6 B.1.1.7   BNT162b2            160  339. 
    # 7 B.1.351   AZD1222              32   44.4
    # 8 B.1.351   BNT162b2            160  183. 
    # 9 B.1.617.2 AZD1222              32   50.9
    # 10 B.1.617.2 BNT162b2            160  153. 
    
    # Get p-values for change in median
    
    pvals <- integer(0)
    for (i in 1:length(strainOrder)){
      currStrain <- strainOrder[i]
      testTable <- relevantData %>% ungroup %>% select(sample_barcode, strain, COVID_vaccineName, ic50) %>% filter(strain == currStrain)
      AZic50 <- testTable %>% filter(COVID_vaccineName == "AZD1222") %>% pull(ic50)
      Pfic50 <- testTable %>% filter(COVID_vaccineName == "BNT162b2") %>% pull(ic50)
      Wresult <- wilcox.test(AZic50, Pfic50)
      pvals[i] <- Wresult$p.value
      print("/////////////////////////////////////")
      print(currStrain)
      CIbootstrap_unpaired(testTable, currStrain, rpt=1000)
      # print(Wresult)
      print("")
    }
    print("/////////////////////////////////////")
    pvals
    p.adjust(pvals, method = "holm")
    p.adjust(pvals, method = "bonferroni")
    
    
    
    

    
    
    

###############################################################################
#   Panel B. BINNED responses per strain following 2nd dose  AZ vs Pfizer     #
###############################################################################
    
    relevantData <- dose2cohort %>% 
      filter(COVID_vaccStatus %in% c(1,2)) %>%
      filter(sampleOrderInVaccStatus == 1) %>%
      pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50") %>% drop_na(ic50) %>%
      mutate(ic50Range = cut(ic50, breaks=c(0,40,256,9999), right=FALSE, label=c("[<40]", "40-256", ">256")))
    relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
    relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)
    
    relevantData$COVID_vaccStatus_pretty <- c("1 Dose", "2 Doses")[relevantData$COVID_vaccStatus]
    
    outplot <- ggplot(relevantData, aes(x=ic50Range))  + geom_bar(aes(fill=strain, alpha=0.75)) +scale_fill_muted() + coord_flip() + 
      facet_grid( strain ~ COVID_vaccineName, switch = "y", scales = "free_x", shrink = FALSE) +
      theme_bw(base_family = "Helvetica Neue Thin") +
      ylab("Participants") +
      xlab( bquote('Virus Neutralisation, '~IC[50]~ '') ) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position="none"   ,
        strip.background = element_blank(),
        strip.placement = "outside"
      )
    # outplot
    PanelB <- outplot
    ggsave("FIGURE-PanelB.jpg", outplot, width=15, height=20, units="cm")

    
######################################################################################
#   STATS - Panel B. OLR: does vaccine mfr influence response indep. of strains?     #
######################################################################################
    
    relevantData <- dose2cohort %>%
      pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50") %>% drop_na(ic50) %>%
      mutate(ic50Range = cut(ic50, breaks=c(0,40,256,9999), right=FALSE, label=c("[<40]", "40-256", ">256")))
    relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
    relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)
    
    olrTable <- relevantData %>% ungroup() %>% select(bc_participant_id , age, COVID_vaccineName, strain, ic50Range)
    nrow(olrTable)

    # Ordered logistical regression model required
    fit_first <- rms::lrm( ic50Range ~ strain*COVID_vaccineName , data=olrTable , x=TRUE, y=TRUE)
    fit_first
    anova(fit_first)
    
    # Hmisc::Ecdf(olrTable~ic50Range, group=group, fun=qlogis, data)
    
    # exact p-values rather than pretty table
    1-pchisq(coef(fit_first)^2/diag(vcov(fit_first)),1)
    
    if (runBootstrapStats == TRUE){
      my.valid <- validate(fit_first, method="boot", B=1000)
      my.calib <- calibrate(fit_first, method="boot", B=1000)
      par(bg="white", las=1)
      plot(my.calib, las=1)
    }
    
    # Output summary table of just above and below limits
    relevantData <- dose2cohort %>%
      pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50") %>% drop_na(ic50) %>%
      mutate(ic50Range = cut(ic50, breaks=c(0,40,9999), right=FALSE, label=c("[<40]", ">40")))
    relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
    relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)
    olrTable <- relevantData %>% ungroup() %>% select(bc_participant_id , age, COVID_vaccineName, strain, ic50Range)
    summaryTable <- olrTable %>% group_by(COVID_vaccineName, strain, ic50Range) %>% summarise(count=n()) 
    
    # COVID_vaccineName strain    ic50Range count
    # <fct>             <fct>     <fct>     <int>
    #   1 AZD1222           Wildtype  >40          63
    # 2 AZD1222           D614G     [<40]         8
    # 3 AZD1222           D614G     >40          55
    # 4 AZD1222           B.1.1.7   [<40]         8
    # 5 AZD1222           B.1.1.7   >40          55
    # 6 AZD1222           B.1.351   [<40]        25
    # 7 AZD1222           B.1.351   >40          38
    # 8 AZD1222           B.1.617.2 [<40]        24
    # 9 AZD1222           B.1.617.2 >40          39
    # 10 BNT162b2          Wildtype  >40         160
    # 11 BNT162b2          D614G     >40         157
    # 12 BNT162b2          B.1.1.7   >40         160
    # 13 BNT162b2          B.1.351   [<40]         9
    # 14 BNT162b2          B.1.351   >40         151
    # 15 BNT162b2          B.1.617.2 [<40]         6
    # 16 BNT162b2          B.1.617.2 >40         151
    
    # B.1.1.7
    above40 <- c(63, 55)
    total <- c(63, 63)
    prop.test(above40, total, conf.level = 0.95)
    
    # B.1.351
    above40 <- c(63, 38)
    total <- c(63, 63)
    prop.test(above40, total, conf.level = 0.95)
    
    # B.1.617.2
    above40 <- c(63, 39)
    total <- c(63, 63)
    prop.test(above40, total, conf.level = 0.95)
    
    # all p vals
    # 5.979e-06
    # 2.815e-10
    # 1.248e-09
    
        

###############################################################################
#  Panel C. AZD1222 responses strat. by prior COVID symptoms after 1st dose   #
###############################################################################
    
    relevantData <- dose1cohort_AZ %>%
      pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50")
    relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
    relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)
    
    relevantData %>% group_by(strain) %>% summarise(counts = n(), qtile = c(0.25, 0.5, 0.75),  value = quantile(ic50, c(0.25, 0.5, 0.75), na.rm=TRUE))
    
    outplot <- ggplot(relevantData, aes(x=strain, y=ic50, color = strain, label = sample_barcode)) + 
      scale_colour_muted() +
      # scale_color_brewer(palette="Set1") + 
      geom_hline(yintercept = referenceIC50,linetype = 3 ) +
      geom_violin(trim=TRUE) + 
      # stat_summary(fun.data = give.n, geom = "text", fun.y = median,
      #              position = position_dodge(width = 0.75)) +
      # scale_y_log10(limits = c(0.1, 2^12))
      scale_y_continuous(trans='log2',
                         breaks=c(5, 10, 64, 256, 1024, 5120), 
                         labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
                         minor_breaks=c(32, 128, 512, 2048)) +
      ylab(bquote('Virus Neutralisation, '~IC[50]~ '')) +
      geom_jitter(shape=20, position=position_jitter(0.2), alpha=0.3) + 
      # facet_wrap(~ COVID_vaccStatus, labeller = label_both, nrow = 1) +
      stat_summary(fun=median, geom = "point", color="black",  shape=5, size=1, stroke=1) + 
      # stat_summary(fun=gm_mean, geom = "point", color="black",  shape=3, size=1, stroke=1) + 
      theme_bw(base_family = "Helvetica Neue Thin") +
      theme(legend.position="none")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + 
      theme(axis.text.y = element_text(size=12))  + 
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y=element_text(size=15),
        axis.title.x=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_blank()
      ) + facet_grid(~COVID_symptoms_pretty, margins = TRUE) + 
      ggtitle("One Dose: Oxford-AstraZeneca (AZD1222)", subtitle="Stratification by prior COVID symptoms")
    
    
    # outplot
    
    PanelC <- outplot
    
    ggsave("FIGURE-PanelC.jpg", outplot, width=40, height=15, units="cm")
    
    # ggplotly(outplot)
    
    ### Stats for Panel C ###
    grouped <- relevantData %>% group_by(strain) %>% group_by(COVID_symptoms_pretty, .add = TRUE)
    summaryTable <- grouped %>% summarise(count=n(), median=median(ic50, na.rm = T), .groups = "keep")
    summaryTable
    # strain    COVID_symptoms_pretty count median
    # <fct>     <chr>                 <int>  <dbl>
    # 1 Wildtype  Prior symptoms: NO       34  106. 
    # 2 Wildtype  Prior symptoms: YES      16 3524. 
    # 3 D614G     Prior symptoms: NO       34   46.3
    # 4 D614G     Prior symptoms: YES      16  552. 
    # 5 B.1.1.7   Prior symptoms: NO       34   10  
    # 6 B.1.1.7   Prior symptoms: YES      16  711. 
    # 7 B.1.351   Prior symptoms: NO       34    5  
    # 8 B.1.351   Prior symptoms: YES      16  407. 
    # 9 B.1.617.2 Prior symptoms: NO       34    5  
    # 10 B.1.617.2 Prior symptoms: YES      16  334. 
    
    # Get p-values for change in median
    
    pvals <- integer(0)
    for (i in 1:length(strainOrder)){
      currStrain <- strainOrder[i]
      testTable <- relevantData %>% ungroup %>% select(sample_barcode, strain, COVID_symptoms, ic50) %>% filter(strain == currStrain)
      NoSymp_ic50 <- testTable %>% filter(COVID_symptoms == 0) %>% pull(ic50)
      YesSymp_ic50 <- testTable %>% filter(COVID_symptoms == 1) %>% pull(ic50)
      Wresult <- wilcox.test(NoSymp_ic50, YesSymp_ic50)
      pvals[i] <- Wresult$p.value
      print("/////////////////////////////////////")
      print(currStrain)
      CIbootstrap_unpaired_symp(testTable, currStrain, rpt=1000)
      print("/////////////////////////////////////")
      print("")
      # print(currStrain)
      # print(Wresult)
    }
    pvals
    p.adjust(pvals, method = "holm")
    p.adjust(pvals, method = "bonferroni")
    
    
###############################################################################
#     Panel D. BINNED responses by COVID symptoms following 1 dose AZ         #
###############################################################################
    
    relevantData <- dose1cohort_AZ %>% 
      filter(COVID_vaccStatus %in% c(1,2)) %>%
      filter(sampleOrderInVaccStatus == 1) %>%
      pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50") %>% drop_na(ic50) %>%
      mutate(ic50Range = cut(ic50, breaks=c(0,40,256,9999), right=FALSE, label=c("[<40]", "40-256", ">256")))
    relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
    relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)
    
    outplot <- ggplot(relevantData, aes(x=ic50Range))  + geom_bar(aes(fill=strain, alpha=0.75)) +scale_fill_muted() + coord_flip() + 
      facet_grid( strain ~ COVID_symptoms_pretty, switch = "y", scales = "free_x", shrink = FALSE) +
      theme_bw(base_family = "Helvetica Neue Thin") +
      ylab("Participants") +
      xlab( bquote('Virus Neutralisation, '~IC[50]~ '') ) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position="none"   ,
        strip.background = element_blank(),
        strip.placement = "outside"
      )
    # outplot
    PanelD <- outplot
    ggsave("FIGURE-PanelD.jpg", outplot, width=15, height=20, units="cm")
    
######################################################################################
#   STATS - Panel D. OLR: does prior COVID symptoms affect NAb indep. of strains?    #
######################################################################################

    
    # NOW CONSIDER WITH SYMPTOMS
    relevantData <- dose1cohort_AZ %>%
      pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50") %>% drop_na(ic50) %>%
      mutate(ic50Range = cut(ic50, breaks=c(0,40,256,9999), right=FALSE, label=c("[<40]", "40-256", ">256")))
    relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
    relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)
    
    olrTable <- relevantData %>% ungroup() %>% select(bc_participant_id , age, COVID_vaccineName, COVID_symptoms, strain, ic50Range)
    nrow(olrTable)
    fit_first <- rms::lrm( ic50Range ~ strain*COVID_symptoms , data=olrTable , x=TRUE, y=TRUE)
    fit_first
    anova(fit_first)
    
    # exact p-values rather than pretty table
    1-pchisq(coef(fit_first)^2/diag(vcov(fit_first)),1)
    
    if (runBootstrapStats == TRUE){
      my.valid <- validate(fit_first, method="boot", B=1000)
      my.calib <- calibrate(fit_first, method="boot", B=1000)
      par(bg="white", las=1)
      plot(my.calib, las=1)
    }
    
    # Output summary table of just above and below limits
    relevantData <- dose1cohort_AZ %>%
      pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50") %>% drop_na(ic50) %>%
      mutate(ic50Range = cut(ic50, breaks=c(0,40,9999), right=FALSE, label=c("[<40]", ">40")))
    relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
    relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)
    olrTable <- relevantData %>% ungroup() %>% select(bc_participant_id , age, COVID_vaccineName, COVID_symptoms, strain, ic50Range)
    summaryTable <- olrTable %>% group_by(COVID_vaccineName, strain, COVID_symptoms, ic50Range) %>% summarise(count=n()) 
  
    # COVID_vaccineName strain    COVID_symptoms ic50Range count
    # <fct>             <fct>              <dbl> <fct>     <int>
    # 1 AZD1222           Wildtype               0 [<40]         3
    # 2 AZD1222           Wildtype               0 >40          31
    # 3 AZD1222           Wildtype               1 >40          16
    # 4 AZD1222           D614G                  0 [<40]        13
    # 5 AZD1222           D614G                  0 >40          20
    # 6 AZD1222           D614G                  1 [<40]         1
    # 7 AZD1222           D614G                  1 >40          14
    # 8 AZD1222           B.1.1.7                0 [<40]        22
    # 9 AZD1222           B.1.1.7                0 >40          12
    # 10 AZD1222           B.1.1.7                1 [<40]         3
    # 11 AZD1222           B.1.1.7                1 >40          13
    # 12 AZD1222           B.1.351                0 [<40]        30
    # 13 AZD1222           B.1.351                0 >40           4
    # 14 AZD1222           B.1.351                1 [<40]         4
    # 15 AZD1222           B.1.351                1 >40          12
    # 16 AZD1222           B.1.617.2              0 [<40]        29
    # 17 AZD1222           B.1.617.2              0 >40           4
    # 18 AZD1222           B.1.617.2              1 [<40]         4
    # 19 AZD1222           B.1.617.2              1 >40          11
    
    # B.1.1.7
    above40 <- c(31, 12)
    total <- c(34, 34)
    prop.test(above40, total, conf.level = 0.95)
    
    # B.1.351
    above40 <- c(31, 4)
    total <- c(34, 34)
    prop.test(above40, total, conf.level = 0.95)
    
    # B.1.617.2
    above40 <- c(31, 5)
    total <- c(34, 34)
    prop.test(above40, total, conf.level = 0.95)
    
    # all p vals
    # 5.979e-06
    # 2.815e-10
    # 1.248e-09
    
    
########################################################################
#                  OUTPUT FIGURE 1 COMBINED                            #
########################################################################

row1 <- plot_grid(PanelA, PanelB, nrow=1, rel_widths = c(3,1.5) )
row2 <- plot_grid(PanelC, PanelD, nrow=1, rel_widths = c(3,1.5) )
block12 <- plot_grid(row1, NULL, row2, rel_heights = c(1, 0.3, 1), ncol=1)
ggsave("FIGURE-Block12_1p25x.jpg", block12, width=25, height=23, units="cm", dpi = 600)
ggsave("FIGURE-Block12_1p25x.svg", block12, width=25, height=23, units="cm")


########################################################################
#    Figure 2 - PanelA - stratification of all data                    #
########################################################################

relevantData <- doseBOTHcohort %>%
  pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50")
relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)

relevantData %>% group_by(strain) %>% summarise(counts = n(), qtile = c(0.25, 0.5, 0.75),  value = quantile(ic50, c(0.25, 0.5, 0.75), na.rm=TRUE))

outplot <- ggplot(relevantData, aes(x=strain, y=ic50, color = strain, label = sample_barcode)) + 
  scale_colour_muted() +
  # scale_color_brewer(palette="Set1") + 
  geom_hline(yintercept = referenceIC50,linetype = 3 ) +
  geom_violin(trim=TRUE) + 
  # stat_summary(fun.data = give.n, geom = "text", fun.y = median,
  #              position = position_dodge(width = 0.75)) +
  # scale_y_log10(limits = c(0.1, 2^12))
  scale_y_continuous(trans='log2',
                     breaks=c(5, 10, 64, 256, 1024, 5120), 
                     labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
                     minor_breaks=c(32, 128, 512, 2048)) +
  ylab(bquote('Virus Neutralisation, '~IC[50]~ '')) +
  geom_jitter(shape=20, position=position_jitter(0.2), alpha=0.3) + 
  # facet_wrap(~ COVID_vaccStatus, labeller = label_both, nrow = 1) +
  stat_summary(fun=median, geom = "point", color="black",  shape=5, size=1, stroke=1) + 
  # stat_summary(fun=gm_mean, geom = "point", color="black",  shape=3, size=1, stroke=1) + 
  theme_bw(base_family = "Helvetica Neue Thin") +
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + 
  theme(axis.text.y = element_text(size=12))  + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y=element_text(size=15),
    axis.title.x=element_blank(),
    strip.text = element_text(size=12),
    strip.background = element_blank(),
    strip.placement = "outside"
  ) + facet_grid(COVID_vaccStatus_pretty~COVID_vaccineName+COVID_symptoms_pretty, margins = FALSE, switch="y") 

# outplot

Fig2PanelA <- outplot

ggsave("FIGURE2-PanelA.jpg", outplot, width=40, height=20, units="cm")


        
###############################################################################
#   Fig3 Resampled - Panel A Vaccine responses per strain following 2nd dose  AZ vs Pfizer    #
###############################################################################
        
        relevantData <- dose2cohort_resam_ALL %>%
          pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50")
        relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
        relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)
        
        relevantData %>% group_by(strain) %>% summarise(counts = n(), qtile = c(0.25, 0.5, 0.75),  value = quantile(ic50, c(0.25, 0.5, 0.75), na.rm=TRUE))
        
        
        outplot <- ggplot(relevantData, aes(x=COVID_vaccineName, y=ic50, color = strain, label = sample_barcode)) + 
          scale_colour_muted() +
          # scale_color_brewer(palette="Set1") + 
          geom_hline(yintercept = referenceIC50,linetype = 3 ) +
          geom_violin(trim=TRUE) + 
          # stat_summary(fun.data = give.n, geom = "text", fun.y = median,
          #              position = position_dodge(width = 0.75)) +
          # scale_y_log10(limits = c(0.1, 2^12))
          scale_y_continuous(trans='log2',
                             breaks=c(5, 10, 64, 256, 1024, 5120), 
                             labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
                             minor_breaks=c(32, 128, 512, 2048)) +
          ylab(bquote('Virus Neutralisation, '~IC[50]~ '')) +
          geom_jitter(shape=20, position=position_jitter(0.2), alpha=0.3) + 
          # facet_wrap(~ COVID_vaccStatus, labeller = label_both, nrow = 1) +
          stat_summary(fun=median, geom = "point", color="black",  shape=5, size=1, stroke=1) + 
          # stat_summary(fun=gm_mean, geom = "point", color="black",  shape=3, size=1, stroke=1) + 
          theme_bw(base_family = "Helvetica Neue Thin") +
          theme(legend.position="none")+
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + 
          theme(axis.text.y = element_text(size=12))  + 
          theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.title.y=element_text(size=15),
            axis.title.x=element_blank(),
            strip.text = element_text(size=12),
            strip.background = element_blank()
          ) + facet_wrap(~strain, nrow=1) +
          ggtitle("Two Doses", subtitle="Cohort-matched Pfizer-BioNTech (BNT162b2) vs. Oxford-AstraZeneca (AZD1222)")
        
        # outplot
        
        Panel3A <- outplot
        
        ggsave("FIGURE3-PanelA.jpg", outplot, width=40, height=15, units="cm")
        
        # ggplotly(outplot)
        
        ### Stats for Fig3 ###
        grouped <- relevantData %>% group_by(strain) %>% group_by(COVID_vaccineName, .add = TRUE)
        summaryTable <- grouped %>% summarise(count=n(), median=median(ic50, na.rm = T), .groups = "keep")
        summaryTable
        # strain    COVID_vaccineName  count median
        # <fct>     <fct>              <int>  <dbl>
        # 1 Wildtype  AZD1222               63  500. 
        # 2 Wildtype  BNT162b2             160  893. 
        # 3 Wildtype  BNT162b2 (matched)    58 1135. 
        # 4 D614G     AZD1222               63  104. 
        # 5 D614G     BNT162b2             160  393. 
        # 6 D614G     BNT162b2 (matched)    58  467. 
        # 7 B.1.1.7   AZD1222               63  112. 
        # 8 B.1.1.7   BNT162b2             160  339. 
        # 9 B.1.1.7   BNT162b2 (matched)    58  415. 
        # 10 B.1.351   AZD1222               63   49.6
        # 11 B.1.351   BNT162b2             160  183. 
        # 12 B.1.351   BNT162b2 (matched)    58  304. 
        # 13 B.1.617.2 AZD1222               63   48.4
        # 14 B.1.617.2 BNT162b2             160  153. 
        # 15 B.1.617.2 BNT162b2 (matched)    58  217. 
        
        # Get p-values for change in median
        
        pvals <- integer(0)
        for (i in 1:length(strainOrder)){
          currStrain <- strainOrder[i]
          testTable <- relevantData %>% ungroup %>% select(sample_barcode, strain, COVID_vaccineName, ic50) %>% filter(strain == currStrain)
          testTable <- testTable %>% filter(COVID_vaccineName != "BNT162b2")
          testTable$COVID_vaccineName <- recode(testTable$COVID_vaccineName, "BNT162b2 (matched)" = "BNT162b2")
          AZic50 <- testTable %>% filter(COVID_vaccineName == "AZD1222") %>% pull(ic50)
          Pfic50 <- testTable %>% filter(COVID_vaccineName == "BNT162b2") %>% pull(ic50)
          Wresult <- wilcox.test(AZic50, Pfic50)
          pvals[i] <- Wresult$p.value
          print("/////////////////////////////////////")
          print(currStrain)
          CIbootstrap_unpaired(testTable, currStrain, rpt=1000)
          # print(Wresult)
          print("")
        }
        print("/////////////////////////////////////")
        pvals
        p.adjust(pvals, method = "holm")
        p.adjust(pvals, method = "bonferroni")
        
        
        ###############################################################################
        #   Panel 3B. BINNED responses per strain following 2nd dose  AZ vs Pfizer cohort matched     #
        ###############################################################################
        
        relevantData <- dose2cohort_resam_ALL %>% 
          filter(COVID_vaccStatus %in% c(1,2)) %>%
          filter(sampleOrderInVaccStatus == 1) %>%
          pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50") %>% drop_na(ic50) %>%
          mutate(ic50Range = cut(ic50, breaks=c(0,40,256,9999), right=FALSE, label=c("[<40]", "40-256", ">256")))
        relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
        relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)
        
        relevantData$COVID_vaccStatus_pretty <- c("1 Dose", "2 Doses")[relevantData$COVID_vaccStatus]
        
        outplot <- ggplot(relevantData, aes(x=ic50Range))  + geom_bar(aes(fill=strain, alpha=0.75)) +scale_fill_muted() + coord_flip() + 
          facet_grid( strain ~ COVID_vaccineName, switch = "y", scales = "free_x", shrink = FALSE) +
          theme_bw(base_family = "Helvetica Neue Thin") +
          ylab("Participants") +
          xlab( bquote('Virus Neutralisation, '~IC[50]~ '') ) +
          theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position="none"   ,
            strip.background = element_blank(),
            strip.placement = "outside"
          )
        # outplot
        Panel3B <- outplot
        ggsave("FIGURE3-PanelB.jpg", outplot, width=15, height=20, units="cm")

        ########################################################################
        #                  OUTPUT FIGURE 3 COMBINED                            #
        ########################################################################
        
        Fig3 <- plot_grid(Panel3A, Panel3B, nrow=1, rel_widths = c(3,1.5) )
        ggsave("FIGURE3-Block12_1p25x.jpg", Fig3, width=25, height=10, units="cm", dpi = 600)
        ggsave("FIGURE3-Block12_1p25x.svg", Fig3, width=25, height=10, units="cm")
        
        
        ######################################################################################
        #   STATS - Panel 3B. OLR: does vaccine type / matched influence response indep. of strains?     #
        ######################################################################################
        
        relevantData <- dose2cohort_resam_ALL %>%
          pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50") %>% drop_na(ic50) %>%
          mutate(ic50Range = cut(ic50, breaks=c(0,40,256,9999), right=FALSE, label=c("[<40]", "40-256", ">256")))
        relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
        relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)
        
        olrTable <- relevantData %>% ungroup() %>% select(bc_participant_id , age, COVID_vaccineName, strain, ic50Range)
        nrow(olrTable)
        
        # Ordered logistical regression model required
        fit_first <- rms::lrm( ic50Range ~ strain*COVID_vaccineName , data=olrTable , x=TRUE, y=TRUE)
        fit_first
        anova(fit_first)
        
        # exact p-values rather than pretty table
        1-pchisq(coef(fit_first)^2/diag(vcov(fit_first)),1)
        
        if (runBootstrapStats == TRUE){
          my.valid <- validate(fit_first, method="boot", B=1000)
          my.calib <- calibrate(fit_first, method="boot", B=1000)
          par(bg="white", las=1)
          plot(my.calib, las=1)
        }
        
        
