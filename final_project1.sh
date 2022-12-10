#!/usr/bin/env bash

install.packages("haven")
library(vegan)
library(haven)
library(ggplot2)


STORY_APS_for_R5 <- read_sav("STORY_APS_for_R5.sav")

adonis2(formula =STORY_APS_for_R5[,20:201]~(SF36_M_Wellbeing),
        data = STORY_APS_for_R5)
adonis2(formula =STORY_APS_for_R5[,20:201]~(PANAS_X_PA),
        data = STORY_APS_for_R5)
adonis2(formula =STORY_APS_for_R5[,20:201]~(PANAS_X_Jov),
        data = STORY_APS_for_R5)
adonis2(formula =STORY_APS_for_R5[,20:201]~(PANAS_X_Serenity),
        data = STORY_APS_for_R5)

nmds<-metaMDS(PA_HI_LO_2[,3:184],trymax = 1000)

plot(nmds)
nmds_points<-as.data.frame(nmds$points)
nmds_points$Participant_ID<-PA_HI_LO_2$Participant_ID
nmds_points_merged<-merge(nmds_points,PA_HI_LO_2,by.x = "Participant_ID",by.y = "Participant_ID")

ggplot(data = nmds_points_merged)+
  aes(x=MDS1,y=MDS2,color=as.factor(PA_cat))+
  geom_point()

STORY_APS_for_R5$Strep_AA_round_off <- round(STORY_APS_for_R5$Strep_AA)
STORY_APS_for_R5$Akkermansia_AA_round_off <- round(STORY_APS_for_R5$Akkermansia_AA)
STORY_APS_for_R5$Dialister_AA_round_off <- round(STORY_APS_for_R5$Dialister_AA)

summary(STORY_APS_for_R5$Strep_AA_round_off)
var(STORY_APS_for_R5$Strep_AA_round_off)

ggplot(STORY_APS_for_R5, aes(x=`Strep_AA_round_off`)) +
  geom_histogram(color="white", fill="#56B4E9") +
  theme_classic()

ggplot(STORY_APS_for_R5, aes(x=`Dialister_AA_round_off`)) +
  geom_histogram(color="white", fill="#56B4E9") +
  theme_classic()

#this needs to be version 7.3-56 to use glm.nb function
install.packages("MASS")
library(MASS)

#this is for Strep and SF36_Wellbeing only
STORY_APS_for_R5.strep.WB <- glm.nb(STORY_APS_for_R5$`Strep_AA_round_off` ~ 
                                      STORY_APS_for_R5$`SF36_M_Wellbeing`)
summary(STORY_APS_for_R5.strep.WB)

#this is for strep and SF36 Wellbeing PLUS covariates
STORY_APS_for_R5.strepCO.WB <- glm.nb(`Strep_AA_round_off` ~ `SF36_M_Wellbeing`
                                      + `MomAge` + `Mom_BMI` + `Grain` + `Meat` + `Fruit_Veg`,
                                      data=STORY_APS_for_R5)
summary(STORY_APS_for_R5.strepCO.WB)

#install this package for zero-inflated negative binomial regression
install.packages("pscl")
library("pscl")

zinb.STORY_APS_for_R5 <- zeroinfl(`Strep_AA_round_off` ~ `SF36_M_Wellbeing`, 
                                  link = "logit", dist = "negbin", data = STORY_APS_for_R5)
summary(zinb.STORY_APS_for_R5)

#this is for Dialister and SF36 Wellbeing only
STORY_APS_for_R5.Dialister.WB <- glm.nb(STORY_APS_for_R5$`Dialister_AA_round_off` ~ 
                                          STORY_APS_for_R5$`SF36_M_Wellbeing`)
summary(STORY_APS_for_R5.Dialister.WB)

#this is for Dialister and SF36 Wellbeing PLUS covariates
zinb.STORY_APS_for_R5_Dial <- zeroinfl(`Dialister_AA_round_off` ~ `SF36_M_Wellbeing`,
                                       + `MomAge` + `Mom_BMI` + `Grain` + `Meat` + `Fruit_Veg`,
                                       link = "logit", dist = "negbin", data = STORY_APS_for_R5)
summary(zinb.STORY_APS_for_R5_Dial)

#this is for Strep and PANAS_X_PA
STORY_APS_for_R5.strep.PA <- glm.nb(STORY_APS_for_R5$`Strep_AA_round_off` ~ 
                                      STORY_APS_for_R5$`PANAS_X_PA`)
summary(STORY_APS_for_R5.strep.PA)

#this is for Dialister and PANAS_X_PA
STORY_APS_for_R5.Dialister.PA <- glm.nb(STORY_APS_for_R5$`Dialister_AA_round_off` ~ 
                                          STORY_APS_for_R5$`PANAS_X_PA`)
summary(STORY_APS_for_R5.Dialister.PA)

#this is for Strep PANAS_X_Jov
STORY_APS_for_R5.strep.Jov <- glm.nb(STORY_APS_for_R5$`Strep_AA_round_off` ~ 
                                       STORY_APS_for_R5$`PANAS_X_Jov`)
summary(STORY_APS_for_R5.strep.Jov)

#this is for Dialister PANAS_X_Jov
STORY_APS_for_R5.Dialister.Jov <- glm.nb(STORY_APS_for_R5$`Dialister_AA_round_off` ~ 
                                           STORY_APS_for_R5$`PANAS_X_Jov`)
summary(STORY_APS_for_R5.Dialister.Jov)

#this is for Strep PANAS_X_Serenity
STORY_APS_for_R5.strep.Serenity <- glm.nb(STORY_APS_for_R5$`Strep_AA_round_off` ~ 
                                            STORY_APS_for_R5$`PANAS_X_Serenity`)
summary(STORY_APS_for_R5.strep.Serenity)

#this is for Dialister PANAS_X_Serenity 
STORY_APS_for_R5.Dialister.Serenity <- glm.nb(STORY_APS_for_R5$`Dialister_AA_round_off` ~ 
                                                STORY_APS_for_R5$`PANAS_X_Serenity`)
summary(STORY_APS_for_R5.Dialister.Serenity)