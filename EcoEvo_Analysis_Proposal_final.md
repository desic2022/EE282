# Project Analysis Proposal

## Happy Gut Happy Mind: Associations Between Gut Microbes and Positive Emotion



Broadly, my goal as a mind-microbiome researcher is to discover low-cost accessible ways in which the gut microbiome and psychological constructs may reciprocally optimize holistic wellness. To date, work assessing mind-microbiome relations has primarily focused on disease and dysfunction (see review by Cryan et.al, 2012) with scant attention given to links to positive psychosocial constructs. The current work will explore whether various types of dispositional (how they usually feel) positive affect (PA) are linked to microbial composition. This research will use an existing dataset that I collected containing all of the microbiome data and metadata needed to explore my research inquiries as described in detail below.
 
Participants self-reported levels of dispositional Serenity, Joviality, and Emotional Well-Being on a 5-point likert scale with ‘1’ indicating very little or no experience of the emotion to ‘5’ indicating intense or extreme experience of that emotion. Participants include 75 healthy women ages 25-55 years old. Additionally, participants were asked to report on diet, age, weight, and height which will be used |as covariates in statistical analyses. Next-generation pyrosequencing of the V1-V3 region of the 16S rRNA gene was used to classify microbial composition from fecal |samples and determine relative and absolute abundances of bacterial genera. This data |is contained in an OTU table available in a BIOM file and text document. 
 
We hypothesize that the endogenous microbiome, both globally and when examining specific taxa, are linked to distinct types of positive affect. I will conduct a series of permutational analysis of variance (PERMANOVA; Anderson, 2005) in R statistical software to assess beta diversity between PA groups. Specifically, I will investigate differences in microbial composition by performing a median split (using SPSS) creating high and low groups for each measure of PA. The PERMANOVA is an appropriate analysis for microbial data because it can assess numerous dependent variables and is robust to non-independent variables and zero-values. If the groups are significantly different, I will probe these findings in R using Random Forest (Liaw & Wiener, 2002). This method uses numerous decision trees to predict and classify bacterial genera based on PA group. Next, analysis of composition of microbiomes (ANCOM; Mandal et. al, 2015) will be used to identify the underlying compositional structure of the data while controlling for covariates (BMI, diet, and age) and the false discovery rate. Global composition will be visualized in |R using a nonparametric multidimensional scaling ordination using Bray-Curtis |dissimilarity to calculate distances between samples on the genus level. Finally, |Spearman's correlations will be conducted in SPSS to reveal possible associations |between Streptococcus, Dialister and each PA type. If significant associations are discovered, linear regressions (adjusting for BMI, diet, and age) will be calculated. 


## References

Anderson, M. J. (2005). Permutational multivariate analysis of variance. Department of Statistics, University of Auckland, Auckland, 26, 32-46.

Cryan,  J. F., & Dinan, T. G. (2012). Mind-altering microorganisms: the impact of the gut microbiota on brain and behaviour. Nature reviews neuroscience, 13(10), 701-712.

Liaw, A., & Wiener, M. (2002). Classification and regression by randomForest. R news, 2(3), 18-22.

Mandal, S., Van Treuren, W., White, R. A., Eggesbø, M., Knight, R., & Peddada, S. D. (2015). Analysis of composition of microbiomes: a novel method for studying microbial composition. Microbial ecology in health and disease, 26(1), 27663.
