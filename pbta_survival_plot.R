# Install
install.packages("survminer")

# Load
library("survminer")

# Load Data
pbta.dat <- read.csv("pbta_clinical_data.csv")
pbta.dat$OS_days <- as.numeric(pbta.dat$OS_days)
pbta.dat$PFS_days <- as.numeric(pbta.dat$PFS_days)

# Drop GNT, Schwannoma, and Teratoma samples
pbta.dat <- pbta.dat[pbta.dat$short_histology!='GNT',]
pbta.dat <- pbta.dat[pbta.dat$short_histology!='Schwannoma',]
pbta.dat <- pbta.dat[pbta.dat$short_histology!='Teratoma',]

# Overall Survival curves
require("survival")
fit <- survfit(Surv(OS_days, OS_status_boolean) ~ short_histology, data = pbta.dat)

# Drawing curves
ggsurvplot(fit, 
           title="Overall Survival (PBTA Cohort)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)

# Progression-free Survival curves
require("survival")
fit <- survfit(Surv(PFS_days) ~ short_histology, data = pbta.dat)

# Drawing curves
ggsurvplot(fit, 
           title="Progression-free Survival (PBTA Cohort)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)


#######################################################################
# GNG Samples
gng.dat <- pbta.dat[pbta.dat$short_histology=='Ganglioglioma',]

# Fit survival curves
require("survival")
#fit <- survfit(Surv(OS_days, OS_event) ~ 1, data = dat)
#fit <- survfit(Surv(OS_days) ~ 1, data = gng.dat)
fit <- survfit(Surv(OS_days) ~ reported_gender, data = gng.dat)

# Drawing curves
ggsurvplot(fit, 
           title="Ganglioglioma Overall Survival (Male vs Female)",
           xlab="Days",
           #ylab="Survival Probability",
           censor=TRUE,
           legend="bottom",
           #palette = "#2E9FDF",
           conf.int = TRUE,
           pval=TRUE
           )

#######################################################################

# CNS Regions Overall Survival curves
require("survival")
fit <- survfit(Surv(OS_days, OS_status_boolean) ~ CNS_region, data = pbta.dat)
ggsurvplot(fit, 
           title="Overall Survival (CNS Region)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)

# Extent of Tumor Resection Overall Survival curves
require("survival")
fit <- survfit(Surv(OS_days, OS_status_boolean) ~ extent_of_tumor_resection, data = pbta.dat)
ggsurvplot(fit, 
           title="Overall Survival (Extent of Tumor Resection)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)

# Race Overall Survival curves
require("survival")
fit <- survfit(Surv(OS_days, OS_status_boolean) ~ race, data = pbta.dat)
ggsurvplot(fit, 
           title="Overall Survival (Race)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)

# Ethnicity Overall Survival curves
require("survival")
fit <- survfit(Surv(OS_days, OS_status_boolean) ~ ethnicity, data = pbta.dat)
ggsurvplot(fit, 
           title="Overall Survival (Ethnicity)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)

# Primary Site Overall Survival curves
require("survival")
fit <- survfit(Surv(OS_days, OS_status_boolean) ~ primary_site, data = pbta.dat)
ggsurvplot(fit, 
           title="Overall Survival (Primary Site)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)

# Tumor Descriptor Overall Survival curves
require("survival")
fit <- survfit(Surv(OS_days, OS_status_boolean) ~ tumor_descriptor, data = pbta.dat)
ggsurvplot(fit, 
           title="Overall Survival (Tumor Descriptor)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)

# Broad Histology Overall Survival curves
require("survival")
fit <- survfit(Surv(OS_days, OS_status_boolean) ~ broad_histology, data = pbta.dat)
ggsurvplot(fit, 
           title="Overall Survival (Broad Histology)",
           xlab="Days",
           legend="bottom",
           censor=TRUE,
           #conf.int = TRUE,
           pval=TRUE
)




