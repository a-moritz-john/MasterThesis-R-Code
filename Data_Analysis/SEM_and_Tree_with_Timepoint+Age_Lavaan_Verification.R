# Load packages
library(haven)
library(lavaan)
library(lavaan)
library(dplyr)
library(tidyr)
library(semtree)
library(ggplot2)
library(reshape2)
library(devtools)

###### Data preparation ######

# Importing Data
long_data <- read.csv("C:/Users/amori/Documents/R/Masterarbeit/data_long.csv")

# Cut data to only include timepoint, age, aaktiv and sdq
long_data <- long_data %>%
  select(timepoint, VAEsdq, aaktiv7, Alter_T1, Alter_T2, Alter_T3)

# Create a new variable: age at the time of measurement my merging Alter_T1, Alter_T2, Alter_T3
long_data <- long_data %>%
  mutate(
    Alter = case_when(
      timepoint == "T1" ~ Alter_T1,
      timepoint == "T2" ~ Alter_T2,
      timepoint == "T3" ~ Alter_T3
    )
  )

# Cut Alter_T1, Alter_T2, Alter_T3 from the dataset
long_data <- long_data %>%
  select(timepoint, VAEsdq, aaktiv7, Alter)

# Remove all rows that have missings on age -> individual did not participate at this measurement point
long_data <- long_data %>%
  filter(!is.na(Alter))
# 15865 observations are left

sum(is.na(long_data$Alter))
# No missing values
sum(is.na(long_data$timepoint))
# No missing values

# Remove all rows that have missing values on aaktiv and sdq -> only one of the variables missing is permissable
long_data <- long_data %>%
  filter(!is.na(VAEsdq) | !is.na(aaktiv7))
# 15844 observations are left

# Standardize sdq and aaktiv
long_data <- long_data %>%
  mutate(
    VAEsdq = as.numeric(scale(VAEsdq)),
    aaktiv7 = as.numeric(scale(aaktiv7))
  )

# Check standardization
mean(long_data$VAEsdq, na.rm = TRUE)
# very close to 0 -> worked
sd(long_data$VAEsdq, na.rm = TRUE)
# 1 -> worked

mean(long_data$aaktiv7, na.rm = TRUE)
# very close to 0 -> worked
sd(long_data$aaktiv7, na.rm = TRUE)
# 1 -> worked

# Change timepoint into a factor
long_data$timepoint <- as.factor(long_data$timepoint)

################################################




################# SEM ################################



# Define lavaan model 
model_lav <- '
  VAEsdq ~~ aaktiv7
  VAEsdq ~~ VAEsdq
  aaktiv7 ~~ aaktiv7
  VAEsdq ~ 1
  aaktiv7 ~ 1
'

# Fit the lavaan model
fit_lav <- sem(model_lav, data = long_data, meanstructure = TRUE, missing = "fiml")

# Check the summary of the model
summary(fit_lav)


############################################






############## SEM Tree #############################

# Load custom version of semtree
if ("semtree" %in% (.packages())) {
  detach("package:semtree", unload = TRUE)
}
load_all("C:/Users/amori/Desktop/Masterarbeit/semtree_Update/semtree")

# Define control settings for the SEM tree
# Use score method 
control <- semtree.control(
  method = "score",
  min.N = 50,
  alpha = 0.05,
  max.depth = 3
)


sem_tree_lav <- semtree(
  model = fit_lav,            # SEM model fitted with openMX
  data = long_data,  # Data used for the tree
  control = control, # Control settings
  predictors = c("timepoint", "Alter")
)

summary(sem_tree_lav)

plot(sem_tree_lav)


