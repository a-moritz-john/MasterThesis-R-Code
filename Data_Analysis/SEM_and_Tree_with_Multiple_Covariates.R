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
long_data <- read.csv("C:/Users/amori/Documents/R/Masterarbeit/data_long_all_cov.csv")

# Remove Trans people 
long_data <- long_data %>%
  filter(!(sexa %in% c(3, 4)))

# Create a new variable: age at the time of measurement my merging Alter_T1, Alter_T2, Alter_T3
long_data <- long_data %>%
  mutate(
    Age = case_when(
      timepoint == "T1" ~ Alter_T1,
      timepoint == "T2" ~ Alter_T2,
      timepoint == "T3" ~ Alter_T3
    )
  )

# Recode Federal states into East and West (States are ordered by their corresponding number in the dataset)

# West Germany: 
## 1. Schleswig-Holstein
## 2. Hamburg
## 3. Niedersachsen
## 4. Bremen
## 5. Nordrhein-Westfalen
## 6. Hessen
## 7. Rheinland-Pfalz
## 8. Baden-Württemberg
## 9. Bayern
## 10. Saarland

# East Germany:
## 11. Berlin (special case, but usually grouped with East)
## 12. Brandenburg
## 13. Mecklenburg-Vorpommern
## 14. Sachsen
## 15. Sachsen-Anhalt
## 16. Thüringen

# Codes for Berlin + new Bundesländer
east_states <- c(11, 12, 13, 14, 15, 16)

long_data$Region <- ifelse(is.na(long_data$bula), NA,
                           ifelse(long_data$bula %in% east_states, "East", "West"))

# Recode Socioeconomicstatus from all 3 Timepoints to be SES at the current measurement
long_data <- long_data %>%
  mutate(
    SES = case_when(
      timepoint == "T1" ~ SDEses,
      timepoint == "T2" ~ SDEsesT2,
      timepoint == "T3" ~ SDEses_T3
    )
  )

# Recode BMI from all 3 TImepoints to be BMI at the current measurement
long_data <- long_data %>%
  mutate(
    BMI = case_when(
      timepoint == "T1" ~ USzbmi,
      timepoint == "T2" ~ MOuszbmi,
      timepoint == "T3" ~ USzbmi_T3
    )
  )

# Cut the dataframe to only include:
# Timepoint, VAEsdq, aaktiv7, Age, sexa, bula, SES, BMI, migrant
long_data <- long_data %>%
  dplyr::select(timepoint, VAEsdq, aaktiv7, Age, sexa, Region, SES, BMI, migrant)

# Remove all rows that have missings on age -> individual did not participate at this measurement point
long_data <- long_data %>%
  filter(!is.na(Age))
# 15860 observations are left

# Remove all rows that have missing values on aaktiv and sdq -> only one of the variables missing is permissable
long_data <- long_data %>%
  filter(!is.na(VAEsdq) | !is.na(aaktiv7))
# 15839 observations are left -> same as other tree

# Count missings on the covariates

# Age
sum(is.na(long_data$Age))
# 0

# Sex
sum(is.na(long_data$sexa))
# 10

# Region
sum(is.na(long_data$Region))
# 4292

# Socioeconomic status
sum(is.na(long_data$SES))
# 2246

# BMI
sum(is.na(long_data$BMI))
# 6645

# Migrant
sum(is.na(long_data$migrant))
# 4335





# Remove Missings on sex
long_data <- long_data %>%
  filter(!is.na(sexa))
# 15829 observations remaining

# Remove Missings on Region
long_data <- long_data %>%
  filter(!is.na(Region))
# 11540 observations remaining

# Remove Missings on SES
long_data <- long_data %>%
  filter(!is.na(SES))
# 9332 observations remaining

# Remove Missings on BMI
long_data <- long_data %>%
  filter(!is.na(BMI))
# 8383 observations remaining

# Remove Missings on migrant
long_data <- long_data %>%
  filter(!is.na(migrant))
# 8347 observations remaining



# Turn Timepoint, sexa, Region, SES and migrant into factors
long_data$timepoint <- as.factor(long_data$timepoint)
long_data$sexa <- as.factor(long_data$sexa)
long_data$Region <- as.factor(long_data$Region)
long_data$SES <- as.factor(long_data$SES)
long_data$migrant <- as.factor(long_data$migrant)

# Rename Sex column
long_data <- long_data %>%
  rename(Sex = sexa)

# Rename sex from 1 to male and 2 to female
long_data <- long_data %>%
  mutate(Sex = factor(Sex,
                      levels = c(1, 2),
                      labels = c("Male", "Female")))


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



################ SEM ################################


# Define OpenMX model
# two manifest varaibles with freely estimated covariances, variances and means
# FIML as estimation method


v1 <- var(long_data$VAEsdq,  na.rm = TRUE)
v2 <- var(long_data$aaktiv7, na.rm = TRUE)

model_mx <- mxModel(
  "covariance_model",
  type = "RAM",
  manifestVars = c("VAEsdq","aaktiv7"),
  # variances: reasonable starts, >0 bounds
  mxPath(from = c("VAEsdq","aaktiv7"), arrows = 2,
         values = c(v1, v2), lbound = 1e-6,
         labels = c("covariance_model_S11", "covariance_model_S22")),
  # covariance: start near 0, keep away from ±1
  mxPath(from = "VAEsdq", to = "aaktiv7", arrows = 2,
         values = 0, lbound = -0.99, ubound = 0.99,
         labels = "covariance_model_S12"),
  # means
  mxPath(from = "one", to = c("VAEsdq","aaktiv7"), arrows = 1,
         labels = c("covariance_model_M11", "covariance_model_M12")),
  mxData(observed = long_data, type = "raw")
) 


# Estimate the model using mxTryHard -> otherwise there is an issue with the starting vaues
fit_mx <- mxTryHard(model_mx)

# Check the summary of the model
summary(fit_mx)


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


sem_tree <- semtree(
  model = fit_mx,            # SEM model fitted with openMX
  data = long_data,  # Data used for the tree
  control = control, # Control settings
  predictors = c("timepoint", "Age", "Sex", "Region", "SES", "BMI", "migrant")
)

summary(sem_tree)

plot(sem_tree)


############# Further examine the splits ################

nl <- getNodeList(sem_tree)

get_top_contrib <- function(node, k = 3) {
  pc <- node$result$par.contrib
  if (is.null(pc)) return(NULL)
  sort(pc, decreasing = TRUE)[1:min(k, length(pc))]
}

fmt_rule_value <- function(val) {
  # val can be numeric, character, factor, or a vector
  if (length(val) > 1) {
    paste(val, collapse = ", ")
  } else {
    as.character(val)
  }
}

top_per_node <- lapply(nl, get_top_contrib)

for (i in seq_along(top_per_node)) {
  if (!is.null(top_per_node[[i]])) {
    nd <- nl[[i]]
    # guard in case any component is missing
    split_var <- if (!is.null(nd$rule$name)) nd$rule$name else NA
    relation  <- if (!is.null(nd$rule$relation)) nd$rule$relation else ""
    split_val <- if (!is.null(nd$rule$value)) fmt_rule_value(nd$rule$value) else NA
    nid       <- if (!is.null(nd$node_id)) nd$node_id else i
    N         <- if (!is.null(nd$N)) nd$N else NA
    pval      <- if (!is.null(nd$p)) nd$p else NA
    
    cat(sprintf(
      "Node %s | N=%s | p=%.3g\nSplit: %s %s %s\n",
      nid, N, pval, split_var, relation, split_val
    ))
    
    print(top_per_node[[i]])
    cat("\n")
  }
}


###########################################
Node 1 | N=8347 | p=1.85e-118
Split: Age >= 10.8045
covariance_model_M12 covariance_model_S22 covariance_model_M11
482.56767             38.02726             28.78137

Node 2 | N=3361 | p=5.72e-37
Split: SES %in% 1, 2
covariance_model_M11 covariance_model_S11 covariance_model_S22
138.793778            42.853160             8.068921

Node 3 | N=850 | p=0.000171
Split: Sex %in% Male
covariance_model_M11 covariance_model_M12 covariance_model_S11
13.404147             6.466386             4.027047

Node 6 | N=2511 | p=1.99e-12
Split: SES %in% 1
covariance_model_M11 covariance_model_S11 covariance_model_S22
37.042359            20.413374             4.553991

Node 9 | N=4986 | p=7.19e-52
Split: SES %in% 1
covariance_model_M11 covariance_model_S11 covariance_model_S22
155.572500           110.310999             1.633222

Node 10 | N=4429 | p=3.38e-27
Split: Sex %in% Male
covariance_model_M12 covariance_model_M11 covariance_model_S11
85.629299            37.176683             7.621822

Node 13 | N=557 | p=0.000152
Split: timepoint %in% T1
covariance_model_S22 covariance_model_M11 covariance_model_M12
14.041176            10.002794             5.697404
###############################

############# Get parameter values for the terminal nodes
leaves <- getLeafs(sem_tree)
length(leaves)

for (node in leaves) {
  cat("\n=== Leaf Node", node$node_id, "===\n")
  print(node$params)
}


##########################
=== Leaf Node 4 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
0.64684272           0.01798078           0.96807771          -0.33182019 
covariance_model_M12 
0.24587450 

=== Leaf Node 5 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
0.7886818            0.0403052            0.8980613           -0.1135843 
covariance_model_M12 
0.4303280 

=== Leaf Node 7 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
0.93373513          -0.08285881           1.02699142           0.11069116 
covariance_model_M12 
0.31053877 

=== Leaf Node 8 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
1.20136324          -0.03484286           1.19402797           0.44057438 
covariance_model_M12 
0.25036058 

=== Leaf Node 11 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
0.88258082          -0.03939268           0.82538270          -0.18519582 
covariance_model_M12 
-0.32661717 

=== Leaf Node 12 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
0.99498141          -0.07051226           0.88843672          -0.01314266 
covariance_model_M12 
-0.06592130 

=== Leaf Node 14 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
1.27727918          -0.03720044           0.72642319           0.31273124 
covariance_model_M12 
-0.30690565 

=== Leaf Node 15 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
1.61902461          -0.06843737           1.15296463           0.47563872 
covariance_model_M12 
-0.11446067 
#################################################################################

################# Get sample size for leave nodes
for (node in leaves) {
  cat("Leaf Node", node$node_id, " | N =", node$N, "\n")
}

################################################################################
Leaf Node 4  | N = 433 
Leaf Node 5  | N = 417 
Leaf Node 7  | N = 2084 
Leaf Node 8  | N = 427 
Leaf Node 11  | N = 2202 
Leaf Node 12  | N = 2227 
Leaf Node 14  | N = 286 
Leaf Node 15  | N = 271 
################################################################################


##### Sanity check ######
# Compare with lavaan model

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

sem_tree_lav <- semtree(
  model = fit_lav,            # SEM model fitted with openMX
  data = long_data,  # Data used for the tree
  control = control, # Control settings
  predictors = c("timepoint", "Age", "sexa", "Region", "SES", "BMI", "migrant")
)

summary(sem_tree_lav)

plot(sem_tree_lav)




#################################################





############### Include focus parameter ####################

fp <- "covariance_model_S12"

sem_tree_fokus <- semtree(
  model = fit_mx,            # SEM model fitted with lavaan
  data = long_data,  # Data used for the tree
  control = control, # Control settings
  constraints = list(focus.parameters=fp),
  predictors = c("timepoint", "Age", "Sex", "Region", "SES", "BMI", "migrant")
)

summary(sem_tree_fokus)

plot(sem_tree_fokus)

