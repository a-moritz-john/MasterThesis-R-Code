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

# Remove Trans people 
long_data <- long_data %>%
  filter(!(sexa %in% c(3, 4)))


# Cut data to only include timepoint, age, aaktiv and sdq
long_data <- long_data %>%
  dplyr:::select(timepoint, VAEsdq, aaktiv7, Alter_T1, Alter_T2, Alter_T3)

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
  dplyr::select(timepoint, VAEsdq, aaktiv7, Alter)

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
# 15839 observations are left

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


# Define OpenMX model
# two manifest varaibles with freely estimated covariances, variances and means
# FIML as estimation method
# model_mx <- mxModel(manifestVars = c("VAEsdq", "aaktiv7"), type = "RAM",
#                mxPath(from = c("VAEsdq", "aaktiv7"), arrows = 2, values = 1, 
#                       labels = c("var_SDQ", "var_aaktiv")),
#                mxPath(from = "VAEsdq", to = "aaktiv7", arrows = 2, values = 1, labels = "cov"),
#                mxPath(from = "one", to = c("VAEsdq", "aaktiv7"),  arrows = 1, 
#                       labels = c("m_VAEsdq", "m_aaktiv7")),
#                mxData(observed = long_data, type = "raw"))


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

## Compute z-values for the estimates ##

# Variance of SDQ
1 / 0.012
# 83.33333

# Covariance of SDQ and PA
- 0.036 / 0.009
# -4.121292

# Variance of PA
1 / 0.011
# 0.90909

# Mean of SDQ
0.001 / 0.009
# 0.1111111

# Mean of PA
0 / 0.008
# 0


# get the corresponding p-value to the z-value
z_value <- -4
pnorm(z_value)
# 3.167124e-05
# 0.00003

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
  predictors = c("timepoint", "Alter")
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

################################################################
Node 1 | N=15839 | p=4.74e-254
Split: Alter >= 10.394250513347
covariance_model_M12 covariance_model_S22 covariance_model_M11 
1065.45991             42.51263             37.83854 

Node 2 | N=5809 | p=1.28e-19
Split: Alter >= 6.65160848733744
covariance_model_S11 covariance_model_M12 covariance_model_S22 
54.061541            34.095933             5.210775 

Node 3 | N=2399 | p=0.000157
Split: timepoint %in% T1, T2
covariance_model_M11 covariance_model_S22 covariance_model_S11 
23.847951             5.266602             3.186764 

Node 6 | N=3410 | p=3.65e-10
Split: timepoint %in% T2
covariance_model_S22 covariance_model_M11 covariance_model_M12 
30.583906            23.786819             5.165537 

Node 9 | N=10030 | p=6.14e-52
Split: Alter >= 13.6112970568104
covariance_model_M12 covariance_model_M11 covariance_model_S11 
116.26414            109.41031             32.43653 

Node 10 | N=3372 | p=9.49e-13
Split: timepoint %in% T1
covariance_model_S22 covariance_model_M11 covariance_model_S11 
33.55790             19.57717             12.14776 

Node 13 | N=6658 | p=2.87e-12
Split: Alter >= 15.763175906913
covariance_model_M11 covariance_model_S11 covariance_model_S12 
35.704884            24.074879             5.630444 
#####################################################################


############# Get parameter values for the terminal nodes
leaves <- getLeafs(sem_tree)
length(leaves)

for (node in leaves) {
  cat("\n=== Leaf Node", node$node_id, "===\n")
  print(node$params)
}

###################################################################
=== Leaf Node 4 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
0.79145873          -0.01112737           0.94984339          -0.12951570 
covariance_model_M12 
0.40928887 

=== Leaf Node 5 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
0.76452514          -0.02134539           1.03497508           0.08925023 
covariance_model_M12 
0.45661804 

=== Leaf Node 7 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
1.07355311          -0.08279814           1.01086998           0.02201621 
covariance_model_M12 
0.26673128 

=== Leaf Node 8 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
0.96454908          -0.02900982           0.80276697           0.14112213 
covariance_model_M12 
0.32405429 

=== Leaf Node 11 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
1.10646286          -0.10897157           0.79144117           0.05123397 
covariance_model_M12 
-0.03790109 

=== Leaf Node 12 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
1.190941541         -0.005000319          1.068236325          0.193306419 
covariance_model_M12 
-0.102296898 

=== Leaf Node 14 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
1.04662210          -0.02237479           0.85887997          -0.05086227 
covariance_model_M12 
-0.26154133 

=== Leaf Node 15 ===
  covariance_model_S11 covariance_model_S12 covariance_model_S22 covariance_model_M11 
0.86606608          -0.08208676           0.93542891          -0.21258930 
covariance_model_M12 
-0.27391189 

###################################################################


################# Get sample size for leave nodes
for (node in leaves) {
  cat("Leaf Node", node$node_id, " | N =", node$N, "\n")
}

###############################################################
Leaf Node 4  | N = 513 
Leaf Node 5  | N = 1886 
Leaf Node 7  | N = 2227 
Leaf Node 8  | N = 1183 
Leaf Node 11  | N = 2425 
Leaf Node 12  | N = 947 
Leaf Node 14  | N = 2313 
Leaf Node 15  | N = 4345 
###############################################################


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
  predictors = c("timepoint", "Alter")
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
  predictors = c("timepoint", "Alter")
)

summary(sem_tree_fokus)

plot(sem_tree_fokus)




#### Tree does not split when using the covariance as a focusparameter ####























# Testing
data <- read_sav("C:/Users/amori/Desktop/Masterarbeit/Daten/A045_MoMo_John.sav")

data <- data %>%
  select(PLSNR, aaktiv7_T1, aaktiv7_T2, aaktiv7T3, VAEsdq, totsdq_e,  VAEsdq_T3, Alter_T1, Alter_T2, Alter_T3)

# Es scheint so als sei Alter doch nicht immer erfasst worden -> oder falsche Eingabe?
