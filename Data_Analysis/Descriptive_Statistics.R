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
library(psych)

###### Data preparation ####### Importing Data
long_data <- read.csv("C:/Users/amori/Documents/R/Masterarbeit/data_long.csv")

# Remove Trans people 
long_data <- long_data %>%
  filter(!(sexa %in% c(3, 4)))

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
# 15860 observations are left

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

############ Sample Size and Participation across Wave
T1 <- long_data[long_data$timepoint == "T1",]
T2 <- long_data[long_data$timepoint == "T2",]
T3 <- long_data[long_data$timepoint == "T3",]

describe(T1$Alter)

describe(T2$Alter)

describe(T3$Alter)

##########################


############### Descriptive Statistics for sample 1 ###############

# Age
describe(long_data$Alter)
# Histogram for Age
ggplot(long_data, aes(x = Alter)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 1) +
  labs(title = "Age at the Time of Measurement",
       x = "Age (years)",
       y = "Count") +
  theme_minimal(base_size = 14)


# Age per Timepoint
ggplot(long_data, aes(x = Alter, fill = timepoint)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.6) +
  facet_wrap(~timepoint, ncol = 1) +
  labs(
    title = "Age Distribution Across Different Timepoints",
    x = "Age (years)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# SDQ
# unstandardized values were used here
describe(long_data$VAEsdq)
ggplot(long_data, aes(x = VAEsdq)) +
  geom_histogram(binwidth = 1, fill = "darkred", color = "white", alpha = 1) +
  labs(title = "SDQ - Raw Scores",
       x = "SDQ scores",
       y = "Count") +
  theme_minimal(base_size = 14)

# aaktiv
describe(long_data$aaktiv7)
ggplot(long_data, aes(x = aaktiv7)) +
  geom_histogram(binwidth = 1, fill = "darkgreen", color = "white", alpha = 1) +
  labs(title = "PA - Raw Scores",
       x = "PA amount",
       y = "Count") +
  theme_minimal(base_size = 14)


################ Descriptive Statistics for sample 2 ###############

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
  select(timepoint, VAEsdq, aaktiv7, Age, sexa, Region, SES, BMI, migrant)

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
# 6647

# Migrant
sum(is.na(long_data$migrant))
# 4335



# Turn Timepoint, sexa, Region, SES and migrant into factors
long_data$timepoint <- as.factor(long_data$timepoint)
long_data$sexa <- as.factor(long_data$sexa)
long_data$Region <- as.factor(long_data$Region)
long_data$SES <- as.factor(long_data$SES)
long_data$migrant <- as.factor(long_data$migrant)


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

# Sex
table(long_data$sexa)


  
long_data <- long_data %>%
  mutate(
    sexa = factor(case_when(
      sexa == 1 ~ "Male",
      sexa == 2 ~ "Female",
      TRUE      ~ NA_character_
    ), levels = c("Male", "Female"))
)


ggplot(long_data, aes(x = sexa)) +
  geom_bar(fill = "darkmagenta") +
  labs(
    title = "Distribution of Sex Categories in the Sample",
    x = "Sex",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

ggplot(long_data, aes(x = timepoint, fill = sexa)) +
  geom_bar(position = "fill") +  # "fill" makes it proportional (0–1)
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Distribution of Sex by Timepoint",
       x = "Timepoint",
       y = "Proportion",
       fill = "Sex") +
  theme_minimal(base_size = 14)

ggplot(long_data, aes(x = factor(timepoint), fill = sexa)) +
  geom_bar(position = "dodge") +  # "dodge" puts categories side by side
  labs(title = "Sex Distribution by Timepoint",
       x = "Timepoint",
       y = "Count",
       fill = "Sex") +
  theme_minimal(base_size = 14)

ggplot(long_data, aes(x = factor(timepoint), fill = sexa)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("Male" = "#1f77b4",         # blue
                               "Female" = "#e377c2"     # pink
  )) +
  labs(title = "Sex Distribution by Timepoint",
       x = "Timepoint",
       y = "Count",
       fill = "Sex") +
  theme_minimal(base_size = 14)

# Region
table(long_data$Region)

ggplot(long_data, aes(x = factor(timepoint), fill = Region)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("East" = "firebrick4",
                               "West" = "dodgerblue4")) +
  labs(title = "Region Distribution by Timepoint",
       x = "Timepoint",
       y = "Count",
       fill = "Region") +
  theme_minimal(base_size = 14)

# SES
table(long_data$SES)

long_data <- long_data %>%
  mutate(
    SES = factor(case_when(
      SES == 1 ~ "Low",
      SES == 2 ~ "Middle",
      SES == 3 ~ "High",
      TRUE      ~ NA_character_
    ), levels = c("Low", "Middle", "High"))
  )

ggplot(long_data, aes(x = factor(timepoint), fill = SES)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("Low" = "orangered3",
                               "Middle" = "goldenrod3",
                               "High" = "chartreuse3")) +
  labs(title = "SES Distribution by Timepoint",
       x = "Timepoint",
       y = "Count",
       fill = "SES") +
  theme_minimal(base_size = 14)

# BMI
describe(long_data$BMI)

ggplot(long_data, aes(x = BMI, fill = timepoint)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.6) +
  facet_wrap(~timepoint, ncol = 1) +
  labs(
    title = "BMI Distribution Across Different Timepoints",
    x = "BMI (z-standardized)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Migrant
table(long_data$migrant)

long_data <- long_data %>%
  mutate(
    migrant = factor(case_when(
      migrant == 1 ~ "Yes",
      migrant == 2 ~ "No",
      TRUE      ~ NA_character_
    ), levels = c("Yes", "No"))
  )

ggplot(long_data, aes(x = factor(timepoint), fill = migrant)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("Yes" = "darkolivegreen4",
                               "No" = "tomato4")) +
  labs(title = "Migrant Distribution by Timepoint",
       x = "Timepoint",
       y = "Count",
       fill = "Migrant") +
  theme_minimal(base_size = 14)





############# Assumptions ###############
library(MVN)

# 1. Mardia’s multivariate normality test
mvn(long_data[, c("VAEsdq", "aaktiv7")], mvn_test = "mardia")

# 2. QQ-plot for Mahalanobis distances
mvn(data = long_data[, c("VAEsdq", "aaktiv7")], mvn_test = "mardia", univariateTest = "SW", multivariatePlot = "qq")

# 3. Scatterplot with ellipse
ggplot(long_data, aes(x = aaktiv7, y = VAEsdq)) +
  geom_point(alpha = 0.3) +
  stat_ellipse(level = 0.95, color = "red", size = 1) +
  labs(
    title = "Bivariate distribution of SDQ and PA with 95% normal ellipse",
    x = "Physical Activity (raw score)",
    y = "SDQ Total Difficulties (raw score)"
  ) +
  theme_minimal(base_size = 14)




######### 3D Plot of the bivariate distribution #############
library(MASS)

# Suppose your variables are x and y in a data frame df
df <- data.frame(long_data$VAEsdq, long_data$aaktiv7)

# Remove rows where either x or y is NA
df_clean <- na.omit(df)

# Extract cleaned vectors
x <- df_clean$long_data.VAEsdq
y <- df_clean$long_data.aaktiv7


# Kernel density estimate on a grid (adjust = bandwidth factor; n = grid size)
dens <- kde2d(x, y, n = 100, h = c(bandwidth.nrd(x), bandwidth.nrd(y))*1.2)

# Pretty colors for height
z <- dens$z
zcols <- colorRampPalette(c("#e6f2ff","#99ccff","#4da3ff","#0066cc","#003d80"))(50)
zfacet <- cut(z, breaks = 50, include.lowest = TRUE)

persp(dens$x, dens$y, z,
      theta = 35, phi = 30, expand = 0.6,
      col = zcols[zfacet],
      shade = 0.3, border = NA,
      xlab = "SDQ", ylab = "PA", zlab = "Density",
      main = "3D Kernel Density Surface")
