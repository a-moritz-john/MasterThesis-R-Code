# Load packages
library(haven)
library(lavaan)
library(lavaan)
library(dplyr)
library(tidyr)
library(semtree)
library(ggplot2)
library(reshape2)

# Importing Data
data <- read_sav("C:/Users/amori/Desktop/Masterarbeit/Daten/A045_MoMo_John.sav")

# Überprüfen, ob es Duplikate in der Spalte PSLNR gibt
# duplicated_rows <- duplicated(data$PLSNR)
# duplicates <- data[duplicated(data$PLSNR), ]
# table(data$PLSNR)[table(data$PLSNR) > 1]
# Keine Duplikate

# Überprüfen, ob es Missings in der Spalte PLSNR gibt
# sum(is.na(data$PLSNR))
# Keine Missings

# Datensatz kürze, sodass nur T1 bis T3 drin sind (dafür wird die Agevariable genutzt).
# Es werden aslo nur Rows behalten, die mind. 1 Wert auf Age1, Age2 and Age3 haben
data_T1_T3 <- data %>%
  filter(if_any(c(Alter_T1, Alter_T2, Alter_T3), ~ !is.na(.)))

# Spalten rasukürzen, die sich auf T4 beziehen
# column_names <- colnames(data_T1_T3)
# T4_column_names <-grep("4", column_names, value = TRUE)
# print(T4_column_names)
# e047z3 nur bezieht sich nicht auf T4
# T4_column_names <- setdiff(T4_column_names, "e047z3")
# Diese sind: "NETTOT4", "momo_nCT4", "Alter_T4", "AltersgruppeT4", "sex_T4", "sexanderesT4", "MOusbmiT4", "MOusbmi_khT4",
# "MOusbmi_iotfT4", "reakxT4", "weitmaxT4", "SP170T4", "PWC170_T4", "weitmaxPerzA_T4", "reakxPerzA_T4",
# "aaktiv7T4","LQkidsB_T4", "LQkidsB_r_T4", "LQkidsB_t_T4", "LQkidsGesund_Feld_T4"
# "LQkidsGesund_GFB_T4", LQkidsGesund_T4", "SD_ses_score_T4", "SD_ses_T4"
# print(T4_column_names)

# Spalten aus T4 rausgekürzt
# data_T1_T3.T4 <- data_T1_T3 %>%
#  select(-all_of(T4_column_names))


# MO_PIDNR_T3 aus dem Datensatz kürzen
# data_T1_T3.T4 <- data_T1_T3.T4 %>% select(-MO_PIDNR_T3)

# Datensatz kürzen, dass nur noch SDQ, aaktiv und die Covaraiten enthalten sind

subset_data <- data_T1_T3 %>%
  select(
    PLSNR, Alter_T1, Alter_T2, Alter_T3, 
    totsdq_e, VAEsdq, VAEsdq_T3,   # SDQ variables
    aaktiv7_T1, aaktiv7_T2, aaktiv7T3,  # Activity variables
    sexa, bula,               # Stable covariates
    SDEses, SDEsesT2, SDEses_T3,
    USzbmi, MOuszbmi, USzbmi_T3,
    migrant
  )

# Rename SDQ columns to reflect the timepoint
subset_data <- subset_data %>%
  rename(
    VAEsdq_T1 = totsdq_e,
    VAEsdq_T2 = VAEsdq,
    aaktiv7_T3 = aaktiv7T3
  )

# Reshape to long format
#data_long <- subset_data %>%
#  pivot_longer(
#    cols = c(VAEsdq_T1, VAEsdq_T2, VAEsdq_T3, aaktiv7_T1, aaktiv7_T2, aaktiv7_T3),
#    names_to = c("variable", "timepoint"), # Split into variable and timepoint
#    names_pattern = "(.*)_T?(\\d)"         # Extract base variable and timepoint
#  ) %>%
#mutate(timepoint = as.numeric(timepoint)) # Convert timepoint to numeric



# Reshape number 2:
long_data <- subset_data %>%
  pivot_longer(
    cols = starts_with("VAEsdq_") | starts_with("aaktiv7_"),
    names_to = c(".value", "timepoint"),  # Splits column names
    names_pattern = "(.*)_(T\\d+)"       # Extracts variable name and timepoint
  )



### Daten wurden in ein Long-format übertragen mit insgesamt 31212 Zeilen (entpsricht 3 Zeilen pro Person: 3x10404, da wir 2 variablen zu 3 zeitpunkten haben).

colnames(long_data)

# Save the long-format data to a CSV file
write.csv(long_data, "data_long_all_cov.csv", row.names = FALSE)

























































# Variablen Gruppieren bezüglich des Timepoints an dem sie gemessen wurden
# col_names_all <- colnames(data_T1_T3.T4)
# print(col_names_all)

# time1_vars <- c("Alter_T1", "AltersgruppeT1", "SDEses_score", "SDEses", "USpwc170", "USpulsmax", "USzbmi", 
#                "USperzbmi", "USbmi", "USbmi_kh", "USbmi_iotf", "totsdq_e", "totsdq_k", "aaktiv7_T1",
#                "USpwc170", "USpulsmax", "weitmax_T1", "weitmaxPerzA_T1", "reakx_T1", "reakxPerzA_T1",
#                "kw100_e", "pw100_e", "sel100_e", "fam100_e", "fre100_e", "sch100_e", "tot100_e", "famil_e",
#                "freun_e", "kw_e", "pw_e", "schsum_e", "selbst_e", "total_e", "kw100_k", "pw100_k", 
#                "sel100_k", "fam100_k", "fre100_k", "sch100_k", "tot100_k", "famil_k", "freun_k", "kw_k",
#                "pw_k", "schsum_k", "selbst_k", "total_k")

#time2_vars <- c("Alter_T2", "AltersgruppeT2", "MOusbmi", "MOusbmi_kh", "MOusbmi_iotf", "MOuszbmi", "MOusperzbmi",
#                "VAEsdq", "VAsdq", "aaktiv7_T2", "ausherzT2", "auswattT2", "weitmax_T2", "weitmaxPerzA_T2", 
#                "reakx_T2", "reakxPerzA_T2", "LQEkidsB", "LQEkidsB_r", "LQEkidsB_t", "LQkidsB", 
#                "LQkidsB_recalc", "LQkidsB_r", "LQkidsB_t", "SDEses_scoreT2", "SDEsesT2")

#time3_vars <- c("Alter_T3", "AltersgruppeT3", "USbmi_T3", "USzbmi_T3", "USperzbmi_T3", "USbmi_iotf_T3", 
#                "USbmiB_kh_T3", "VAEsdq_T3", "VAsdq_T3", "aaktiv7T3", "LQEzufrB1_T3", "LQEzufrB3_T3",
#                "ausherzT3", "auswattT3", "weitmaxT3", "weitmaxPerzA_T3", "reakxT3", "reakxPerzA_T3",
#                "LQkidsB_r_T3", "LQkidsB_t_T3", "LQEkidsB_T3", "LQEkidsB_r_T3", "LQEkidsB_t_T3",
#                "SDEses_score_T3", "SDEses_T3")

#stable_vars <- c("sex", "sexa", "bula", "SDEcasminmz", "SDEcasminm", "SDEcasminvz", "SDEcasminv",
#                 "SDEisced97mz", "SDEisced97meu", "SDEisced97vz", "SDEisced97veu", "migrant", "e005z1",
#                 "e007m", "e007v", "e008m", "e008v", "e010m1", "e017az", "e017bz", "e047z3", "e089m",
#                 "e089v", "e090m", "e090v", "k008", "k008m", "k008v", "e022", "still",
#                 "händig_T1", "Mbmi", "Vbmi", "mbmi_k", "vbmi_k" )


# Variablen, die für alle Zeitpunkte erhoben wurde:
# "sex", "sexa", "bula"

# Variablen, die als Covariaten zu T1 erhoben wurden, aber für alle Zeitpunkte gelten:
# "SDEcasminmz", "SDEcasminm", "SDEcasminvz"     "SDEcasminv"      "SDEisced97mz"    
# "SDEisced97meu", "SDEisced97vz", "SDEisced97veu", "migrant", "e005z1", "e007m", "e007v",
# "e008m", "e008v", "e010m1", "e017az", "e017bz", "e047z3", "e089m", "e089v", "e090m", "e090v",
# "k008", "k008m", "k008v", "e022", "still", "händig_T1", "Mbmi", "Vbmi", "mbmi_k", "vbmi_k"


# Variablen, die ausgelassen wurden: "PLSNR" -> ID-Variable

