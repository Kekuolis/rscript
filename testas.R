# instaliuojamos bibliotekos, jei dar nėra instaliuotos
# install.packages(c("psych", "caret", "corrplot"))

library(psych)
library(caret)
library(dplyr)
library(ggplot2)
library(corrplot)

# įkeliamas duomenų failas
ekg_data <- read.csv("./EKG_pupsniu_analize.csv", na.strings=c(""))

# Pirmasis stulpelis suskaidomas į daug stulpelių
ekg_data_split <- data.frame(do.call('rbind', strsplit(as.character(ekg_data[,1]), ';')))

# Priskiriami stulpelių pavadinimai pagal požymių reikšmes
colnames(ekg_data_split) <- c('RR_l_0', 'RR_l_0/RR_l_1', 'RR_l_1', 'RR_l_1/RR_l_2', 'RR_l_2', 'RR_l_2/RR_l_3', 
                              'RR_l_3', 'RR_l_3/RR_l_4', 'RR_r_0', 'RR_r_0/RR_r_1', 'RR_r_1', 'RR_r_1/RR_r_2', 
                              'RR_r_2', 'RR_r_2/RR_r_3', 'RR_r_3', 'RR_r_3/RR_r_4', 'seq_size', 'signal_mean', 
                              'signal_std', 'wl_side', 'wr_side', 'P_val', 'Q_val', 'R_val', 'S_val', 'T_val', 
                              'P_pos', 'Q_pos', 'R_pos', 'S_pos', 'T_pos', 'label')

# konvertuojami stulpeliai į skaičiaus formatus
ekg_data_split[] <- lapply(ekg_data_split, function(x) as.numeric(as.character(x)))

# Padalinti duomenis į tris dalis, po 1000 kiekviename, arba mažiau paskutiniame rinkinyje
data_splits <- split(ekg_data_split, rep(1:3, each=1000, length.out=nrow(ekg_data_split)))
# Now you have three datasets: set1, set2, and set3

# Analizės funkcija, kuri bus naudojama trims rinkiniams
analyze_subset <- function(ekg_data_split, split_name) {
  
  ### 1. Aprašomoji statistika
  selected_features <- ekg_data_split[, c('RR_l_0', 'RR_r_0', 'RR_l_1/RR_l_2', 'signal_mean', 'signal_std', 'R_val')]
  
  summary(selected_features)
  
  desc_stats <- describe(selected_features)
  variance <- apply(selected_features, 2, var, na.rm = TRUE)
  desc_stats$variance <- variance
  
  print(desc_stats)
  
  # dispersija
  print(desc_stats$variance)
  
  ### 2. Užpildomos praleistos reikšmės naudojant mediana
  ekg_data_split <- ekg_data_split %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
  
  ### 3. Taškų atsiskyrėlių radimas ir pašalinimas iš duomenų aibės
  find_outliers <- function(data_column) {
    Q1 <- quantile(data_column, 0.25, na.rm = TRUE)
    Q3 <- quantile(data_column, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    return(data_column < lower_bound | data_column > upper_bound)
  }
  
  ekg_data_split$IsOutlier <- apply(selected_features, 1, function(row) {
    return(any(find_outliers(row)))
  })
  
  ekg_data_clean <- ekg_data_split[!ekg_data_split$IsOutlier, ]
  
  # Palyginimas prieš taškų atsiskyrėlių pašalinimą ir po to
  before_outliers <- describe(selected_features)  
  after_outliers <- describe(ekg_data_clean[, c('RR_l_0', 'RR_r_0', 'RR_l_1/RR_l_2', 'signal_mean', 'signal_std', 'R_val')])
  
  print("before outliers")
  summary(selected_features)
  print(before_outliers)
  
  print("after outliers")
  summary(ekg_data_clean[, c('RR_l_0', 'RR_r_0', 'RR_l_1/RR_l_2', 'signal_mean', 'signal_std', 'R_val')])
  print(after_outliers)
  
  ### 4. Normalizacija
  ekg_data_standardized <- scale(ekg_data_split[, c('RR_l_0', 'RR_r_0', 'RR_l_1/RR_l_2', 'signal_mean', 'signal_std', 'R_val')])
  
  # Dynamic plot titles based on the split_name
  par(mfrow = c(2, 3))  # Arrange plots in 2 rows, 3 columns
  boxplot(ekg_data_split$RR_l_0, main = paste("RR_l_0 with outliers -", split_name), xlab = "RR_l_0")
  boxplot(ekg_data_split$RR_r_0, main = paste("RR_r_0 with outliers -", split_name), xlab = "RR_r_0")
  boxplot(ekg_data_split$`RR_l_1/RR_l_2`, main = paste("RR_l_1/RR_l_2 with outliers -", split_name), xlab = "RR_l_1/RR_l_2")
  boxplot(ekg_data_split$signal_mean, main = paste("Signal Mean with outliers -", split_name), xlab = "Signal Mean")
  boxplot(ekg_data_split$signal_std, main = paste("Signal Std with outliers -", split_name), xlab = "Signal Std")
  boxplot(ekg_data_split$R_val, main = paste("R_val with outliers -", split_name), xlab = "R_val")

  # Dynamic plot titles based on the split_name
  par(mfrow = c(2, 3))  # Arrange plots in 2 rows, 3 columns

  # Replace boxplots with histograms
  hist(ekg_data_split$RR_l_0, main = paste("Histogram of RR_l_0 with outliers -", split_name), xlab = "RR_l_0", col = "lightblue", breaks = 30)
  hist(ekg_data_split$RR_r_0, main = paste("Histogram of RR_r_0 with outliers -", split_name), xlab = "RR_r_0", col = "lightgreen", breaks = 30)
  hist(ekg_data_split$`RR_l_1/RR_l_2`, main = paste("Histogram of RR_l_1/RR_l_2 with outliers -", split_name), xlab = "RR_l_1/RR_l_2", col = "lightcoral", breaks = 30)
  hist(ekg_data_split$signal_mean, main = paste("Histogram of Signal Mean with outliers -", split_name), xlab = "Signal Mean", col = "lightyellow", breaks = 30)
  hist(ekg_data_split$signal_std, main = paste("Histogram of Signal Std with outliers -", split_name), xlab = "Signal Std", col = "lightcyan", breaks = 30)
  hist(ekg_data_split$R_val, main = paste("Histogram of R_val with outliers -", split_name), xlab = "R_val", col = "lightpink", breaks = 30)
  
  par(mfrow = c(2, 3))  # Arrange plots in 2 rows, 3 columns

  # Replace boxplots with histograms
  boxplot(ekg_data_split$RR_l_0, main = paste("RR_l_0 norm outliers -", split_name), xlab = "RR_l_0")
  boxplot(ekg_data_split$RR_r_0, main = paste("RR_r_0 norm outliers -", split_name), xlab = "RR_r_0")
  boxplot(ekg_data_split$`RR_l_1/RR_l_2`, main = paste("RR_l_1/RR_l_2 norm outliers -", split_name), xlab = "RR_l_1/RR_l_2")
  boxplot(ekg_data_split$signal_mean, main = paste("Signal Mean norm  outliers -", split_name), xlab = "Signal Mean")
  boxplot(ekg_data_split$signal_std, main = paste("Signal Std norm outliers -", split_name), xlab = "Signal Std")
  boxplot(ekg_data_split$R_val, main = paste("R_val norm outliers -", split_name), xlab = "R_val")

  ### 5. Min-max normalizacija
  minmax_preprocess <- preProcess(ekg_data_clean[, c('RR_l_0', 'RR_r_0', 'RR_l_1/RR_l_2', 'signal_mean', 'signal_std', 'R_val')], method = "range")
  ekg_data_minmax <- predict(minmax_preprocess, ekg_data_clean)
  
  ### 6. Požymių tarpusavio priklausomybė
  # Koreliacijos matrica
  correlation_matrix <- cor(ekg_data_clean[, c('RR_l_0', 'RR_r_0', 'RR_l_1/RR_l_2', 'signal_mean', 'signal_std', 'R_val')], use = "complete.obs")
  
  print(correlation_matrix)
  
  # Atvaizduojama koreliacijos matrica
  corrplot(correlation_matrix, method = "circle")
  
  return(ekg_data_clean)
}

# Taikyti funkciją kiekvienam rinkiniui ir perduoti atitinkamus pavadinimus
split_names <- c("1", "2", "3")

cleaned_sets <- lapply(1:length(data_splits), function(i) {
  analyze_subset(data_splits[[i]], split_names[i])
})

# Išsaugoti duomenis po apdorojimo
for (i in 1:length(cleaned_sets)) {
  write.csv(cleaned_sets[[i]], paste0("processed_ekg_data_set_", i, ".csv"), row.names=FALSE)
}
