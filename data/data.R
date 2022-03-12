#############################################################
## Here we work on the data provided by Ipeadata
## gslindenmeyer@gmail.com - Nov 2021
#############################################################

##############################
# 0 - Load libraries
##############################
rm(list = ls())
##############################
library(ipeadatar)
library(readxl)
library(dplyr)
library(lubridate)
library(urca)
library(tidyr)
library(tseries)
##############################
# 1 - Source file
##############################
dataPath <- "C:/Users/guili/Documents/GitHub/Research_UFRGS/data/"
myfiles <- list.files(paste0(dataPath, "code_data"))
dataset <- read_excel("data/code_master2.xlsx",
  col_types = c(
    "text", "text", "text",
    "date", "date", "numeric", "text"
  )
)
source(paste0(dataPath, "helper_functions.R"))
##############################
# 2 - Start my code
##############################

## Acquiring the metadata from each CODE
metadados <- metadata(dataset$codigo)

## Acquiring the values from each CODE
data <- ipeadata(metadados$code)

## Converting data to wide format
df <- data %>%
  pivot_wider(names_from = "code") %>%
  select(-c(uname, tcode))

## Sorting by date ascending
df_sorted <- df[order(as.Date(df$date, format = "%Y/%m/%d")), ]

## Filtering by dates
df_filtered <- subset(df_sorted, date >= as.Date("1996-01-01") & date < as.Date("2020-01-01"))



## We want to know how to make the data set stationary using ADF Test ##

non_stationary <- df_filtered

## Identifying the columns
tipo <- c(4)
for (i in 2:ncol(non_stationary)) {
  print(i)
  if (sum(non_stationary[, i] < 0, na.rm = T) > 0) { ## if there is any negative value

    if (sum(non_stationary[, i] == 0, na.rm = T) > 0) { ## and if there is any zero value

      tipo <- append(tipo, 3) # negative and zero
    } else {
      if (sum(non_stationary[, i] < 0, na.rm = T) == length(non_stationary[, i])) {
        tipo <- append(tipo,5) # only negative values
      } else {
      tipo <- append(tipo, 1) }# some negative values
    }
  } else { ## no negative values

    if (sum(non_stationary[, i] == 0, na.rm = T) > 0) { ## no negative values but a zero
      tipo <- append(tipo, 2)
    } else {
      tipo <- append(tipo, 0) ## no negative values and no zero
    }
  }
}

## 0 - no neg nor zero, 1 - any neg values, 2 - only pos but zero, 3 - neg values and zero.
test <- c(4)
for (i in 2:(ncol(non_stationary))) {
  print(i)
  X <- na.exclude(as.matrix(non_stationary[, i]))

  k <- 0

  j <- 0

  status <- "non-stationary"

  while (status == "non-stationary") {
    adf_test <- summary(ur.df(X, "none", lags = 12))
    if (k == 2) {
      dale <- colnames(non_stationary)[i]
    }
    # adf.test(X)[4] <= 0.05

    if (adf_test@teststat[1] <= adf_test@cval[1, 2] && kpss.test(X, null = "T")$p.value >= 0.05) {
      status <- "stationary"

      if (j == 1) {
        k <- 0.5
      }
    } else {
      if (j == 0) {
        X <- preparacao(X, i)

        j <- 1
      } else {
        k <- k + 1

        if (tipo[i] == 0 | tipo[i] == 2 | tipo[i] == 3 | tipo[i]==1) {
          X <- diff(X)
        }
        if (tipo[i] == 5) {
          X <- cresc_discreto(X)
        }

        j <- j + 1
      }
    }


    #  print(k)
  }

  test <- append(test, k)
}


## Assessing transformation value
transformation <- c(-1)

for (i in 2:ncol(non_stationary)) {
  if (tipo[i] == 0 && test[i] == 0) {
    transformation[i] <- 4
  } else if (tipo[i] == 0 && test[i] == 1) {
    transformation[i] <- 5
  } else if (tipo[i] == 0 && test[i] == 2) {
    transformation[i] <- 6
  } else if (tipo[i] == 1 && test[i] == 0) {
    transformation[i] <- 1
  } else if (tipo[i] == 1 && test[i] == 1) {
    transformation[i] <- 2
  } else if (tipo[i] == 2 && test[i] == 0) {
    transformation[i] <- 1
  } else if (tipo[i] == 3 && test[i] == 0) {
    transformation[i] <- 1
  } else if (tipo[i] == 3 & test[i] == 1) {
    transformation[i] <- 2
  } else if (tipo[i] == 2 & test[i] == 1) {
    transformation[i] <- 2
  } else if (tipo[i] == 5 && test[i] == 1) {
    transformation[i] <- 7}
}


month <- t(as.matrix(month(df_filtered$date)))
year <- t(as.matrix(year(df_filtered$date)))
reference <- available_subjects()
metadados <- metadados %>% inner_join(reference)
is_na_2019 <- as.vector(which(sapply(as.data.frame(is.na(df_filtered)), sum) != 0))

dataset <- list(
  "data" = as.matrix(df_filtered[2:ncol(df_filtered)]),
  "transformation" = t(as.matrix(transformation[2:ncol(df_filtered)])),
  month = month, year = year, "metadados" = metadados,
  "na_2019" = is_na_2019
)
dataset$names <- names(dataset$data)
names(dataset$data) <- NULL

save(dataset, file = paste0(dataPath, "dataset.RData"))

# library("writexl")
write_xlsx(dataset$metadados, path = paste0(dataPath, "metadados.xlsx"))


