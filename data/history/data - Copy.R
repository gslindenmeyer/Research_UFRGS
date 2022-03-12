#############################################################
## Here we work on the data provided by Ipeadata
## gslindenmeyer@gmail.com - Nov 2021
#############################################################

##############################
# 0 - Load librairies
############################## library(readr)
rm(list = ls())
library(dplyr)
library(tibble)
library(hflights)
library(readr)
library(urca)
library(readxl)
library(tidyr)
##############################
# 1 - Source file
##############################
dataPath <- "C:/Users/guili/Documents/GitHub/Research_UFRGS/data/"
myfiles <- list.files(paste0(dataPath, "code_data"))
code_master <- read_csv(paste0(dataPath, "code_master.csv"), show_col_types = FALSE)
dataset <- read_excel("data/code_master2.xlsx", 
                           col_types = c("text", "text", "text", 
                                         "date", "date", "numeric"))

source(paste0(dataPath, "helper_functions.R"))
##############################
# 2 - Start my code
##############################



master_data <- list() # creating an empty list to put all data into
## Saving all the csv's into the same list
for (k in 1:length(myfiles)) {
  master_data[[k]] <- read_csv(paste("data/code_data/", myfiles[k], sep = ""), show_col_types = F)
}
save(master_data, file = paste0(dataPath, "master_data.RData"))


filtered_data <- list() # creating an empty list to put all data into
i <- 1
## Filtering the list to get only data that include 1995 and 2020
for (k in 1:length(myfiles)) {
  df <- tibble(as.data.frame(master_data[k]))
  test <- df %>% filter(DATE >= "1995-01-01" & DATE < "2021-01-01")
  if (c(test[1, 2] == 1995 && test[1, 3] == 1) && c(test[1, 4] == 1 && length(which(test$YEAR == 2020)) == 12)) {
    filtered_data[[i]] <- test
    i <- i + 1
  }
  print(k / length(myfiles))
}
save(filtered_data, file = paste0(dataPath, "filtered_data.RData"))


## If the previous data already have been run, then we only need to load the data
load(paste0(dataPath, "filtered_data.RData"))

## Which of the variables do not have 312 rows
index <- which(data_checker(filtered_data) != 312)
filtered_data_nodup <- filtered_data
filtered_data_nodup[index] <- NULL

## We want to remove NA

index <- which(na_checker(filtered_data_nodup) != 0)
filtered_data_nodup_nona <- filtered_data_nodup
filtered_data_nodup_nona[index] <- NULL

## We want to remove series that contains zeroes
index <- which(zero_checker(filtered_data_nodup_nona) != 0)
filtered_data_nodup_nona_no0 <- filtered_data_nodup_nona
filtered_data_nodup_nona_no0[index] <- NULL


save(filtered_data_nodup_nona_no0, file = paste0(dataPath, "filtered_data_3.RData"))

##
df <- filtered_data_nodup_nona_no0

### Now we want to build a proper data set
names <- getting_names(df)
values <- getting_values(df)
time <- df[[1]][1]


df1 <- data.frame(matrix(unlist(values), ncol = length(values), byrow = F))
colnames(df1) <- names
df1
save(df1, file = paste0(dataPath, "df1.RData"))


## We want to make the data set stationary using ADF Test####

non_stationary <- df1

## Deteccao de tipo de dado

tipo <- c()
for (i in 1:ncol(non_stationary)) {
  print(i)
  if (sum(non_stationary[, i] < 0) > 0) { ## if there is any negative value

    if (sum(non_stationary[, i] == 0) > 0) { ## and if there is any zero value

      tipo <- append(tipo, 3) # negative and zero
    } else {
      tipo <- append(tipo, 1) # only negative values
    }
  } else { ## no negative values

    if (sum(non_stationary[, i] == 0) > 0) { ## no negative values but a zero
      tipo <- append(tipo, 2)
    } else {
      tipo <- append(tipo, 0) ## no negative values and no zero
    }
  }
}
## 0 - no neg nor zero, 1 - any neg values, 2 - only pos but zero, 3 - neg values and zero.
test <- c()
for (i in 1:ncol(non_stationary)) {
  print(i)
  X <- non_stationary[, i]
  #X = preparacao(X,i)
  k <- 0

  j <- 0

  status <- "non-stationary"

  while (status == "non-stationary") {
    adf_test <- summary(ur.df(X, "trend"))

    # adf.test(X)[4] <= 0.05

    if (adf_test@teststat[1] <= adf_test@cval[1, 1]) {
      status <- "stationary"

      if (j == 1) {
        k <- 0.5
      }
    } else {
      if (j == 0) {
       X = preparacao(X,i)
        j <- 1
      } else {
        k <- k + 1

        if (tipo[i] == 0 | tipo[i] == 2) {
          X <- diff(X)
        }
        if (tipo[i] == 1 | tipo[i] == 3) {
          X <- cresc_discreto(X)
        }

        j <- j + 1
      }
    }


    print(k)
  }

  test <- append(test, k)
}

stationary_df = non_stationary
stationary_df = stationary_df[-1,]
for (i in 1:ncol(non_stationary)) {
  print(i)
  X = non_stationary[,i]
  
  if (test[i] == 0) {
    
    #X = X[-1:-2]
    X = X[-1]
  }
  
  if (test[i] == 0.5) {
    
    X = preparacao(X, i)
    
    #X = X[-1:-2]
    X = X[-1]
    
  }
  
  if (test[i] == 1) {
    
    X = preparacao(X, i)
    
    if (tipo[i] == 0 | tipo[i] == 2) {
      X = diff(X)
    }
    if (tipo[i] == 1 | tipo[i] == 3) {
      X = cresc_discreto(X)
    }
    
    #    X = X[-1]
  }
  
  if (test[i] == 2) {
    
    X = preparacao(X, i)
    
    if (tipo[i] == 0 | tipo[i] == 2) {
      X = diff(diff(X))
    }
    if (tipo[i] == 1 | tipo[i] == 3) {
      X = crescimento_discreto(cresc_discreto(X))
    }
    
  }
  
  stationary_df[,i] = X
  
}

#############################
library(ipeadatar)

## Acquiring the metadata from each CODE
metadados = metadata(dataset$Codigo)

## Acquiring the values from each CODE
data = ipeadata(metadados$code)

## Converting data to wide format
df = data %>% pivot_wider(names_from="code") %>% select(-c(uname,tcode))

## Sorting by date ascending
df_sorted = df[order(as.Date(df$date, format="%Y/%m/%d")),]

## Filtering by dates
df_filtered =  subset(df_sorted, date >= as.Date("1996-01-01") & date < as.Date("2020-01-01"))


