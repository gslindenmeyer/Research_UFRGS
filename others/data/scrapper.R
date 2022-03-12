install.packages("devtools")
devtools::install_github("gomesleduardo/ipeadatar")
library(ipeadatar)
library(knitr);library(tidyr);library(dplyr);library(DT);library(magrittr)
PIB_search<-as.data.frame(search_series(terms = c("PIB - paridade", "Mensal"), fields = c("name"),language = c("br")))
PIB_search %<>% dplyr::slice(1:500L)
datatable(PIB_search)

?search_series
