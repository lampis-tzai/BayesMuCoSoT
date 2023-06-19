setwd("C:/Users/ltzai/Desktop/PhD/BayesMuCoSoT")
library(readxl)
library(dplyr)
library(MASS)
library(matrixcalc)
library(Matrix)
library(matlib)

set.seed(2)
adoq_data <- read_excel("adoq colonnes.xls")
adoq_data = as.data.frame(adoq_data)


#coefficients
deg2rad <- function(deg) {deg * pi/180}
coef_data = data.frame(Surface = adoq_data$Surface)

for(h in 1:4){
  ampl = paste0('Ampl',h)
  phase = paste0('Phase',h)
  a_h = adoq_data[,ampl]*cos(deg2rad(adoq_data[,phase]))
  b_h = adoq_data[,ampl]*sin(deg2rad(adoq_data[,phase]))
  har_coef = data.frame(a_h,b_h)
  colnames(har_coef) = c(paste0('a_',h),paste0('b_',h))
  coef_data = cbind(coef_data,har_coef)
}

adoq_data = cbind(adoq_data[,1:4],scale(coef_data))

library(fastDummies)


adoq_data = dummy_cols(adoq_data, select_columns = c("Lettre"),remove_first_dummy = TRUE)
adoq_data['Lettre_1']=1

writer_data = adoq_data[(adoq_data$N==1),]
background_data = adoq_data[(adoq_data$N!=1),]
#background_data = adoq_data[(adoq_data$N %in% c(6,7,8)),]

questioned_data = data.frame()
suspect_data = data.frame()
for (c in 1:4){
  writer_data_c = writer_data[(writer_data$Lettre==c),]
  random_percentage = runif(1,0.35,0.65)
  smp_size <- floor(random_percentage * nrow(writer_data_c))
  suspect_ind <- sample(seq_len(nrow(writer_data_c)), size = smp_size)

  questioned_data = rbind(questioned_data,writer_data_c[suspect_ind, ])
  suspect_data = rbind(suspect_data, writer_data_c[-suspect_ind, ])
}

y = c("Surface","a_1","b_1","a_2","b_2","a_3","b_3","a_4","b_4")
x = c( "Lettre_1","Lettre_2","Lettre_3","Lettre_4")
background_data_id = "N"

known_data = suspect_data

devtools::load_all()
BayesMuCoSoT_fit(y,x,questioned_data,known_data,background_data,background_data_id)
