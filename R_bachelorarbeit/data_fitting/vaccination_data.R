
library("ggpubr")
#vaccination_data = read.csv2('C:/Users/ge69fup/Documents/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/data_fitting/merged_impfdaten.csv', header=TRUE)
vaccination_data = read.csv2('D:/Dokumente/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/data_fitting/merged_impfdaten.csv',header=TRUE);

#hist(vaccination_data$Anzahl)
#hist(vaccination_data$Inzidenz)
#ggscatter(vaccination_data, x = "Anzahl", y = "Inzidenz", add = "reg.line", conf.int = TRUE,cor.coef = TRUE, cor.method = "pearson",xlab = "Impfrate", ylab = "Inzidenz in %");
#impfrate = vaccination_data$Wert / 100;


lul = vaccination_data$Inzidenz;

for(i in 1:length(lul)){
  if (is.na(lul[i])){
    cat('NA');
  }
}

#correlation.pearson = cor(impfrate,inzidenz,method='pearson');
#correlation.kendall = cor(impfrate,inzidenz,method='kendall');
#correlation.spearmn = cor(impfrate,inzidenz,method='spearman');

#plot(inzidenz,impfrate);

#res1 <- cor.test(impfrate,inzidenz,method='pearson');

#res2 <- cor.test(impfrate,inzidenz,method='kendall');

#res3 <- cor.test(impfrate,inzidenz,method='spearman');
