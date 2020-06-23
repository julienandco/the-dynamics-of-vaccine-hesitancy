#
#
# Estimate the parameters for the my model.
#
#  (a) For each election, estimate
#     - Zealot model (beta-distrib)
#     - All reinfrocement paras equal
#     - All parameters can be chosen independently.
# (b) Compare the models by the logl-likelihood-ration test
# (c) Compare empirical and theoretical distributions 
#     by the Kolmogorovv-Smirnov tests.
# 

#
# in case: define work directory
#setwd("C:/Users/ge69fup/Documents/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/data_fitting")
setwd("D:/Dokumente/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/data_fitting");

post <- function(nme){
  # remove blanks
  xx = strsplit(nme, " ");
  nme1 = paste(xx[[1]], collapse="");
  cat(nme1, "\n");
  postscript(file=nme1, pointsize=32, paper="special",
             width=8, height=9, horizontal=FALSE);
}

###############
# read data
#
#################
#vaccination_data = read.csv2('C:/Users/ge69fup/Documents/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/data_fitting/merged_impfdaten.csv', header=TRUE)
vaccination_data = read.csv2('D:/Dokumente/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/data_fitting/merged_impfdaten.csv',header=TRUE);


###############
# read tools
#
#################
source("my_parameter_estimation_reinforcement.R");


######################################################
# initialisation
######################################################

theta_hut.init = 0.5;
psi_hut.init  = 0.5;
ksi_hut.init = 0.2;
s_hut.init   = 100; 

##change this
dummy.i = 0;

info = data.frame(vaccination_data$Wert[1:400] / 100,vaccination_data$Inzidenz[1:400]);
colnames(info) = c("Wert", "Inzidenz");

#remove all values below 0.5 as we esteem them to be fixed voters
info = info[info$Wert >= 0.5,];

#renormalise vaccination values, such that 0.5 -> 0 and 1 -> 1
info$Wert = info$Wert * 2 - 1;

#make sure every datapoint gets an incidence value
#(datapoints with NA get replaced by dummy.i)
info$Inzidenz = replace(info$Inzidenz, dummy.i);

#check whether we kept our functionality
info$Inzidenz = info$Inzidenz * 0;
#myDataEsti = c();
#impfer = vaccination_data$Wert[1:400]; #last few ones are NA
#myDataEsti = impfer/100; 


#remove all values below 0.5 as we esteem them to be fixed voters
#myDataEsti = Filter(remove_first_half, myDataEsti);

#renormalise it, such that 0.5 -> 0 and 1 -> 1
#myDataEsti = myDataEsti * 2 - 1;


##change this
#i.value = numeric(length(myDataEsti));


#myDataInzi = c();
#inzidenz = vaccination_data$Inzidenz;

#make sure every datapoint gets an incidence value
#(datapoints with NA get replaced by dummy.i)
#myDataInzi = replace(inzidenz,dummy.i);
#myDataInzi = get_rid_of_zeroes(inzidenz);

#what do we want to do?
analyse = TRUE;
#analyse = FALSE;
produce.table = FALSE;
produce.figures = FALSE;



######################################################
# data fit
######################################################


if (analyse){
  # check optima
  s_hut.max=2500;
  
  res.tab = c();
  
  
  {
    #mmean = mean(myDataEsti);
    mmean = mean(info$Wert);
    
    
    hist(info$Wert, freq = FALSE, nclass=30, 
         main="full reinforcement",
         xlim=c(0,1),xlab = "percentage of pro-vaxxers");
    
    
    ###################################################################
    # first run: estimate the full reinforcement model
    ###################################################################
    cat("first run", "\n");
    ##para.ref sind diese Hut parameter in der reihenfolge: nu_hut,theta_hut,psi_hut,ksi_hut,s_hut
    para.ref = c(mmean, theta_hut.init, psi_hut.init,ksi_hut.init,s_hut.init);   # define init para
    lll.last = lll(para.ref);
    cat("lll.last: ", lll.last, "\n");
    curve(g(x), add=TRUE, col="blue", lwd=2);
    
    unrestricted.model     = TRUE;      # we aim at the full model
    unrestrict.theta_hut   = TRUE;
    opti.cyclic(para.ref);
    cat(lll.last, "\n");
    para.unrest = para.last;        # store the result
    lll.unrest  = lll.last;
    theta.res   = theta_hut*s_hut;
    
    # Kolmogorov-Smirnov test
    res.ks = ks.test(info$Wert, function(x){pReinforce(x)}); 
    
    # produce a figure with histogram and estimated distribution
    party.x = "pro-vaxx";
    
    
    ###################################################################
    # second run: reinforcement model, force equal reinforcement parameters 
    ###################################################################
    cat("second run","\n");
    para.ref = c(mmean, theta_hut.init, psi_hut.init,ksi_hut.init,s_hut.init);  # define init para
    lll.last = lll(para.ref);
    cat(lll.last, "\n");
    hist(info$Wert, freq = FALSE, nclass=30, 
         main="reinforcement with equal reinf. params",
         xlim=c(0,1),xlab = "percentage of pro-vaxxers");
    curve(g(x), add=TRUE, col="blue", lwd=2);
    
    unrestricted.model     = TRUE;      # we aim at the full model
    unrestrict.theta_hut   = FALSE;     # we want to keep equal parameters for reinforcement
    opti.cyclic(para.ref);
    cat(lll.last, "\n");
    para.halfRestrict = para.last;        # store the result
    lll.halfRestrict  = lll.last;
    theta.res   = theta_hut*s_hut;
    
    # Kolmogorov-Smirnov test
    res.halfRestrict.ks = ks.test(info$Wert, function(x){pReinforce(x)}); 
    
    
    ###################################################################
    # third run: estimate the zealot model (beta-distrib)
    ###################################################################
    cat("third run", "\n");
    unrestricted.model     = FALSE;           # we fix all reinforcement-paras
    unrestrict.theta_hut   = FALSE;
    para.ref = c(mmean, theta_hut.init, psi_hut.init,ksi_hut.init,s_hut.init);
    hist(info$Wert, freq = FALSE, nclass=30, 
         main="zealot model",
         xlim=c(0,1),xlab="percentage of pro-vaxxers");
    curve(g(x), add=TRUE, col="blue", lwd=2);
    
    opti.cyclic(para.ref);
    cat("lll.last: ", lll.last, "\n");
    para.restrict = para.last;            # store result
    lll.restrict  = lll.last;
    
    # kolmogorov-smirnov-test
    res.restric.ks = ks.test(info$Wert, function(x){pReinforce(x)});
    
    #in line muss noch das ergebnis für A und ksi_hut rein...
    line = c("run", party.x, 
             theta.res,
             para.unrest,
             lll.unrest,
             para.halfRestrict,
             lll.halfRestrict,
             para.restrict,
             lll.restrict, 
             pchisq(2*(lll.unrest-lll.restrict),df=2, 
                    lower.tail=FALSE), 
             pchisq(2*(lll.unrest-lll.halfRestrict),df=2, 
                    lower.tail=FALSE), 
             res.ks$p.value,
             res.halfRestrict.ks$p.value,
             res.restric.ks$p.value
    );
    
    res.tab = rbind(res.tab, line);
  }
  # names orient themselves ar the supplement II of the paper
  col.names = c(
    "run", "opinion", 
    "Theta1PlusTheta2.unr",  
    "nu.unr", "theta.hat.unr", "psi.unr", "ksi.unr","s.unr",
    "lll.unr",
    "nu.halfr", "theta.hat.halfr", "psi.halfr", "ksi.halfr","s.halfr",
    "lll.halfr",
    "nu.restr", "theta.hat.restr", "psi.restr", "ksi.restr","s.restr",
    "lll.restr",
    "lll.unr.rest", "lll.unr.halfr",
    "ks.unres", "ks.halfRestr", "ks.restr");
  dimnames(res.tab)[[2]] =col.names;
  
  save(file="datAnaMyModel_V1.rSave", res.tab);
  
}



if (produce.table){
  # produce a table
  load(file="datAnaMyModel_V1.rSave");
  sink(file="datMyModel.tex");
  cat(dimnames(res.tab)[[2]][c(1,2,4,5,6,7,8)]); cat(" theta1 ");cat(" theta2 "); cat(dimnames(res.tab)[[2]][c(22,24,26)]);
  cat("\n");
  nn = dim(res.tab)[1];
  for (j in 1:nn){
    cat(res.tab[j,1], " & ", res.tab[j,2], " & ");
    cat(res.tab[j,4], " & ", res.tab[j,5], " & ");
    cat(res.tab[j,6], " & ", res.tab[j,7], " & ", res.tab[j,8], "&");
    # theta_2 = h.s*h.theta*(1-h.psi)
    cat(as.double(res.tab[j,5])*as.double(res.tab[j,8])*as.double(res.tab[j,6]), " & ");
    cat(as.double(res.tab[j,5])*as.double(res.tab[j,8])*(1-as.double(res.tab[j,6])), " & ");
    cat(as.double(res.tab[j,22]), " & ", 
        as.double(res.tab[j,24]), " & ", 
        as.double(res.tab[j,26]), 
        "\\\\\n");
  }
  
  cat("Test hat.psi=0.5 versus free model (where is the reinforcement?)\n");
  cat(as.double(res.tab[j,23]), " \n");
  
  sink();
  
}



if (produce.figures){
  # produce figures
  impfer = vaccination_data$Wert[1:400];
  myDataEsti = impfer/100;       
  mmean <<- mean(info$Wert);
  
  
  post(paste("My_model",as.character(j),".eps",sep=""));
  hist(myDataEsti, freq = FALSE, 
       main=paste("Vaccinational behaviour"),
       xlim=c(0,1), xlab="amount of pro-vaxxers x", nclass=30);
  
  para.last = as.double(res.tab[j, 4:8]);
  
  lll.last = lll(para.last);
  cat(lll.last, "\n");  mmean <<- mean(info$Wert);
  
  curve(g(x), add=TRUE,  lwd=2);
  
  para.last = as.double(res.tab[j, 16:20]);
  lll.last = lll(para.last);
  cat(lll.last, "\n");
  curve(g(x), add=TRUE,  lwd=2, lty=2);
  dev.off();
  
  
  post(paste("My_model_2",as.character(j),".eps",sep=""));
  hist(myDataEsti, freq = FALSE, 
       main=paste("Vaccinational behaviour"),
       xlim=c(0,1), xlab="amount of pro-vaxxers x", nclass=30);
  
  para.last = as.double(res.tab[j, 4:8]);
  
  lll.last = lll(para.last);
  cat(lll.last, "\n");  mmean <<- mean(info$Wert);
  
  curve(g(x), add=TRUE,  lwd=2);
  
  para.last = as.double(res.tab[j, 10:14]);
  lll.last = lll(para.last);
  cat(lll.last, "\n");
  curve(g(x), add=TRUE,  lwd=2, col = "green");
  
  para.last = as.double(res.tab[j, 16:20]);
  lll.last = lll(para.last);
  cat(lll.last, "\n");
  curve(g(x), add=TRUE,  lwd=2, lty=2);
  dev.off();
  
  
}
