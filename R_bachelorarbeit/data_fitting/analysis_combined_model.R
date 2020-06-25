#
# Estimate the parameters for the combined model.
#
#  (a) For each election, estimate
#     - Zealot model (beta-distrib)
#     - All reinforcement paras equal
#     - All parameters can be chosen independently.
# (b) Compare the models by the logl-likelihood-ration test
# (c) Compare empirical and theoretical distributions 
#     by the Kolmogorov-Smirnov tests.
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
source("parameter_estimation_reinforcement.R");


######################################################
# initialisation (apply changes ONLY here)
######################################################

#initial parameters
theta_hut.init = 0.5;
psi_hut.init  = 0.5;
ksi_hut.init = 0.5;
s_hut.init   = 100; 


#switches for work todo
analyse = TRUE;
produce.table = FALSE;
produce.figures = FALSE;


######################################################
# prepare the data
######################################################

myDataEsti = c();
impfer = vaccination_data$Wert[1:400]; #last few ones are NA
myDataEsti = impfer/100; 

#remove all values below 0.5 as we esteem them to be fixed voters
myDataEsti=myDataEsti[myDataEsti>0.5];

#renormalise it, such that 0.5 -> 0 and 1 -> 1
myDataEsti = myDataEsti * 2 - 1;

incidences = vaccination_data$Inzidenz[1:400];
inc_no_NA = get_rid_of(incidences);
#get the mean value of all incidences w/o the NA ones
dummy.i = mean(inc_no_NA);

#replace all NA values in incidences by the dummy.i
incidences = replace(incidences,dummy.i);


#set the incidence that will be used to the mean of all observed values
i.param = mean(incidences);


######################################################
# data fit
######################################################


if (analyse){
  # check optima
  s_hut.max=2500;
  
  res.tab = c();
  
  
  {
    mmean = mean(myDataEsti);
    
    
    hist(myDataEsti, freq = FALSE, nclass=30, 
         main="full reinforcement",
         xlim=c(0,1));
    
    
    ###################################################################
    # first run: estimate the full reinforcement model
    ###################################################################
    cat("first run", "\n");
    
    #para.ref in following order: nu_hut,theta_hut,psi_hut,ksi_hut,s_hut
    para.ref = c(mean(myDataEsti), theta_hut.init, psi_hut.init,ksi_hut.init,s_hut.init);   # define init para
    lll.last = lll(para.ref);
    cat("lll.last: ", lll.last, "\n");
    
    unrestricted.model     = TRUE;      # we aim at the full model
    unrestrict.theta_hut   = TRUE;
    opti.cyclic(para.ref);
    
    cat(lll.last, "\n");
    para.unrest = para.last;        # store the result
    lll.unrest  = lll.last;
    theta.res   = theta_hut*s_hut;
    
    # Kolmogorov-Smirnov test
    res.ks = ks.test(myDataEsti, function(x){pReinforce(x)}); 
    
    
    ###################################################################
    # second run: reinforcement model, force equal reinforcement parameters 
    ###################################################################
    
    cat("second run","\n");
    para.ref = c(mean(myDataEsti), theta_hut.init, psi_hut.init,ksi_hut.init,s_hut.init);  # define init para
    lll.last = lll(para.ref);
    cat(lll.last, "\n");
    hist(myDataEsti, freq = FALSE, nclass=30, 
         main="reinforcement with equal reinf. params",
         xlim=c(0,1));
    
    unrestricted.model     = TRUE;      # we aim at the full model
    unrestrict.theta_hut   = FALSE;     # we want to keep equal parameters for reinforcement
    opti.cyclic(para.ref);
    cat(lll.last, "\n");
    para.halfRestrict = para.last;        # store the result
    lll.halfRestrict  = lll.last;
    theta.res   = theta_hut*s_hut;
    
    # Kolmogorov-Smirnov test
    res.halfRestrict.ks = ks.test(myDataEsti, function(x){pReinforce(x)}); 
    
    
    ###################################################################
    # third run: estimate the zealot model (beta-distrib)
    ###################################################################
    
    cat("third run", "\n");
    unrestricted.model     = FALSE;           # we fix all reinforcement-paras
    unrestrict.theta_hut   = FALSE;
    
    theta_hut.zealot = 0;
    
    para.ref = c(mean(myDataEsti), theta_hut.zealot, psi_hut.init,ksi_hut.init,s_hut.init);
    
    hist(myDataEsti, freq = FALSE, nclass=30, 
         main="zealot model",
         xlim=c(0,1));
    
    opti.cyclic(para.ref);
    cat("lll.last: ", lll.last, "\n");
    para.restrict = para.last;            # store result
    lll.restrict  = lll.last;
    
    # kolmogorov-smirnov-test
    res.restric.ks = ks.test(myDataEsti, function(x){pReinforce(x)});
    
    line = c("run", "pro-vaxx", 
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
  load(file="datAnaCombinedModel_V1.rSave");
  sink(file="datCombinedModel.tex");
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
  mmean <<- mean(myDataEsti);
  
  post(paste("Combined_model",as.character(j),".eps",sep=""));
  hist(myDataEsti, freq = FALSE, 
       main=paste("Vaccinational behaviour"),
       xlim=c(0,1), xlab="amount of pro-vaxxers x", nclass=30);
  
  para.last = as.double(res.tab[j, 4:8]);
  
  lll.last = lll(para.last);
  cat(lll.last, "\n");  mmean <<- mean(myDataEsti);
  
  curve(g(x), add=TRUE,  lwd=2);
  
  para.last = as.double(res.tab[j, 16:20]);
  lll.last = lll(para.last);
  cat(lll.last, "\n");
  curve(g(x), add=TRUE,  lwd=2, lty=2);
  dev.off();
  
  
  post(paste("Combined_model_2",as.character(j),".eps",sep=""));
  hist(myDataEsti, freq = FALSE, 
       main=paste("Vaccinational behaviour"),
       xlim=c(0,1), xlab="amount of pro-vaxxers x", nclass=30);
  
  para.last = as.double(res.tab[j, 4:8]);
  
  lll.last = lll(para.last);
  cat(lll.last, "\n");  mmean <<- mean(myDataEsti);
  
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