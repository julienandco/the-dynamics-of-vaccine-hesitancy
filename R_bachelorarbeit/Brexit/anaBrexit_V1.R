#
#
# Estimate the parameters for the Brexit.
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
#setwd("C:/Users/ge69fup/Documents/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/Brexit")
setwd("D:/Dokumente/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/Brexit");

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
load("electBrexitV4.rSave");


###############
# read tools
#
#################
source("estiReinforceV0.txt");

######################################
#
#  do it

if (1==1){
  # check optima
  skal.max=2500;
  
  res.tab = c();
  
  
  {
    e =    electGBbrexit[[1]];
    myDataEsti = c();
    ep = e$PartyVotes;
    myDataEsti = ep[,1]/(ep[,1]+ep[,2]);     
    
    mmean = mean(myDataEsti);
    
    
    hist(myDataEsti, freq = FALSE, nclass=30, 
         main=as.character(e$year),
         xlim=c(0,1));
    ###################################################################
    # first run: estimate the full reinforcement model
    ###################################################################
    
    para.ref = c(mean(myDataEsti), 0.5, 0.5,100);   # define init para
    lll.last = lll(para.ref);
  }}
    cat(lll.last, "\n");
    'curve(g(x), add=TRUE, col="blue", lwd=2);
    
    unrestricted.model = TRUE;      # we aim at the full model
    unrestrict.theta   = TRUE;
    opti.cylcic(para.ref);
    cat(lll.last, "\n");
    para.unrest = para.last;        # store the result
    lll.unrest  = lll.last;
    theta.res   = theta*skal;
    
    # Kolmogorov-Smirnov test
    res.ks = ks.test(myDataEsti, function(x){pReinforce(x)}); 
    
    # produce a figure with histogram and estimated distribution
    party.x = "remain";
    
    ###################################################################
    # second run: reinforcement model, force equal reinforcement parameters 
    ###################################################################
    para.ref = c(mean(myDataEsti), 0.5, 0.5,100);   # define init para
    lll.last = lll(para.ref);
    cat(lll.last, "\n");
    hist(myDataEsti, freq = FALSE, nclass=30, 
         main=as.character(e$year),
         xlim=c(0,1));
    curve(g(x), add=TRUE, col="green", lwd=2);
    
    unrestricted.model = TRUE;      # we aim at the full model
    unrestrict.theta   = FALSE;     # we want to keep equal parameters for reinforcement
    opti.cylcic(para.ref);
    cat(lll.last, "\n");
    para.halfRestrict = para.last;        # store the result
    lll.halfRestrict  = lll.last;
    theta.res   = theta*skal;
    
    # Kolmogorov-Smirnov test
    res.halfRestrict.ks = ks.test(myDataEsti, function(x){pReinforce(x)}); 
    
    
    ###################################################################
    # third run: estimate the zealot model (beta-distrib)
    ###################################################################
    
    unrestricted.model = FALSE;           # we fix all reinfrcement-paras
    unrestrict.theta   = FALSE;
    para.ref = c(mean(myDataEsti), 0.0, 0.5,100);
    hist(myDataEsti, freq = FALSE, nclass=30, 
         main=as.character(e$year),
         xlim=c(0,1));
    curve(g(x), add=TRUE, col="red", lwd=2);
    
    opti.cylcic(para.ref);
    cat(lll.last, "\n");
    para.restrict = para.last;            # store result
    lll.restrict  = lll.last;
    
    # komogorov-smirnov-test
    res.restric.ks = ks.test(myDataEsti, function(x){pReinforce(x)});
    
    line = c(e$year, party.x, 
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
  # names orient themseves ar the supplement II of the paper
  col.names = c(
    "year", "party", 
    "Theta1PlusTheta2.unr",  
    "nu.unr", "theta.hat.unr", "psi.unr", "s.unr",
    "lll.unr",
    "nu.halfr", "theta.hat.halfr", "psi.halfr", "s.halfr",
    "lll.halfr",
    "nu.restr", "theta.hat.restr", "psi.restr", "s.restr",
    "lll.restr",
    "lll.unr.rest", "lll.unr.halfr",
    "ks.unres", "ks.halfRestr", "ks.restr");
  dimnames(res.tab)[[2]] =col.names;
  
  save(file="datAnaBrexit_V1.rSave", res.tab);
  
}
'



if (0==1){
  # produce a table
  load(file="datAnaBrexit_V1.rSave");
  sink(file="datBrexit.tex");
  cat(dimnames(res.tab)[[2]][c(1,2,4,5,6,7)]); cat(" theta1 ");cat(" theta2 "); cat(dimnames(res.tab)[[2]][c(19,21,23)]);
  cat("\n");
  nn = dim(res.tab)[1];
  for (i in 1:nn){
    cat(res.tab[i,1], " & ", res.tab[i,2], " & ");
    cat(res.tab[i,4], " & ", res.tab[i,5], " & ");
    cat(res.tab[i,6], " & ", res.tab[i,7], " & ");
    # theta_2 = h.s*h.theta*(1-h.psi)
    cat(as.double(res.tab[i,5])*as.double(res.tab[i,7])*as.double(res.tab[i,6]), " & ");
    cat(as.double(res.tab[i,5])*as.double(res.tab[i,7])*(1-as.double(res.tab[i,6])), " & ");
    cat(as.double(res.tab[i,19]), " & ", 
        as.double(res.tab[i,21]), " & ", 
        as.double(res.tab[i,23]), 
        "\\\\\n");
  }
  
  cat("Test hat.psi=0.5 versus free model (where is the reinforcement?)\n");
  cat(as.double(res.tab[i,20]), " \n");
  
  sink();
  
}



if (0==1){
  # produce figures
  e =    electGBbrexit[[1]];
  ep = e$PartyVotes;
  myDataEsti = ep[,1]/(ep[,1]+ep[,2]);     
  mmean <<- mean(myDataEsti);
  
  
  post(paste("GBbrexit",as.character(i),".eps",sep=""));
  hist(myDataEsti, freq = FALSE, 
       main=paste("Brexit"),
       xlim=c(0,1), xlab="vote share ``remain''", nclass=30);
  
  para.last = as.double(res.tab[i, 4:7]);
  
  lll.last = lll(para.last);
  cat(lll.last, "\n");  mmean <<- mean(myDataEsti);
  
  curve(g(x), add=TRUE,  lwd=2);
  
  para.last = as.double(res.tab[i, 14:17]);
  lll.last = lll(para.last);
  cat(lll.last, "\n");
  curve(g(x), add=TRUE,  lwd=2, lty=2);
  dev.off();
  
  
  post(paste("GBbrexit2",as.character(i),".eps",sep=""));
  hist(myDataEsti, freq = FALSE, 
       main=paste("Brexit"),
       xlim=c(0,1), xlab="vote share ``remain''", nclass=30);
  
  para.last = as.double(res.tab[i, 4:7]);
  
  lll.last = lll(para.last);
  cat(lll.last, "\n");  mmean <<- mean(myDataEsti);
  
  curve(g(x), add=TRUE,  lwd=2);
  
  para.last = as.double(res.tab[i, 9:12]);
  lll.last = lll(para.last);
  cat(lll.last, "\n");
  curve(g(x), add=TRUE,  lwd=2, col = "green");
  
  para.last = as.double(res.tab[i, 14:17]);
  lll.last = lll(para.last);
  cat(lll.last, "\n");
  curve(g(x), add=TRUE,  lwd=2, lty=2);
  dev.off();
  
  
}
