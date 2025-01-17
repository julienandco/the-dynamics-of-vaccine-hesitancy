
#
# Estimate the paras for the reinforcement model
#####################################################################
#
#    estimate the paras of the reinforcement model - general functions
#
#####################################################################

#
# myDataEsti: vector with values in (0,1); the data it uses for the estimation
#

##�berall wo verteilung aufgerufen wird, muss mein param rein
#
# in case: define work directory

#setwd("C:/Users/ge69fup/Documents/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/data_fitting")
#setwd("D:/Dokumente/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/data_fitting");

replace <- function(my_vec,replacement){
  for (i in 1:length(my_vec)){
    if (is.na(my_vec[i])){
      my_vec[i] = replacement;
    }
  }
  return (my_vec);
}

post <- function(nme){
  # remove blanks
  xx = strsplit(nme, " ");
  nme1 = paste(xx[[1]], collapse="");
  cat(nme1, "\n");
  postscript(file=nme1, pointsize=28, paper="special",
             width=8, height=9, horizontal=FALSE);
}

##################
#
# define distrib
#
##################

s_hut.max = 2000;

f.norm <- function(){
  # an additive constant to avoid large numbers
  # in the exponent; this number is chosen in dependence
  # on the mean value of the data at hand.
  x = mmean;
  i = mean(info$Inzidenz);
  
  return(-1*(s_hut*theta_hut*(0.5*x**2-phi_hut*x)
             +log(x)*((1-theta_hut)*nu_hut*(1-ksi_hut)+nu_hut*ksi_hut*i)*s_hut
             +log(1-x)*((1-theta_hut)*(1-nu_hut)*(1-ksi_hut)+nu_hut*ksi_hut*(1-i))*s_hut)
  );
}

f <- function(x){
  # N1 = (1-theta_hut)*(1-ksi_hut)*nu_hut*scal+1, N2 = (1-theta_hut)*(1-nu_hut)*(1-ksi_hut)*scal +1
  # theta_hut \in (0,1), nu_hut \in (0,1), ksi_hut \in (0,1) scal >0
  return(
    exp(  s_hut*theta_hut*(0.5*x**2-phi_hut*x)
          +log(x)*((1-theta_hut)*nu_hut*(1-ksi_hut)+nu_hut*ksi_hut*i.param)*s_hut
          +log(1-x)*((1-theta_hut)*(1-nu_hut)*(1-ksi_hut)+nu_hut*ksi_hut*(1-i.param))*s_hut +f.norm()
    )
  );
}

get.cc <- function(){
  # get the normalisation constant
  CC=(integrate(f, lower=10**-6,upper=1-10**-6)$value)**(-1);
  return(CC);
}


#g<-function(x){
  # density (note: CC is computed using the constant f.norm().
  # That constant cancels at the end of the day.
#  return(
#    ccs*
#      exp(s_hut*theta_hut*(0.5*x**2-phi_hut*x)
#          +log(x)*((1-theta_hut)*nu_hut*(1-ksi_hut)+nu_hut*ksi_hut*i.param)*s_hut
#          +log(1-x)*((1-theta_hut)*(1-nu_hut)*(1-ksi_hut)+nu_hut*ksi_hut*(1-i.param))*s_hut +f.norm()
#      )
#  );
#}

get.g <- function(j){
  g <- function(x){
    return(
          ccs[j]*
            exp(s_hut*theta_hut*(0.5*x**2-phi_hut*x)
                +log(x)*((1-theta_hut)*nu_hut*(1-ksi_hut)+nu_hut*ksi_hut*info$Inzidenz[j])*s_hut
                +log(1-x)*((1-theta_hut)*(1-nu_hut)*(1-ksi_hut)+nu_hut*ksi_hut*(1-info$Inzidenz[j]))*s_hut +f.norm()
            )
        );
  }
  
  return(g);
}

lll.dat <- function(x,nu_hut,theta_hut,phi_hut,ksi_hut,i,s_hut, ccs){
  res = numeric(length(x));
  for (j in 1:length(x)){
    res[j] =  s_hut*theta_hut*(0.5*(x[j])**2-phi_hut*x[j])
    +log(x[j])*((1-theta_hut)*nu_hut*(1-ksi_hut)+nu_hut*ksi_hut*i[j])*s_hut
    +log(1-x[j])*((1-theta_hut)*(1-nu_hut)*(1-ksi_hut)+nu_hut*ksi_hut*(1-i[j]))*s_hut
    +log(ccs[j])+f.norm()
  }
  
  return(res);
}

pReinforce.loc <- function(x) {
  # parameters given by global parameters;
  # particularly, CC is defined.
  res = integrate(g,lower = 0, upper = x);
  return(res$value);
}

pReinforce <- function(x){
  return(sapply(x, pReinforce.loc));
}

#######################
# estimate paras:
# nu_hut, N2, phi_hut, theta_hut2, A
########################
theta_hut.max = 1800;
theta_hut.max = 1900;
tryCatch.W.E <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}

##einstiegsfunktion
# compute log likeli
lll <- function(para){
  # nu_hut = (1-theta_hut)*nu_hut+1, N2 = (1-theta_hut)*(1-nu_hut)*scal +1
  # theta_hut \in (0,1), nu_hut \in (0,1), scal >0
  ppara     <<- para;
  nu_hut    <<- para[1];
  theta_hut <<- para[2];
  phi_hut   <<- para[3];
  ksi_hut   <<- para[4];
  s_hut     <<- min(s_hut.max,abs(para[5]));
  OK = TRUE;
  
  #compute all C (for every datapoint)
  ccs = numeric(length(info$Wert));
  ccs <<- ccs;
  
  for(j in 1:length(ccs)){
    #get current incidence
    i.param <<- info[j,2]
    aa = tryCatch.W.E(get.cc());
    aa <<- aa;
    if (!(is.double(aa$value))>0) return(-10000);
    ccs[j]<<- aa$value;
  }
  
  temp <<- lll.dat(info$Wert, nu_hut,theta_hut,phi_hut,ksi_hut,info$Inzidenz,s_hut,ccs);
  return(sum(temp));    
}


#
# for the optimization: vary only one of the parameters
#
search.p1 <- function(px){
  # para.last gives the framework; we modify parameter 1 only
  pxx = para.last; pxx[1] = px;
  return( lll(pxx) );
}
search.p2 <- function(px){
  # para.last gives the framework; we modify parameter 2 only
  pxx = para.last; pxx[2] = px;
  return( lll(pxx) );
}
search.p3 <- function(px){
  # para.last gives the framework; we modify parameter 3 only
  pxx = para.last; pxx[3] = px;
  return( lll(pxx) );
}
search.p4 <- function(px){
  # para.last gives the framework; we modify parameter 4 only
  pxx = para.last; pxx[4] = px;
  return( lll(pxx) );
}



############################################
#
# optimization
#
############################################

opti.cyclic <- function(para.init){
  # optimize cyclically the parameters.
  #
  # we have different modes 
  # unrestricted.model == FALSE: fix all paras expect of s_hut, and nu_hut => we can take
  #                              the reinforcement to zero and fit a beta distribution.
  # unrestrict.theta_hut   == TRUE:  allow theta_hut to vary.
  #                              (if FALSE: we can fix theta_hut=0.5, and in this
  #                               try to find out where the reinforcement takes place)
  
  mmean <<- mean(info$Wert);
  
  para.last <- para.init; 
  lll.last  <- lll(para.ref);
  last.s_hut = -1; no.s_hut.const = 0;last.lll=-1e10;
  fertig = FALSE;
  
  i = 0;
  while ((i<10000)&(fertig==FALSE)){
    i = i+1;
    cat("para.last:", para.last,"\n");
    lll.x = lll.last;
    para.last <<- para.last;
    res1 = optimize(search.p1, interval=c(0,1), maximum=TRUE);
    para.loc1 = para.last; para.loc1[1]=res1$maximum;
    lll.lok  = lll(para.loc1);
    if (lll.lok>lll.last){
      para.last <- para.loc1;
      lll.last  <- lll.lok;
    }
    
    ##optimize ksi
    para.last <<- para.last;
    res4 = optimize(search.p4, interval=c(0,1), maximum=TRUE);
    para.loc4 = para.last; para.loc4[4]=res4$maximum;
    lll.lok  = lll(para.loc4);
    if (lll.lok>lll.last){
      para.last <- para.loc4;
      lll.last  <- lll.lok;
    }
    
    if (unrestricted.model){
      para.last <<- para.last;
      res2 = optimize(search.p2, interval=c(0,1), maximum=TRUE);
      para.loc2 = para.last; para.loc2[2]=res2$maximum;
      lll.lok  = lll(para.loc2);
      if (lll.lok>lll.last){
        para.last <- para.loc2;
        lll.last  <- lll.lok;
      }
      
      if (unrestrict.theta_hut) {
        para.last <<- para.last;
        res3 = optimize(search.p3, interval=c(0,1), maximum=TRUE);
        para.loc3 = para.last; para.loc3[3]=res3$maximum;
        lll.lok  = lll(para.loc3);
        if (lll.lok>lll.last){
          para.last <- para.loc3;
          lll.last  <- lll.lok;
        }
      }
    }
    
    lll.1     = lll(para.last);
    Delta = 0.01;
    if (para.last[5]>10)  Delta=0.1;
    if (para.last[5]>50)  Delta=0.5;
    if (para.last[5]>100) Delta=1;
    
    lll.p2    = lll(para.last+c(0,0, 0,0, Delta));
    lll.m2    = lll(para.last+c(0,0,0,0, -Delta));
    if (lll.p2>lll.1){
      paral.loc3 = para.last+c(0,0,0,0, Delta);
      para.last  = para.last+c(0,0,0,0, Delta);
      last.lll = lll.p2;
    } else {
      if (lll.m2>lll.1){
        paral.loc3 = para.last+c(0,0,0,0,-Delta);
        para.last  <- para.last+c(0,0,0,0,-Delta);
        last.lll   <- lll.m2
      } else {
        paral.loc3 = para.last+c(0,0,0,0,0);
        para.last  <- para.last+c(0,0,0,0,0);
        last.lll   <- lll.1
      }
    }
    if (last.s_hut!=para.last[3]){
      last.s_hut = para.last[3];
      no.s_hut.const = 0;
      lll.x =last.lll;
    } else {
      no.s_hut.const = no.s_hut.const+1;
      loc.lll = lll(para.last);
      if (last.lll>lll.x+1e-6) {
        no.s_hut.const = 0;
      }
      if(no.s_hut.const>100) fertig=TRUE;
    }
    
    para.last <<- para.last;
    lll.last  <<- last.lll;
    cat(i," ", lll(para.last), "\n");     
    #curve(g(x), add=TRUE, col="blue");
    
    for(j in 1:length(info$Wert)){
      g = get.g(j);
      
      #did we finish optimising?
      if(fertig){
        curve(g(x), add=TRUE, col="red", lwd=2);
        #save the optimal curve
        g.res[current.g] = g(x);
      } else{
        curve(g(x), add=TRUE, col="green", lwd=2);
      }
    }
    
    
    if(fertig){
      cat("plotted final optimisation results of round ", i, " in red colour.", "\n");
    }else{
      cat("plotted opti results of round ", i, "\n");
    }
    
  }
}