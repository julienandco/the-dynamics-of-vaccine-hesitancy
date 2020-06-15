#
# SIRS model with reinforcement
#
# Aim:
# find a Hopf bifurcation
#
#setwd("C:/Users/ge69fup/Documents/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/Brexit")
setwd("D:/Dokumente/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit");




# parameter
bbeta = 5    # contact rate
alpha = 0.45  # recovery rate
B     = 2    # birthrate - deterimes thetime scale -- fixed.

# we calibrate the system s.t. x=1/2 is a stat states
s.star = (B+alpha)/bbeta;
i.star = (1-0.5)*B/(alpha+B)-B/bbeta;
cat("i.star = ", i.star, "\n");


theta1 = 1; # reinforcement parameter of vacc pro (x)
theta2 = theta1; # reinforcement parameter of vacc con (1-x)
n1     = 10;    # zelots pro vacc
n2     = 10;    # zealots contra vacc
a      = 0.1;    # influence I
Omega  = 0.0;    # influence vaccination oponent
c      = 200# 1;   # time scale reinforce

# adapt n1, n2 s.t. a*i+n1 = a*(1-i)+n2
# choose n1 minial, s.t. we have a non-negative n2 (=0)
n1 = max(c(a*(1-2*i.star), 0));
n2 = a*(2*i.star-1)+n1;
n.star = a*i.star+n1
theta1 = (1-2*n.star)/(1+2*n.star);
theta2 = theta1;
cat("a = ", a, " n1 = ", n1, "n2 = ", n2, " n.star = ", n.star, 
    " theta.pich = ", (1-2*n.star)/(1+2*n.star), "\n"); 

#
#
theo.stat.point = c(s.star, i.star, 0.5);



#init
s = 0.8;
i = 0.2;
r = 0;
x = 0.20;
#x = 0.50;

state = c(s,i,x);   # no r-component
# rhs of ODE
rhs <- function(state){
  s = state[1]; i = state[2]; 
  x = state[3];
  s1 = -bbeta*s*i+(1-x)*B-B*s;
  i1 = bbeta*s*i - alpha*i-B*i;
  x1 = -c*x*theta2*(1-x+n2+a*(1-i))/((x+n1+a*i)+theta2*(1-x+n2+a*(1-i)));
  x1 = x1+c*(1-x)*theta1*(x+n1+a*i)/(theta1*(x+n1+a*i)+(1-x+n2+a*(1-i)));
  return(c(s1,i1,x1));
}

alu <- function(x, i){
  return(
    -x*theta2*(1-x+n2+a*(1-i))/((x+n1+a*i)+theta2*(1-x+n2+a*(1-i)))
    +(1-x)*theta1*(x+n1+a*i)/(theta1*(x+n1+a*i)+(1-x+n2+a*(1-i)))
  );
}

get.jacobian <- function(state){
  # compute jacobian of the vector field at point "state"
  #     ( (f_1)_x, (f_2)_y, (f_3)_z )
  # J = ( (f_2)_x, (f_2)_y, (f_3)_z )
  #     ( (f_3)_x, (f_2)_y, (f_3)_z )
  J = matrix(NA, 3, 3);
  hh = 1e-3;
  
  rhsp  = rhs(state+c(hh,0,0));
  rhsm  = rhs(state+c(-hh,0,0));
  grad  = (rhsp-rhsm)/(2*hh);
  J[1,] = grad;
  rhsp  = rhs(state+c(0,hh,0));
  rhsm  = rhs(state+c(0,-hh,0));
  grad  = (rhsp-rhsm)/(2*hh);
  J[2,] = grad;
  rhsp  = rhs(state+c(0,0,hh));
  rhsm  = rhs(state+c(0,0,-hh));
  grad  = (rhsp-rhsm)/(2*hh);
  J[3,] = grad;
  
  return(J);
}

simul.plot<-function(horizont){
  # simulate
  res <<- numeric();
  
  h = 0.001; tt = 0; h = 0.01;
  while (tt<horizont){
    state <<- state+h*rhs(state);
    tt    = tt + h;
    res   <<- rbind(res, c(tt, state));
    #cat(tt, " ");
  }
  
  plot(res[,1], res[,2], t="l", xlab = "tt", ylab="I", ylim=c(0,1))
  lines(res[,1], res[,3], t="l", col="blue")
  lines(res[,1], res[,4], t="l", col="orange")
  legend(0,1,legend=c("s","i","x"),col=c("black","blue","orange"), lty=1:2, cex=0.6)
}


################################################################################

curve(alu(x,i.star), from=0, to=1, ylab="antivaxx(x,i)");
curve(alu(x,i.star+0.01), from=0, to=1, add=TRUE, col="blue");
curve(alu(x,i.star+0.1), from=0, to=1, add=TRUE, col="blue");
curve(alu(x,i.star-0.01), from=0, to=1, add=TRUE, col="red");
curve(alu(x,i.star-0.1), from=0, to=1, add=TRUE, col="red");
abline(h=0);

input <- scan();
simul.plot(input)

