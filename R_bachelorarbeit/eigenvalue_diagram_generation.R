setwd("C:/Users/ge69fup/Documents/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/Brexit")
#setwd("D:/Dokumente/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit");


# parameter
bbeta = 5    # contact rate
alpha = 0.45  # recovery rate
B     = 2    # birthrate - deterimes thetime scale -- fixed.

theta1 = 1; # reinforcement parameter of vacc pro (x)
theta2 = theta1; # reinforcement parameter of vacc con (1-x)
n1     = 10;    # zelots pro vacc
n2     = 10;    # zealots contra vacc
a      = 0.1;    # influence I
Omega  = 0.0;    # influence vaccination oponent
c      = 200# 1;   # time scale reinforce

# adapt n1, n2 s.t. a*i+n1 = a*(1-i)+n2
# choose n1 minial, s.t. we have a non-negative n2 (=0)
adapt_params <- function(a,i.star, givenN1, givenN2, givenTheta1){
  if(missing(givenN1)){
    givenN1 = max(c(a*(1-2*i.star), 0));
  }
  if(missing(givenN2)){
    givenN2 = a*(2*i.star-1)+givenN1;
  }
  
  n.star = a*i.star+givenN1
  
  if(missing(givenTheta1)){
    givenTheta1 = (1-2*n.star)/(1+2*n.star);
  }
  
  givenTheta2 = givenTheta1

  return(c(givenN1,givenN2,givenTheta1,givenTheta2))
}

#my system
rhs.param <- function(params, state){
  s = state[1]; i = state[2]; 
  x = state[3];
  alpha = params[1];
  bbeta= params[2];
  c= params[3];
  a= params[4];
  B= params[5];
  n1= params[6];
  n2= params[7];
  theta1= params[8];
  theta2= params[9];
  
  s1 = -bbeta*s*i+(1-x)*B-B*s;
  i1 = bbeta*s*i - alpha*i-B*i;
  x1 = -c*x*theta2*(1-x+n2+a*(1-i))/((x+n1+a*i)+theta2*(1-x+n2+a*(1-i)));
  x1 = x1+c*(1-x)*theta1*(x+n1+a*i)/(theta1*(x+n1+a*i)+(1-x+n2+a*(1-i)));
  return(c(s1,i1,x1));
}

get.jacobian <- function(state){
  # compute jacobian of the vector fieeld at point "state"
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

check_for_imaginary_ev <- function(evals){
  t = 1;
  while(t <= length(evals)){
    if(Re(evals[t]) == 0 && Im(evals[t]) != 0){
      return (TRUE);
    }
    t = t+1;
  }
  
  return(FALSE);
}

#init
s  = 0.8;
i = 0.2;
x = 0.2;


#alpha, beta, mu, a, b, nx, ny, thetax, thetay
params_names = c('alpha', 'bbeta', 'c', 'a', 'B', 'n1', 'n2', 'theta1', 'theta2');
params = c(alpha,bbeta,c,a,B,n1,n2,theta1,theta2);

#my_fav_param = c(index, value); (index in params)
my_fav_param = c(9,0);

state = c(s,i,x);
horizon = 100;

h = 0.01; tt = 0;
#start
cat('start');
J = get.jacobian(state)
evals = eigen(J)$values;
title = paste('Spektrum in Abhängigkeit von ',params_names[my_fav_param[1]]);
plot(main=title,Re(evals), Im(evals), xlab = 'Re', ylab = 'Im');
res <<- numeric();

while(tt < horizon){
  params[my_fav_param[1]] = params[my_fav_param[1]] + h
  state = state + h*rhs.param(params,state);
  J = get.jacobian(state)
  evals = eigen(J)$values;
  points(Re(evals), Im(evals), cex = 0.01);
  
  if(check_for_imaginary_ev(evals)){
    cat('Purely imaginary eigenvalue at tt = ', tt);
    cat('Re:', Re(evals));
    cat('Im:', Im(evals));
  }
  tt = tt+h;
  res   <<- rbind(res, c(tt, state));
}
