setwd("C:/Users/ge69fup/Documents/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit/Brexit")
#setwd("D:/Dokumente/Uni/TUM/Mathe_B_Sc/SS_20/Bachelorarbeit/bachelorarbeit-repo/R_bachelorarbeit");


"
#my system (no params neded)
rhs <- function(state){
  s = state[1]; i = state[2]; 
  x = state[3];
  s1 = -bbeta*s*i+(1-x)*B-B*s;
  i1 = bbeta*s*i - alpha*i-B*i;
  x1 = -c*x*theta2*(1-x+n2+a*(1-i))/((x+n1+a*i)+theta2*(1-x+n2+a*(1-i)));
  x1 = x1+c*(1-x)*theta1*(x+n1+a*i)/(theta1*(x+n1+a*i)+(1-x+n2+a*(1-i)));
  return(c(s1,i1,x1));
}"
#my system (have to pass params)
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

get.jacobian <- function(state,params){
  # compute jacobian of the vector fieeld at point "state"
  #     ( (f_1)_x, (f_2)_y, (f_3)_z )
  # J = ( (f_2)_x, (f_2)_y, (f_3)_z )
  #     ( (f_3)_x, (f_2)_y, (f_3)_z )
  J = matrix(NA, 3, 3);
  hh = 1e-3;
  
  rhsp  = rhs.param(params,state+c(hh,0,0));
  rhsm  = rhs.param(params,state+c(-hh,0,0));
  grad  = (rhsp-rhsm)/(2*hh);
  J[1,] = grad;
  rhsp  = rhs.param(params,state+c(0,hh,0));
  rhsm  = rhs.param(params,state+c(0,-hh,0));
  grad  = (rhsp-rhsm)/(2*hh);
  J[2,] = grad;
  rhsp  = rhs.param(params,state+c(0,0,hh));
  rhsm  = rhs.param(params,state+c(0,0,-hh));
  grad  = (rhsp-rhsm)/(2*hh);
  J[3,] = grad;
  
  return(J);
}

check_for_imaginary_ev <- function(evals, tol){
  t = 1;
  while(t <= length(evals)){
    if(abs(Re(evals[t])) <= tol && Im(evals[t]) != 0){
      return (TRUE);
    }
    t = t+1;
  }
  
  return(FALSE);
}

get.stat.state <- function(init.state,my.tol,max.iter, params){
  # forward simulation until max iter, 
  # or ||rhs|| <- my.tol
  
  my.tol2 = my.tol**2;
  my.state = init.state;
  h        = 0.01;
  for (i in 1:max.iter){
    loc.rhs = rhs.param(params, my.state);
    if (sum(loc.rhs**2)<my.tol2){
      return(my.state);
    }
    my.state = my.state+h*loc.rhs;
  }
  cat("# get.stat.state DID NOT CONVERGE!!!\n");
  cat("# ",loc.rhs, "\n");
  cat("# get.stat.state DID NOT CONVERGE!!!\n");
  
  return(my.state);
}


#init
s  = 0.8;
i = 0.2;
x = 0.2;

# parameter
bbeta = 5    # contact rate 5
alpha = 0.45  # recovery rate 0.45
B     = 2    # birthrate - deterimes thetime scale -- fixed. 2

theta1 = 1; # reinforcement parameter of vacc pro (x)
theta2 = theta1; # reinforcement parameter of vacc con (1-x)
n1     = 10;    # zelots pro vacc
n2     = 10;    # zealots contra vacc
a      = 0.1;    # influence I: 0.1
c      = 0.3# 1;   # time scale reinforce

# we calibrate the system s.t. x=1/2 is a stat states
s.star = (B+alpha)/bbeta;
i.star = (1-0.5)*B/(alpha+B)-B/bbeta;
cat("i.star = ", i.star, "\n");

n1 = max(c(a*(1-2*i.star), 0));
n2 = a*(2*i.star-1)+n1;
n.star = a*i.star+n1
theta1 = (1-2*n.star)/(1+2*n.star);
theta2 = theta1;

#alpha, beta, mu, a, b, nx, ny, thetax, thetay
params_names = c('alpha', 'bbeta', 'c', 'a', 'B', 'n1', 'n2', 'theta1', 'theta2');
params = c(alpha,bbeta,c,a,B,n1,n2,theta1,theta2);

#my_fav_param = c(index, value); (index in params)

my_fav_param = c(5,0);

params[my_fav_param[1]] = my_fav_param[2];

state = c(s,i,x);
horizon = 100;

h = 0.01; tt = 0;
#start
J = get.jacobian(state,params)
evals = eigen(J)$values;
title = paste('Spektrum in Abhängigkeit von',params_names[my_fav_param[1]]);
#plot(main=title,Re(evals), Im(evals), xlab = 'Re', ylab = 'Im', cex = 0.01);
#res <<- numeric();
res <<- data.frame('tt', params_names[my_fav_param[1]], 's*','i*','x*','eig_val_1','eig_val_2','eig_val_3');
plot(main=title,0,0,xlab='Re',ylab='Im',cex=0.01);
tol = 1e-6;
max_iter = 500000;

while(tt < horizon){
  cat('tt', tt);
  steady = get.stat.state(state, tol, max_iter,params);
  #todo: wie bekomme ich diese Translation durch Parameter hin?
  #J = get.jacobian(steady,params) + 0.25*diag(3);
  
  J = get.jacobian(steady,params);
  evals = eigen(J)$values;
  points(Re(evals), Im(evals), cex = 0.01);
  
  if(check_for_imaginary_ev(evals,tol)){
    cat('Purely imaginary eigenvalue at steady state = ', steady, 'tt=', tt);
    cat('Re:', Re(evals));
    cat('Im:', Im(evals));
  }
  tt = tt+h;
  params[my_fav_param[1]] = params[my_fav_param[1]] + h
  res   <<- rbind(res, c(tt,params[my_fav_param[1]], steady, evals));
}