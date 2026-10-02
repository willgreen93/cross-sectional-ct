functions {
  real vlfunc(real a1, real b1, real l1, real t1, real tmax1, real vlmax1, int k){
    real vl;
    if(k==1) vl = (t1<tmax1)? (vlmax1+exp(a1)*(t1-tmax1)):(vlmax1-exp(b1)*(t1-tmax1));
    if(k==2) vl = (t1<tmax1)? (vlmax1+exp(a1)*(t1-tmax1)):(vlmax1-exp(b1)/exp(a1)*(t1-tmax1));
    if(k==3) vl = (t1<tmax1)? (vlmax1*exp(a1)+exp(a1)*(t1-tmax1)):(vlmax1*exp(a1)-exp(b1)/exp(a1)*(t1-tmax1));
    if(k==4) vl = (t1<tmax1)? (vlmax1+vlmax1/exp(a1)*(t1-tmax1)):(vlmax1-exp(b1)/vlmax1*(t1-tmax1));
    
    if(k==6) vl = (t1<tmax1)? (vlmax1+(0.5+exp(a1))*(t1-tmax1)):(vlmax1-(0.25+exp(b1))*(t1-tmax1));
    if(k==7) vl = (t1 < tmax1 - 0.5*exp(l1))? (vlmax1+(0.5+exp(a1))*(t1-tmax1+0.5*exp(l1)))
                : (t1 < tmax1 + 0.5*exp(l1))?  vlmax1
                :                              (vlmax1-(0.25+exp(b1))*(t1-tmax1-0.5*exp(l1)));
    
    //if(k==4) vl = (t1<tmax1)? (vlmax1+exp(a1)*(t1-tmax1)): (t1<tmax1+exp(l1))? vlmax1: vlmax1-exp(b1)*(t1-tmax1-exp(l1));
    //if(k==5) vl = (t1<tmax1)? (vlmax1+exp(a1)*(t1-tmax1)): (t1<tmax1+exp(l1))? vlmax1: vlmax1-exp(b1)/exp(a1)*(t1-tmax1-exp(l1));
    //if(k==6) vl = (t1<tmax1)? (vlmax1*exp(a1)+exp(a1)*(t1-tmax1)): (t1<tmax1+exp(l1))? vlmax1*exp(a1): vlmax1*exp(a1)-exp(b1)/exp(a1)*(t1-tmax1-exp(l1));
     
		return(vl);
  }
}

data { 
  int N;             // number of data points
  int<lower = 0> K;  // number of individuals
  real vl_min[N];       // baseline viral load
  //real tmax_bar;
  
  int k;
  //vector[N] study;

  real t[N];         // time points of observations 
  int individual[N]; // identity of individuals
  vector[N] vl;      // viral load reads
}

parameters{
  //non-hierarchical parmaeters
  
  //hierarchical parameters
  
  real<lower=0> sigma;
  real a_bar;
  real<lower=0> a_sigma;
  real b_bar;
  real<lower=0> b_sigma;
  real l1;
  real vl_max_bar;
  real<lower=0> vl_max_sigma;
  
  real tmax[K];
  
  real a_raw[K];
  real b_raw[K];
  real vl_max_raw[K];
  //real tmax_raw[K];
  
  real<lower=0, upper=1> l_fn; // false negative probability
}

transformed parameters{
  real a[K];
  real b[K];
  real vl_max[K];
  vector[N] log_lik;
  //real tmax[K];
  
  real ll_total;
  
  for(i in 1:K){
    a[i] = a_bar + a_raw[i]*a_sigma;
    b[i] = b_bar + b_raw[i]*b_sigma;
    vl_max[i] = vl_max_bar + vl_max_raw[i]*vl_max_sigma;
    //tmax[i] = tmax_bar + tmax_raw[i]*tmax_sigma;
  }
  
  //non-hierarchical likelihood
  for(n in 1:N){
    if(vl[n] > vl_min[n])   log_lik[n] = log(1-l_fn) + normal_lpdf(vl[n] | vlfunc(a[individual[n]], b[individual[n]], l1, t[n], tmax[individual[n]], vl_max[individual[n]], k), sigma);
    else                    log_lik[n] = log(l_fn + (1-l_fn) * exp(normal_lcdf(vl[n] | vlfunc(a[individual[n]], b[individual[n]], l1, t[n], tmax[individual[n]], vl_max[individual[n]], k), sigma)));
  }
    
  ll_total = sum(log_lik);
  
}


model{
  {
  target += ll_total;
  }
  //non-hierarchical priors
  
  // hyper-priors
  a_bar ~ normal(-1,1);
  a_sigma ~ normal(1,1);
  b_bar ~ normal(-1,1);
  b_sigma ~ normal(1,1);
  l1 ~ normal(-0.5,1);
  sigma ~ normal(0,2);

  tmax ~ normal(0,1);
  vl_max_bar ~ normal(8,2);
  vl_max_sigma ~ normal(2,2);
  
  a_raw ~ std_normal();
  b_raw ~ std_normal();
  //tmax_raw ~ std_normal();
  vl_max_raw ~ std_normal();
}
