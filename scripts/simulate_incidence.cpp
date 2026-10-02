#include <Rcpp.h>
//[[Rcpp::depends(dqrng, BH, sitmo)]]
#include <xoshiro.h>
#include <dqrng_distribution.h>
//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::plugins(openmp)]]
#include <omp.h>

using namespace Rcpp;

//[[Rcpp::depends(dqrng)]]

/* Capacities of the per-thread frequency array ft. Previously this was a bare
 * double ft[32][41][51] while t_max was CLAMPED to 55, so t_n could reach 56 and the write
 * loop ran five elements past the end of a fixed-size stack array -- silently. The same
 * defect was found and fixed in simulate_onset.cpp; these two kernels share ancestry.
 *
 *   MAX_CORES  threads; ncores is clamped to this.
 *   MAX_CT     viral load categories; ct_n = 2*(vl_max-vli)+1, at most 31 given the
 *              clamps vl_max <= 15 and vli >= 0.
 *   MAX_T      timesteps. For the incidence form the trajectory must span
 *              max_test_day + symp_delay_lim, so this needs more headroom than the
 *              threshold form, not less.
 */
#define MAX_CORES 32
#define MAX_CT    41
#define MAX_T     64

//[[Rcpp::export]]
NumericMatrix ct_func2_cpp_inc(const NumericVector& a, const NumericVector& b, const NumericVector& c, const NumericVector& l, 
                            const NumericVector& p, const int nthreads, const NumericVector& pseeds, const int vg, const int vk, const int vt, double stoch) 
{
  dqrng::normal_distribution distA(a[0],a[1]); 
  dqrng::normal_distribution distB(b[0],b[1]);
  dqrng::normal_distribution distC(c[0],c[1]);
  dqrng::normal_distribution distL2(l[0],l[1]);
  dqrng::uniform_distribution distT(-stoch, stoch);
  
  //double lower;
  //double upper;
  //if(stoch<0.1){
  //  lower=0;
  //  upper=stoch;
  //} 
  //else{
  //  lower=-stoch;
  //  upper=stoch;
  //} 
  
  /* Set-up and some bound checking */
  double flat = l[0];

  /* Growth starts at the fitted immune-activation threshold vl_thresh, passed as p[4], rather
   * than at the limit of detection. That makes the trajectory clock here measure time since
   * crossing vl_thresh -- the same origin simulate_onset.cpp uses -- so the convolution in
   * p_array_func indexes incidence on threshold crossing and the delay to symptom onset keeps
   * exactly the meaning it has under the threshold form. Anchoring at the limit of detection
   * instead would make that delay a property of the assay rather than of the host.
   *
   * This kernel is simulate_onset.cpp with tsv set to zero: there the onset delay is drawn
   * per individual inside the C++, here it is convolved in R, which is what makes the
   * composition depend on calendar time. */
  double vthresh = p[4];

  double vli = p[1];
  if(vli<0) {
    Rprintf("Max CT is %lg, now reduced to 40\n",vli);
    vli=0;
  }
  double vl_max = p[2];
  if(vl_max>15) {
    Rprintf("Max VL is %lg, now decreased to 1\n",vl_max);
    vl_max=15;
  }
  if(vl_max<=vli) {
    Rprintf("Max VL %lg is <= min VL, now adjusted\n",vl_max);
    vl_max=15;
    vli=0;
  }
  int ct_n=(int) 2*(vl_max-vli)+1;
  if(ct_n>MAX_CT)
    Rcpp::stop("ct_n is %d but the frequency array holds only %d viral load categories", ct_n, MAX_CT);

  double t_max = p[0];
  int t_n=(int) (t_max+1);
  /* Fail loudly rather than clamping. A silent clamp both wrote past the end of ft and
     returned a matrix whose row count no longer matched what R expects when it assigns
     dimnames, which is far harder to diagnose than a stop. */
  if(t_n>MAX_T)
    Rcpp::stop("t_n is %d but the frequency array holds only %d timesteps; raise MAX_T", t_n, MAX_T);

  double a_lower = 0.5;
  double b_lower = 0.25;

  int ncores=nthreads; /* need a copy to be able to change it */
  if(ncores>MAX_CORES) {
    Rprintf("ncores is %i, now reduced to %i\n",ncores,MAX_CORES);
    ncores=MAX_CORES;
  }
  int test_pop = (int) p[3];
  if(test_pop<100) {
    Rprintf("Popiulation size is %i, now increased to 100\n",test_pop);
    test_pop=100;		
  }
  int n_s = test_pop/ncores; /* pop to simulate per thread */
  
  /* initialise thread specific random number generators with the seeds provided
   - neeed to *8 to avoid false-sharing (cache line conflicts */
  dqrng::xoroshiro128plus* rng = new dqrng::xoroshiro128plus[ncores*8];
  for (int i = 0; i < ncores*8; i++) rng[i] = dqrng::xoroshiro128plus(pseeds(i/8));
  
  omp_set_num_threads(ncores); /* set number of threads */
  
  /* static thread specific frequency matrix - need to zero */
  double ft[MAX_CORES][MAX_CT][MAX_T];
  for(int tn=0;tn<ncores;tn++)
    for(int i=0;i<ct_n;i++)
      for(int j=0;j<t_n;j++)
        ft[tn][i][j]=0;
  
  /* The one parallel loop which does all the work */
#pragma omp parallel for schedule(static,1)
  for(int tn=0;tn<ncores;tn++) {
    int l=tn*n_s;
    int m=l+n_s;
    for(int k=l; k<m; k++){
      double asv;
      double bsv;
      double lsv;
      double tstoch = distT(rng[tn*8]);
      
      if(flat==0) lsv=0;
      else lsv = exp(distL2(rng[tn*8]));
      
      double csv;// = distC(rng[tn*8]);
      double cv;
      
      if(vg==0){
        csv = distC(rng[tn*8]);
        do {
          asv = exp(distA(rng[tn*8]));
        } while (asv<(csv-vli)/14);
        do {
          bsv = exp(distB(rng[tn*8]));
        } while (bsv<(csv-vli)/28);
      }
      else if(vg==1){
        csv = distC(rng[tn*8]);
        do {
          asv = exp(distA(rng[tn*8]));
          bsv = exp(distB(rng[tn*8]))/asv;
        } while (asv<(csv-vli)/14 && bsv<(csv-vli)/28);
      }
      else if(vg==2){
        do {
          asv = exp(distA(rng[tn*8]));
          csv = distC(rng[tn*8])*asv;
          bsv = exp(distB(rng[tn*8]))/asv;
        } while (asv<(csv-vli)/14 && bsv<(csv-vli)/28);
      }
      else if(vg==3){
        csv = distC(rng[tn*8]);
        asv = exp(distA(rng[tn*8]));
        bsv = exp(distB(rng[tn*8]));
      }
      else if(vg==4){
        csv = distC(rng[tn*8]);
        asv = exp(distA(rng[tn*8]));
        bsv = exp(distB(rng[tn*8]))/asv;
      }
      else if(vg==5){
        asv = exp(distA(rng[tn*8]));
        csv = distC(rng[tn*8])*asv;
        bsv = exp(distB(rng[tn*8]))/asv;
      }
      else if(vg==6){
        csv = distC(rng[tn*8]);
        asv = a_lower+exp(distA(rng[tn*8]));
        bsv = b_lower+exp(distB(rng[tn*8]));
      }

      
      if(vk==1){
        /* time from crossing vl_thresh to peak */
        double t1=(csv-vthresh)/asv;
        for(int j=0;j<t_n;j++) { // 7 should be 0
          int t=(double) j;
          /* A peak below the threshold means this individual never triggers symptoms and is
             never observed; the threshold kernel guards the same case the same way. */
          if(vthresh > csv) cv = vli;
          else {
            if(j <= t1+tstoch) cv = vthresh+asv*(t-tstoch);
            else if(j <= t1+lsv+tstoch) cv = csv;
            else cv = csv-bsv*(t-t1-lsv-tstoch);
          }

          /* Clamped at the limit of detection, which is separate from the threshold: with
             vl_thresh fitted near zero the early trajectory sits below vli and is simply
             unobservable, exactly as in the threshold form. */
          if(cv<vli)
            cv=vli; //+1;
          else if(cv>vl_max)
            cv=vl_max;
          ft[tn][(int) (2*vl_max-2*cv)][j]++;
        }
      }
      
    }
  }
  
  /* ensure total simulated pop accounts for any rounding after division by ncores */
  double test_pop_d=(double) (n_s*ncores);
  
  NumericMatrix freq(t_n,ct_n); /* final freq matrix in Rcpp matrix format */
  
  /* sum thread specific frequencies and divide by total pop */
  for(int j=0;j<t_n;j++)
    for(int i=0;i<ct_n;i++) {
      double s=0;
      for(int tn=0;tn<ncores;tn++)
        s+=ft[tn][i][j];
      freq(j,i)=s/test_pop_d;
    }

  /* the per-thread generators were heap allocated above and were previously never released,
     leaking on every call -- which in an MCMC means every iteration. The bound checks that
     can throw are all placed before that allocation, so a plain delete here is sufficient. */
  delete[] rng;

  return freq;

}


