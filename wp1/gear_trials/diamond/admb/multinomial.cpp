#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <multinomial.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  n.allocate("n");
  m.allocate("m");
  p.allocate("p");
  Y.allocate(1,n,1,m,"Y");
  X.allocate(1,n,1,p,"X");
 nchol = (m-1)*(m)/2; // m-1 
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  beta0.allocate(1,p,1,m-1,"beta0");
  L.allocate(1,m-1,1,m-1,"L");
  #ifndef NO_AD_INITIALIZE
    L.initialize();
  #endif
  a.allocate(1,nchol,"a");
  nll.allocate("nll");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  nll =0.0;
  //----------------
  // RANDOM EFFECTS
  //----------------
  int k=0;
  L(1,1) = 1.0;
  for(int i=2;i<=(m-1);i++){
    for(int j=1;j<=i;j++){
      L(i,j) = a(k++);
    }
  }
  dvar_matrix Sigma = L * trans(L);
  //---------------
  // FIXED EFFECTS 
  //---------------
  // to add extra column of zeros for to beta for baseline
  dvar_matrix betaconvert(1,p,1,p+1);
  for (int i=1;i<=p;i++){
      betaconvert(i,i+1) = 1.0;
  }
  dvar_matrix beta = beta0 * betaconvert;
  // linear predictor
  dvar_matrix eta = X*beta;
  // get rowsums for probabilities
  dvar_matrix ones(1,m,1,1);
  ones = ones + 1.0;
  dvar_matrix rowsums = exp(eta) * ones;
  dvar_matrix P(1,n,1,m);
  for (int i=1;i<=n;i++){
    for (int j=1;j<=m;j++){
      P(i,j) = exp(eta(i,j))/rowsums(i,1); 
    }
  }
  nll = 0.0;
  for (int i=1;i<=n;i++){
      nll += nllMultiNomial(Y(i), P(i));
  }
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::report(const dvector& gradients){}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
