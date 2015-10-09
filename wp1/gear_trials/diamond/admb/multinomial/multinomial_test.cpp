#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <multinomial_test.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  n.allocate("n");
  m.allocate("m");
  p.allocate("p");
  ngp.allocate("ngp");
  Y.allocate(1,n,1,m,"Y");
  X.allocate(1,n,1,p,"X");
  gp.allocate(1,n,"gp");
 nchol = (m-1)*(m)/2; // m-1 
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  beta0.allocate(1,p,1,m-1,"beta0");
  nll.allocate("nll");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  nll =0.0;
  // to add extra column of zeros for to beta for baseline
  dvar_matrix betaconvert(1,p,1,p+1);
  for (int i=1;i<=p;i++){
      betaconvert(i,i+1) = 1.0;
  }
  //cout << beta0 << endl << endl;
  //cout << betaconvert << endl;
  dvar_matrix beta = beta0 * betaconvert;
  //cout << beta << endl << endl;
  //dvar_matrix eta = X*beta;
  dvar_matrix eta(1,n,1,m);
  //dvar_matrix eta = X * beta;//+ u(gp(i)); 
  for (int i=1;i<=n;i++){
    for (int j=1;j<=m;j++){
      eta(i,j) = sum(elem_prod(X(i),trans(beta)(j)));
    }
  }
  // get rowsums
  dvar_matrix ones(1,m,1,1);
  ones = ones + 1.0;
  //cout << ones << endl << endl;  
  dvar_matrix rowsums = exp(eta) * ones;
  //cout << rowsums << endl << endl;  
  dvar_matrix P(1,n,1,m);
  for (int i=1;i<=n;i++){
    for (int j=1;j<=m;j++){
      P(i,j) = exp(eta(i,j))/rowsums(i,1); 
    }
  }
  //cout << P << endl << endl;  
  //cout << Y(1) << endl << endl;      
  //dvar_matrix P = elem_div(exp(eta), (1 + ));
  nll = 0.0;
  //
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
