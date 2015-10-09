#include <admodel.h>
#include <contrib.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <multinomial.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
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
  model_parameters_ptr=this;
  initializationfunction();
  beta0.allocate(1,p,1,m-1,"beta0");
  L.allocate(1,m-1,1,m-1,"L");
  #ifndef NO_AD_INITIALIZE
    L.initialize();
  #endif
  a.allocate(1,nchol,3,"a");
  u0.allocate(1,ngp,1,m-1,2,"u0");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  nll.allocate("nll");  /* ADOBJECTIVEFUNCTION */
}
void model_parameters::userfunction(void)
{
  nll =0.0;
  nll = 0.0;
  //----------------
  // RANDOM EFFECTS
  //----------------
  int k=1;
  for(int i=1;i<=(m-1);i++){
    for(int j=1;j<=i;j++){
      L(i,j) = a(k);
      k++;
    }
  }
  dvar_matrix Sigma = L * trans(L);
  //cout << Sigma << endl << endl;
  // to add extra column of zeros for baseline random effect
  dvar_matrix uconvert(1,m-1,1,m);
  for (int i=1;i<=m-1;i++){
      uconvert(i,i+1) = 1.0;
  }
  //cout << "unconvert" << endl << uconvert << endl << endl;
  //cout << "u0" << endl << u0 << endl << endl;
  dvar_matrix u = u0 * uconvert;
  //cout << "u" << endl << u << endl << endl;
  // likelihood
  for (int i=1;i<=ngp;i++){
    //nll += 0.5 * ln_det(Sigma) + 0.5 * u0(i) * solve(Sigma,u0(i));
    nll += 0.5 * (log(2*M_PI) + ln_det(Sigma) + u0(i) * solve(Sigma,u0(i)));
  }
  //---------------
  // FIXED EFFECTS 
  //---------------
  // to add extra column of zeros for to beta for baseline
  dvar_matrix betaconvert(1,m-1,1,m);
  for (int i=1;i<=(m-1);i++){
      betaconvert(i,i+1) = 1.0;
  }
  //cout << beta0 << endl << endl;
  //cout << betaconvert << endl;
  dvar_matrix beta =  beta0 * betaconvert;
  //cout << "beta" << endl << beta << endl << endl;
  //dvar_matrix eta = X * beta;//+ u(gp(i)); 
  dvar_matrix eta(1,n,1,m);
  for (int i=1;i<=n;i++){
    for (int j=1;j<=m;j++){
      eta(i,j) = sum(elem_prod(X(i),trans(beta)(j))) + u(gp(i),j);
    }
  }
  //cout << eta << endl << endl;
  //cout << u << endl << endl;
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
  //cout << "P" << endl << trans(P)(1) << endl << endl;  
  //cout << Y(1) << endl << endl;      
  //dvar_matrix P = elem_div(exp(eta), (1 + ));
  //
  for (int i=1;i<=n;i++){
      nll += nllMultiNomial(Y(i), P(i));
  }
}
  long int arrmblsize=0;

int main(int argc,char * argv[])
{
  ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
      if (!arrmblsize) arrmblsize=150000;
    df1b2variable::noallocate=1;
    df1b2_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;

    function_minimizer::random_effects_flag=1;
    df1b2variable::noallocate=0;
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

void df1b2_parameters::user_function(void)
{
  nll =0.0;
  nll = 0.0;
  //----------------
  // RANDOM EFFECTS
  //----------------
  int k=1;
  for(int i=1;i<=(m-1);i++){
    for(int j=1;j<=i;j++){
      L(i,j) = a(k);
      k++;
    }
  }
  df1b2matrix Sigma = L * trans(L);
  //cout << Sigma << endl << endl;
  // to add extra column of zeros for baseline random effect
  df1b2matrix uconvert(1,m-1,1,m);
  for (int i=1;i<=m-1;i++){
      uconvert(i,i+1) = 1.0;
  }
  //cout << "unconvert" << endl << uconvert << endl << endl;
  //cout << "u0" << endl << u0 << endl << endl;
  df1b2matrix u = u0 * uconvert;
  //cout << "u" << endl << u << endl << endl;
  // likelihood
  for (int i=1;i<=ngp;i++){
    //nll += 0.5 * ln_det(Sigma) + 0.5 * u0(i) * solve(Sigma,u0(i));
    nll += 0.5 * (log(2*M_PI) + ln_det(Sigma) + u0(i) * solve(Sigma,u0(i)));
  }
  //---------------
  // FIXED EFFECTS 
  //---------------
  // to add extra column of zeros for to beta for baseline
  df1b2matrix betaconvert(1,m-1,1,m);
  for (int i=1;i<=(m-1);i++){
      betaconvert(i,i+1) = 1.0;
  }
  //cout << beta0 << endl << endl;
  //cout << betaconvert << endl;
  df1b2matrix beta =  beta0 * betaconvert;
  //cout << "beta" << endl << beta << endl << endl;
  //df1b2matrix eta = X * beta;//+ u(gp(i)); 
  df1b2matrix eta(1,n,1,m);
  for (int i=1;i<=n;i++){
    for (int j=1;j<=m;j++){
      eta(i,j) = sum(elem_prod(X(i),trans(beta)(j))) + u(gp(i),j);
    }
  }
  //cout << eta << endl << endl;
  //cout << u << endl << endl;
  // get rowsums
  df1b2matrix ones(1,m,1,1);
  ones = ones + 1.0;
  //cout << ones << endl << endl;  
  df1b2matrix rowsums = exp(eta) * ones;
  //cout << rowsums << endl << endl;  
  df1b2matrix P(1,n,1,m);
  for (int i=1;i<=n;i++){
    for (int j=1;j<=m;j++){
      P(i,j) = exp(eta(i,j))/rowsums(i,1); 
    }
  }
  //cout << "P" << endl << trans(P)(1) << endl << endl;  
  //cout << Y(1) << endl << endl;      
  //df1b2matrix P = elem_div(exp(eta), (1 + ));
  //
  for (int i=1;i<=n;i++){
      nll += nllMultiNomial(Y(i), P(i));
  }
}
   
void df1b2_pre_parameters::setup_quadprior_calcs(void) 
{ 
  df1b2_gradlist::set_no_derivatives(); 
  quadratic_prior::in_qp_calculations=1; 
}  
  
void df1b2_pre_parameters::begin_df1b2_funnel(void) 
{ 
  (*re_objective_function_value::pobjfun)=0; 
  other_separable_stuff_begin(); 
  f1b2gradlist->reset();  
  if (!quadratic_prior::in_qp_calculations) 
  { 
    df1b2_gradlist::set_yes_derivatives();  
  } 
  funnel_init_var::allocate_all();  
}  
 
void df1b2_pre_parameters::end_df1b2_funnel(void) 
{  
  lapprox->do_separable_stuff(); 
  other_separable_stuff_end(); 
} 
  
void model_parameters::begin_df1b2_funnel(void) 
{ 
  if (lapprox)  
  {  
    {  
      begin_funnel_stuff();  
    }  
  }  
}  
 
void model_parameters::end_df1b2_funnel(void) 
{  
  if (lapprox)  
  {  
    end_df1b2_funnel_stuff();  
  }  
} 

void df1b2_parameters::allocate(void) 
{
  beta0.allocate(1,p,1,m-1,"beta0");
  L.allocate(1,m-1,1,m-1,"L");
  #ifndef NO_AD_INITIALIZE
    L.initialize();
  #endif
  a.allocate(1,nchol,3,"a");
  u0.allocate(1,ngp,1,m-1,2,"u0");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  nll.allocate("nll");  /* ADOBJECTIVEFUNCTION */
}
