#include <admodel.h>
#include <contrib.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <gdbprintlib.cpp>

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
  a.allocate(1,nchol,-10.0,10.0,3,"a");
  u0.allocate(1,ngp,1,m-1,2,"u0");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  nll.allocate("nll");  /* ADOBJECTIVEFUNCTION */
}
void model_parameters::userfunction(void)
{
  nll =0.0;
  nll = 0.0;
  for(int i=1;i<=n;i++){
    //cout << u0(gp(i)) << endl << endl;
    MultiNomialSep(i, extract_column(beta0,1), u0(gp(i)), a);
  }
}

void SEPFUN1  model_parameters::MultiNomialSep(int i, const dvar_vector& beta0, const dvar_vector& u0, const dvar_vector& a)
{
  begin_df1b2_funnel();
  //cout << i << endl << endl;
  //----------------
  // RANDOM EFFECTS
  //----------------
  dvar_matrix L(1,m-1,1,m-1);   // Cholesky factor
  int k=1;
  for(int ii=1;ii<=(m-1);ii++){
    for(int j=1;j<=ii;j++){
      L(ii,j) = a(k);
      k++;
    }
  }
  dvar_matrix Sigma = L * trans(L);
  //cout << Sigma << endl << endl;
  // likelihood for the random effects - MVN
  nll += 0.5 * (log(2*M_PI) + ln_det(Sigma) + u0 * solve(Sigma,u0));
  //---------------
  // FIXED EFFECTS 
  //---------------
  // should be able to do this with matrix multiplication
  dvar_matrix eta0(1,1,1,m-1);
  for (int j=1;j<=m-1;j++){
    eta0(1,j) = sum(elem_prod(X(i),trans(beta0)(j))) + u0(j);
  }
  //dvar_matrix eta0 = X * beta0;
  // add a column of zeros to eta
  dvar_matrix etaconvert(1,m-1,1,m);
  for (int ii=1;ii<=m-1;ii++){
      etaconvert(ii,ii+1) = 1.0;
  }
  dvar_matrix eta =  eta0 * etaconvert;
  // get rowsums
  dvar_matrix ones(1,m,1,1);
  ones = ones + 1.0;
  dvar_matrix rowsums = exp(eta) * ones;
  // Probability matrix
  dvar_matrix P(1,1,1,m);
  for (int j=1;j<=m;j++){
    P(1,j) = exp(eta(1,j))/rowsums(1,1); 
  }
  //cout << "P" << endl << trans(P)(1) << endl << endl;  
  //cout << Y(1) << endl << endl;      
  //dvar_matrix P = elem_div(exp(eta), (1 + ));
  //
  nll += nllMultiNomial(Y(i), P(1));
  end_df1b2_funnel();
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
    initial_df1b2params::separable_flag=1;
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
  for(int i=1;i<=n;i++){
    //cout << u0(gp(i)) << endl << endl;
    MultiNomialSep(i, extract_column(beta0,1), u0(gp(i)), a);
  }
}

void   df1b2_pre_parameters::MultiNomialSep(int i, const funnel_init_df1b2vector& beta0, const funnel_init_df1b2vector& u0, const funnel_init_df1b2vector& a)
{
  begin_df1b2_funnel();
  //cout << i << endl << endl;
  //----------------
  // RANDOM EFFECTS
  //----------------
  df1b2matrix L(1,m-1,1,m-1);   // Cholesky factor
  int k=1;
  for(int ii=1;ii<=(m-1);ii++){
    for(int j=1;j<=ii;j++){
      L(ii,j) = a(k);
      k++;
    }
  }
  df1b2matrix Sigma = L * trans(L);
  //cout << Sigma << endl << endl;
  // likelihood for the random effects - MVN
  nll += 0.5 * (log(2*M_PI) + ln_det(Sigma) + u0 * solve(Sigma,u0));
  //---------------
  // FIXED EFFECTS 
  //---------------
  // should be able to do this with matrix multiplication
  df1b2matrix eta0(1,1,1,m-1);
  for (int j=1;j<=m-1;j++){
    eta0(1,j) = sum(elem_prod(X(i),trans(beta0)(j))) + u0(j);
  }
  //df1b2matrix eta0 = X * beta0;
  // add a column of zeros to eta
  df1b2matrix etaconvert(1,m-1,1,m);
  for (int ii=1;ii<=m-1;ii++){
      etaconvert(ii,ii+1) = 1.0;
  }
  df1b2matrix eta =  eta0 * etaconvert;
  // get rowsums
  df1b2matrix ones(1,m,1,1);
  ones = ones + 1.0;
  df1b2matrix rowsums = exp(eta) * ones;
  // Probability matrix
  df1b2matrix P(1,1,1,m);
  for (int j=1;j<=m;j++){
    P(1,j) = exp(eta(1,j))/rowsums(1,1); 
  }
  //cout << "P" << endl << trans(P)(1) << endl << endl;  
  //cout << Y(1) << endl << endl;      
  //df1b2matrix P = elem_div(exp(eta), (1 + ));
  //
  nll += nllMultiNomial(Y(i), P(1));
  end_df1b2_funnel();
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
  a.allocate(1,nchol,-10.0,10.0,3,"a");
  u0.allocate(1,ngp,1,m-1,2,"u0");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  nll.allocate("nll");  /* ADOBJECTIVEFUNCTION */
}
