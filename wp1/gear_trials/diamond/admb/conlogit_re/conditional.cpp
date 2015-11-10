#include <admodel.h>
#include <contrib.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <conditional.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
 ofstream ofs("data.log");
  n.allocate("n");
 ofs << "n" << endl << n <<  endl;
  m.allocate("m");
 ofs << "m" << endl << m <<  endl;
  p.allocate("p");
 ofs << "p" << endl << p <<  endl;
  q.allocate("q");
 ofs << "q" << endl << q <<  endl;
  ngp.allocate("ngp");
 ofs << "ngp" << endl << ngp <<  endl;
  Y.allocate(1,n,1,m,"Y");
 ofs << "Y" << endl << Y <<  endl;
  X.allocate(1,n,1,p,"X");
 ofs << "X" << endl << X <<  endl;
  Xcond.allocate(1,n,1,q,"Xcond");
 ofs << "Xcond" << endl << Xcond <<  endl;
  gp.allocate(1,n,"gp");
 ofs << "gp" << endl << gp <<  endl;
  Offset.allocate(1,n,1,m,"Offset");
 ofs << "Offset" << endl << Offset <<  endl;
 nchol = (m-1)*m/2;
  npred.allocate("npred");
 ofs << "npred" << endl << npred <<  endl;
  Xpred.allocate(1,npred,1,p,"Xpred");
 ofs << "Xpred" << endl << Xpred <<  endl;
  Xcondpred.allocate(1,npred,1,q,"Xcondpred");
 ofs << "Xcondpred" << endl << Xcondpred <<  endl;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  model_parameters_ptr=this;
  initializationfunction();
  beta0.allocate(1,p,1,m-1,"beta0");
  betacond.allocate("betacond");
  L.allocate(1,m-1,1,m-1,"L");
  #ifndef NO_AD_INITIALIZE
    L.initialize();
  #endif
  a.allocate(1,nchol,-5.0,5.0,3,"a");
  u0.allocate(1,ngp,1,m-1,2,"u0");
  etapred.allocate(1,npred,1,m,"etapred");
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
      //L(i,i) = a(k);
      k++;
    }
  }
  dvar_matrix Sigma = L * trans(L);
  //cout << Sigma << endl << endl;
  // likelihood for the random effects - MVN 
  for (int i=1;i<=ngp;i++){
    nll += 0.5 * (log(2*M_PI) + ln_det(Sigma) + u0(i) * solve(Sigma,u0(i)));
  }
  //---------------
  // FIXED EFFECTS 
  //---------------
  // should be able to do this with matrix multiplication
  dvar_matrix eta0(1,n,1,m-1);
  for (int i=1;i<=n;i++){
    for (int j=1;j<=m-1;j++){
      eta0(i,j) = sum(elem_prod(X(i),trans(beta0)(j))) + u0(gp(i),j);
    }
  }
  //dvar_matrix eta0 = X * beta0;
  // add a column of zeros to eta
  dvar_matrix etaconvert(1,m-1,1,m);
  for (int i=1;i<=m-1;i++){
      etaconvert(i,i+1) = 1.0;
  }
  dvar_matrix eta =  eta0 * etaconvert + Offset + betacond * Xcond;
  // get rowsums
  dvar_matrix ones(1,m,1,1);
  ones = ones + 1.0;
  dvar_matrix rowsums = exp(eta) * ones;
  // Probability matrix
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
      //nll += nllMultiNomial(Y(i), P(i));
      nll += -1. * (Y(i) * log(P(i)) + gammln(sum(Y(i)) + 1.) - sum(gammln(Y(i) + 1.)));
  }
  // sd report
  if (sd_phase())
  {
    dvar_matrix etapred0(1,npred,1,m-1);
    for (int i=1;i<=npred;i++){
      for (int j=1;j<=m-1;j++){
        etapred0(i,j) = sum(elem_prod(Xpred(i),trans(beta0)(j)));
      }
    }
    etapred =  etapred0 * etaconvert + betacond * Xcondpred;
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
      //L(i,i) = a(k);
      k++;
    }
  }
  df1b2matrix Sigma = L * trans(L);
  //cout << Sigma << endl << endl;
  // likelihood for the random effects - MVN 
  for (int i=1;i<=ngp;i++){
    nll += 0.5 * (log(2*M_PI) + ln_det(Sigma) + u0(i) * solve(Sigma,u0(i)));
  }
  //---------------
  // FIXED EFFECTS 
  //---------------
  // should be able to do this with matrix multiplication
  df1b2matrix eta0(1,n,1,m-1);
  for (int i=1;i<=n;i++){
    for (int j=1;j<=m-1;j++){
      eta0(i,j) = sum(elem_prod(X(i),trans(beta0)(j))) + u0(gp(i),j);
    }
  }
  //df1b2matrix eta0 = X * beta0;
  // add a column of zeros to eta
  df1b2matrix etaconvert(1,m-1,1,m);
  for (int i=1;i<=m-1;i++){
      etaconvert(i,i+1) = 1.0;
  }
  df1b2matrix eta =  eta0 * etaconvert + Offset + betacond * Xcond;
  // get rowsums
  df1b2matrix ones(1,m,1,1);
  ones = ones + 1.0;
  df1b2matrix rowsums = exp(eta) * ones;
  // Probability matrix
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
      //nll += nllMultiNomial(Y(i), P(i));
      nll += -1. * (Y(i) * log(P(i)) + gammln(sum(Y(i)) + 1.) - sum(gammln(Y(i) + 1.)));
  }
  // sd report
  if (sd_phase())
  {
    df1b2matrix etapred0(1,npred,1,m-1);
    for (int i=1;i<=npred;i++){
      for (int j=1;j<=m-1;j++){
        etapred0(i,j) = sum(elem_prod(Xpred(i),trans(beta0)(j)));
      }
    }
    etapred =  etapred0 * etaconvert + betacond * Xcondpred;
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
  betacond.allocate("betacond");
  L.allocate(1,m-1,1,m-1,"L");
  #ifndef NO_AD_INITIALIZE
    L.initialize();
  #endif
  a.allocate(1,nchol,-5.0,5.0,3,"a");
  u0.allocate(1,ngp,1,m-1,2,"u0");
  etapred.allocate(1,npred,1,m,"etapred");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  nll.allocate("nll");  /* ADOBJECTIVEFUNCTION */
}
