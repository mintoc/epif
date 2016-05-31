//-------------------------------------------------------------
// ADMB-RE code for implementing multinomial with random effects implemented in: 
// A general catch comparison method for multi-gear trials: application to a quad-rig trawling fishery for Nephrops
// Notes:
// If getting aborted error: ret != -1
// clean temporary files from the directory
// Run with ./multinomialre -l1 50000000 -l2 200000000 -l3 50000000
//-------------------------------------------------------------

DATA_SECTION
  !! ofstream ofs("data.log"); // make sure data inputted correctly by examining data log file
  init_int n // number of rows
  !! ofs << "n" << endl << n <<  endl;
  init_int m // number of cod-ends/categories of response
  !! ofs << "m" << endl << m <<  endl;
  init_int p // dimension of the model matrix
  !! ofs << "p" << endl << p <<  endl;
  init_int q // dimension of the conditional variable - here same as m
  !! ofs << "q" << endl << q <<  endl;
  init_int ngp // number of groups for random effects
  !! ofs << "ngp" << endl << ngp <<  endl;
  init_matrix Y(1,n,1,m) // response counts
  !! ofs << "Y" << endl << Y <<  endl;
  init_matrix X(1,n,1,p) // model matrix
  !! ofs << "X" << endl << X <<  endl;
  init_matrix Xcond(1,n,1,q) // conditional variable matrix
  !! ofs << "Xcond" << endl << Xcond <<  endl;
  init_vector gp(1,n) // vector of groups
  !! ofs << "gp" << endl << gp <<  endl;
  init_matrix Offset(1,n,1,m) // offset matrix
  !! ofs << "Offset" << endl << Offset <<  endl;
  int nchol // number of parameters in the Cholesky decomposition of the random effect covariance matrix
  !! nchol = (m-1)*m/2;
  // prediction matrices
  init_int npred
  !! ofs << "npred" << endl << npred <<  endl;
  init_matrix Xpred(1,npred,1,p)
  !! ofs << "Xpred" << endl << Xpred <<  endl;
  init_matrix Xcondpred(1,npred,1,q)
  !! ofs << "Xcondpred" << endl << Xcondpred <<  endl;

PARAMETER_SECTION
  // fixed effect coefficients
  init_matrix beta0(1,p,1,m-1);
  // conditional logit coefficient
  init_number betacond; // only one implemented here
  // random effects vcov
  matrix L(1,m-1,1,m-1) // Cholesky factor
  init_bounded_vector a(1,nchol,-5.0,5.0,3)   // Free parameters in C


  // random effects
  random_effects_matrix u0(1,ngp,1,m-1,2)
  // linear predictor predictions
  sdreport_matrix etapred(1,npred,1,m);
  objective_function_value nll;

PROCEDURE_SECTION
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
  // log-likelihood for the random effects - MVN 
  for (int i=1;i<=ngp;i++){
    nll += 0.5 * (log(2*M_PI) + ln_det(Sigma) + u0(i) * solve(Sigma,u0(i)));
  }
  //---------------
  // FIXED EFFECTS 
  //---------------
  // matrix multiplication not working in my version of admb-re so done a long way
  dvar_matrix eta0(1,n,1,m-1);
  for (int i=1;i<=n;i++){
    for (int j=1;j<=m-1;j++){
      eta0(i,j) = sum(elem_prod(X(i),trans(beta0)(j))) + u0(gp(i),j);
    }
  }
  // add a column of zeros to eta
  dvar_matrix etaconvert(1,m-1,1,m);
  for (int i=1;i<=m-1;i++){
      etaconvert(i,i+1) = 1.0;
  }
  // linear predictor
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
  // log-likelihood for the observations
  for (int i=1;i<=n;i++){
      nll += -1. * (Y(i) * log(P(i)) + gammln(sum(Y(i)) + 1.) - sum(gammln(Y(i) + 1.)));
  }
  // sd report of linear predictor predictions
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
