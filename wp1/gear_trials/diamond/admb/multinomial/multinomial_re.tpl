DATA_SECTION
  init_int n
  init_int m
  init_int p
  init_int ngp // number of groups for random effects
  init_matrix Y(1,n,1,m)
  init_matrix X(1,n,1,p)
  init_vector gp(1,n)
  int nchol
  !! nchol = (m-1)*(m)/2; // m-1 
  
//PRELIMINARY_CALCS_SECTION
  // number of Cholesky decomp parameters
  //int nchol = (m-1)*(m-2)/2;

PARAMETER_SECTION
  // fixed effect coefficients
  init_matrix beta0(1,p,1,m-1);
  init_vector a(1,nchol)   // Free parameters in C
  random_effects_matrix u0(1,ngp,1,m-1)
  // random effects vcov
  matrix L(1,m-1,1,m-1)   // Cholesky factor
  objective_function_value nll;

PROCEDURE_SECTION
  nll = 0.0;
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
  // to add extra column of zeros for baseline random effect
  dvar_matrix uconvert(1,m-1,1,m);
  for (int i=1;i<=m-1;i++){
      uconvert(i,i+1) = 1.0;
  }
  cout << Sigma << endl << endl;
  dvar_matrix u = u0 * uconvert;
  // likelihood - using u0 here
  for (int i=1;i<=n;i++){
    nll -= 0.5 * ln_det(Sigma) + 0.5 * u0(i) * solve(Sigma,u0(i));
  }
  //---------------
  // FIXED EFFECTS 
  //---------------
  // to add extra column of zeros for to beta for baseline
  dvar_matrix betaconvert(1,p,1,p+1);
  for (int i=1;i<=p;i++){
      betaconvert(i,i+1) = 1.0;
  }
  dvar_matrix beta = beta0 * betaconvert;
  //------------------
  // LINEAR PREDICTOR 
  //------------------
  // problem here, matrix mult not implemented for random effects?
  // dvar_matrix eta = X*beta;
  // looping through instead
  dvar_matrix eta(1,n,1,m);
  for (int i=1;i<=n;i++){
    eta(i) = X(i)*beta + u(gp(i)); 
  }
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
  //
  for (int i=1;i<=n;i++){
      nll += nllMultiNomial(Y(i), P(i));
  }
