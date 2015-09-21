DATA_SECTION
  init_int n
  init_int m
  init_int p
  init_int ngp // number of groups for random effects
  init_matrix Y(1,n,1,m)
  init_matrix X(1,n,1,p)
  int nchol
  !! nchol = (m-1)*(m)/2; // m-1 
  
//PRELIMINARY_CALCS_SECTION
  // number of Cholesky decomp parameters
  //int nchol = (m-1)*(m-2)/2;

PARAMETER_SECTION
  // fixed effect coefficients
  init_matrix beta0(1,p,1,m-1);
  // random effects vcov
  matrix L(1,m-1,1,m-1)   // Cholesky factor
  init_vector a(1,nchol)   // Free parameters in C
  random_effects_matrix u(1,ngp,1,m-1)
  objective_function_value nll;

PROCEDURE_SECTION
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
