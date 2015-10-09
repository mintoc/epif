
DATA_SECTION
  init_int n
  init_int m
  init_int p
  init_int ngp // number of groups for random effects
  init_matrix Y(1,n,1,m)
  init_matrix X(1,n,1,p)
  init_vector gp(1,n)
  //!! cout << "gp " << gp<< endl << endl;
  int nchol
  !! nchol = (m-1)*(m)/2; // m-1 

PARAMETER_SECTION
  // fixed effect coefficients
  init_matrix beta0(1,p,1,m-1);
  // random effects vcov
  init_bounded_vector a(1,nchol,-10.0,10.0,3)   // Free parameters in C
  random_effects_matrix u0(1,ngp,1,m-1,2)
  objective_function_value nll;

PROCEDURE_SECTION
  nll = 0.0;
  for(int i=1;i<=n;i++){
    //cout << u0(gp(i)) << endl << endl;
    MultiNomialSep(i, beta0, u0(gp(i)), a);
  }

SEPARABLE_FUNCTION void MultiNomialSep(int i, const dvar_matrix& beta0, const dvar_vector& u0, const dvar_vector& a)
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
