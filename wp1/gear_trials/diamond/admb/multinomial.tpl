
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
  matrix L(1,m-1,1,m-1)   // Cholesky factor
  init_vector a(1,nchol,3)   // Free parameters in C
  random_effects_matrix u0(1,ngp,1,m-1,2)
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
