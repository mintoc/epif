DATA_SECTION
  init_int n
  init_int m
  init_int p
  init_matrix Y(1,n,1,m)
  init_matrix X(1,n,1,p)

PARAMETER_SECTION
  init_matrix beta0(1,p,1,m-1);
  objective_function_value nll;

PROCEDURE_SECTION
  // to add extra column of zeros for to beta for baseline
  dvar_matrix betaconvert(1,p,1,p+1);
  for (int i=1;i<=p;i++){
      betaconvert(i,i+1) = 1.0;
  }
  //cout << beta0 << endl << endl;
  //cout << betaconvert << endl;
  dvar_matrix beta = beta0 * betaconvert;
  //cout << beta << endl << endl;
  dvar_matrix eta = X*beta;
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
