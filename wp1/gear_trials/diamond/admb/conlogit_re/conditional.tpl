// note if getting aborted error: ret != -1
// clean temporary files from the directory
// also run with ./conditional -est -l1 50000000 -l2 200000000 -l3 50000000

DATA_SECTION
  !! ofstream ofs("data.log");
  init_int n
  !! ofs << "n" << endl << n <<  endl;
  init_int m
  !! ofs << "m" << endl << m <<  endl;
  init_int p
  !! ofs << "p" << endl << p <<  endl;
  init_int q
  !! ofs << "q" << endl << q <<  endl;
  init_int ngp // number of groups for random effects
  !! ofs << "ngp" << endl << ngp <<  endl;
  init_matrix Y(1,n,1,m)
  !! ofs << "Y" << endl << Y <<  endl;
  init_matrix X(1,n,1,p)
  !! ofs << "X" << endl << X <<  endl;
  init_matrix Xcond(1,n,1,q)
  !! ofs << "Xcond" << endl << Xcond <<  endl;
  init_vector gp(1,n)
  !! ofs << "gp" << endl << gp <<  endl;
  init_matrix Offset(1,n,1,m)
  !! ofs << "Offset" << endl << Offset <<  endl;
  int nchol
  !! nchol = (m-1)*m/2; 

PARAMETER_SECTION
  // fixed effect coefficients
  init_matrix beta0(1,p,1,m-1);
  // conditional logit coefficient
  init_number betacond(1); // only one implemented here
  // random effects vcov
  matrix L(1,m-1,1,m-1) // Cholesky factor
  init_bounded_vector a(1,nchol,-5.0,5.0,3)   // Free parameters in C
  //init_bounded_vector a(1,m-1,-5.0,5.0,3)   // Free parameters in C
  //init_bounded_vector a(1,m-1,-5.0,5.0,-1)   // Free parameters in C
  random_effects_matrix u0(1,ngp,1,m-1,2)
  //random_effects_matrix u0(1,ngp,1,m-1,-1)
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
