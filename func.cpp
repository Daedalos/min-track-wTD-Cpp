// EDITED BY URI TO INCLUDE DISCRETE MAP/JACOBIAN
// discF, discDF


#define NX 6   	     // dim of state variable + number of parameters
#define ND 5         // dim of state variable
using namespace alglib;

void func_origin(real_1d_array &x, real_1d_array &func);
void func_DF(real_1d_array &x, real_2d_array &Jac);

void func_origin(real_1d_array &x, real_1d_array &func){
	//lorenz96
	//const ae_int_t D = x.length();
	//const double v = 8.17;
	int i;
	for(i=0;i<ND;i++){
		func[i] = x[(i-1+ND)%ND]*(x[(i+1)%ND]-x[(i-2+ND)%ND])-x[i]+x[ND];
	}

}

void func_DF(real_1d_array &x, real_2d_array &Jac){
	//lorenz96 Jacobian
	//const ae_int_t D = x.length();
	//const double v = 8.17;
	int i,j;
	for(i=0;i<ND;i++){
		for(j=0;j<ND;j++)
			Jac[i][j]=0;
		Jac[i][i] = -1;
		Jac[i][(i+1)%ND] = x[(i-1+ND)%ND];
		Jac[i][(i-2+ND)%ND] = - x[(i-1+ND)%ND];
		Jac[i][(i-1+ND)%ND] = x[(i+1)%ND] - x[(i-2+ND)%ND];
		Jac[i][ND]=1;
	}
	for(i=0;i<NX;i++) 
	  Jac[ND][i]=0;

}

// produces the vector-func map f(x) for given state vector (x)
// Currently uses RK2
// x should be state vector plus at one time point, f is output vector. same length DX
void discF(real_1d_array x, real_1d_array &f){
  real_1d_array k1, k2, xtmp;
  k1.setlength(NX);
  k2.setlength(NX);
  xtmp.setlength(NX);

  func_origin(x, k1);

  int i;
  for(i=0; i<ND; i++){
    xtmp[i] = x[i] + 0.5*DT*k1[i];
  }

  func_origin(xtmp, k2);

  for(i=0; i<ND; i++){
    f[i] = x[i] + DT*k2[i] ; 
  }
  
}

int delta(int i, int j){
  if(i==j){
    return 1;
  }
  else{
    return 0;
  }
}

void discDF(real_1d_array x, real_2d_array &Jac){

  int i, j, k;
  real_1d_array F, xmid;
  real_2d_array JacMid, Jac0;

  F.setlength(NX);
  func_origin(x,F);

  Jac0.setlength(NX,NX);
  func_DF(x,Jac0);

  JacMid.setlength(NX,NX);

  xmid.setlength(NX);
  for(k=0; k<ND; k++){
    xmid[k] = x[k] + 0.5*DT*F[k];	
  }
  func_DF(xmid,JacMid);
          
  for(i=0; i<NX; i++){
    for(j=0; j<NX; j++){
      Jac[i][j] = delta(i,j) + DT*JacMid[i][j]*( delta(i,j) + 0.5*DT*Jac0[i][j]);
    }
  }

  
}
