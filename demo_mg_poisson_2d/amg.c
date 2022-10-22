/* Classic algebraic multigrid (C-AMG)
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "sparse.h"
#include "solver.h"

icsr_t icsr_transpose(icsr_t S);
csr_t csr_transpose(csr_t P);
void icsr_close(icsr_t *S);
void csr_close(csr_t *P);

void gauss_jordon(int n, double *a, int m, double *b);
void solve_gmres(int n, double *x, double *b, op_t Aop, int niter, double tol, int m, int verb);

double dotprod(int n, double *a, double *b);

//This routine will set: is_cg[i]=1, C-point; is_cg[i]=0, F-point
icsr_t mark_strong_connection(double theta, csr_t A)
{
  int i, j, k;
  double threshold;
  icsr_t S;

  S.nrow = A.nrow;
  S.ncol = A.ncol;
  S.nnz = A.nnz;
  S.row_ptr = alloc1int(S.nrow+1);
  S.col_ind = alloc1int(S.nnz);
  S.val = alloc1int(S.nnz);

  memcpy(S.row_ptr, A.row_ptr, (S.nrow+1)*sizeof(int));
  memcpy(S.col_ind, A.col_ind, S.nnz*sizeof(int));
  for(i=0; i<A.nrow; i++){
    //first pass in row-i, find the threshold
    threshold = 0;//set it to be the first value in row i
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];
      if(j!=i) threshold = MAX(threshold, fabs(A.val[k]));//max_{k!=i}{-a_ik}
    }
    threshold *= theta;//threshold=theta* max_{k!=i} {-a_ik}

    //second pass in row-i, set up the connection matrix M:
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];
      S.val[k] = (-A.val[k]>=threshold)?1:0;//1 for strong dependence; 0 for weak dependence
    }    
  }

  return S;
}

//coloring points into C/F; called C/F splitting
//see Wagner 1998 p50-51, algorithm 4.6.1 Ruge-Stuben algorithm RS_Coarsen1
//this is to satisfy property 4.6.2
void cf_splitting_first_pass(icsr_t S, icsr_t St, int *is_cgp)
{
  int i, j, k, k1, k2, maxlam;
  int *lambda, nvisited;

  lambda = alloc1int(St.nrow);//lambda_i, cardinanlity of S_i^t, importance measure to be C-point
  for(i=0; i<St.nrow; i++) {
    is_cgp[i] = -1;//set all point unvisited with flag=-1
    lambda[i] = 0;
    for(k=St.row_ptr[i]; k<St.row_ptr[i+1]; k++){
      j = St.col_ind[k];
      if(St.val[k]) lambda[i]++;//z_i=#{S_i^t}
    }
  }

  //In this pass, is_cgp[:] will mark C-point (coarse grid) as 1, and F-point (fine grid) as 0
  //is_cgp[i]=-1, unvisited; is_cgp[i]=1, visited, C-point; is_cgp[i]=0, visited, F-point
  nvisited = 0;//all n points are unvisited
  while(nvisited<S.nrow){
    //find lambda[i]=max{lambda[:]}, everytime O(n) complexity, not optimal
    i = 0;
    maxlam = lambda[i];
    for(j=0; j<S.nrow; j++) {
      if(maxlam<lambda[j]){
	maxlam = lambda[j];
	i = j;
      }
    }
    is_cgp[i] = 1;//point i goes to coarse grid, set flag=1
    nvisited++;//point i visited, number of unvisited points decrease 1
    lambda[i] = -1;    

    for(k1=St.row_ptr[i]; k1<St.row_ptr[i+1]; k1++){
      j = St.col_ind[k1];
      if(St.val[k1] && is_cgp[j]==-1){//unvisited points with flag=-1, (j in St_i\cap U)
	is_cgp[j] = 0;//point j should go to fine grid, set flag=0 (F=F\cup{j}, U=U\{j})
	nvisited++;//point j has been visited
	lambda[j] = -1;
      
	for(k2=S.row_ptr[j]; k2<S.row_ptr[j+1]; k2++){
	  k = S.col_ind[k2];
	  if(S.val[k2] && is_cgp[k]==-1) lambda[k]++;//cardinality increased by 1, k in S_j\cap U
	}
      }//end if
    }//end for k1

    for(k1=S.row_ptr[i]; k1<S.row_ptr[i]; k1++){
      j = S.col_ind[k1];
      if(S.val[k1] && is_cgp[j]==-1) lambda[j]--;//j in Si\cap U
    }
  }//end while

  free1int(lambda);
}

void cf_splitting_refine_pass(icsr_t S, icsr_t St, int *is_cgp)
{
  int i, j, k, k1, k2, k3;
  int connect;
  
  for(i=0; i<S.nrow; i++){
    if(!is_cgp[i]){//i in F
      for(k1=S.row_ptr[i]; k1<S.row_ptr[i+1]; k1++){
	j = S.col_ind[k1];
	if(!is_cgp[j] && S.val[k1]){// j in F, S_ij=1
	  connect = 0;
	  for(k2=S.row_ptr[i]; k2<S.row_ptr[i+1]; k2++){
	    k = S.col_ind[k2];
	    if(S.val[k2] && is_cgp[k]){//S_ik=1, k in C
	      for(k3=St.row_ptr[k]; k3<St.row_ptr[k+1]; k3++){
	  	if(j==St.col_ind[k3] && St.val[k3]) {//St_kj=S_jk=1
	  	  connect = 1;//found k is strongly connected to both i and j, k in C
	  	  break;
	  	}//end if
	      }//end for k3
	      if(connect) break;//found k is strongly connected to both i and j, k in C
	    }
	  }

	  //if no common C-point found strongly connected to F-point i and F-point j
	  if(!connect) {
	    is_cgp[i] = 1;//then, i goes from F to C; j remains in F
	    continue;//check i+1
	  }
	}
      }//end for k1
    }//end if
  }//end for i

}


//build interpolation matrix P
csr_t build_P(csr_t A, double theta, int interp)
{
  int i, j, k, m, k1, k2, k3, nnz, ntotal, ncoarse;
  int *is_cgp, *id_cgp, *in_Ci, *in_Dsi;
  double a_mj, sum_amk, di, dj;
  double s1p, s2p, s12p, s1n, s2n, s12n;
  icsr_t S, St;
  csr_t P;

  ntotal = A.ncol;//should be the same as A.nrow
  is_cgp = alloc1int(ntotal);//flag of coarse grid point, 1=C-point, 0=F-point
  id_cgp = alloc1int(ntotal);//coarse grid id

  S = mark_strong_connection(theta, A);
  St = icsr_transpose(S);
  cf_splitting_first_pass(S, St, is_cgp);
  cf_splitting_refine_pass(S, St, is_cgp);

  ncoarse = 0;//count the number of coarse grid point
  for(i=0; i<ntotal; i++){
    if(is_cgp[i]) {
      id_cgp[i] = ncoarse;//point i in C, mark its coarse-grid-index
      ncoarse++;//count total number of C-point
    }else{
      id_cgp[i] = -1;//point i in F, initialize with -1 (invalid index)
    }
  }
  
  P.ncol = ncoarse;
  P.nrow = ntotal;
  nnz = 0;
  for(i=0; i<ntotal; i++){
    if(is_cgp[i]){//i in C: (Pe)_i = e_i, P_ii=1
      nnz++;
    }else{//i in F: (Pe)_i=\sum_j w_ij e_j, if j in Ci, i.e., P_ij=w_ij
      for(k1=S.row_ptr[i]; k1<S.row_ptr[i+1]; k1++){
	j = S.col_ind[k1];//a_ij
	if(is_cgp[j] && S.val[k1]) {//j in Ci
	  nnz++;
	}
      }//end for k1      
    }//end if
  }//end for i
  P.nnz = nnz;
  P.row_ptr = alloc1int(P.nrow+1);
  P.col_ind = alloc1int(P.nnz);
  P.val = alloc1double(P.nnz);

  if(interp==1){//1=standard interpolation
    /*--------------------------------------------------------------------*/
    //compute w_ij, see eqn 8.12 in Briggs 2000
    in_Ci = alloc1int(A.ncol);
    in_Dsi = alloc1int(A.ncol);
    memset(in_Ci, 0, A.ncol*sizeof(int));
    memset(in_Dsi, 0, A.ncol*sizeof(int));
    nnz = 0;
    P.row_ptr[0] = 0;
    for(i=0; i<ntotal; i++){
      if(is_cgp[i]){//i in C, (Pe)_i = e_i
	P.col_ind[nnz] = id_cgp[i];//point i has index id_cgp[i] on coarse grid
	P.val[nnz] = 1.;//P_ii=1
	nnz++;
      }else{//i in F: (Pe)_i=\sum_j w_ij e_j, if j in Ci, i.e., P_ij=w_ij
	di = 0;
	for(k1=A.row_ptr[i]; k1<A.row_ptr[i+1]; k1++){
	  j = A.col_ind[k1];//a_ij
	  if(is_cgp[j] && S.val[k1]) {//j in Ci
	    in_Ci[j] = 1;
	  }
	  if(!is_cgp[j] && S.val[k1]){//j in Dsi (m in F)
	    in_Dsi[j] = 1;
	  }
	  if(j==i) di += A.val[k1];//di += a_ii
	  else if(!S.val[k1]) di += A.val[k1];//di += a_ij, j in Dwi
	}

	for(k1=A.row_ptr[i]; k1<A.row_ptr[i+1]; k1++){
	  j = A.col_ind[k1];//a_ij
	  if(in_Ci[j]) {//j in Ci
	    dj = A.val[k1];//dj = a_ij, j in Ci

	    for(k2=A.row_ptr[i]; k2<A.row_ptr[i+1]; k2++){
	      m = A.col_ind[k2];//a_im=A.val[k2]
	      if(in_Dsi[m]){//m in Dsi (m in F)
		sum_amk = 0.;
		a_mj = 0.;
		for(k3=A.row_ptr[m]; k3<A.row_ptr[m+1]; k3++){
		  k = A.col_ind[k3];
		  if(in_Ci[k]) sum_amk += A.val[k3];//sum_{k in Ci} a_mk
		  if(k==j) a_mj = A.val[k3];//a_mj=A.val[k3];
		}//end for k3
		dj += A.val[k2]*a_mj/sum_amk;//a_im=A.val[k2]
	      }//end if
	    }//end for k2

	    P.col_ind[nnz] = id_cgp[j];
	    P.val[nnz] = -dj/di;
	    //printf("wij=%e wij_den=%e\n", P.val[nnz], wij_den);
	    nnz++;
	  }//end if
	}//end for k1
      
	for(k1=A.row_ptr[i]; k1<A.row_ptr[i+1]; k1++){
	  j = A.col_ind[k1];//a_ij
	  //reset flag back to 0 for next row
	  if(in_Ci[j]) in_Ci[j] = 0;//j in Ci
	  if(in_Dsi[j])in_Dsi[j] = 0;//j in Dsi
	}//end for k1

      }//end if
      P.row_ptr[i+1] = nnz;
    }//end for i

    free1int(in_Ci);
    free1int(in_Dsi);

  }else if(interp==2){//2=direct interpolation
    nnz = 0;
    P.row_ptr[0] = 0;
    for(i=0; i<ntotal; i++){
      if(is_cgp[i]){//i in C, (Pe)_i = e_i
	P.col_ind[nnz] = id_cgp[i];//point i has index id_cgp[i] on coarse grid
	P.val[nnz] = 1.;//P_ii=1
	nnz++;
      }else{//i in F: (Pe)_i=\sum_j w_ij e_j, if j in Ci, i.e., P_ij=w_ij
	s1p = 0;
	s1n = 0;
	s2p = 0;
	s2n = 0;
	di = 0;
	for(k1=A.row_ptr[i]; k1<A.row_ptr[i+1]; k1++){
	  j = A.col_ind[k1];//a_ij
	  if(A.val[k1]>0) s1p += A.val[k1];
	  else s1n += A.val[k1];
	  if(is_cgp[j] && S.val[k1]) {//j in Ci
	    if(A.val[k1]>0) s2p += A.val[k1];//sum over postive terms
	    else s2n += A.val[k1];//sum over negative terms
	  }
	  if(j==i) di = A.val[k1];//di = a_ii
	}
	s12p = s1p/s2p;
	s12n = s1n/s2n;
	for(k1=A.row_ptr[i]; k1<A.row_ptr[i+1]; k1++){
	  j = A.col_ind[k1];//a_ij
	  if(is_cgp[j] && S.val[k1]) {//j in Ci
	    P.col_ind[nnz] = id_cgp[j];
	    if(A.val[k1]>0) P.val[nnz] = -s12p*A.val[k1]/di;//positive term
	    else P.val[nnz] = -s12n*A.val[k1]/di;//negative term
	    //printf("wij=%e wij_den=%e\n", P.val[nnz], wij_den);
	    nnz++;
	  }
	}//end for k1

      }//end if
      P.row_ptr[i+1] = nnz;
    }//end for i

  }
  
  icsr_close(&S);
  icsr_close(&St);
  free1int(is_cgp);
  free1int(id_cgp);
  
  return P;
}


//CSR matrix multiplication: T=RAP, T_mn += R_mi*A_ij P_jn
csr_t build_RAP(csr_t R, csr_t A, csr_t P)
{
  int i, j, m, n, k;
  int k1, k2, k3, *pos;
  double Rmi, Rmi_Aij;
  csr_t T;

  T.nrow = R.nrow;
  T.ncol = P.ncol;
  T.row_ptr = alloc1int(T.nrow+1);
  pos = alloc1int(T.ncol);
  
  for(i=0; i<T.ncol; i++) pos[i] = -1;
  k = 0;
  for(m=0; m<R.nrow; m++){
    T.row_ptr[m] = k;//row_ptr[m] stores total number of nonzeros in previous (m-1) rows
    for(k1=R.row_ptr[m]; k1<R.row_ptr[m+1]; k1++){//R_mi
      i = R.col_ind[k1];
      for(k2=A.row_ptr[i]; k2<A.row_ptr[i+1]; k2++){//A_ij
	j = A.col_ind[k2];
	for(k3=P.row_ptr[j]; k3<P.row_ptr[j+1]; k3++){//P_jn
	  n = P.col_ind[k3];
	  if(pos[n]<T.row_ptr[m]){//this is a new column because the pos index is not up to date
	    pos[n] = k;//T_mn is the k-th element of matrix T
	    k++;//increase global element counter for matrix T
	  }
	}//end for k3
      }//end for k2
    }//end for k1
  }//end for m
  T.row_ptr[T.nrow] = k;
  T.nnz = k;
  T.col_ind = alloc1int(T.nnz);
  T.val = alloc1double(T.nnz);

  for(i=0; i<T.ncol; i++) pos[i] = -1;
  k = 0;
  for(m=0; m<R.nrow; m++){
    T.row_ptr[m] = k;//row_ptr[m] stores total number of nonzeros in previous (m-1) rows
    for(k1=R.row_ptr[m]; k1<R.row_ptr[m+1]; k1++){//R_mi
      i = R.col_ind[k1];
      Rmi = R.val[k1];
      for(k2=A.row_ptr[i]; k2<A.row_ptr[i+1]; k2++){//A_ij
	j = A.col_ind[k2];
	Rmi_Aij = Rmi*A.val[k2];
	for(k3=P.row_ptr[j]; k3<P.row_ptr[j+1]; k3++){//P_jn
	  n = P.col_ind[k3];	  
	  if(pos[n]<T.row_ptr[m]){//this is a new column because the pos index is not up to date
	    T.col_ind[k] = n;//record column index for k-th element of T
	    pos[n] = k;//T_mn is the k-th element of matrix T
	    k++;//increase global counter for the elements of matrix T
	  }
	  T.val[pos[n]] += Rmi_Aij*P.val[k3];//T_mn += R_mi*A_ij P_jn
	}//end for k3
      }//end for k2
    }//end for k1
  }//end for m
  
  free1int(pos);

  return T;
}

void gauss_seidel(int n, csr_t A, double *x, double *b)
{
  int i, j, k;
  double aii;
  
  for(i=0; i<A.nrow; i++){
    x[i] = b[i];
    aii = 1;
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];
      if(j!=i) x[i] -= A.val[k]*x[j];// x[i] -= a_ij*x_j, (j!=i)
      else aii = A.val[k];
    }
    x[i] /= aii;
  }
}

typedef struct node{
  int lev;
  csr_t A;
  csr_t P;
  csr_t R;
  struct node *next;
} amg_t;//AMG is a linked list struct
amg_t *amg;

int maxlevel;
float maxratio;//maximum ratio measuring the sparsity of the matrix
int ndirect;//dimension for direct solve
int amg_niter;//number of amg cycles
int nu1;//number of pre-smoothing
int nu2;//number of post-smoothing


//AMG setup phase
void amg_init(csr_t A)
{
  int interp, lev = 0;
  float theta, ratio;
  amg_t *current, *newnode;
  
  if(!getparfloat("theta", &theta)) theta = 0.25;/* dimension in x */
  if(!getparfloat("maxratio", &maxratio)) maxratio = 0.8;/* dimension in x */
  if(!getparint("ndirect", &ndirect)) ndirect = 100;/* dimension for direct solve */
  if(!getparint("maxlevel", &maxlevel)) maxlevel = 10;/* dimension for direct solve */
  if(!getparint("amg_niter", &amg_niter)) amg_niter = 10;/* dimension for direct solve */
  if(!getparint("nu1", &nu1)) nu1 = 2;/* number of pre-smoothing */
  if(!getparint("nu2", &nu2)) nu2 = 2;/* number of post-smoothing */
  if(!getparint("interp", &interp)) interp = 1;/*1=standard interpolation, 2=direct interpolation */

  amg = (amg_t*)malloc(sizeof(amg_t));
  amg->lev = lev++;
  amg->A = A;
  amg->next = NULL;
  current = amg;
  ratio = (float)(current->A.nnz)/(current->A.nrow*current->A.ncol);
  printf("level=%d, nrow=%d, nnz=%d, nnz_ratio=%.2e\n",
	 current->lev, current->A.nrow, current->A.nnz, ratio);

  while(current->A.nrow > ndirect && lev<maxlevel && ratio<maxratio){
    newnode = (amg_t*)malloc(sizeof(amg_t));
    newnode->lev = lev++;
    current->P = build_P(current->A, theta, interp);
    current->R = csr_transpose(current->P);//R=Pt
    newnode->A = build_RAP(current->R, current->A, current->P);//A_{l+1}=R*A_l*P
    newnode->next = NULL;
    
    current->next = newnode;//append to the end of current linked list
    current = current->next;
    ratio = (float)(current->A.nnz)/(current->A.nrow*current->A.ncol);//nnz ratio/percentage
    printf("level=%d, nrow=%d, nnz=%d, nnz_ratio=%.2e\n",
	   current->lev, current->A.nrow, current->A.nnz, ratio);
  }
}

void amg_close()
{
  amg_t *current = amg;
  while(current->next!= NULL){
    csr_close(&current->P);
    csr_close(&current->R);
    csr_close(&current->next->A);
    current = current->next;
    printf("level=%d\n", current->lev);
  }
}

void amg_residual(csr_t A, double *x, double *b, double *r)
{
  int i, j, k;
  double Ax;

  for(i=0; i<A.nrow; i++){
    Ax = 0;
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];
      Ax += A.val[k]*x[j];
    }
    r[i] = b[i] - Ax;
  }
}

//y=A*x
void mat_vec_prod(csr_t A, double *x, double *y)
{
  int i, j, k;

  for(i=0; i<A.nrow; i++){
    y[i] = 0.;
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];
      y[i] += A.val[k]*x[j];
    }
  }
}

csr_t Amat;
void Aop_init(csr_t Amat_)
{
  Amat = Amat_;
}

void Aop_apply(int n, double *x, double *y)
{
  mat_vec_prod(Amat, x, y);
}

void amg_cycle(amg_t *amg, int n, double *x, double *b)
{
  int i, j, k, nc;
  double *r, *xc, *bc, **Ac;
  amg_t *amgnode = amg;

  if(amgnode->next==NULL){
    nc = amgnode->A.nrow;
    if(nc<ndirect){//direct solve at coarsest grid when dimension is small
      Ac = alloc2double(nc, nc);
      memset(&Ac[0][0], 0, nc*nc*sizeof(double));
      for(i=0; i<amgnode->A.nrow; i++){
	for(k=amgnode->A.row_ptr[i]; k<amgnode->A.row_ptr[i+1]; k++){
	  j = amgnode->A.col_ind[k];
	  Ac[i][j] = amgnode->A.val[k];
	}
	x[i] = b[i];
      }
      gauss_jordon(nc, &Ac[0][0], 1, x);//solution will be in x[]
      free2double(Ac);
    }else{//iterative solve at coarsest grid when dimension is large
      Aop_init(amgnode->A);
      solve_gmres(nc, x, b, Aop_apply, 1, 1e-6, 10, 0);
    }
  }else{
    nc = amgnode->next->A.nrow;
    r = alloc1double(n);
    xc = alloc1double(nc);
    bc = alloc1double(nc);

    for(i=0; i<nu1; i++) gauss_seidel(n, amgnode->A, x, b);//pre-smoothing
    amg_residual(amgnode->A, x, b, r);//r=b-Ax
    mat_vec_prod(amgnode->R, r, bc);//restriction, bc=R*r
    amg_cycle(amgnode->next, nc, xc, bc);

    mat_vec_prod(amgnode->P, xc, r);//prolongation, r=P*xc
    for(i=0; i<n; i++) x[i] += r[i];//correction
    for(i=0; i<nu2; i++) gauss_seidel(n, amgnode->A, x, b);//post-smoothing
    
    free1double(r);
    free1double(xc);
    free1double(bc);
  }
}


//x=invA*b
void amg_apply(int n, double *b, double *x)
{
  int iter;
  double rs0, rs, *r;

  r = alloc1double(n);
  
  amg_residual(amg->A, x, b, r);
  rs = dotprod(n, r, r);
  rs = sqrt(rs);
  rs0 = rs;
  for(iter=0; iter<amg_niter; iter++) {
    printf("AMG, iter=%d |r|=%e\n", iter, rs);
    amg_cycle(amg, n, x, b);

    amg_residual(amg->A, x, b, r);
    rs = dotprod(n, r, r);
    rs = sqrt(rs);
    if(rs<1e-6*rs0||rs<1e-15) break;
  }

  free1double(r);
}
