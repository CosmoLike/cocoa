// slik_support.c
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
void comp_indexes(int** pfis, int **pfls, int p, int q, int r, int s, int iperms, int ich) {
  int perms_0, perms_1, perms_2;
  int fulperms[4];
  int terms[12];
  int ichm[18] = {0,1,2,0,2,1,1,0,2,1,2,0,2,0,1,2,1,0};
  int odd[3] = {3,5,7};
  int even[3] = {2,4,6};
  int qrs[3];
  int freqinds[8];
  int lterms[8];
  int i;
  int *fis,*fls;

  fis = malloc(sizeof(int)*8);
  fls = malloc(sizeof(int)*8);

  *pfls = fls;
  *pfis = fis;
  
  qrs[0] = q;
  qrs[1] = r;
  qrs[2] = s;

  perms_0 =  (-(iperms&4)>>1)+1;
  perms_1 =  -(iperms&2)+1;
  perms_2 =  (-(iperms&1)<<1)+1;

  fulperms[0] = 1;
  fulperms[1] = perms_0;
  fulperms[2] = perms_1;
  fulperms[3] = perms_2;

  terms[0*3+0] = p;
  terms[0*3+1] = 0;
  terms[0*3+2] = 1;

  terms[1*3+0] = qrs[ichm[ich*3+0]];
  terms[1*3+1] = even[ichm[ich*3+0]];
  terms[1*3+2] = odd[ichm[ich*3+0]];

  terms[2*3+0] = qrs[ichm[ich*3+1]];
  terms[2*3+1] = even[ichm[ich*3+1]];
  terms[2*3+2] = odd[ichm[ich*3+1]];

  terms[3*3+0] = qrs[ichm[ich*3+2]];
  terms[3*3+1] = even[ichm[ich*3+2]];
  terms[3*3+2] = odd[ichm[ich*3+2]];

  for(i=0;i<4;i++) {
    lterms[i*2+0] = terms[i*3+0];
    lterms[i*2+1] = terms[i*3+0];
  
    if (fulperms[i]>0) {
      freqinds[i*2+0] = terms[i*3+1];
      freqinds[i*2+1] = terms[i*3+2];
    } else {
      freqinds[i*2+0] = terms[i*3+2];
      freqinds[i*2+1] = terms[i*3+1];
    } 
  }
  
  fis[0] = freqinds[7];
  fls[0] = lterms[7];
  for(i=1;i<8;i++) {
    fis[i] = freqinds[i-1];
    fls[i] = lterms[i-1];
  }

}

int mind0(int L,int i0,int i1,int i2,int i3,int i4,int i5,int outlen) {
  int o2;
  o2 = outlen*outlen;
  return L*64*o2 + i0*32*o2 + i1*16*o2 + i2*8*o2 + i3*4*o2 + i4*2*o2 + i5*o2;
}

void symkappa4_12_fullinner_loop(double *k4, double* m, int lmin, int lmax, int outlen) {
  int a,b,c,d;
  int p,q,r,s;

  for(p=lmin;p<lmax+1;p++) {
    for(q=p;q<lmax+1;q++) {
      for(r=q;r<lmax+1;r++) {
        for(s=r;r<lmax+1;r++) {
          fprintf(stderr, "%s\n", "rien");

        }
      }
    }
  }
}
void symkappa4_12_inner_loop(double *z, int p, int q, int r, int s, double* m, int lmin, int outlen) {
  int ich,ip,a,b;
  int *fis, *fls;
  int inds1, inds1_0, inds1_1, inds1_2, inds1_3;
  int inds2, inds2_0, inds2_1, inds2_2, inds2_3;
  int mp10 ,mp20 ,mp1_p,mp1_pstart,mp1_q,mp1_qstart,mp2_p,mp2_pstart,mp2_q,mp2_qstart;
  double zval;
  
  //fprintf(stderr, "%d\n",__LINE__ );
  memset(z,0,sizeof(double)*3*3*3*3);
  //fprintf(stderr, "%d\n",__LINE__ );
  for(ich=0;ich<6;ich++) {
    for(ip=0;ip<8;ip++) {
      comp_indexes(&fis,&fls,p,q,r,s,ip,ich);
      for(inds1=0;inds1<16;inds1++) {
        inds1_0 = (inds1&8)>>3;
        inds1_1 = (inds1&4)>>2;
        inds1_2 = (inds1&2)>>1;
        inds1_3 = inds1&1;
        mp10 = mind0(fls[1]-lmin,fis[1]%2,fis[2]%2,inds1_0,inds1_1,fis[0]%2,fis[3]%2,outlen);
        mp20 = mind0(fls[5]-lmin,fis[5]%2,fis[6]%2,inds1_2,inds1_3,fis[4]%2,fis[7]%2,outlen);
        mp1_p = fls[0];
        mp1_pstart=(mp1_p*mp1_p-lmin*lmin)*2;
        mp1_q = fls[3];
        mp1_qstart=(mp1_q*mp1_q-lmin*lmin)*2;
        mp2_p = fls[4];
        mp2_pstart=(mp2_p*mp2_p-lmin*lmin)*2;
        mp2_q = fls[7];
        mp2_qstart=(mp2_q*mp2_q-lmin*lmin)*2;
        //fprintf(stderr,"%d %d %d %d %d %d\n",mp10,mp20,mp1_p,mp1_q,mp2_p,mp2_q);
        for(inds2=0;inds2<16;inds2++) {
          inds2_0 = (inds2&8)>>3;
          inds2_1 = (inds2&4)>>2;
          inds2_2 = (inds2&2)>>1;
          inds2_3 = inds2&1;
          zval=0;
          for(a=0;a<2*mp1_p+1;a++) {
            for(b=0;b<2*mp1_q+1;b++) {
              zval += m[mp10 + (mp1_pstart+inds2_0 + 2*a) * outlen + mp1_qstart+inds2_1 + 2*b] * m[mp20 + (mp2_pstart+inds2_2 + 2*b) * outlen + mp2_qstart+inds2_3 + 2*a];
              //fprintf(stderr,"%g %d %d %d %d\n",zval,a,b,mp10 + (mp1_pstart+inds2_0 + 2*a) * outlen + mp1_qstart+inds2_1 + 2*b,mp20 + (mp2_pstart+inds2_2 + 2*b) * outlen + mp2_qstart+inds2_3 + 2*a);
            }
          }
          //return;
          //fprintf(stderr,"%g\n",zval);
          z[ (inds1_0+inds1_1) *27+ (inds2_0+inds2_1) *9+ (inds1_2+inds1_3) *3+ (inds2_2+inds2_3)] += zval;
        }
      }
      free(fis);
      free(fls);            
    }
  //fprintf(stderr, "%d\n",__LINE__ );
  }

}


void dec_kap3_inner(double* dk3, double *ci,double *k3, int n) {
  int p,q,k,s,j,i,n2;
  double *tmp1, *tmp2;

  n2=n*n;

  tmp1 = malloc(sizeof(double)*n*n*n);
  tmp2 = malloc(sizeof(double)*n*n*n);
  memset(tmp1,0,sizeof(double)*n*n*n);
  memset(tmp2,0,sizeof(double)*n*n*n);
  for(p=0;p<n;p++) {
    for(q=0;q<n;q++) {
      for(k=0;k<n;k++) {
        for(s=0;s<n;s++) {
          tmp1[p*n2+q*n+k]+=ci[k*n+s]*k3[p*n2+q*n+s];
        }
      }
    }
  }
  for(p=0;p<n;p++) {
    for(j=0;j<n;j++) {
      for(k=0;k<n;k++) {
        for(s=0;s<n;s++) {
          tmp2[p*n2+j*n+k]+=ci[j*n+s]*tmp1[p*n2+s*n+k];
        }
      }
    }
  }
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      for(k=0;k<n;k++) {
        for(s=0;s<n;s++) {
          dk3[i*n2+j*n+k]+=ci[i*n+s]*tmp2[s*n2+j*n+k];
        }
      }
    }
  }
  free(tmp1);
  free(tmp2);
  
  return;
}

void dec_derivkap2_inner(double* dk2, double*  ci,double * k2, int n) {
  int n2,p,j,x,s,i; 
  double *tmp;

  n2 = n*n;
  tmp = malloc(sizeof(double)*n*n*n);
  memset(tmp,0,sizeof(double)*n*n*n);
  
  for(p=0;p<n;p++) {
    for(j=0;j<n;j++) {
      for(x=0;x<n;x++) {
        for(s=0;s<n;s++) {
          tmp[p*n2+j*n+x] += ci[j*n+s]*k2[p*n2+s*n+x];
        }
      }
    }
  }

  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      for(x=0;x<n;x++) {
        for(s=0;s<n;s++) {
          dk2[i*n2+j*n+x] += ci[i*n+s]*tmp[s*n2+j*n+x];
        }
      }
    }
  }
  free(tmp);
  return;
}


double pairdot(double* A, int As, int AT, double * B, int Bs, int BT, double* C, int Cs, int CT, int Ap, int Aop, int Aoq, int Bp, int Bop, int Boq, int Cp, int Cop, int Coq, int lmin) {
  int i,k,j,lmin2;
  int Apstart, Aqstart, Bpstart, Bqstart, Cpstart, Cqstart;
  int Aq,Cq,Bq, Ax, Ay, Bx, By, Cx, Cy;
  double res, AB_ij;

  Cq = Ap;
  Aq = Bp; 
  Bq = Cp; 

  lmin2 = lmin*lmin;
  Apstart = (Ap*Ap - lmin2)*2;
  Bpstart = (Bp*Bp - lmin2)*2;
  Cpstart = (Cp*Cp - lmin2)*2;
  Aqstart = (Aq*Aq - lmin2)*2;
  Bqstart = (Bq*Bq - lmin2)*2;
  Cqstart = (Cq*Cq - lmin2)*2;
  
  Ay = As;
  Ax = 1;
  if (AT==1) {
    Ay = 1;
    Ax = As;
  }
  By = Bs;
  Bx = 1;
  if (BT==1) {
    By = 1;
    Bx = Bs;
  }
  Cy = Cs;
  Cx = 1;
  if (CT==1) {
    Cy = 1;
    Cx = Cs;
  }

  res = 0;
  for(i=0;i<2*Ap+1;i++) {
    for(j=0;j<2*Cp+1;j++) {
      AB_ij = 0;
      for(k=0;k<2*Bp+1;k++){
        AB_ij += A[(Apstart+Aop+2*i) *Ay+ (Aqstart+Aoq+2*k) *Ax] * B[(Bpstart+Bop+2*k) *By+ (Bqstart+Boq+2*j)  *Bx];
      }
      res += AB_ij * C[(Cpstart+Cop+2*j) *Cy+ (Cqstart+Coq+2*i) *Cx];
    }
  }
  return res;
}

 void d_symcov12_inner(double *dcov, double* mat11,double* mat12,double* mat22,double* qty11,double* qty22,double *bl, int shape1, int lmin,int lmax) {
  int ni,i_0,i_1,i_2,i_3,i_4,i_5,l_0,l_1,l_2,ll,ll2,i,j,k;
  double *zz1;

  // don't use yet, not yet fully debugged.
  zz1 = malloc(sizeof(double)*3*3*3);
  ll = lmax+1-lmin;
  ll=ll*3;
  ll2 = ll*ll;

  for(l_0=lmin;l_0<lmax+1;l_0++) {
    for(l_1=lmin;l_1<lmax+1;l_1++) {
      for(l_2=lmin;l_2<lmax+1;l_2++) {

        memset(zz1,0,sizeof(double)*3*3*3);

        for(ni=0;i<64;i++) {
          i_0 = (ni&32)>>5;
          i_1 = (ni&16)>>4;
          i_2 = (ni&8)>>3;
          i_3 = (ni&4)>>2;
          i_4 = (ni&2)>>1;
          i_5 = ni&1;
          
             
          zz1[(i_0+i_1) * 9+ (i_2+i_3) *3 + i_4+i_5] += (
            + pairdot(qty22,shape1,0,qty11,shape1,0,mat12,shape1,1,l_0,i_1,i_4,l_2,i_5,i_2,l_1,i_3,i_0,lmin)
            + pairdot(mat12,shape1,1,qty22,shape1,0,qty11,shape1,0,l_0,i_1,i_2,l_1,i_3,i_5,l_2,i_4,i_0,lmin)
            + pairdot(qty22,shape1,0,qty22,shape1,0,mat11,shape1,0,l_0,i_1,i_4,l_2,i_5,i_3,l_1,i_2,i_0,lmin)
            + pairdot(mat22,shape1,0,qty11,shape1,0,qty11,shape1,0,l_0,i_1,i_3,l_1,i_2,i_5,l_2,i_4,i_0,lmin)
          )*bl[l_2]*bl[l_2];
        }
        for(i=0;i<3;i++) {
          for(j=0;j<3;j++) {
            for(k=0;k<3;k++) {
              dcov[(3*(l_0-lmin)+i) *ll2+ (3*(l_1-lmin)+j) * ll+ 3*(l_2-lmin)+k]=zz1[i*9+j*3+k];      
            }
          }
        }  
      }
    }
  }
  free(zz1);
}

void bsort4(int* arg4,int *in4 ) {
  int i,stop,tmp;
  arg4[0]=0;
  arg4[1]=1;
  arg4[2]=2;
  arg4[3]=3;

  stop=0;
  while(stop==0) {
    stop=1;
    for(i=0;i<3;i++) {
      if (in4[arg4[i]]>in4[arg4[i+1]]) {
        stop=0;
        tmp = arg4[i];
        arg4[i] = arg4[i+1];
        arg4[i+1] = tmp;
      }
    }
  }
}

void symkappa4_12_finalize(double* k4, int lmin,int lmax) {
  int p,q,r,s,pp,qp,rp,sp;
  int i,j,k,l,ip,jp,kp,lp;
  int ll,ll2,ll3;

  int in4[4],arg4[4],ind[4],rnd[4];

  ll = (lmax+1-lmin)*3;
  ll2 = ll*ll;
  ll3 = ll*ll2;

  for (p=lmin;p<lmax+1;p++) {
    in4[0] = p;
    for (q=lmin;q<lmax+1;q++) {
      in4[1] = q;
      for (r=lmin;r<lmax+1;r++) {
        in4[2] = r;
        for (s=lmin;s<lmax+1;s++) {
          in4[3] = s;
          bsort4(arg4,in4);
          pp = in4[arg4[0]];
          qp = in4[arg4[1]];
          rp = in4[arg4[2]];
          sp = in4[arg4[3]];

          for(i=0;i<3;i++) {
            ind[0]=i;
            for(j=0;j<3;j++) {
              ind[1]=j;
              for(k=0;k<3;k++) {
                ind[2]=k;
                for(l=0;l<3;l++) {
                  ind[3]=l;
                  ip=ind[arg4[0]];
                  jp=ind[arg4[1]];
                  kp=ind[arg4[2]];
                  lp=ind[arg4[3]];
                  k4[(3*(p-lmin)+i) *ll3+ (3*(q-lmin)+j) *ll2+ (3*(r-lmin)+k) *ll+ (3*(s-lmin)+l)] = k4[(3*(pp-lmin)+ip) *ll3+ (3*(qp-lmin)+jp) *ll2+ (3*(rp-lmin)+kp) *ll+ (3*(sp-lmin)+lp)];
                } 
              }
            }  
          }
        }
      }
    }  
  }
}