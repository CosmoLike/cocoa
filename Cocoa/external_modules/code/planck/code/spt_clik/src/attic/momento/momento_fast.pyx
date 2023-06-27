cimport numpy as nm
import numpy as nm
nm.import_array()
cimport stdlib as stdlib
cimport stdio as stdio
import itertools
cimport cython
from cython.parallel import prange

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef nm.ndarray[nm.double_t, ndim=2] matpart(int p,int q,nm.ndarray[nm.double_t, ndim=2] mat,int lmin, int lmax):
  pstart=p**2-lmin**2
  qstart=q**2-lmin**2
  return mat[2*pstart:2*(pstart+2*p+1),2*qstart:2*(qstart+2*q+1)]

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double _pairtrace(nm.ndarray[nm.double_t,ndim=2] m1, nm.ndarray[nm.double_t,ndim=2] m2):
  #print m1.shape[0],m1.shape[1],m2.shape[0],m2.shape[1]
  return nm.sum(m1*(m2.transpose()))

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef double pairtrace(nm.ndarray[nm.double_t,ndim=2] m1, nm.ndarray[nm.double_t,ndim=2] m2):
  cdef int i,j,n,m
  cdef double res
  n = m1.shape[0]
  m = m1.shape[1]
  #print m1.shape[0],m1.shape[1],m2.shape[0],m2.shape[1],n,m
  res=0
  for i in range(n):
    for j in range(m):
      res += m1[i,j]*m2[j,i]
  #r = _pairtrace(m1,m2)
  #print res,r
  return res


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def symcov12_short(nm.ndarray[nm.double_t, ndim=2] mat11,nm.ndarray[nm.double_t, ndim=2] mat12,nm.ndarray[nm.double_t, ndim=2] mat22,sih):
  cdef int lmin,lmax
  cdef int l0,l1,i0,i1,i2,i3,i4
  cdef nm.ndarray[nm.double_t, ndim=2] m23,m41,m24,m31
  #print "fast"
  lmin,lmax = sih.lmin,sih.lmax
  #mat21=nm.transpose(mat12).copy()
  mat21 = mat12.T
  cov=nm.zeros((3*(lmax-lmin+1),3*(lmax-lmin+1)))
  #for l in itertools.product(xrange(lmin,lmax+1),repeat=2):
  for l0 in range(lmin,lmax+1):
    for l1 in range(lmin,lmax+1):
      zz1=nm.zeros((3,3))
      m23=matpart(l0,l1,mat21,lmin,lmax)
      m41=matpart(l1,l0,mat21,lmin,lmax)
      m24=matpart(l0,l1,mat22,lmin,lmax)
      m31=matpart(l1,l0,mat11,lmin,lmax)
      #for i in itertools.product([0,1],repeat=4):
      for i0 in range(2):
        for i1 in range(2):
          for i2 in range(2):
            for i3 in range(2):
              zz1[i0+i1,i2+i3] += (
                + pairtrace(m23[i1::2,i2::2],m41[i3::2,i0::2])
                + pairtrace(m24[i1::2,i3::2],m31[i2::2,i0::2])
                )
      cov[3*(l0-lmin):3*(l0+1-lmin),3*(l1-lmin):3*(l1+1-lmin)]=zz1
  return cov

import time

cdef extern void symkappa4_12_inner_loop(double *z, int p, int q, int r, int s, double* m, int lmin, int outlen) nogil
cdef extern void comp_indexes(int** pfis, int **pfls, int p, int q, int r, int s, int iperm, int ich) nogil
cdef extern int mind0(int L,int i0,int i1,int i2,int i3,int i4,int i5,int outlen) nogil
cdef extern void dec_kap3_inner(double* dk3, double *ci,double *k3, int n) nogil
cdef extern void dec_derivkap2_inner(double* dk2, double*  ci,double * k2, int n) nogil
cdef extern double pairdot(double* A, int As, int AT, double * B, int Bs, int BT, double* C, int Cs, int CT, int Ap, int Aop, int Aoq, int Bp, int Bop, int Boq, int Cp, int Cop, int Coq, int lmin) nogil
cdef extern void d_symcov12_inner(double *dcov, double* mat11,double* mat12,double* mat22,double* qty11,double* qty22,double *bl, int shape1, int lmin,int lmax) nogil
cdef extern void symkappa4_12_finalize(double* k4, int lmin,int lmax) nogil

#cdef int mind0(int L,int i0,int i1,int i2,int i3,int i4,int i5,int outlen) nogil:
#  cdef int o2 = outlen*outlen
#  return L*64*o2 + i0*32*o2 + i1*16*o2 + i2*8*o2 + i3*4*o2 + i4*2*o2 + i5*o2

import sys
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def symkappa4_12(nm.ndarray mat11,nm.ndarray mat12,nm.ndarray mat22,sih):
  cdef int lmin,lmax,inds1_0, inds1_1,inds1_2,inds1_3,inds1,inds2_0, inds2_1,inds2_2,inds2_3 ,inds2
  cdef int perms,perms_0,perms_1,perms_2
  cdef int p,q,r,s,pos,ich,a,b,c,d,outlen
  cdef double zval
  cdef nm.ndarray[nm.double_t,ndim=4] z=nm.zeros((3,3,3,3))
  #cdef nm.ndarray[nm.double_t,ndim=2] mp1full,mp2full,mp1,mp2
  cdef int *fis, *fls
  cdef int mp10 ,mp20 ,mp1_p,mp1_pstart,mp1_q,mp1_qstart,mp2_p,mp2_pstart,mp2_q,mp2_qstart
  cdef double *zarr
  cdef int i,i_0,i_1,i_2,i_3,i_4,i_5,r_o
  cdef nm.ndarray[nm.double_t,ndim=4] mats=nm.zeros((2,2,mat11.shape[0],mat11.shape[1]))
  
  lmin,lmax = sih.lmin,sih.lmax
  
  mats[0,0,:,:]=mat11
  mats[0,1,:,:]=mat12
  mats[1,0,:,:]=mat12.T
  mats[1,1,:,:]=mat22

  count=0
  
  cdef nm.ndarray[nm.double_t, ndim=4] k4=nm.zeros((3*(lmax-lmin+1),3*(lmax-lmin+1),3*(lmax-lmin+1),3*(lmax-lmin+1)))
  #cdef nm.ndarray[nm.double_t, ndim=4] k4b

  outlen=2*((lmax+1)**2-lmin**2)
  # l, inside freqs, inside t/e, outside freqs, outside t/e by ::2 or 1::2  
  #cdef nm.ndarray[nm.double_t,ndim=4] m=nm.zeros((lmax-lmin+1,2*2*2*2*2*2,outlen,outlen))
  #cdef nm.ndarray[nm.double_t, ndim=1] m=nm.zeros((lmax-lmin+1)*64*outlen*outlen)
  cdef double *m
  m = <double*> stdlib.malloc(sizeof(double)*(lmax-lmin+1)*64*outlen*outlen)
  cdef nm.ndarray[nm.double_t,ndim=2] mdot
  
  t1 = time.time()
  for p in range(lmin,lmax+1):
    for i in range(64):
      i_0,i_1,i_2,i_3,i_4,i_5 = (i&32)>>5,(i&16)>>4,(i&8)>>3,(i&4)>>2,(i&2)>>1,i&1
      m1=mats[i_4,i_0,:,i_2+2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2]
      m2=mats[i_1,i_5,i_3+2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:]
      mdot = nm.dot(m1,m2)
      r_o = mind0(p-lmin,i_0 , i_1 , i_2 , i_3 , i_4 , i_5,outlen)
      for a in range(outlen):
        for b in range(outlen):
          m[r_o+a*outlen+b] = mdot[a,b]
  
  #print time.time()-t1  
  
  t1 = time.time()
  #for p in prange(lmin,lmax+1,nogil=True):
  zarr = <double*>stdlib.malloc(sizeof(double)*3*3*3*3)
  for p in range(lmin,lmax+1):
    #print p
    
    for q in range(p,lmax+1):
      for r in range(q,lmax+1):
        for s in range(r,lmax+1):
          symkappa4_12_inner_loop(zarr,p,q,r,s,m,lmin,outlen)
          for a in range(3):
            for b in range(3):
              for c in range(3):
                for d in range(3):
                  k4[3*(p-lmin)+a,3*(q-lmin)+b,3*(r-lmin)+c,3*(s-lmin)+d] = zarr[a*27+b*9+c*3+d]
  stdlib.free(zarr)
  stdlib.free(m)
  #print time.time()-t1  

  #k4b = k4.copy()

  #k4.tofile("momento_fast_k4.dat") 
  ##t1 = time.time()
  ##for ells in (itertools.product(range(lmin,lmax+1),repeat=4)):
  ##  p,q,r,s=list(ells)
  ##  sargs=nm.argsort(ells)
  ##  psort,qsort,rsort,ssort=sorted(list(ells))
  ##  zz1 = k4[3*(psort-lmin):3*(psort+1-lmin),3*(qsort-lmin):3*(qsort+1-lmin),3*(rsort-lmin):3*(rsort+1-lmin),3*(ssort-lmin):3*(ssort+1-lmin)].copy()
  ##  zz2 = nm.zeros_like(zz1)
  ##  for inds in itertools.product(range(0,3),repeat=4):
  ##    tinds=list(inds)
  ##    #sinds=tinds
  ##    sinds2=list()
  ##    for i in sargs:
  ##        sinds2.append(tinds[i])
  ##    #sinds=array(tinds)[sargs]
  ##    zz2[tuple(tinds)] = zz1[tuple(sinds2)]
  ##  k4[3*(p-lmin):3*(p+1-lmin),3*(q-lmin):3*(q+1-lmin),3*(r-lmin):3*(r+1-lmin),3*(s-lmin):3*(s+1-lmin)] = zz2
  ##
  ##print time.time()-t1  
  ##k4.tofile("momento_fast_k4.dat") 
  
  t1 = time.time()
  symkappa4_12_finalize(<double*>k4.data,lmin,lmax)
  #print time.time()-t1  
  #k4b.tofile("momento_fast_k4b.dat") 
  
  return k4
  
def d_symcov12(nm.ndarray[nm.double_t,ndim=2] mat11,nm.ndarray[nm.double_t,ndim=2] mat12,nm.ndarray[nm.double_t,ndim=2] mat22,nm.ndarray[nm.double_t,ndim=2] qty11,nm.ndarray[nm.double_t,ndim=2] qty22,sih):
  cdef int lmin,lmax,i,i_0,i_1,i_2,i_3,i_4,i_5,l_0,l_1,l_2
  cdef nm.ndarray[nm.double_t,ndim=2] mat21
  lmin,lmax = sih.lmin,sih.lmax
  mat21=mat12.T
  cdef nm.ndarray[nm.double_t,ndim=3] dcov=nm.zeros((3*(lmax-lmin+1),3*(lmax-lmin+1),3*(lmax-lmin+1)))
  cdef nm.ndarray[nm.double_t,ndim=1] bl = sih.b_l
  #d_symcov12_inner(<double*>dcov.data, <double*> mat11.data,<double*> mat12.data,<double*> mat22.data,<double*> qty11.data,<double*> qty22.data,<double*> bl.data, mat11.shape[1], sih.lmin,sih.lmax)
  #return dcov
  #for l in itertools.product(xrange(lmin,lmax+1),repeat=3):
        #l_0,l_1,l_2 = l
  for l_0 in range(lmin,lmax+1):
    for l_1 in range(lmin,lmax+1):
      for l_2 in range(lmin,lmax+1):

        zz1=nm.zeros((3,3,3))
        #for i in itertools.product([0,1],repeat=6): 
        for i in range(64):
          i_0,i_1,i_2,i_3,i_4,i_5 = (i&32)>>5,(i&16)>>4,(i&8)>>3,(i&4)>>2,(i&2)>>1,i&1
          #i_0,i_1,i_2,i_3,i_4,i_5 = i   
          zz1[i_0+i_1,i_2+i_3,i_4+i_5] += (
            + pairdot(<double*>qty22.data,qty22.shape[1],0,<double*>qty11.data,qty11.shape[1],0,<double*>mat21.data,mat21.shape[1],1,l_0,i_1,i_4,l_2,i_5,i_2,l_1,i_3,i_0,lmin)
            + pairdot(<double*>mat21.data,mat21.shape[1],1,<double*>qty22.data,qty22.shape[1],0,<double*>qty11.data,qty11.shape[1],0,l_0,i_1,i_2,l_1,i_3,i_5,l_2,i_4,i_0,lmin)
            + pairdot(<double*>qty22.data,qty22.shape[1],0,<double*>qty22.data,qty22.shape[1],0,<double*>mat11.data,mat11.shape[1],0,l_0,i_1,i_4,l_2,i_5,i_3,l_1,i_2,i_0,lmin)
            + pairdot(<double*>mat22.data,mat22.shape[1],0,<double*>qty11.data,qty11.shape[1],0,<double*>qty11.data,qty11.shape[1],0,l_0,i_1,i_3,l_1,i_2,i_5,l_2,i_4,i_0,lmin)
          )*sih.b_l[l_2]**2
          #zz1[i[0]+i[1],i[2]+i[3],i[4]+i[5]] += (
          #  + _pairdot(qty22,qty11,mat21,l[0],i[1],i[4],l[2],i[5],i[2],l[1],i[3],i[0],lmin,0,0,1)
          #  + _pairdot(mat21,qty22,qty11,l[0],i[1],i[2],l[1],i[3],i[5],l[2],i[4],i[0],lmin,1,0,0)
          #  + _pairdot(qty22,qty22,mat11,l[0],i[1],i[4],l[2],i[5],i[3],l[1],i[2],i[0],lmin,0,0,0)
          #  + _pairdot(mat22,qty11,qty11,l[0],i[1],i[3],l[1],i[2],i[5],l[2],i[4],i[0],lmin,0,0,0)
          #)*sih.b_l[l[2]]**2

        dcov[3*(l_0-lmin):3*(l_0+1-lmin),3*(l_1-lmin):3*(l_1+1-lmin),3*(l_2-lmin):3*(l_2+1-lmin)]=zz1
  #dcov.tofile("dcov_fast.dat")
  return dcov


cdef double _pairdot(nm.ndarray[nm.double_t,ndim=2] A, nm.ndarray[nm.double_t,ndim=2] B, nm.ndarray[nm.double_t,ndim=2] C,int Ap, int Aop, int Aoq, int Bp, int Bop, int Boq, int Cp, int Cop, int Coq, int lmin,int AT,int BT,int CT):
  cdef int i,k,j,lmin2
  cdef int Apstart, Aqstart, Bpstart, Bqstart, Cpstart, Cqstart
  cdef int Aq,Cq,Bq
  cdef double res, AB_ij,r2

  Cq = Ap
  Aq = Bp 
  Bq = Cp 

  lmin2 = lmin*lmin
  Apstart = (Ap*Ap - lmin2)*2
  Bpstart = (Bp*Bp - lmin2)*2
  Cpstart = (Cp*Cp - lmin2)*2
  Aqstart = (Aq*Aq - lmin2)*2
  Bqstart = (Bq*Bq - lmin2)*2
  Cqstart = (Cq*Cq - lmin2)*2

  r2 = pairdot(<double*>A.data,A.shape[1],AT,<double*>B.data,B.shape[1],BT,<double*>C.data,C.shape[1],CT,Ap,Aop,  Aoq,  Bp, Bop, Boq, Cp, Cop, Coq,lmin)
  return r2
  

def dec_derivkap2(nm.ndarray[nm.double_t, ndim=2]  ci,nm.ndarray[nm.double_t, ndim=3] k2):
  cdef int n,p,j,x,s,i  
  cdef nm.ndarray[nm.double_t,ndim=3] dk2

  n=k2.shape[0]
  dk2=nm.zeros_like(k2)
  dec_derivkap2_inner(<double*>dk2.data,<double*>ci.data,<double*>k2.data,n);
  return dk2

  dk2=nm.zeros_like(k2)
  for p in range(n):
      for j in range(n):
          for x in range(n):
              for s in range(n):
                  dk2[p,j,x]+=ci[j,s]*k2[p,s,x]

  k2t=dk2.copy()
  dk2=nm.zeros_like(k2)
  for i in range(n):
      for j in range(n):
          for x in range(n):
              for s in range(n):
                  dk2[i,j,x]+=ci[i,s]*k2t[s,j,x]

  return dk2


def dec_kap3(nm.ndarray[nm.double_t, ndim=2] ci,nm.ndarray[nm.double_t, ndim=3] k3):
  cdef int n,p,q,k,s,j,i
  cdef nm.ndarray[nm.double_t, ndim=3] dk3
  n=k3.shape[0]
  #print k3.shape[0],k3.shape[1],k3.shape[2]

  dk3=nm.zeros_like(k3)
  dec_kap3_inner(<double*> dk3.data, <double*> ci.data ,<double*> k3.data,n)
  return dk3

    
def symkappa3_12(mat11,mat12,mat22,sih):
    lmin,lmax = sih.lmin,sih.lmax
    mat21=nm.transpose(mat12).copy()
    k3=nm.zeros((3*(lmax-lmin+1),3*(lmax-lmin+1),3*(lmax-lmin+1)))
    outlen=2*((lmax+1)**2-lmin**2)
    m_11_21=nm.zeros((2,2,outlen,outlen))
    m_11_22=nm.zeros((2,2,outlen,outlen))
    m_21_21=nm.zeros((2,2,outlen,outlen))
    m_21_22=nm.zeros((2,2,outlen,outlen))

    t1 = time.time()
    for p in range(lmin,lmax+1):
        m_11_21[0,0,:,:]=nm.dot(mat11[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_11_21[0,1,:,:]=nm.dot(mat11[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])
        m_11_21[1,0,:,:]=nm.dot(mat11[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_11_21[1,1,:,:]=nm.dot(mat11[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])

        m_11_22[0,0,:,:]=nm.dot(mat11[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_11_22[0,1,:,:]=nm.dot(mat11[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])
        m_11_22[1,0,:,:]=nm.dot(mat11[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_11_22[1,1,:,:]=nm.dot(mat11[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])

        m_21_21[0,0,:,:]=nm.dot(mat21[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_21_21[0,1,:,:]=nm.dot(mat21[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])
        m_21_21[1,0,:,:]=nm.dot(mat21[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_21_21[1,1,:,:]=nm.dot(mat21[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])

        m_21_22[0,0,:,:]=nm.dot(mat21[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_21_22[0,1,:,:]=nm.dot(mat21[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])
        m_21_22[1,0,:,:]=nm.dot(mat21[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_21_22[1,1,:,:]=nm.dot(mat21[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])

        for q in range(p,lmax+1):
            for r in range(q,lmax+1):
                z=nm.zeros((3,3,3))
                for inds in itertools.product([0,1],repeat=6):
                    biginds=tuple(nm.array(inds[::2])+nm.array(inds[1::2]))
                    z[biginds] += (
                        + pairtrace(matpart(r,q,m_21_21[inds[0],inds[1],:,:],lmin,lmax)[inds[5]::2,inds[2]::2],matpart(q,r,mat21,lmin,lmax)[inds[3]::2,inds[4]::2])
                        + pairtrace(matpart(r,q,m_11_21[inds[0],inds[1],:,:],lmin,lmax)[inds[4]::2,inds[2]::2],matpart(q,r,mat22,lmin,lmax)[inds[3]::2,inds[5]::2])
                        + pairtrace(matpart(r,q,m_21_22[inds[0],inds[1],:,:],lmin,lmax)[inds[5]::2,inds[3]::2],matpart(q,r,mat11,lmin,lmax)[inds[2]::2,inds[4]::2])
                        + pairtrace(matpart(r,q,m_11_22[inds[0],inds[1],:,:],lmin,lmax)[inds[4]::2,inds[3]::2],matpart(q,r,mat12,lmin,lmax)[inds[2]::2,inds[5]::2])
                        + pairtrace(matpart(q,r,m_21_21[inds[0],inds[1],:,:],lmin,lmax)[inds[3]::2,inds[4]::2],matpart(r,q,mat21,lmin,lmax)[inds[5]::2,inds[2]::2])
                        + pairtrace(matpart(q,r,m_11_21[inds[0],inds[1],:,:],lmin,lmax)[inds[2]::2,inds[4]::2],matpart(r,q,mat22,lmin,lmax)[inds[5]::2,inds[3]::2])
                        + pairtrace(matpart(q,r,m_21_22[inds[0],inds[1],:,:],lmin,lmax)[inds[3]::2,inds[5]::2],matpart(r,q,mat11,lmin,lmax)[inds[4]::2,inds[2]::2])
                        + pairtrace(matpart(q,r,m_11_22[inds[0],inds[1],:,:],lmin,lmax)[inds[2]::2,inds[5]::2],matpart(r,q,mat12,lmin,lmax)[inds[4]::2,inds[3]::2])
                        )
                k3[3*(p-lmin):3*(p+1-lmin),3*(q-lmin):3*(q+1-lmin),3*(r-lmin):3*(r+1-lmin)]=z     
            #print p,q,r
    t2 = time.time()
    #print "s",t2-t1
    for ells in (itertools.product(range(lmin,lmax+1),repeat=3)):
        p,q,r=list(ells)
        sargs=nm.argsort(ells)
        psort,qsort,rsort=sorted(list(ells))
        zz1 = k3[3*(psort-lmin):3*(psort+1-lmin),3*(qsort-lmin):3*(qsort+1-lmin),3*(rsort-lmin):3*(rsort+1-lmin)].copy()
        zz2 = nm.zeros_like(zz1)
        for inds in itertools.product(range(0,3),repeat=3):
            tinds=list(inds)
            sinds=nm.array(tinds)[sargs]
            zz2[tuple(tinds)] = zz1[tuple(sinds)]
        k3[3*(p-lmin):3*(p+1-lmin),3*(q-lmin):3*(q+1-lmin),3*(r-lmin):3*(r+1-lmin)] = zz2
    #print "s",time.time()-t2
    #k3.tofile("sk3_slow.dat")
    return k3
