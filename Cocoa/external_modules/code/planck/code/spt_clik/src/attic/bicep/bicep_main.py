#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import bicep_util as bu
import numpy as nm
import sys
import os

dirn = sys.stdin.readline()
os.chdir(dirn.strip())

C_l, C_l_hat, N_l, C_fl, M_inv, bpwf_l, bpwf_Cs_l =  bu.init("bicep2", "EB")

while True:
  vo = sys.stdin.readline()
  if vo.strip()=="stop":
    print>>sys.stderr, "bye"
    break
  cls = nm.array([float(vo)]+[float(sys.stdin.readline()) for i in range(1201)])
  cl = nm.zeros((601,9))
  llp1 = nm.arange(601)*nm.arange(1,602)/2./nm.pi
  cl[:,2] = cls[:601]*llp1
  cl[:,3] = cls[601:]*llp1

  expv = bu.calc_expvals(nm.arange(601),cl,bpwf_l,bpwf_Cs_l)
  
  C_l[:,0,0] = expv[:,2]
  C_l[:,0,1] = expv[:,5]
  C_l[:,1,0] = expv[:,5]
  C_l[:,1,1] = expv[:,3]

  C_l = C_l + N_l
  
  logL = bu.evaluateLikelihood(C_l,C_l_hat,C_fl,M_inv)
  print>>sys.stdout, "READY"
  print logL
  sys.stdout.flush()
