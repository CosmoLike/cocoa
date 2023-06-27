#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path
import numpy as nm
import sys
import os


import cosmoslik as K
import os.path as osp

#the order of passed in Cls / returns lmaxs. Karim, set accordingly. 
order=['TT','EE','BB','TE','TB','EB']

class mspec_clik(object):

  def __init__(self,mspec_folder):
    """
    Initialize the mspec likelihood. 
    mspec_folder points to the top folder for the given run.  
    """
    self.script = K.load_script(osp.join(mspec_folder,'script.py'),cls_set_externally=True)

  def get_nuisance_params(self):
    """
    Returns the names of needed nuisance parameter as a list of strings.
    """
    return list(self.script.get_sampled().keys())

  def get_lmaxs(self):
    """
    Returns a list of lmaxs in the order given by the order parameter.
    """
    return [max([-1]+[lmax for ((x1,_),(x2,_)),(_,lmax) in list(self.script.params.mspec.use.items()) if x1+x2==x]) for x in order]

  def get_lnl(self,dls,nuisance_params):
    """
    Return negative log likelihood. 
    dls is an array of CMB l(l+1)/2/pi Cl's in units of muK^2, in the order given by order. 
    0th index of the array is l=0. lmaxs need to be at least those given get get_lmaxs()
    nuisance_params is a list of the values of the nuisance parameters, in the order given by get_nuisance_params.
    """
    params = dict(list(zip(self.script.get_sampled(),nuisance_params)))
    params['cmb_result'] = {'cl_%s'%x:dl for x,dl in zip(order,dls)}
    return self.script.evaluate(**params)[0]

def main(argv):
  
  #print >>sys.stderr,"ivi"
  dirn = sys.stdin.readline().strip()
  #print >>sys.stderr, dirn
  #print >>sys.stderr,sys.path
  msp = mspec_clik(dirn)

  nuis = msp.get_nuisance_params()
  lmax = msp.get_lmaxs()
  #print >>sys.stderr,lmax
  #for lm in lmax:
  #  print lm
  #  sys.stdout.flush()
  #  resp = sys.stdin.readline()
  #  if resp.strip()!='ok':
  #    print >>sys.stderr,"bad !"
  #    sys.exit(-1)

  #print len(nuis)
  #sys.stdout.flush()
  #if resp.strip()!='ok':
  #  print >>sys.stderr,"bad !"
  #  sys.exit(-1)
  #for n in nuis:
  #  print n
  #  sys.stdout.flush()
  #  if resp.strip()!='ok':
  #    print >>sys.stderr,"bad !"
  #    sys.exit(-1)
  #
  lm = max(lmax)
  lmax = [lm if ll!=-1 else -1 for ll in lmax]
  llp1 = [nm.arange(ll+1)*(nm.arange(1,ll+2))/2./nm.pi for ll in lmax]
  #print >>sys.stderr,lmax
  for i in range(6):
    if lmax[i]!=-1:
      llp1[i][0]=1
    
  ltot = nm.sum(lmax)+6
  ntot = len(nuis)
  #print >>sys.stderr,ntot,ltot,ntot+ltot
  print("READY")
  sys.stdout.flush()
  
  while(True):
    vo = sys.stdin.readline()
    if vo.strip()=="stop":
      print("bye", file=sys.stderr)
      break
    cls = nm.array([float(vo)]+[float(sys.stdin.readline()) for i in range(ltot+ntot-1)])
    dls = []
    cnt = 0

    for i in range(6):
      dls += [cls[cnt:cnt+lmax[i]+1]*llp1[i]]
      cnt += lmax[i]+1

    nu = cls[cnt:]

    lkl = msp.get_lnl(dls,nu)
    print("READY")
    print(-lkl)
    sys.stdout.flush()
    #if resp.strip()!='ok':
    #  print >>sys.stderr,"bad !"
    #  sys.exit(-1)





import sys
if __name__=="__main__":
 main(sys.argv)