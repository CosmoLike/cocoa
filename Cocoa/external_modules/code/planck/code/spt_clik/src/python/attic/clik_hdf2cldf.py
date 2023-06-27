#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import clik.cldf as cldf

cldf.hdf2cldf(sys.argv[1],sys.argv[2])

