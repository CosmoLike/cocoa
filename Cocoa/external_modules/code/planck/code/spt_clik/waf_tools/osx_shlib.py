from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
def configure(ctx):
  import sys
  ctx.env.shsuffix = "so"
  if sys.platform.lower()=="darwin":
    ctx.env.shsuffix = "dylib"