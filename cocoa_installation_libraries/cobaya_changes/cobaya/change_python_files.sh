#!/bin/bash
# ---------------------------------------------------
# ---------------------------------------------------
# ---------------------------------------------------
PLANCK_CLIK='conventions.py'

sed --in-place --regexp-extended "s@from cobaya import __version__@from cobaya.__init__ import __version__@g" 'conventions.py'
sed --in-place --regexp-extended "s@from cobaya import __obsolete__@from cobaya.__init__ import __obsolete__@g" 'tools.py'


