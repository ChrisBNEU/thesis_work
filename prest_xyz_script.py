import sys
import os
prest_path = "/work/westgroup/ChrisB/_00_local_packages/prest/"
if prest_path not in sys.path:
    sys.path.append(prest_path)
    
meoh_unc_path = '/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-Py/rmgpy'
if meoh_unc_path in sys.path:
    sys.path.remove(meoh_unc_path)
    
import rmgpy
rmgpy.__path__
