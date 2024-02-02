import sys
import chimera
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
import AddH
import AddCharge
import WriteMol2

lig_file = sys.argv[3]
out_lig_file = sys.argv[4]
#open ligand file for reading
lig = chimera.openModels.open(lig_file)

#chimera.update.checkForChanges()

rc("del H")
rc("addh")
rc("addcharge all method gas")
rc("write format mol2 #0 %s"%(out_lig_file))

rc('close all')

print ('completed')

exit()
