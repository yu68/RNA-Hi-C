import sys
from RNAstructure import dot2block

stems = dot2block(sys.argv[1],"test")

print stems
