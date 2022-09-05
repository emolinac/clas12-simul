import sys
import random

# Extension of the cryotarget
# z_limits[0] : 1cmlD2
# z_limits[1] : 2cmlD2
# z_limits[2] : 3cmlD2
z_limits = [[-4.2,4.5],[-9.2,9.5],[-14.2,14.5]]

variation = int(sys.argv[1])
rdm       = float(sys.argv[2])
#print ("Variations used ",variation)
#print ("Vertex = ",z_limits[variation-1][0] + rdm*(z_limits[variation-1][1] - z_limits[variation-1][0]))

print (z_limits[variation-1][0] + rdm*(z_limits[variation-1][1] - z_limits[variation-1][0]))
