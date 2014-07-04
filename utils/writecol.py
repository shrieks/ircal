import numpy as np

output = np.column_stack((arrA.flatten(),arrB.flatten(),arrC.flatten()))
np.savetxt('emissivity-py.outt',output)
