import numpy as np
import h5py as h5
import sys
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))
print(sys.argv[1])
nu = sys.argv[1];
filename = "data.h5.old"
new_file = h5.File("data_run"+str(nu)+".h5","w")
with h5.File(filename, "r") as f:
    print("Opening:",filename)
    dataset = "solution/run_"+str(nu)
    print("copying:",dataset)
    data = f[dataset]
    dset_b = f[dataset]["boundaries"]
    dset_i = f[dataset]["interior"]
    gr = new_file.create_group("solution")
    grgr = gr.create_group("run_"+str(nu))
    grgr.create_dataset("boundaries",data=dset_b[:])
    grgr.create_dataset("interior",data=dset_i[:])


