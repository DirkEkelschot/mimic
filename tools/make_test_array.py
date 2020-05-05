import numpy as np
import h5py as h5
import sys

arr = np.zeros((8,8));
#Element 0
arr[0,0] = 0;
arr[0,1] = 1;
arr[0,2] = 6;
arr[0,3] = 5;
arr[0,4] = 14+0;
arr[0,5] = 14+1;
arr[0,6] = 14+6;
arr[0,7] = 14+5;
#Element 1
arr[1,0] = 1;
arr[1,1] = 2;
arr[1,2] = 7;
arr[1,3] = 6;
arr[1,4] = 14+1;
arr[1,5] = 14+2;
arr[1,6] = 14+7;
arr[1,7] = 14+6;
#Element 2
arr[2,0] = 2;
arr[2,1] = 3;
arr[2,2] = 8;
arr[2,3] = 7;
arr[2,4] = 14+2;
arr[2,5] = 14+3;
arr[2,6] = 14+8;
arr[2,7] = 14+7;
#Element 3
arr[3,0] = 3;
arr[3,1] = 4;
arr[3,2] = 9;
arr[3,3] = 8;
arr[3,4] = 14+3;
arr[3,5] = 14+4;
arr[3,6] = 14+9;
arr[3,7] = 14+8;
#Element 4
arr[4,0] = 5;
arr[4,1] = 6;
arr[4,2] = 11;
arr[4,3] = 10;
arr[4,4] = 14+5;
arr[4,5] = 14+6;
arr[4,6] = 14+11;
arr[4,7] = 14+10;
#Element 5
arr[5,0] = 6;
arr[5,1] = 7;
arr[5,2] = 12;
arr[5,3] = 11;
arr[5,4] = 14+6;
arr[5,5] = 14+7;
arr[5,6] = 14+12;
arr[5,7] = 14+11;
#Element 6
arr[6,0] = 7;
arr[6,1] = 8;
arr[6,2] = 13;
arr[6,3] = 12;
arr[6,4] = 14+7;
arr[6,5] = 14+8;
arr[6,6] = 14+13;
arr[6,7] = 14+12;
#Element 7
arr[7,0] = 8;
arr[7,1] = 9;
arr[7,2] = 14;
arr[7,3] = 13;
arr[7,4] = 14+8;
arr[7,5] = 14+9;
arr[7,6] = 14+14;
arr[7,7] = 14+13;

filename = "test_grid.h5"
new_file = h5.File(filename,"w")
new_file.create_dataset("ien",data=arr)
