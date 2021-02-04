import numpy as np
import h5py as h5
import sys

arr = np.zeros((9,8));
#Element 0
arr[0,0] = 0;
arr[0,1] = 1;
arr[0,2] = 4;
arr[0,3] = 5;
arr[0,4] = 15+0;
arr[0,5] = 15+1;
arr[0,6] = 15+4;
arr[0,7] = 15+5;
#Element 1
arr[1,0] = 1;
arr[1,1] = 2;
arr[1,2] = 5;
arr[1,3] = 6;
arr[1,4] = 15+1;
arr[1,5] = 15+2;
arr[1,6] = 15+5;
arr[1,7] = 15+6;
#Element 2
arr[2,0] = 2;
arr[2,1] = 3;
arr[2,2] = 6;
arr[2,3] = 7;
arr[2,4] = 15+2;
arr[2,5] = 15+3;
arr[2,6] = 15+6;
arr[2,7] = 15+7;
#Element 3
arr[3,0] = 4;
arr[3,1] = 5;
arr[3,2] = 8;
arr[3,3] = 9;
arr[3,4] = 15+4;
arr[3,5] = 15+5;
arr[3,6] = 15+8;
arr[3,7] = 15+9;
#Element 4
arr[4,0] = 5;
arr[4,1] = 6;
arr[4,2] = 9;
arr[4,3] = 10;
arr[4,4] = 15+5;
arr[4,5] = 15+6;
arr[4,6] = 15+9;
arr[4,7] = 15+10;
#Element 5
arr[5,0] = 6;
arr[5,1] = 7;
arr[5,2] = 10;
arr[5,3] = 11;
arr[5,4] = 15+6;
arr[5,5] = 15+7;
arr[5,6] = 15+10;
arr[5,7] = 15+11;
#Element 6
arr[6,0] = 8;
arr[6,1] = 9;
arr[6,2] = 12;
arr[6,3] = 13;
arr[6,4] = 15+8;
arr[6,5] = 15+9;
arr[6,6] = 15+12;
arr[6,7] = 15+13;
#Element 7
arr[7,0] = 9;
arr[7,1] = 10;
arr[7,2] = 13;
arr[7,3] = 14;
arr[7,4] = 15+9;
arr[7,5] = 15+10;
arr[7,6] = 15+13;
arr[7,7] = 15+14;

#Element 8
arr[8,0] = 10;
arr[8,1] = 11;
arr[8,2] = 14;
arr[8,3] = 15;
arr[8,4] = 15+10;
arr[8,5] = 15+11;
arr[8,6] = 15+14;
arr[8,7] = 15+15;
filename = "test_grid.h5"
new_file = h5.File(filename,"w")
new_file.create_dataset("ien",data=arr)
