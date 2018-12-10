import h5py
import numpy as np


data_i = np.zeros((3,3,3))
data_j = np.zeros((3,3,3))
data_k = np.zeros((3,3,3))
data_0 = np.zeros((3,3,3))
x = np.array(range(3), dtype=float)

for i in range(3):
    for j in range(3):
        for k in range(3):
            data_i[i,j,k] = i
            data_j[i,j,k] = j
            data_k[i,j,k] = k


with h5py.File("read_test.h5", "w") as f:
    f.create_dataset("x", data=x)
    f.create_dataset("y", data=x/2)
    f.create_dataset("z", data=x+1)
    f.create_dataset("magnetic_field_x", data=data_i)
    f.create_dataset("magnetic_field_y", data=data_j)
    f.create_dataset("magnetic_field_z", data=data_k)
    f.create_dataset("electric_field_x", data=data_0+2)
    f.create_dataset("electric_field_y", data=data_0+1)
    f.create_dataset("electric_field_z", data=data_0)
