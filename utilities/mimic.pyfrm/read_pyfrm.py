import h5py

def print_hdf5_item(name, obj):
    if isinstance(obj, h5py.Group):
        print(f'Group: {name}')
    elif isinstance(obj, h5py.Dataset):
        print(f'Dataset: {name}, shape: {obj.shape}, dtype: {obj.dtype}')

with h5py.File('mesh/hemisphere.pyfrm', 'r') as file:
    file.visititems(print_hdf5_item)