import numpy as np
import meshio

# Read SU2 file
mesh = meshio.read("mimic_hemi.su2")

# Convert all integer arrays to standard Int32 type
for key in mesh.point_data:
    if np.issubdtype(mesh.point_data[key].dtype, np.integer):
        mesh.point_data[key] = mesh.point_data[key].astype(np.int32)

for key in mesh.cell_data:
    for i, data in enumerate(mesh.cell_data[key]):
        if np.issubdtype(data.dtype, np.integer):
            mesh.cell_data[key][i] = data.astype(np.int32)

# Write corrected VTK file
mesh.write("output_v2.vtu", file_format="vtu")  # Explicitly specify VTK format
