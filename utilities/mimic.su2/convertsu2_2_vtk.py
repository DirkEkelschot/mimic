import meshio


mesh = meshio.read("mimic_hemi.su2")

print("Points:", mesh.points.shape)
#print("Cells:", {cell_type: len(data) for cell_type, data in mesh.cells})
print("Point data:", mesh.point_data.keys())
print("Cell data:", mesh.cell_data.keys())

mesh.write("output.vtk")
