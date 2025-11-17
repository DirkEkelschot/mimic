import meshio

# Read TetGen output files
mesh = meshio.read("flow_domain.1.ele")

# Write to VTK
meshio.write("flow_domain.vtk", mesh)

print(f"Converted to VTK with {mesh.points.shape[0]} points and {mesh.cells[0].data.shape[0]} cells")