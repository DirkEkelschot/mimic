
#============================BOUNDARY LAYER MESH GENERATION=====================================
Run c_stl2step.py geometry.stl -> to convert stl to step file. outputs geometry.step

./baseMeshGenerator geometry.step output.vtk 5 0.01 1.2 -> outputs the boundary layer mesh but also a file called output_outer.vtk which contains the outer triangles of the boundary layer mesh.

#============================VOLUME  MESH GENERATION============================================
To generate a volume mesh with output_outer.vtk as input
run python3 vol_mesh_working.py
/Users/dirk/Software/TetGen/build/tetgen -pq1.4/15YAa5.0V flow_domain.poly
python3 convert_ele_2_vtk.py