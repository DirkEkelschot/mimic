#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <gp_Pnt.hxx>
#include <STEPControl_Writer.hxx>
#include <StlAPI_Writer.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>

// Simple structure for a tetrahedron
struct Tetrahedron {
    int v[4]; // Four vertex indices
    
    Tetrahedron(int v0, int v1, int v2, int v3) {
        v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
    }
};

// Check if a point is inside the sphere
bool isInsideSphere(const gp_Pnt& p, double radius) {
    double dist = std::sqrt(p.X() * p.X() + p.Y() * p.Y() + p.Z() * p.Z());
    return dist < radius;
}

// Check if a point is inside the box
bool isInsideBox(const gp_Pnt& p, double size) {
    return (std::abs(p.X()) <= size/2.0 &&
            std::abs(p.Y()) <= size/2.0 &&
            std::abs(p.Z()) <= size/2.0);
}

// Simple tetrahedral mesh generator
void generateTetrahedralMesh(
    std::vector<gp_Pnt>& vertices,
    std::vector<Tetrahedron>& tetrahedra,
    double boxSize,
    double sphereRadius,
    int divisions)
{
    std::cout << "Generating volumetric tetrahedral mesh..." << std::endl;
    
    // Create a regular grid of points
    double step = boxSize / divisions;
    std::map<std::tuple<int,int,int>, int> gridToIndex;
    
    // Generate grid points (only those outside sphere but inside box)
    for (int i = 0; i <= divisions; i++) {
        for (int j = 0; j <= divisions; j++) {
            for (int k = 0; k <= divisions; k++) {
                double x = -boxSize/2.0 + i * step;
                double y = -boxSize/2.0 + j * step;
                double z = -boxSize/2.0 + k * step;
                
                gp_Pnt p(x, y, z);
                
                // Only add points outside the sphere
                if (!isInsideSphere(p, sphereRadius)) {
                    int idx = vertices.size();
                    vertices.push_back(p);
                    gridToIndex[std::make_tuple(i, j, k)] = idx;
                }
            }
        }
    }
    
    std::cout << "Generated " << vertices.size() << " grid vertices" << std::endl;
    
    // Generate tetrahedra from grid cubes
    // Each cube is divided into 6 tetrahedra
    for (int i = 0; i < divisions; i++) {
        for (int j = 0; j < divisions; j++) {
            for (int k = 0; k < divisions; k++) {
                // Get the 8 corners of the cube
                std::vector<int> corners;
                for (int di = 0; di <= 1; di++) {
                    for (int dj = 0; dj <= 1; dj++) {
                        for (int dk = 0; dk <= 1; dk++) {
                            auto key = std::make_tuple(i+di, j+dj, k+dk);
                            if (gridToIndex.find(key) != gridToIndex.end()) {
                                corners.push_back(gridToIndex[key]);
                            } else {
                                corners.push_back(-1); // Point inside sphere
                            }
                        }
                    }
                }
                
                // Skip if any corner is inside the sphere
                bool hasInvalid = false;
                for (int c : corners) {
                    if (c == -1) {
                        hasInvalid = true;
                        break;
                    }
                }
                if (hasInvalid) continue;
                
                // Cube corners indices:
                // 0: (0,0,0), 1: (1,0,0), 2: (0,1,0), 3: (1,1,0)
                // 4: (0,0,1), 5: (1,0,1), 6: (0,1,1), 7: (1,1,1)
                int v0 = corners[0], v1 = corners[1], v2 = corners[2], v3 = corners[3];
                int v4 = corners[4], v5 = corners[5], v6 = corners[6], v7 = corners[7];
                
                // Divide cube into 6 tetrahedra (using a consistent pattern)
                tetrahedra.push_back(Tetrahedron(v0, v1, v2, v5));
                tetrahedra.push_back(Tetrahedron(v0, v2, v4, v5));
                tetrahedra.push_back(Tetrahedron(v2, v4, v5, v6));
                tetrahedra.push_back(Tetrahedron(v2, v3, v5, v7));
                tetrahedra.push_back(Tetrahedron(v2, v5, v6, v7));
                tetrahedra.push_back(Tetrahedron(v1, v2, v3, v5));
            }
        }
    }
    
    std::cout << "Generated " << tetrahedra.size() << " tetrahedra" << std::endl;
}

int main() {
    std::cout << "Generating tetrahedral tessellation around sphere..." << std::endl;
    
    double boxSize = 5.0;
    double sphereRadius = 0.5;
    int divisions = 20; // Adjust for mesh density (higher = finer mesh)
    
    std::vector<gp_Pnt> vertices;
    std::vector<Tetrahedron> tetrahedra;
    
    // Generate the tetrahedral mesh
    generateTetrahedralMesh(vertices, tetrahedra, boxSize, sphereRadius, divisions);
    
    // Export to VTK format for visualization
    std::ofstream vtkFile("sphere_mesh.vtk");
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Tetrahedral mesh around sphere\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";
    vtkFile << "POINTS " << vertices.size() << " float\n";
    
    for (const auto& vertex : vertices) {
        vtkFile << vertex.X() << " " << vertex.Y() << " " << vertex.Z() << "\n";
    }
    
    vtkFile << "\nCELLS " << tetrahedra.size() << " " << tetrahedra.size() * 5 << "\n";
    for (const auto& tet : tetrahedra) {
        vtkFile << "4 " << tet.v[0] << " " << tet.v[1] << " " << tet.v[2] << " " << tet.v[3] << "\n";
    }
    
    vtkFile << "\nCELL_TYPES " << tetrahedra.size() << "\n";
    for (size_t i = 0; i < tetrahedra.size(); i++) {
        vtkFile << "10\n"; // VTK_TETRA (tetrahedron)
    }
    
    vtkFile.close();
    std::cout << "Tetrahedral mesh exported to sphere_mesh.vtk" << std::endl;
    
    // Also export surface mesh for comparison
    TopoDS_Shape box = BRepPrimAPI_MakeBox(gp_Pnt(-boxSize/2, -boxSize/2, -boxSize/2),
                                            boxSize, boxSize, boxSize).Shape();
    TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), sphereRadius).Shape();
    TopoDS_Shape meshRegion = BRepAlgoAPI_Cut(box, sphere).Shape();
    
    StlAPI_Writer stlWriter;
    stlWriter.Write(meshRegion, "sphere_surface.stl");
    std::cout << "Surface mesh exported to sphere_surface.stl" << std::endl;
    
    std::cout << "\nDone! Generated files:" << std::endl;
    std::cout << "  - sphere_mesh.vtk       (tetrahedral volume mesh for ParaView)" << std::endl;
    std::cout << "  - sphere_surface.stl    (surface mesh for reference)" << std::endl;
    std::cout << "\nTo visualize in ParaView:" << std::endl;
    std::cout << "  1. Open sphere_mesh.vtk" << std::endl;
    std::cout << "  2. Apply 'Extract Surface' filter to see outer surface" << std::endl;
    std::cout << "  3. Or use 'Slice' filter to see interior tetrahedra" << std::endl;
    std::cout << "  4. Adjust 'divisions' variable (line 136) for mesh density" << std::endl;
    
    return 0;
}
