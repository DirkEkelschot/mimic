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

// TetGen library
#include "tetgen.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <cmath>

int main() {
    std::cout << "=== Body-Fitted Tetrahedral Mesh Generator ===" << std::endl;
    std::cout << "Using TetGen library directly" << std::endl << std::endl;
    
    double boxSize = 5.0;
    double sphereRadius = 0.5;
    double meshDensity = 0.2;  // Surface mesh density
    
    // Create geometry using OpenCASCADE
    std::cout << "Creating geometry..." << std::endl;
    gp_Pnt boxCorner(-boxSize/2, -boxSize/2, -boxSize/2);
    TopoDS_Shape box = BRepPrimAPI_MakeBox(boxCorner, boxSize, boxSize, boxSize).Shape();
    
    gp_Pnt sphereCenter(0.0, 0.0, 0.0);
    TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(sphereCenter, sphereRadius).Shape();
    
    TopoDS_Shape meshRegion = BRepAlgoAPI_Cut(box, sphere).Shape();
    
    std::cout << "Generating surface mesh..." << std::endl;
    
    // Generate surface mesh with OpenCASCADE
    BRepMesh_IncrementalMesh surfaceMesh(meshRegion, meshDensity, Standard_False, 0.5);
    
    // Extract surface triangulation
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::map<std::array<double, 3>, int> vertexMap;
    
    for (TopExp_Explorer explorer(meshRegion, TopAbs_FACE); explorer.More(); explorer.Next()) {
        TopoDS_Face face = TopoDS::Face(explorer.Current());
        TopLoc_Location location;
        Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, location);
        
        if (!triangulation.IsNull()) {
            // Extract vertices
            for (int i = 1; i <= triangulation->NbNodes(); i++) {
                gp_Pnt point = triangulation->Node(i).Transformed(location.Transformation());
                
                std::array<double, 3> key = {
                    std::round(point.X() * 1e6) / 1e6,
                    std::round(point.Y() * 1e6) / 1e6,
                    std::round(point.Z() * 1e6) / 1e6
                };
                
                if (vertexMap.find(key) == vertexMap.end()) {
                    vertexMap[key] = vertices.size();
                    vertices.push_back({point.X(), point.Y(), point.Z()});
                }
            }
            
            // Extract triangles
            for (int i = 1; i <= triangulation->NbTriangles(); i++) {
                int n1, n2, n3;
                triangulation->Triangle(i).Get(n1, n2, n3);
                
                gp_Pnt p1 = triangulation->Node(n1).Transformed(location.Transformation());
                gp_Pnt p2 = triangulation->Node(n2).Transformed(location.Transformation());
                gp_Pnt p3 = triangulation->Node(n3).Transformed(location.Transformation());
                
                auto getVertexIndex = [&](const gp_Pnt& p) {
                    std::array<double, 3> key = {
                        std::round(p.X() * 1e6) / 1e6,
                        std::round(p.Y() * 1e6) / 1e6,
                        std::round(p.Z() * 1e6) / 1e6
                    };
                    return vertexMap[key];
                };
                
                int idx1 = getVertexIndex(p1);
                int idx2 = getVertexIndex(p2);
                int idx3 = getVertexIndex(p3);
                
                triangles.push_back({idx1, idx2, idx3});
            }
        }
    }
    
    std::cout << "Extracted " << vertices.size() << " surface vertices and "
              << triangles.size() << " surface triangles" << std::endl;
    
    // Prepare TetGen input
    std::cout << "\nPreparing TetGen input..." << std::endl;
    
    tetgenio in, out;
    
    // Set vertices
    in.numberofpoints = vertices.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    
    for (int i = 0; i < vertices.size(); i++) {
        in.pointlist[i * 3 + 0] = vertices[i][0];
        in.pointlist[i * 3 + 1] = vertices[i][1];
        in.pointlist[i * 3 + 2] = vertices[i][2];
    }
    
    // Set facets (surface triangles)
    in.numberoffacets = triangles.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    
    for (int i = 0; i < triangles.size(); i++) {
        tetgenio::facet *f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[1];
        f->numberofholes = 0;
        f->holelist = NULL;
        
        tetgenio::polygon *p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[3];
        p->vertexlist[0] = triangles[i][0];
        p->vertexlist[1] = triangles[i][1];
        p->vertexlist[2] = triangles[i][2];
    }
    
    // CRITICAL: Define hole point inside the sphere
    // This tells TetGen to NOT generate tetrahedra inside the sphere
    in.numberofholes = 1;
    in.holelist = new REAL[in.numberofholes * 3];
    in.holelist[0] = 0.0;  // x coordinate inside sphere
    in.holelist[1] = 0.0;  // y coordinate inside sphere
    in.holelist[2] = 0.0;  // z coordinate inside sphere
    
    std::cout << "Added hole point at (0, 0, 0) to create void inside sphere" << std::endl;
    
    // Call TetGen
    std::cout << "Running TetGen..." << std::endl;
    std::cout << "Parameters: -pq1.414a0.1" << std::endl;
    std::cout << "  -p : Tetrahedralize a piecewise linear complex" << std::endl;
    std::cout << "  -q1.414 : Quality mesh (radius-edge ratio)" << std::endl;
    std::cout << "  -a0.1 : Maximum tetrahedron volume" << std::endl;
    std::cout << std::endl;
    
    try {
        // TetGen switches:
        // p - Tetrahedralizes a piecewise linear complex (PLC)
        // q - Quality mesh generation with specified quality bound
        // a - Applies maximum tetrahedron volume constraint
        // Q - Quiet mode
        tetrahedralize((char*)"pq1.414a0.1Q", &in, &out);
        
        std::cout << "TetGen completed successfully!" << std::endl;
        std::cout << "Generated " << out.numberofpoints << " vertices and "
                  << out.numberoftetrahedra << " tetrahedra" << std::endl;
        
    } catch (int e) {
        std::cerr << "TetGen error code: " << e << std::endl;
        return 1;
    }
    
    // Export to VTK format
    std::cout << "\nExporting to VTK format..." << std::endl;
    
    std::ofstream vtkFile("sphere_mesh_bodyfitted.vtk");
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Body-fitted tetrahedral mesh generated by TetGen\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write vertices
    vtkFile << "POINTS " << out.numberofpoints << " float\n";
    for (int i = 0; i < out.numberofpoints; i++) {
        vtkFile << out.pointlist[i * 3 + 0] << " "
                << out.pointlist[i * 3 + 1] << " "
                << out.pointlist[i * 3 + 2] << "\n";
    }
    
    // Write tetrahedra
    vtkFile << "\nCELLS " << out.numberoftetrahedra << " "
            << out.numberoftetrahedra * 5 << "\n";
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        vtkFile << "4 "
                << out.tetrahedronlist[i * 4 + 0] << " "
                << out.tetrahedronlist[i * 4 + 1] << " "
                << out.tetrahedronlist[i * 4 + 2] << " "
                << out.tetrahedronlist[i * 4 + 3] << "\n";
    }
    
    // Write cell types
    vtkFile << "\nCELL_TYPES " << out.numberoftetrahedra << "\n";
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        vtkFile << "10\n";  // VTK_TETRA
    }
    
    vtkFile.close();
    
    std::cout << "VTK file written: sphere_mesh_bodyfitted.vtk" << std::endl;
    
    // Also export the surface for reference
    StlAPI_Writer stlWriter;
    stlWriter.Write(meshRegion, "sphere_surface.stl");
    std::cout << "Surface STL written: sphere_surface.stl" << std::endl;
    
    // Print statistics
    std::cout << "\n=== Mesh Statistics ===" << std::endl;
    std::cout << "Input surface vertices:  " << vertices.size() << std::endl;
    std::cout << "Input surface triangles: " << triangles.size() << std::endl;
    std::cout << "Output volume vertices:  " << out.numberofpoints << std::endl;
    std::cout << "Output tetrahedra:       " << out.numberoftetrahedra << std::endl;
    
    if (out.numberoftetrahedronattributes > 0) {
        std::cout << "Tetrahedron attributes:  " << out.numberoftetrahedronattributes << std::endl;
    }
    
    std::cout << "\n=== Files Generated ===" << std::endl;
    std::cout << "  sphere_mesh_bodyfitted.vtk - Body-fitted tetrahedral mesh" << std::endl;
    std::cout << "  sphere_surface.stl         - Surface reference mesh" << std::endl;
    
    std::cout << "\n=== Visualization ===" << std::endl;
    std::cout << "Open in ParaView:" << std::endl;
    std::cout << "  paraview sphere_mesh_bodyfitted.vtk" << std::endl;
    std::cout << "\nRecommended filters:" << std::endl;
    std::cout << "  - Extract Surface: See the outer boundary" << std::endl;
    std::cout << "  - Clip: Cut through the mesh" << std::endl;
    std::cout << "  - Slice: View cross-sections" << std::endl;
    std::cout << "  - Shrink: Separate tetrahedra visually" << std::endl;
    
    std::cout << "\n=== Adjusting Mesh Density ===" << std::endl;
    std::cout << "Edit the code to change:" << std::endl;
    std::cout << "  meshDensity (line 28): Surface triangle size (smaller = finer)" << std::endl;
    std::cout << "  TetGen 'a' parameter (line 126): Max tet volume (smaller = finer)" << std::endl;
    std::cout << "  Example: \"pq1.414a0.01\" for finer mesh" << std::endl;
    std::cout << "  Example: \"pq1.414a0.5\" for coarser mesh" << std::endl;
    
    return 0;
}
