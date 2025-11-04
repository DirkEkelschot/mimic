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
#include <iostream>
#include <fstream>
#include <vector>
#include <set>

int main() {
    std::cout << "Generating tetrahedral tessellation around sphere..." << std::endl;
    
    // Create a box of 5x5x5 centered at origin
    gp_Pnt boxCorner(-2.5, -2.5, -2.5);
    TopoDS_Shape box = BRepPrimAPI_MakeBox(boxCorner, 5.0, 5.0, 5.0).Shape();
    
    // Create a sphere of radius 0.5 at the origin
    gp_Pnt sphereCenter(0.0, 0.0, 0.0);
    TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(sphereCenter, 0.5).Shape();
    
    // Subtract the sphere from the box to create the void region
    TopoDS_Shape meshRegion = BRepAlgoAPI_Cut(box, sphere).Shape();
    
    std::cout << "Box and sphere created. Generating mesh..." << std::endl;
    
    // Generate the mesh with specified deflection
    // Smaller deflection = finer mesh
    double linearDeflection = 0.2;
    double angularDeflection = 0.5;
    BRepMesh_IncrementalMesh mesh(meshRegion, linearDeflection, Standard_False, angularDeflection);
    
    std::cout << "Mesh generated. Extracting triangulation data..." << std::endl;
    
    // Extract triangulation data
    std::vector<gp_Pnt> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::set<std::array<double, 3>> uniqueVertices;
    
    // Explore all faces and extract triangulation
    for (TopExp_Explorer explorer(meshRegion, TopAbs_FACE); explorer.More(); explorer.Next()) {
        TopoDS_Face face = TopoDS::Face(explorer.Current());
        TopLoc_Location location;
        Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, location);
        
        if (!triangulation.IsNull()) {
            int numTriangles = triangulation->NbTriangles();
            int numNodes = triangulation->NbNodes();
            int vertexOffset = vertices.size();
            
            // Extract vertices
            for (int i = 1; i <= numNodes; i++) {
                gp_Pnt point = triangulation->Node(i).Transformed(location.Transformation());
                vertices.push_back(point);
            }
            
            // Extract triangles
            for (int i = 1; i <= numTriangles; i++) {
                int n1, n2, n3;
                triangulation->Triangle(i).Get(n1, n2, n3);
                
                // Adjust indices based on vertex offset
                triangles.push_back({vertexOffset + n1 - 1,
                                   vertexOffset + n2 - 1,
                                   vertexOffset + n3 - 1});
            }
        }
    }
    
    std::cout << "Extracted " << vertices.size() << " vertices and "
              << triangles.size() << " triangles" << std::endl;
    
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
    
    vtkFile << "\nCELLS " << triangles.size() << " " << triangles.size() * 4 << "\n";
    for (const auto& tri : triangles) {
        vtkFile << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }
    
    vtkFile << "\nCELL_TYPES " << triangles.size() << "\n";
    for (size_t i = 0; i < triangles.size(); i++) {
        vtkFile << "5\n"; // VTK_TRIANGLE
    }
    
    vtkFile.close();
    std::cout << "Mesh exported to sphere_mesh.vtk" << std::endl;
    
    // Also export as STL for 3D viewing
    StlAPI_Writer stlWriter;
    stlWriter.Write(meshRegion, "sphere_mesh.stl");
    std::cout << "Mesh also exported to sphere_mesh.stl" << std::endl;
    
    // Export as STEP for CAD compatibility
    STEPControl_Writer stepWriter;
    stepWriter.Transfer(meshRegion, STEPControl_AsIs);
    stepWriter.Write("sphere_mesh.step");
    std::cout << "Geometry exported to sphere_mesh.step" << std::endl;
    
    std::cout << "\nDone! Generated files:" << std::endl;
    std::cout << "  - sphere_mesh.vtk  (for ParaView/VisIt)" << std::endl;
    std::cout << "  - sphere_mesh.stl  (for 3D viewers)" << std::endl;
    std::cout << "  - sphere_mesh.step (for CAD software)" << std::endl;
    
    return 0;
}
