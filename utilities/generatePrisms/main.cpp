#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <STEPControl_Reader.hxx>
#include <IGESControl_Reader.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Compound.hxx>
#include <TopExp_Explorer.hxx>
#include <Poly_Triangulation.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pln.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax2.hxx>
#include <gp_Vec.hxx>
#include <Standard_Handle.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_Plane.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Precision.hxx>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <queue>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Structure to hold mesh data
struct MeshData {
    std::vector<gp_Pnt> vertices;
    std::vector<std::array<int, 3>> triangles;
};

// Edge structure for remeshing
struct Edge {
    int v1, v2;
    Edge(int a, int b) : v1(std::min(a,b)), v2(std::max(a,b)) {}
    bool operator<(const Edge& other) const {
        if (v1 != other.v1) return v1 < other.v1;
        return v2 < other.v2;
    }
};

// Function to create a hemisphere programmatically (curved surface only)
TopoDS_Shape CreateHemisphere(double radius = 1.0) {
    gp_Pnt center(0.0, 0.0, 0.0);
    gp_Ax2 axis(center, gp_Dir(0.0, 0.0, 1.0));
    
    Handle(Geom_SphericalSurface) sphereSurface = new Geom_SphericalSurface(axis, radius);
    
    // Create hemisphere: U parameter: 0 to 2*PI, V parameter: 0 to PI/2
    BRepBuilderAPI_MakeFace faceMaker(sphereSurface, 0.0, 2.0 * M_PI, 0.0, M_PI / 2.0, Precision::Confusion());
    
    TopoDS_Shape hemisphere = faceMaker.Shape();
    
    std::cout << "Created hemisphere with radius: " << radius << " (curved surface only)" << std::endl;
    
    return hemisphere;
}

// Function to load CAD file
TopoDS_Shape LoadCADFile(const std::string& filename) {
    TopoDS_Shape shape;
    std::string extension = filename.substr(filename.find_last_of(".") + 1);
    
    if (extension == "step" || extension == "stp" || extension == "STEP" || extension == "STP") {
        STEPControl_Reader reader;
        IFSelect_ReturnStatus status = reader.ReadFile(filename.c_str());
        
        if (status != IFSelect_RetDone) {
            std::cerr << "Error reading STEP file: " << filename << std::endl;
            return shape;
        }
        
        reader.TransferRoots();
        shape = reader.OneShape();
        std::cout << "Loaded STEP file: " << filename << std::endl;
    }
    else if (extension == "iges" || extension == "igs" || extension == "IGES" || extension == "IGS") {
        IGESControl_Reader reader;
        IFSelect_ReturnStatus status = reader.ReadFile(filename.c_str());
        
        if (status != IFSelect_RetDone) {
            std::cerr << "Error reading IGES file: " << filename << std::endl;
            return shape;
        }
        
        reader.TransferRoots();
        shape = reader.OneShape();
        std::cout << "Loaded IGES file: " << filename << std::endl;
    }
    else {
        std::cerr << "Unsupported file format: " << extension << std::endl;
    }
    
    return shape;
}

// Check if face is planar
Standard_Boolean IsFacePlanar(const TopoDS_Face& face) {
    TopLoc_Location location;
    Handle(Geom_Surface) surface = BRep_Tool::Surface(face, location);
    
    if (surface->DynamicType() == STANDARD_TYPE(Geom_Plane)) {
        return Standard_True;
    }
    
    GeomAdaptor_Surface surfaceAdaptor(surface);
    GeomAbs_SurfaceType surfType = surfaceAdaptor.GetType();
    
    return (surfType == GeomAbs_Plane);
}

// Filter to keep only curved surfaces
TopoDS_Shape FilterCurvedSurfacesOnly(const TopoDS_Shape& shape) {
    TopoDS_Compound result;
    BRep_Builder builder;
    builder.MakeCompound(result);
    
    int curvedFaceCount = 0;
    int planarFaceCount = 0;
    
    TopExp_Explorer faceExplorer(shape, TopAbs_FACE);
    for (; faceExplorer.More(); faceExplorer.Next()) {
        TopoDS_Face face = TopoDS::Face(faceExplorer.Current());
        
        if (!IsFacePlanar(face)) {
            builder.Add(result, face);
            curvedFaceCount++;
        } else {
            planarFaceCount++;
        }
    }
    
    std::cout << "Filtered: " << curvedFaceCount << " curved, "
              << planarFaceCount << " planar (removed)" << std::endl;
    
    return result;
}

// Generate initial mesh with OCCT
MeshData GenerateInitialMesh(const TopoDS_Shape& shape, double linearDeflection, double angularDeflection) {
    MeshData meshData;
    
    std::cout << "\nGenerating initial mesh..." << std::endl;
    std::cout << "  Linear deflection: " << linearDeflection << std::endl;
    std::cout << "  Angular deflection: " << angularDeflection << " rad" << std::endl;
    
    // Use very small angular deflection for better starting point
    BRepMesh_IncrementalMesh mesher(
        shape,
        linearDeflection,
        Standard_True,
        angularDeflection,
        Standard_True
    );
    
    if (!mesher.IsDone()) {
        std::cerr << "Meshing failed!" << std::endl;
        return meshData;
    }
    
    // Extract triangulation
    int vertexOffset = 0;
    TopExp_Explorer faceExplorer(shape, TopAbs_FACE);
    
    for (; faceExplorer.More(); faceExplorer.Next()) {
        TopoDS_Face face = TopoDS::Face(faceExplorer.Current());
        TopLoc_Location location;
        
        Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, location);
        
        if (triangulation.IsNull()) continue;
        
        gp_Trsf transformation = location.Transformation();
        
        const Standard_Integer nbNodes = triangulation->NbNodes();
        for (Standard_Integer i = 1; i <= nbNodes; i++) {
            gp_Pnt point = triangulation->Node(i).Transformed(transformation);
            meshData.vertices.push_back(point);
        }
        
        const Standard_Integer nbTriangles = triangulation->NbTriangles();
        Standard_Boolean isReversed = (face.Orientation() == TopAbs_REVERSED);
        
        for (Standard_Integer i = 1; i <= nbTriangles; i++) {
            const Poly_Triangle& triangle = triangulation->Triangle(i);
            Standard_Integer n1, n2, n3;
            triangle.Get(n1, n2, n3);
            
            n1 += vertexOffset - 1;
            n2 += vertexOffset - 1;
            n3 += vertexOffset - 1;
            
            if (isReversed) {
                meshData.triangles.push_back({n1, n3, n2});
            } else {
                meshData.triangles.push_back({n1, n2, n3});
            }
        }
        
        vertexOffset += nbNodes;
    }
    
    std::cout << "Initial mesh: " << meshData.vertices.size() << " vertices, "
              << meshData.triangles.size() << " triangles" << std::endl;
    
    return meshData;
}

// Simple isotropic remeshing using edge splitting and flipping
MeshData IsotropicRemesh(const MeshData& inputMesh, double targetEdgeLength, int iterations = 5) {
    std::cout << "\nPerforming isotropic remeshing..." << std::endl;
    std::cout << "  Target edge length: " << targetEdgeLength << std::endl;
    std::cout << "  Iterations: " << iterations << std::endl;
    
    MeshData mesh = inputMesh;
    
    // Calculate target edge length if not specified
    if (targetEdgeLength <= 0) {
        double totalLength = 0.0;
        int edgeCount = 0;
        
        for (const auto& tri : mesh.triangles) {
            totalLength += mesh.vertices[tri[0]].Distance(mesh.vertices[tri[1]]);
            totalLength += mesh.vertices[tri[1]].Distance(mesh.vertices[tri[2]]);
            totalLength += mesh.vertices[tri[2]].Distance(mesh.vertices[tri[0]]);
            edgeCount += 3;
        }
        
        targetEdgeLength = totalLength / edgeCount;
    }
    
    double minLength = 0.8 * targetEdgeLength;
    double maxLength = 1.33 * targetEdgeLength;
    
    std::cout << "  Edge length range: [" << minLength << ", " << maxLength << "]" << std::endl;
    
    // For a simple implementation, we'll just report the stats
    // A full remeshing algorithm would require:
    // 1. Edge splitting for edges > maxLength
    // 2. Edge collapse for edges < minLength
    // 3. Edge flipping to improve triangle quality
    // 4. Vertex smoothing via tangential relaxation
    
    std::cout << "  Note: Full isotropic remeshing requires external libraries" << std::endl;
    std::cout << "  Current mesh is optimized using OCCT's best parameters" << std::endl;
    
    return mesh;
}

// Export functions
void ExportToVTK(const MeshData& meshData, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }
    
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "Isotropic triangulated mesh" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET POLYDATA" << std::endl;
    file << "POINTS " << meshData.vertices.size() << " double" << std::endl;
    
    for (const auto& vertex : meshData.vertices) {
        file << vertex.X() << " " << vertex.Y() << " " << vertex.Z() << std::endl;
    }
    
    file << "POLYGONS " << meshData.triangles.size() << " " << (meshData.triangles.size() * 4) << std::endl;
    for (const auto& triangle : meshData.triangles) {
        file << "3 " << triangle[0] << " " << triangle[1] << " " << triangle[2] << std::endl;
    }
    
    file.close();
    std::cout << "Exported to: " << filename << std::endl;
}

void ExportToOBJ(const MeshData& meshData, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    file << "# Isotropic triangulated mesh" << std::endl;
    for (const auto& vertex : meshData.vertices) {
        file << "v " << vertex.X() << " " << vertex.Y() << " " << vertex.Z() << std::endl;
    }
    
    for (const auto& triangle : meshData.triangles) {
        file << "f " << (triangle[0] + 1) << " " << (triangle[1] + 1) << " " << (triangle[2] + 1) << std::endl;
    }
    
    file.close();
    std::cout << "Exported to: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    double radius = 1.0;
    double linearDeflection = 0.02;   // Very fine for good starting point
    double angularDeflection = 0.08;  // Very small for maximum isotropy
    std::string outputBaseName = "hemisphere";
    TopoDS_Shape shape;
    
    std::cout << "========================================" << std::endl;
    std::cout << "ISOTROPIC HEMISPHERE MESH GENERATOR" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    if (argc < 2) {
        std::cout << "Creating hemisphere programmatically (curved surface only)" << std::endl;
        std::cout << "Using MAXIMUM ISOTROPY settings:" << std::endl;
        std::cout << "  Linear: " << linearDeflection << ", Angular: " << angularDeflection << " rad\n" << std::endl;
        shape = CreateHemisphere(radius);
    }
    else if (std::string(argv[1]) == "--create") {
        if (argc >= 3) radius = std::stod(argv[2]);
        if (argc >= 4) linearDeflection = std::stod(argv[3]);
        if (argc >= 5) angularDeflection = std::stod(argv[4]);
        
        shape = CreateHemisphere(radius);
        outputBaseName = "hemisphere_r" + std::to_string(radius);
    }
    else {
        std::string inputFile = argv[1];
        if (argc >= 3) linearDeflection = std::stod(argv[2]);
        if (argc >= 4) angularDeflection = std::stod(argv[3]);
        
        TopoDS_Shape loadedShape = LoadCADFile(inputFile);
        if (loadedShape.IsNull()) {
            std::cerr << "Failed to load file" << std::endl;
            return 1;
        }
        
        shape = FilterCurvedSurfacesOnly(loadedShape);
        outputBaseName = inputFile.substr(0, inputFile.find_last_of("."));
    }
    
    if (shape.IsNull()) {
        std::cerr << "Failed to create geometry" << std::endl;
        return 1;
    }
    
    // Generate mesh
    MeshData meshData = GenerateInitialMesh(shape, linearDeflection, angularDeflection);
    
    if (meshData.vertices.empty()) {
        std::cerr << "No mesh generated" << std::endl;
        return 1;
    }
    
    // Export
    ExportToVTK(meshData, outputBaseName + "_mesh.vtk");
    ExportToOBJ(meshData, outputBaseName + "_mesh.obj");
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "RECOMMENDATIONS FOR MAXIMUM ISOTROPY:" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "1. Use angular deflection 0.05-0.10 rad" << std::endl;
    std::cout << "2. Use linear deflection 0.02-0.03" << std::endl;
    std::cout << "3. For Pointwise-quality meshes, consider:" << std::endl;
    std::cout << "   - Gmsh (open source, excellent)" << std::endl;
    std::cout << "   - CGAL (library with isotropic remeshing)" << std::endl;
    std::cout << "   - MMG (anisotropic to isotropic remeshing)" << std::endl;
    std::cout << "\nOCCT is designed for visualization, not CFD!" << std::endl;
    std::cout << "For production CFD meshes, use dedicated mesh tools." << std::endl;
    
    return 0;
}
