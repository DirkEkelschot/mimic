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
#include <map>
#include <cmath>
#include <algorithm>
#include <array>

// Simple structure for a tetrahedron
struct Tetrahedron {
    int v[4]; // Four vertex indices
    
    Tetrahedron(int v0, int v1, int v2, int v3) {
        v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
    }
};

// Structure for circumsphere
struct Circumsphere {
    gp_Pnt center;
    double radiusSquared;
};

// Calculate circumsphere of a tetrahedron
Circumsphere getCircumsphere(const gp_Pnt& p0, const gp_Pnt& p1,
                              const gp_Pnt& p2, const gp_Pnt& p3) {
    // Using matrix determinant method
    double x0 = p0.X(), y0 = p0.Y(), z0 = p0.Z();
    double x1 = p1.X(), y1 = p1.Y(), z1 = p1.Z();
    double x2 = p2.X(), y2 = p2.Y(), z2 = p2.Z();
    double x3 = p3.X(), y3 = p3.Y(), z3 = p3.Z();
    
    double a11 = x1 - x0, a12 = y1 - y0, a13 = z1 - z0;
    double a21 = x2 - x0, a22 = y2 - y0, a23 = z2 - z0;
    double a31 = x3 - x0, a32 = y3 - y0, a33 = z3 - z0;
    
    double b1 = (x1*x1 + y1*y1 + z1*z1 - x0*x0 - y0*y0 - z0*z0) / 2.0;
    double b2 = (x2*x2 + y2*y2 + z2*z2 - x0*x0 - y0*y0 - z0*z0) / 2.0;
    double b3 = (x3*x3 + y3*y3 + z3*z3 - x0*x0 - y0*y0 - z0*z0) / 2.0;
    
    double det = a11 * (a22 * a33 - a23 * a32) -
                 a12 * (a21 * a33 - a23 * a31) +
                 a13 * (a21 * a32 - a22 * a31);
    
    if (std::abs(det) < 1e-10) {
        // Degenerate tetrahedron
        Circumsphere cs;
        cs.center = gp_Pnt((x0+x1+x2+x3)/4, (y0+y1+y2+y3)/4, (z0+z1+z2+z3)/4);
        cs.radiusSquared = 1e20;
        return cs;
    }
    
    double cx = (b1 * (a22 * a33 - a23 * a32) -
                 a12 * (b2 * a33 - b3 * a32) +
                 a13 * (b2 * a32 - b3 * a22)) / det;
    double cy = (a11 * (b2 * a33 - b3 * a32) -
                 b1 * (a21 * a33 - a23 * a31) +
                 a13 * (a21 * b3 - a31 * b2)) / det;
    double cz = (a11 * (a22 * b3 - b2 * a32) -
                 a12 * (a21 * b3 - b2 * a31) +
                 b1 * (a21 * a32 - a22 * a31)) / det;
    
    Circumsphere cs;
    cs.center = gp_Pnt(cx, cy, cz);
    cs.radiusSquared = (cx-x0)*(cx-x0) + (cy-y0)*(cy-y0) + (cz-z0)*(cz-z0);
    
    return cs;
}

// Check if point is inside circumsphere
bool isInsideCircumsphere(const gp_Pnt& p, const Circumsphere& cs) {
    double dx = p.X() - cs.center.X();
    double dy = p.Y() - cs.center.Y();
    double dz = p.Z() - cs.center.Z();
    double distSquared = dx*dx + dy*dy + dz*dz;
    return distSquared < cs.radiusSquared - 1e-10;
}

// Bowyer-Watson Delaunay tetrahedralization
void delaunayTetrahedralization(
    std::vector<gp_Pnt>& vertices,
    std::vector<Tetrahedron>& tetrahedra)
{
    std::cout << "Performing Delaunay tetrahedralization on " << vertices.size() << " points..." << std::endl;
    
    if (vertices.size() < 4) {
        std::cout << "Not enough vertices for tetrahedralization" << std::endl;
        return;
    }
    
    // Find bounding box
    double minX = vertices[0].X(), maxX = vertices[0].X();
    double minY = vertices[0].Y(), maxY = vertices[0].Y();
    double minZ = vertices[0].Z(), maxZ = vertices[0].Z();
    
    for (const auto& v : vertices) {
        minX = std::min(minX, v.X()); maxX = std::max(maxX, v.X());
        minY = std::min(minY, v.Y()); maxY = std::max(maxY, v.Y());
        minZ = std::min(minZ, v.Z()); maxZ = std::max(maxZ, v.Z());
    }
    
    // Create super-tetrahedron that contains all points
    double dx = maxX - minX, dy = maxY - minY, dz = maxZ - minZ;
    double scale = 10.0 * std::max({dx, dy, dz});
    
    int nOriginal = vertices.size();
    vertices.push_back(gp_Pnt((minX + maxX) / 2.0, (minY + maxY) / 2.0 - scale, (minZ + maxZ) / 2.0));
    vertices.push_back(gp_Pnt((minX + maxX) / 2.0 + scale, (minY + maxY) / 2.0 + scale, (minZ + maxZ) / 2.0));
    vertices.push_back(gp_Pnt((minX + maxX) / 2.0 - scale, (minY + maxY) / 2.0 + scale, (minZ + maxZ) / 2.0));
    vertices.push_back(gp_Pnt((minX + maxX) / 2.0, (minY + maxY) / 2.0, (minZ + maxZ) / 2.0 + scale));
    
    // Initialize with super-tetrahedron
    tetrahedra.clear();
    tetrahedra.push_back(Tetrahedron(nOriginal, nOriginal+1, nOriginal+2, nOriginal+3));
    
    std::cout << "Inserting points..." << std::endl;
    
    // Insert each point
    for (int i = 0; i < nOriginal; i++) {
        if (i % 100 == 0) {
            std::cout << "  Inserted " << i << "/" << nOriginal << " points" << std::endl;
        }
        
        const gp_Pnt& point = vertices[i];
        std::vector<Tetrahedron> badTets;
        std::set<std::array<int, 3>> boundary;
        
        // Find all tetrahedra whose circumsphere contains the point
        for (const auto& tet : tetrahedra) {
            Circumsphere cs = getCircumsphere(
                vertices[tet.v[0]], vertices[tet.v[1]],
                vertices[tet.v[2]], vertices[tet.v[3]]
            );
            
            if (isInsideCircumsphere(point, cs)) {
                badTets.push_back(tet);
                
                // Add faces to boundary (each face appears twice, once from each adjacent tet)
                std::array<std::array<int, 3>, 4> faces = {{
                    {tet.v[0], tet.v[1], tet.v[2]},
                    {tet.v[0], tet.v[1], tet.v[3]},
                    {tet.v[0], tet.v[2], tet.v[3]},
                    {tet.v[1], tet.v[2], tet.v[3]}
                }};
                
                for (auto face : faces) {
                    std::sort(face.begin(), face.end());
                    if (boundary.count(face)) {
                        boundary.erase(face);
                    } else {
                        boundary.insert(face);
                    }
                }
            }
        }
        
        // Remove bad tetrahedra
        tetrahedra.erase(
            std::remove_if(tetrahedra.begin(), tetrahedra.end(),
                [&badTets](const Tetrahedron& tet) {
                    return std::find_if(badTets.begin(), badTets.end(),
                        [&tet](const Tetrahedron& bad) {
                            return tet.v[0] == bad.v[0] && tet.v[1] == bad.v[1] &&
                                   tet.v[2] == bad.v[2] && tet.v[3] == bad.v[3];
                        }) != badTets.end();
                }),
            tetrahedra.end()
        );
        
        // Create new tetrahedra from boundary faces to new point
        for (const auto& face : boundary) {
            tetrahedra.push_back(Tetrahedron(face[0], face[1], face[2], i));
        }
    }
    
    std::cout << "Removing super-tetrahedron vertices..." << std::endl;
    
    // Remove tetrahedra that use super-tetrahedron vertices
    tetrahedra.erase(
        std::remove_if(tetrahedra.begin(), tetrahedra.end(),
            [nOriginal](const Tetrahedron& tet) {
                return tet.v[0] >= nOriginal || tet.v[1] >= nOriginal ||
                       tet.v[2] >= nOriginal || tet.v[3] >= nOriginal;
            }),
        tetrahedra.end()
    );
    
    // Remove super-tetrahedron vertices
    vertices.resize(nOriginal);
    
    std::cout << "Delaunay tetrahedralization complete: " << tetrahedra.size() << " tetrahedra" << std::endl;
}

int main() {
    std::cout << "Generating body-fitted tetrahedral mesh around sphere..." << std::endl;
    
    double boxSize = 5.0;
    double sphereRadius = 0.5;
    double meshDensity = 0.15; // Smaller = finer mesh
    
    // Create geometry
    gp_Pnt boxCorner(-boxSize/2, -boxSize/2, -boxSize/2);
    TopoDS_Shape box = BRepPrimAPI_MakeBox(boxCorner, boxSize, boxSize, boxSize).Shape();
    
    gp_Pnt sphereCenter(0.0, 0.0, 0.0);
    TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(sphereCenter, sphereRadius).Shape();
    
    TopoDS_Shape meshRegion = BRepAlgoAPI_Cut(box, sphere).Shape();
    
    std::cout << "Generating surface mesh..." << std::endl;
    
    // Generate fine surface mesh
    BRepMesh_IncrementalMesh surfaceMesh(meshRegion, meshDensity, Standard_False, 0.5);
    
    // Extract surface vertices
    std::vector<gp_Pnt> vertices;
    std::map<std::array<double, 3>, int> vertexMap;
    
    for (TopExp_Explorer explorer(meshRegion, TopAbs_FACE); explorer.More(); explorer.Next()) {
        TopoDS_Face face = TopoDS::Face(explorer.Current());
        TopLoc_Location location;
        Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, location);
        
        if (!triangulation.IsNull()) {
            for (int i = 1; i <= triangulation->NbNodes(); i++) {
                gp_Pnt point = triangulation->Node(i).Transformed(location.Transformation());
                
                // Round to avoid duplicates
                std::array<double, 3> key = {
                    std::round(point.X() * 1e6) / 1e6,
                    std::round(point.Y() * 1e6) / 1e6,
                    std::round(point.Z() * 1e6) / 1e6
                };
                
                if (vertexMap.find(key) == vertexMap.end()) {
                    vertexMap[key] = vertices.size();
                    vertices.push_back(point);
                }
            }
        }
    }
    
    std::cout << "Extracted " << vertices.size() << " surface vertices" << std::endl;
    
    // Add interior points for better mesh quality
    std::cout << "Adding interior points..." << std::endl;
    int interiorDivisions = 8;
    double step = boxSize / interiorDivisions;
    
    for (int i = 1; i < interiorDivisions; i++) {
        for (int j = 1; j < interiorDivisions; j++) {
            for (int k = 1; k < interiorDivisions; k++) {
                double x = -boxSize/2.0 + i * step;
                double y = -boxSize/2.0 + j * step;
                double z = -boxSize/2.0 + k * step;
                
                gp_Pnt p(x, y, z);
                double dist = std::sqrt(x*x + y*y + z*z);
                
                // Only add if outside sphere with some margin
                if (dist > sphereRadius + meshDensity) {
                    vertices.push_back(p);
                }
            }
        }
    }
    
    std::cout << "Total vertices (surface + interior): " << vertices.size() << std::endl;
    
    // Perform Delaunay tetrahedralization
    std::vector<Tetrahedron> tetrahedra;
    delaunayTetrahedralization(vertices, tetrahedra);
    
    // Export to VTK
    std::cout << "Exporting to VTK..." << std::endl;
    std::ofstream vtkFile("sphere_mesh_bodyfitted.vtk");
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Body-fitted tetrahedral mesh around sphere\n";
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
        vtkFile << "10\n"; // VTK_TETRA
    }
    
    vtkFile.close();
    
    // Also export surface for reference
    StlAPI_Writer stlWriter;
    stlWriter.Write(meshRegion, "sphere_surface.stl");
    
    std::cout << "\n=== Mesh Generation Complete ===" << std::endl;
    std::cout << "Vertices: " << vertices.size() << std::endl;
    std::cout << "Tetrahedra: " << tetrahedra.size() << std::endl;
    std::cout << "\nGenerated files:" << std::endl;
    std::cout << "  - sphere_mesh_bodyfitted.vtk  (body-fitted tetrahedral mesh)" << std::endl;
    std::cout << "  - sphere_surface.stl          (surface reference)" << std::endl;
    std::cout << "\nVisualization in ParaView:" << std::endl;
    std::cout << "  1. Open sphere_mesh_bodyfitted.vtk" << std::endl;
    std::cout << "  2. Apply 'Clip' or 'Slice' filter to see interior" << std::endl;
    std::cout << "  3. Apply 'Extract Surface' to see outer boundary" << std::endl;
    std::cout << "  4. Color by 'vtkGhostType' to verify mesh quality" << std::endl;
    std::cout << "\nAdjust 'meshDensity' (line 185) for mesh fineness (smaller = finer)" << std::endl;
    
    return 0;
}
