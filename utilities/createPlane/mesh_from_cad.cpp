#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <gp_Pnt.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>

// STEP and IGES readers
#include <STEPControl_Reader.hxx>
#include <IGESControl_Reader.hxx>
#include <StlAPI_Writer.hxx>

// TetGen library
#include "tetgen.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <cmath>
#include <algorithm>
#include <string>

int main(int argc, char* argv[]) {
    std::cout << "=== Tetrahedral Mesh Around CAD Geometry ===" << std::endl;
    std::cout << "Generates volume mesh AROUND imported CAD (external flow domain)" << std::endl << std::endl;
    
    // Parse arguments
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <input.step|.iges> [box_scale] [mesh_density] [tetgen_params]" << std::endl;
        std::cout << "\nParameters:" << std::endl;
        std::cout << "  box_scale     : Bounding box size multiplier (default: 3.0)" << std::endl;
        std::cout << "                  3.0 = box is 3x larger than CAD in each dimension" << std::endl;
        std::cout << "  mesh_density  : Surface mesh density (default: 0.2, smaller = finer)" << std::endl;
        std::cout << "  tetgen_params : TetGen quality parameters (default: pq1.414a0.1)" << std::endl;
        std::cout << "\nExamples:" << std::endl;
        std::cout << "  " << argv[0] << " airplane.step" << std::endl;
        std::cout << "  " << argv[0] << " car.step 5.0 0.15" << std::endl;
        std::cout << "  " << argv[0] << " turbine.step 4.0 0.1 pq1.2a0.05" << std::endl;
        std::cout << "\nThis creates a box around your CAD and meshes the SPACE BETWEEN" << std::endl;
        std::cout << "the box and your CAD geometry (like air/water flow around an object)." << std::endl;
        return 1;
    }
    
    std::string inputFile = argv[1];
    double boxScale = (argc >= 3) ? std::atof(argv[2]) : 3.0;
    double meshDensity = (argc >= 4) ? std::atof(argv[3]) : 0.2;
    std::string tetgenParams = (argc >= 5) ? argv[4] : "pq1.414a0.1";
    
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Input CAD file:   " << inputFile << std::endl;
    std::cout << "  Bounding box:     " << boxScale << "x CAD size" << std::endl;
    std::cout << "  Surface density:  " << meshDensity << std::endl;
    std::cout << "  TetGen params:    " << tetgenParams << std::endl << std::endl;
    
    // ========================================
    // Step 1: Read CAD file
    // ========================================
    std::cout << "Step 1: Reading CAD file..." << std::endl;
    
    TopoDS_Shape cadGeometry;
    std::string extension = inputFile.substr(inputFile.find_last_of(".") + 1);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    
    if (extension == "step" || extension == "stp") {
        STEPControl_Reader reader;
        IFSelect_ReturnStatus status = reader.ReadFile(inputFile.c_str());
        if (status != IFSelect_RetDone) {
            std::cerr << "ERROR: Could not read STEP file" << std::endl;
            return 1;
        }
        reader.TransferRoots();
        cadGeometry = reader.OneShape();
    } else if (extension == "iges" || extension == "igs") {
        IGESControl_Reader reader;
        int status = reader.ReadFile(inputFile.c_str());
        if (status != IFSelect_RetDone) {
            std::cerr << "ERROR: Could not read IGES file" << std::endl;
            return 1;
        }
        reader.TransferRoots();
        cadGeometry = reader.OneShape();
    } else {
        std::cerr << "ERROR: Unsupported file format. Use .step, .stp, .iges, or .igs" << std::endl;
        return 1;
    }
    
    if (cadGeometry.IsNull()) {
        std::cerr << "ERROR: No geometry found in CAD file" << std::endl;
        return 1;
    }
    
    std::cout << "  ✓ CAD geometry loaded successfully" << std::endl;
    
    // ========================================
    // Step 2: Compute bounding box of CAD
    // ========================================
    std::cout << "\nStep 2: Computing bounding box..." << std::endl;
    
    Bnd_Box cadBBox;
    BRepBndLib::Add(cadGeometry, cadBBox);
    
    if (cadBBox.IsVoid()) {
        std::cerr << "ERROR: Could not compute bounding box" << std::endl;
        return 1;
    }
    
    double xmin, ymin, zmin, xmax, ymax, zmax;
    cadBBox.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    
    // Calculate CAD dimensions
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    double dz = zmax - zmin;
    double cadCenter_x = (xmin + xmax) / 2.0;
    double cadCenter_y = (ymin + ymax) / 2.0;
    double cadCenter_z = (zmin + zmax) / 2.0;
    
    std::cout << "  CAD bounding box:" << std::endl;
    std::cout << "    X: [" << xmin << ", " << xmax << "] (width: " << dx << ")" << std::endl;
    std::cout << "    Y: [" << ymin << ", " << ymax << "] (height: " << dy << ")" << std::endl;
    std::cout << "    Z: [" << zmin << ", " << zmax << "] (depth: " << dz << ")" << std::endl;
    std::cout << "    Center: (" << cadCenter_x << ", " << cadCenter_y << ", " << cadCenter_z << ")" << std::endl;
    
    // ========================================
    // Step 3: Create outer bounding box
    // ========================================
    std::cout << "\nStep 3: Creating outer bounding box..." << std::endl;
    
    // Expand box around CAD
    double boxHalfWidth_x = (dx / 2.0) * boxScale;
    double boxHalfWidth_y = (dy / 2.0) * boxScale;
    double boxHalfWidth_z = (dz / 2.0) * boxScale;
    
    gp_Pnt boxCorner1(cadCenter_x - boxHalfWidth_x,
                      cadCenter_y - boxHalfWidth_y,
                      cadCenter_z - boxHalfWidth_z);
    gp_Pnt boxCorner2(cadCenter_x + boxHalfWidth_x,
                      cadCenter_y + boxHalfWidth_y,
                      cadCenter_z + boxHalfWidth_z);
    
    TopoDS_Shape outerBox = BRepPrimAPI_MakeBox(boxCorner1, boxCorner2).Shape();
    
    std::cout << "  Outer box dimensions:" << std::endl;
    std::cout << "    Width:  " << 2 * boxHalfWidth_x << " (CAD × " << boxScale << ")" << std::endl;
    std::cout << "    Height: " << 2 * boxHalfWidth_y << " (CAD × " << boxScale << ")" << std::endl;
    std::cout << "    Depth:  " << 2 * boxHalfWidth_z << " (CAD × " << boxScale << ")" << std::endl;
    
    // ========================================
    // Step 4: Subtract CAD from box (Boolean operation)
    // ========================================
    std::cout << "\nStep 4: Creating mesh domain (box - CAD)..." << std::endl;
    
    TopoDS_Shape meshDomain = BRepAlgoAPI_Cut(outerBox, cadGeometry).Shape();
    
    if (meshDomain.IsNull()) {
        std::cerr << "ERROR: Boolean operation failed" << std::endl;
        return 1;
    }
    
    std::cout << "  ✓ Mesh domain created (volume around CAD)" << std::endl;
    
    // ========================================
    // Step 5: Find hole point (inside CAD)
    // ========================================
    std::cout << "\nStep 5: Finding hole point inside CAD..." << std::endl;
    
    // Use CAD centroid as hole point
    gp_Pnt holePoint;
    
    // Try to get mass properties for accurate centroid
    try {
        GProp_GProps props;
        BRepGProp::VolumeProperties(cadGeometry, props);
        holePoint = props.CentreOfMass();
        std::cout << "  Using CAD centroid: (" << holePoint.X() << ", "
                  << holePoint.Y() << ", " << holePoint.Z() << ")" << std::endl;
    } catch (...) {
        // Fallback to bounding box center
        holePoint = gp_Pnt(cadCenter_x, cadCenter_y, cadCenter_z);
        std::cout << "  Using bounding box center: (" << holePoint.X() << ", "
                  << holePoint.Y() << ", " << holePoint.Z() << ")" << std::endl;
    }
    
    // ========================================
    // Step 6: Generate surface mesh
    // ========================================
    std::cout << "\nStep 6: Generating surface triangulation..." << std::endl;
    
    BRepMesh_IncrementalMesh surfaceMesh(meshDomain, meshDensity, Standard_False, 0.5);
    
    // Extract surface vertices and triangles
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::map<std::array<double, 3>, int> vertexMap;
    
    for (TopExp_Explorer explorer(meshDomain, TopAbs_FACE); explorer.More(); explorer.Next()) {
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
                
                triangles.push_back({getVertexIndex(p1), getVertexIndex(p2), getVertexIndex(p3)});
            }
        }
    }
    
    std::cout << "  Extracted " << vertices.size() << " surface vertices" << std::endl;
    std::cout << "  Extracted " << triangles.size() << " surface triangles" << std::endl;
    
    // ========================================
    // Step 7: Prepare TetGen input
    // ========================================
    std::cout << "\nStep 7: Preparing TetGen input..." << std::endl;
    
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
    
    // CRITICAL: Set hole point (tells TetGen to NOT mesh inside CAD)
    in.numberofholes = 1;
    in.holelist = new REAL[3];
    in.holelist[0] = holePoint.X();
    in.holelist[1] = holePoint.Y();
    in.holelist[2] = holePoint.Z();
    
    std::cout << "  ✓ TetGen input prepared" << std::endl;
    std::cout << "  ✓ Hole point set (mesh will avoid CAD interior)" << std::endl;
    
    // ========================================
    // Step 8: Run TetGen
    // ========================================
    std::cout << "\nStep 8: Running TetGen..." << std::endl;
    std::cout << "  Parameters: " << tetgenParams << std::endl;
    
    try {
        std::string params = tetgenParams;
        if (params.find('Q') == std::string::npos && params.find('V') == std::string::npos) {
            params += "Q"; // Add quiet mode
        }
        
        tetrahedralize((char*)params.c_str(), &in, &out);
        
        std::cout << "  ✓ TetGen completed successfully!" << std::endl;
        std::cout << "  Generated " << out.numberofpoints << " volume vertices" << std::endl;
        std::cout << "  Generated " << out.numberoftetrahedra << " tetrahedra" << std::endl;
        
    } catch (int e) {
        std::cerr << "ERROR: TetGen failed with error code " << e << std::endl;
        return 1;
    }
    
    // ========================================
    // Step 9: Export to VTK
    // ========================================
    std::cout << "\nStep 9: Exporting results..." << std::endl;
    
    std::string baseName = inputFile.substr(0, inputFile.find_last_of("."));
    std::string vtkFile = baseName + "_volume_mesh.vtk";
    
    std::ofstream vtk(vtkFile);
    vtk << "# vtk DataFile Version 3.0\n";
    vtk << "Tetrahedral mesh around " << inputFile << "\n";
    vtk << "ASCII\n";
    vtk << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write vertices
    vtk << "POINTS " << out.numberofpoints << " float\n";
    for (int i = 0; i < out.numberofpoints; i++) {
        vtk << out.pointlist[i * 3 + 0] << " "
            << out.pointlist[i * 3 + 1] << " "
            << out.pointlist[i * 3 + 2] << "\n";
    }
    
    // Write tetrahedra
    vtk << "\nCELLS " << out.numberoftetrahedra << " "
        << out.numberoftetrahedra * 5 << "\n";
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        vtk << "4 "
            << out.tetrahedronlist[i * 4 + 0] << " "
            << out.tetrahedronlist[i * 4 + 1] << " "
            << out.tetrahedronlist[i * 4 + 2] << " "
            << out.tetrahedronlist[i * 4 + 3] << "\n";
    }
    
    // Write cell types
    vtk << "\nCELL_TYPES " << out.numberoftetrahedra << "\n";
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        vtk << "10\n";  // VTK_TETRA
    }
    
    vtk.close();
    
    // Export surface STL
    std::string stlFile = baseName + "_surface.stl";
    StlAPI_Writer stlWriter;
    stlWriter.Write(meshDomain, stlFile.c_str());
    
    std::cout << "  ✓ " << vtkFile << std::endl;
    std::cout << "  ✓ " << stlFile << std::endl;
    
    // ========================================
    // Summary
    // ========================================
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "SUCCESS! Mesh generation complete" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    std::cout << "\nInput:" << std::endl;
    std::cout << "  CAD file:         " << inputFile << std::endl;
    std::cout << "  CAD dimensions:   " << dx << " × " << dy << " × " << dz << std::endl;
    
    std::cout << "\nMesh Domain:" << std::endl;
    std::cout << "  Outer box:        " << 2*boxHalfWidth_x << " × "
              << 2*boxHalfWidth_y << " × " << 2*boxHalfWidth_z << std::endl;
    std::cout << "  Volume meshes:    AROUND CAD (external domain)" << std::endl;
    
    std::cout << "\nMesh Statistics:" << std::endl;
    std::cout << "  Surface vertices: " << vertices.size() << std::endl;
    std::cout << "  Surface triangles:" << triangles.size() << std::endl;
    std::cout << "  Volume vertices:  " << out.numberofpoints << std::endl;
    std::cout << "  Tetrahedra:       " << out.numberoftetrahedra << std::endl;
    
    std::cout << "\nOutput Files:" << std::endl;
    std::cout << "  " << vtkFile << " - Volume tetrahedral mesh" << std::endl;
    std::cout << "  " << stlFile << " - Surface mesh" << std::endl;
    
    std::cout << "\nVisualization:" << std::endl;
    std::cout << "  paraview " << vtkFile << std::endl;
    std::cout << "\nRecommended ParaView filters:" << std::endl;
    std::cout << "  • Extract Surface - See outer/inner boundaries" << std::endl;
    std::cout << "  • Clip - Cut through to see tetrahedra" << std::endl;
    std::cout << "  • Slice - View cross-sections" << std::endl;
    std::cout << std::endl;
    
    return 0;
}
