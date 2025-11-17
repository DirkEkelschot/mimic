#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include "tetgen.h"

struct Point3D {
    double x, y, z;
    Point3D(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}
};

struct Triangle {
    int v1, v2, v3;
    Triangle(int a, int b, int c) : v1(a), v2(b), v3(c) {}
};

struct Tetrahedron {
    int v1, v2, v3, v4;
    Tetrahedron(int a, int b, int c, int d) : v1(a), v2(b), v3(c), v4(d) {}
};

// Read VTK file with boundary faces
bool readVTKFile(const std::string& filename, std::vector<Point3D>& points, std::vector<Triangle>& triangles) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }
    
    std::string line;
    int numPoints = 0;
    int numCells = 0;
    bool readingPoints = false;
    bool readingCells = false;
    int pointsRead = 0;
    int cellsRead = 0;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string firstWord;
        iss >> firstWord;
        
        // Parse POINTS line
        if (firstWord == "POINTS") {
            iss >> numPoints;
            points.reserve(numPoints);
            readingPoints = true;
            readingCells = false;
            continue;
        }
        
        // Parse CELLS line
        if (firstWord == "CELLS") {
            iss >> numCells;
            readingPoints = false;
            readingCells = true;
            continue;
        }
        
        // Parse CELL_TYPES line - skip it
        if (firstWord == "CELL_TYPES") {
            readingCells = false;
            continue;
        }
        
        // Read point coordinates
        if (readingPoints && pointsRead < numPoints) {
            double x, y, z;
            std::istringstream pointStream(line);
            if (pointStream >> x >> y >> z) {
                points.push_back(Point3D(x, y, z));
                pointsRead++;
            }
        }
        
        // Read cell connectivity (triangles have format: 3 v1 v2 v3)
        if (readingCells && cellsRead < numCells) {
            int numVerts;
            std::istringstream cellStream(line);
            if (cellStream >> numVerts) {
                if (numVerts == 3) {  // Triangle
                    int v1, v2, v3;
                    if (cellStream >> v1 >> v2 >> v3) {
                        triangles.push_back(Triangle(v1, v2, v3));
                        cellsRead++;
                    }
                } else {
                    // Skip non-triangle cells
                    cellsRead++;
                }
            }
        }
    }
    
    file.close();
    std::cout << "Read " << points.size() << " points and " << triangles.size() << " triangles from VTK file" << std::endl;
    return true;
}

// Compute bounding box and return center and maximum diameter
void computeBoundingBox(const std::vector<Point3D>& points, 
                        double& xmin, double& xmax,
                        double& ymin, double& ymax,
                        double& zmin, double& zmax,
                        double& maxDiameter) {
    if (points.empty()) {
        std::cerr << "Error: No points to compute bounding box" << std::endl;
        return;
    }
    
    xmin = xmax = points[0].x;
    ymin = ymax = points[0].y;
    zmin = zmax = points[0].z;
    
    for (const auto& p : points) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
        zmin = std::min(zmin, p.z);
        zmax = std::max(zmax, p.z);
    }
    
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    double dz = zmax - zmin;
    maxDiameter = std::max({dx, dy, dz});
    
    std::cout << "Bounding box: [" << xmin << ", " << xmax << "] x [" 
              << ymin << ", " << ymax << "] x [" << zmin << ", " << zmax << "]" << std::endl;
    std::cout << "Maximum diameter: " << maxDiameter << std::endl;
}

// Create box vertices and faces for TetGen
void createBoxFaces(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                    std::vector<Point3D>& boxPoints, std::vector<std::vector<int> >& boxFaces) {
    // Create 8 box vertices
    boxPoints.clear();
    boxPoints.push_back(Point3D(xmin, ymin, zmin));  // 0
    boxPoints.push_back(Point3D(xmax, ymin, zmin));  // 1
    boxPoints.push_back(Point3D(xmax, ymax, zmin));  // 2
    boxPoints.push_back(Point3D(xmin, ymax, zmin));  // 3
    boxPoints.push_back(Point3D(xmin, ymin, zmax));  // 4
    boxPoints.push_back(Point3D(xmax, ymin, zmax));  // 5
    boxPoints.push_back(Point3D(xmax, ymax, zmax));  // 6
    boxPoints.push_back(Point3D(xmin, ymax, zmax));  // 7
    
    // Create 6 faces (each face has 2 triangles)
    boxFaces.clear();
    
    // Bottom face (z = zmin): 0-1-2, 0-2-3
    boxFaces.push_back({0, 1, 2});
    boxFaces.push_back({0, 2, 3});
    
    // Top face (z = zmax): 4-7-6, 4-6-5
    boxFaces.push_back({4, 7, 6});
    boxFaces.push_back({4, 6, 5});
    
    // Front face (y = ymin): 0-4-5, 0-5-1
    boxFaces.push_back({0, 4, 5});
    boxFaces.push_back({0, 5, 1});
    
    // Back face (y = ymax): 3-2-6, 3-6-7
    boxFaces.push_back({3, 2, 6});
    boxFaces.push_back({3, 6, 7});
    
    // Left face (x = xmin): 0-3-7, 0-7-4
    boxFaces.push_back({0, 3, 7});
    boxFaces.push_back({0, 7, 4});
    
    // Right face (x = xmax): 1-5-6, 1-6-2
    boxFaces.push_back({1, 5, 6});
    boxFaces.push_back({1, 6, 2});
}

// Create TetGen input from boundary faces and box
void createTetGenInput(const std::vector<Point3D>& boundaryPoints,
                       const std::vector<Triangle>& boundaryTriangles,
                       const std::vector<Point3D>& boxPoints,
                       const std::vector<std::vector<int> >& boxFaces,
                       tetgenio& in) {
    
    int numBoundaryPoints = boundaryPoints.size();
    int numBoxPoints = boxPoints.size();
    int totalPoints = numBoundaryPoints + numBoxPoints;
    
    // Allocate point list
    in.numberofpoints = totalPoints;
    in.pointlist = new REAL[totalPoints * 3];
    
    // Copy boundary points
    for (int i = 0; i < numBoundaryPoints; i++) {
        in.pointlist[i * 3 + 0] = boundaryPoints[i].x;
        in.pointlist[i * 3 + 1] = boundaryPoints[i].y;
        in.pointlist[i * 3 + 2] = boundaryPoints[i].z;
    }
    
    // Copy box points
    for (int i = 0; i < numBoxPoints; i++) {
        int idx = numBoundaryPoints + i;
        in.pointlist[idx * 3 + 0] = boxPoints[i].x;
        in.pointlist[idx * 3 + 1] = boxPoints[i].y;
        in.pointlist[idx * 3 + 2] = boxPoints[i].z;
    }
    
    // Create facets: boundary triangles + box faces
    int numBoundaryFacets = boundaryTriangles.size();
    int numBoxFacets = boxFaces.size();
    int totalFacets = numBoundaryFacets + numBoxFacets;
    
    in.numberoffacets = totalFacets;
    in.facetlist = new tetgenio::facet[totalFacets];
    in.facetmarkerlist = new int[totalFacets];
    
    // Add boundary triangles (marker = 1)
    for (int i = 0; i < numBoundaryFacets; i++) {
        tetgenio::facet* f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[1];
        f->numberofholes = 0;
        f->holelist = NULL;
        
        tetgenio::polygon* p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[3];
        p->vertexlist[0] = boundaryTriangles[i].v1;
        p->vertexlist[1] = boundaryTriangles[i].v2;
        p->vertexlist[2] = boundaryTriangles[i].v3;
        
        in.facetmarkerlist[i] = 1;  // Boundary marker
    }
    
    // Add box faces (marker = 2)
    for (int i = 0; i < numBoxFacets; i++) {
        int facetIdx = numBoundaryFacets + i;
        tetgenio::facet* f = &in.facetlist[facetIdx];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[1];
        f->numberofholes = 0;
        f->holelist = NULL;
        
        tetgenio::polygon* p = &f->polygonlist[0];
        p->numberofvertices = boxFaces[i].size();
        p->vertexlist = new int[p->numberofvertices];
        for (size_t j = 0; j < boxFaces[i].size(); j++) {
            // Offset by number of boundary points
            p->vertexlist[j] = numBoundaryPoints + boxFaces[i][j];
        }
        
        in.facetmarkerlist[facetIdx] = 2;  // Box marker
    }
    
    std::cout << "Created TetGen input: " << totalPoints << " points, " 
              << totalFacets << " facets (" << numBoundaryFacets 
              << " boundary + " << numBoxFacets << " box)" << std::endl;
}

// Read VTK file with tetrahedra (volume mesh)
bool readVTKVolumeMesh(const std::string& filename, 
                       std::vector<Point3D>& points, 
                       std::vector<Tetrahedron>& tetrahedra) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }
    
    std::string line;
    int numPoints = 0;
    int numCells = 0;
    bool readingPoints = false;
    bool readingCells = false;
    int pointsRead = 0;
    int cellsRead = 0;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string firstWord;
        iss >> firstWord;
        
        // Parse POINTS line
        if (firstWord == "POINTS") {
            iss >> numPoints;
            points.reserve(numPoints);
            readingPoints = true;
            readingCells = false;
            continue;
        }
        
        // Parse CELLS line
        if (firstWord == "CELLS") {
            iss >> numCells;
            readingPoints = false;
            readingCells = true;
            continue;
        }
        
        // Parse CELL_TYPES line - skip it
        if (firstWord == "CELL_TYPES") {
            readingCells = false;
            continue;
        }
        
        // Read point coordinates
        if (readingPoints && pointsRead < numPoints) {
            double x, y, z;
            std::istringstream pointStream(line);
            if (pointStream >> x >> y >> z) {
                points.push_back(Point3D(x, y, z));
                pointsRead++;
            }
        }
        
        // Read cell connectivity (tetrahedra have format: 4 v1 v2 v3 v4)
        if (readingCells && cellsRead < numCells) {
            int numVerts;
            std::istringstream cellStream(line);
            if (cellStream >> numVerts) {
                if (numVerts == 4) {  // Tetrahedron
                    int v1, v2, v3, v4;
                    if (cellStream >> v1 >> v2 >> v3 >> v4) {
                        tetrahedra.push_back(Tetrahedron(v1, v2, v3, v4));
                        cellsRead++;
                    }
                } else {
                    // Skip non-tetrahedron cells
                    cellsRead++;
                }
            }
        }
    }
    
    file.close();
    std::cout << "Read " << points.size() << " points and " << tetrahedra.size() << " tetrahedra from VTK file" << std::endl;
    return true;
}

// Compute distance between two points
double pointDistance(const Point3D& p1, const Point3D& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// Merge two meshes by combining points and remapping cell connectivity
// Points within tolerance are considered duplicates
void mergeMeshes(const std::vector<Point3D>& points1,
                 const std::vector<Tetrahedron>& cells1,
                 const std::vector<Point3D>& points2,
                 const std::vector<Tetrahedron>& cells2,
                 double tolerance,
                 std::vector<Point3D>& mergedPoints,
                 std::vector<Tetrahedron>& mergedCells,
                 std::vector<int>& pointMap1,
                 std::vector<int>& pointMap2) {
    
    mergedPoints.clear();
    mergedCells.clear();
    pointMap1.resize(points1.size(), -1);
    pointMap2.resize(points2.size(), -1);
    
    // First, add all points from mesh1
    for (size_t i = 0; i < points1.size(); i++) {
        pointMap1[i] = mergedPoints.size();
        mergedPoints.push_back(points1[i]);
    }
    
    // Add points from mesh2, checking for duplicates
    for (size_t i = 0; i < points2.size(); i++) {
        bool found = false;
        // Check against all existing points (from both meshes)
        for (size_t j = 0; j < mergedPoints.size(); j++) {
            if (pointDistance(points2[i], mergedPoints[j]) < tolerance) {
                pointMap2[i] = j;
                found = true;
                break;
            }
        }
        if (!found) {
            pointMap2[i] = mergedPoints.size();
            mergedPoints.push_back(points2[i]);
        }
    }
    
    // Add cells from mesh1 (remap indices)
    for (const auto& tet : cells1) {
        mergedCells.push_back(Tetrahedron(
            pointMap1[tet.v1],
            pointMap1[tet.v2],
            pointMap1[tet.v3],
            pointMap1[tet.v4]
        ));
    }
    
    // Add cells from mesh2 (remap indices)
    for (const auto& tet : cells2) {
        mergedCells.push_back(Tetrahedron(
            pointMap2[tet.v1],
            pointMap2[tet.v2],
            pointMap2[tet.v3],
            pointMap2[tet.v4]
        ));
    }
    
    std::cout << "Merged meshes: " << mergedPoints.size() << " unique points, " 
              << mergedCells.size() << " tetrahedra" << std::endl;
    std::cout << "  Mesh 1: " << points1.size() << " points, " << cells1.size() << " cells" << std::endl;
    std::cout << "  Mesh 2: " << points2.size() << " points, " << cells2.size() << " cells" << std::endl;
    std::cout << "  Removed " << (points1.size() + points2.size() - mergedPoints.size()) 
              << " duplicate points" << std::endl;
}

// Write merged mesh to VTK format
void writeMergedVTK(const std::vector<Point3D>& points,
                    const std::vector<Tetrahedron>& tetrahedra,
                    const std::string& filename) {
    std::ofstream vtkFile(filename);
    if (!vtkFile.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << " for writing" << std::endl;
        return;
    }
    
    // Write VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Merged volume mesh\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write points
    vtkFile << "POINTS " << points.size() << " float\n";
    for (const auto& p : points) {
        vtkFile << p.x << " " << p.y << " " << p.z << "\n";
    }
    
    // Write cells (tetrahedra)
    int totalSize = tetrahedra.size() * 5;  // 4 vertices + 1 count
    vtkFile << "\nCELLS " << tetrahedra.size() << " " << totalSize << "\n";
    for (const auto& tet : tetrahedra) {
        vtkFile << "4 "  // 4 vertices for tetrahedron
                << tet.v1 << " "
                << tet.v2 << " "
                << tet.v3 << " "
                << tet.v4 << "\n";
    }
    
    // Write cell types (10 = VTK_TETRA)
    vtkFile << "\nCELL_TYPES " << tetrahedra.size() << "\n";
    for (size_t i = 0; i < tetrahedra.size(); i++) {
        vtkFile << "10\n";  // VTK_TETRA = 10
    }
    
    vtkFile.close();
    std::cout << "Written merged mesh to " << filename << std::endl;
}

// Write TetGen output to VTK format
void writeVTKOutput(tetgenio& out, const std::string& filename) {
    std::ofstream vtkFile(filename);
    if (!vtkFile.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << " for writing" << std::endl;
        return;
    }
    
    // Write VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Volume mesh from TetGen\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write points
    vtkFile << "POINTS " << out.numberofpoints << " float\n";
    for (int i = 0; i < out.numberofpoints; i++) {
        vtkFile << out.pointlist[i*3] << " "
                << out.pointlist[i*3+1] << " "
                << out.pointlist[i*3+2] << "\n";
    }
    
    // Write cells (tetrahedra)
    // VTK format: CELLS numCells totalSize
    // For each cell: numVertices v1 v2 v3 v4
    int totalSize = out.numberoftetrahedra * 5;  // 4 vertices + 1 count
    vtkFile << "\nCELLS " << out.numberoftetrahedra << " " << totalSize << "\n";
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        vtkFile << "4 "  // 4 vertices for tetrahedron
                << out.tetrahedronlist[i*4] << " "
                << out.tetrahedronlist[i*4+1] << " "
                << out.tetrahedronlist[i*4+2] << " "
                << out.tetrahedronlist[i*4+3] << "\n";
    }
    
    // Write cell types (10 = VTK_TETRA)
    vtkFile << "\nCELL_TYPES " << out.numberoftetrahedra << "\n";
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        vtkFile << "10\n";  // VTK_TETRA = 10
    }
    
    vtkFile.close();
    std::cout << "Written VTK output to " << filename << std::endl;
}

int main(int argc, char** argv) {
    // Check if we're in merge mode (4 arguments: volume.vtk bl_mesh.vtk output.vtk [tolerance])
    if (argc >= 4) {
        std::string volumeFile = argv[1];
        std::string blMeshFile = argv[2];
        std::string outputFile = argv[3];
        double tolerance = 1e-6;
        if (argc > 4) {
            tolerance = std::atof(argv[4]);
        }
        
        std::cout << "=== Mesh Merger ===" << std::endl;
        std::cout << "Volume mesh: " << volumeFile << std::endl;
        std::cout << "Boundary layer mesh: " << blMeshFile << std::endl;
        std::cout << "Output merged mesh: " << outputFile << std::endl;
        std::cout << "Merge tolerance: " << tolerance << std::endl;
        std::cout << std::endl;
        
        // Read both meshes
        std::vector<Point3D> volumePoints, blPoints;
        std::vector<Tetrahedron> volumeCells, blCells;
        
        std::cout << "Reading volume mesh..." << std::endl;
        if (!readVTKVolumeMesh(volumeFile, volumePoints, volumeCells)) {
            return 1;
        }
        
        std::cout << "\nReading boundary layer mesh..." << std::endl;
        if (!readVTKVolumeMesh(blMeshFile, blPoints, blCells)) {
            return 1;
        }
        
        // Merge meshes
        std::cout << "\nMerging meshes..." << std::endl;
        std::vector<Point3D> mergedPoints;
        std::vector<Tetrahedron> mergedCells;
        std::vector<int> pointMap1, pointMap2;
        
        mergeMeshes(volumePoints, volumeCells, blPoints, blCells, 
                   tolerance, mergedPoints, mergedCells, pointMap1, pointMap2);
        
        // Write merged mesh
        std::cout << "\nWriting merged mesh..." << std::endl;
        writeMergedVTK(mergedPoints, mergedCells, outputFile);
        
        std::cout << "\nMesh merging completed successfully!" << std::endl;
        return 0;
    }
    
    // Original mode: Generate volume mesh from boundary faces
    if (argc < 3) {
        std::cerr << "Usage (Volume Generation): " << argv[0] << " <boundaryface.vtk> <output.vtk> [quality]" << std::endl;
        std::cerr << "  quality: TetGen quality parameter (default: 1.2)" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Usage (Mesh Merging): " << argv[0] << " <volume.vtk> <bl_mesh.vtk> <output.vtk> [tolerance]" << std::endl;
        std::cerr << "  tolerance: Point merge tolerance (default: 1e-6)" << std::endl;
        return 1;
    }
    
    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    double quality = 1.2;
    if (argc > 3) {
        quality = std::atof(argv[3]);
    }
    
    std::cout << "=== Volume Mesh Generator using TetGen ===" << std::endl;
    std::cout << "Input boundary faces: " << inputFile << std::endl;
    std::cout << "Output mesh: " << outputFile << std::endl;
    std::cout << "Quality parameter: " << quality << std::endl;
    std::cout << std::endl;
    
    // Read boundary faces from VTK file
    std::vector<Point3D> boundaryPoints;
    std::vector<Triangle> boundaryTriangles;
    
    std::cout << "Reading boundary faces..." << std::endl;
    if (!readVTKFile(inputFile, boundaryPoints, boundaryTriangles)) {
        return 1;
    }
    
    if (boundaryPoints.empty() || boundaryTriangles.empty()) {
        std::cerr << "Error: No points or triangles found in input file" << std::endl;
        return 1;
    }
    
    // Compute bounding box
    std::cout << "\nComputing bounding box..." << std::endl;
    double xmin, xmax, ymin, ymax, zmin, zmax, maxDiameter;
    computeBoundingBox(boundaryPoints, xmin, xmax, ymin, ymax, zmin, zmax, maxDiameter);
    
    // Calculate box dimensions (10x the maximum diameter)
    double boxSize = 10.0 * maxDiameter;
    double centerX = (xmin + xmax) / 2.0;
    double centerY = (ymin + ymax) / 2.0;
    double centerZ = (zmin + zmax) / 2.0;
    
    double boxXmin = centerX - boxSize / 2.0;
    double boxXmax = centerX + boxSize / 2.0;
    double boxYmin = centerY - boxSize / 2.0;
    double boxYmax = centerY + boxSize / 2.0;
    double boxZmin = centerZ - boxSize / 2.0;
    double boxZmax = centerZ + boxSize / 2.0;
    
    std::cout << "\nCreating bounding box (10x maximum diameter)..." << std::endl;
    std::cout << "Box size: " << boxSize << std::endl;
    std::cout << "Box bounds: [" << boxXmin << ", " << boxXmax << "] x [" 
              << boxYmin << ", " << boxYmax << "] x [" 
              << boxZmin << ", " << boxZmax << "]" << std::endl;
    
    // Create box faces
    std::vector<Point3D> boxPoints;
    std::vector<std::vector<int> > boxFaces;
    createBoxFaces(boxXmin, boxXmax, boxYmin, boxYmax, boxZmin, boxZmax, 
                   boxPoints, boxFaces);
    
    // Create TetGen input
    std::cout << "\nPreparing TetGen input..." << std::endl;
    tetgenio in, out;
    createTetGenInput(boundaryPoints, boundaryTriangles, boxPoints, boxFaces, in);
    
    // Run TetGen
    std::cout << "\nRunning TetGen..." << std::endl;
    try {
        // TetGen switches:
        // p - Tetrahedralize a piecewise linear complex (PLC)
        // q - Quality mesh generation with specified quality
        // A - Assign attributes to regions
        // V - Verbose output
        std::ostringstream switches;
        switches << "pq" << quality << "AV";
        
        std::cout << "TetGen switches: " << switches.str() << std::endl;
        tetrahedralize(const_cast<char*>(switches.str().c_str()), &in, &out);
        
        std::cout << "\n=== TetGen Results ===" << std::endl;
        std::cout << "Generated " << out.numberofpoints << " nodes" << std::endl;
        std::cout << "Generated " << out.numberoftetrahedra << " tetrahedra" << std::endl;
        std::cout << "Generated " << out.numberoftrifaces << " boundary faces" << std::endl;
        
        // Write output
        std::cout << "\nWriting output..." << std::endl;
        writeVTKOutput(out, outputFile);
        
        std::cout << "\nMesh generation completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during tetrahedralization: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
