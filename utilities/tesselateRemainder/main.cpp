#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include "tetgen.h"

struct Point3D {
    double x, y, z;
    Point3D(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}
};

struct Triangle {
    int v1, v2, v3;
    Triangle(int a, int b, int c) : v1(a), v2(b), v3(c) {}
};

bool readVTKFile(const std::string& filename, std::vector<Point3D>& points, std::vector<Triangle>& triangles) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }
    
    std::string line;
    int numPoints = 0;
    bool readingPoints = false;
    bool readingTriangles = false;
    int pointsRead = 0;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string firstWord;
        iss >> firstWord;
        
        // Parse POINTS line
        if (firstWord == "POINTS") {
            iss >> numPoints;
            points.reserve(numPoints);
            readingPoints = true;
            continue;
        }
        
        // Parse POLYGONS line
        if (firstWord == "POLYGONS") {
            readingPoints = false;
            readingTriangles = true;
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
        
        // Read triangle connectivity
        if (readingTriangles && firstWord == "3") {
            int v1, v2, v3;
            iss >> v1 >> v2 >> v3;
            triangles.push_back(Triangle(v1, v2, v3));
        }
    }
    
    file.close();
    std::cout << "Read " << points.size() << " points and " << triangles.size() << " triangles from VTK file" << std::endl;
    return true;
}

void createBoxAndCombineMesh(const std::vector<Point3D>& hemiPoints,
                              const std::vector<Triangle>& hemiTriangles,
                              tetgenio& in) {
    
    // Box dimensions
    double xmin = 0.0, xmax = 5.0;
    double ymin = 0.0, ymax = 5.0;
    double zmin = 0.0, zmax = 5.0;
    
    // First, add all hemisphere points
    int hemiPointCount = hemiPoints.size();
    
    // Define box vertices (8 corners)
    // We need to add box vertices that are NOT on x=0 plane
    std::vector<Point3D> boxVertices;
    
    // Back face (x = xmax)
    boxVertices.push_back(Point3D(xmax, ymin, zmin));  // 0
    boxVertices.push_back(Point3D(xmax, ymax, zmin));  // 1
    boxVertices.push_back(Point3D(xmax, ymax, zmax));  // 2
    boxVertices.push_back(Point3D(xmax, ymin, zmax));  // 3
    
    // Total points
    int totalPoints = hemiPointCount + boxVertices.size();
    
    in.numberofpoints = totalPoints;
    in.pointlist = new REAL[totalPoints * 3];
    
    // Copy hemisphere points
    for (int i = 0; i < hemiPointCount; i++) {
        in.pointlist[i * 3 + 0] = hemiPoints[i].x;
        in.pointlist[i * 3 + 1] = hemiPoints[i].y;
        in.pointlist[i * 3 + 2] = hemiPoints[i].z;
    }
    
    // Copy box vertices
    for (size_t i = 0; i < boxVertices.size(); i++) {
        int idx = hemiPointCount + i;
        in.pointlist[idx * 3 + 0] = boxVertices[i].x;
        in.pointlist[idx * 3 + 1] = boxVertices[i].y;
        in.pointlist[idx * 3 + 2] = boxVertices[i].z;
    }
    
    // Create facets (surfaces)
    // Facets: hemisphere triangles + box faces
    
    std::vector<std::vector<int>> allFacets;
    
    // Add hemisphere triangles
    for (const auto& tri : hemiTriangles) {
        allFacets.push_back({tri.v1, tri.v2, tri.v3});
    }
    
    // Find hemisphere boundary points on x=0 plane
    std::vector<int> boundaryPoints;
    for (int i = 0; i < hemiPointCount; i++) {
        if (std::abs(hemiPoints[i].x) < 1e-6) {
            boundaryPoints.push_back(i);
        }
    }
    
    std::cout << "Found " << boundaryPoints.size() << " boundary points on x=0 plane" << std::endl;
    
    // Box faces (connecting to hemisphere boundary)
    int b0 = hemiPointCount + 0;  // Back face vertices
    int b1 = hemiPointCount + 1;
    int b2 = hemiPointCount + 2;
    int b3 = hemiPointCount + 3;
    
    // Back face (x = xmax): 2 triangles
    allFacets.push_back({b0, b1, b2});
    allFacets.push_back({b0, b2, b3});
    
    // For the side faces, we need to connect the hemisphere boundary to the box edges
    // This is complex - we'll create simplified connecting faces
    
    // Bottom face (z = zmin) - simplified
    // Connect hemisphere boundary points at z~zmin to box bottom edge
    std::vector<int> bottomBoundary;
    for (int idx : boundaryPoints) {
        if (hemiPoints[idx].z < -0.3) {  // Bottom hemisphere points
            bottomBoundary.push_back(idx);
        }
    }
    
    // Create triangulation for bottom face
    if (bottomBoundary.size() > 0) {
        // Simple fan triangulation from first point
        for (size_t i = 1; i < bottomBoundary.size() - 1; i++) {
            allFacets.push_back({bottomBoundary[0], bottomBoundary[i], bottomBoundary[i+1]});
        }
        // Connect to box
        if (bottomBoundary.size() > 0) {
            allFacets.push_back({bottomBoundary[bottomBoundary.size()-1], b0, bottomBoundary[0]});
            allFacets.push_back({b0, b1, bottomBoundary[bottomBoundary.size()-1]});
        }
    }
    
    // Top face (z = zmax) - similar approach
    std::vector<int> topBoundary;
    for (int idx : boundaryPoints) {
        if (hemiPoints[idx].z > 0.3) {
            topBoundary.push_back(idx);
        }
    }
    
    if (topBoundary.size() > 0) {
        for (size_t i = 1; i < topBoundary.size() - 1; i++) {
            allFacets.push_back({topBoundary[0], topBoundary[i+1], topBoundary[i]});
        }
        if (topBoundary.size() > 0) {
            allFacets.push_back({topBoundary[0], b3, topBoundary[topBoundary.size()-1]});
            allFacets.push_back({b3, b2, topBoundary[0]});
        }
    }
    
    // Front (y = ymin) and back (y = ymax) faces - simplified
    // These connect hemisphere boundary to box edges
    
    // Set up facets in TetGen format
    in.numberoffacets = allFacets.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    
    for (int i = 0; i < in.numberoffacets; i++) {
        tetgenio::facet* f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[1];
        f->numberofholes = 0;
        f->holelist = NULL;
        
        tetgenio::polygon* p = &f->polygonlist[0];
        p->numberofvertices = allFacets[i].size();
        p->vertexlist = new int[p->numberofvertices];
        
        for (int j = 0; j < p->numberofvertices; j++) {
            p->vertexlist[j] = allFacets[i][j];
        }
        
        in.facetmarkerlist[i] = i < (int)hemiTriangles.size() ? 1 : 2;
    }
    
    std::cout << "Created " << in.numberoffacets << " facets" << std::endl;
}

void writeTetgenOutput(tetgenio& out, const std::string& prefix) {
    // Write nodes
    std::string nodeFile = prefix + ".node";
    std::ofstream nf(nodeFile);
    nf << out.numberofpoints << " 3 0 0\n";
    for (int i = 0; i < out.numberofpoints; i++) {
        nf << i << " "
           << out.pointlist[i*3] << " "
           << out.pointlist[i*3+1] << " "
           << out.pointlist[i*3+2] << "\n";
    }
    nf.close();
    std::cout << "Written " << nodeFile << std::endl;
    
    // Write elements (tetrahedra)
    std::string eleFile = prefix + ".ele";
    std::ofstream ef(eleFile);
    ef << out.numberoftetrahedra << " 4 0\n";
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        ef << i << " "
           << out.tetrahedronlist[i*4] << " "
           << out.tetrahedronlist[i*4+1] << " "
           << out.tetrahedronlist[i*4+2] << " "
           << out.tetrahedronlist[i*4+3] << "\n";
    }
    ef.close();
    std::cout << "Written " << eleFile << std::endl;
    
    // Write faces
    std::string faceFile = prefix + ".face";
    std::ofstream ff(faceFile);
    ff << out.numberoftrifaces << " 0\n";
    for (int i = 0; i < out.numberoftrifaces; i++) {
        ff << i << " "
           << out.trifacelist[i*3] << " "
           << out.trifacelist[i*3+1] << " "
           << out.trifacelist[i*3+2] << "\n";
    }
    ff.close();
    std::cout << "Written " << faceFile << std::endl;
}

int main(int argc, char** argv) {
    std::string vtkFile = "tess.vtk";
    if (argc > 1) {
        vtkFile = argv[1];
    }
    
    std::cout << "=== TetGen Mesh Generator ===" << std::endl;
    std::cout << "Reading hemisphere surface from: " << vtkFile << std::endl;
    
    // Read VTK file
    std::vector<Point3D> hemiPoints;
    std::vector<Triangle> hemiTriangles;
    
    if (!readVTKFile(vtkFile, hemiPoints, hemiTriangles)) {
        return 1;
    }
    
    // Create TetGen input
    tetgenio in, out;
    
    std::cout << "\nCreating combined geometry (hemisphere + box)..." << std::endl;
    createBoxAndCombineMesh(hemiPoints, hemiTriangles, in);
    
    // Run TetGen
    std::cout << "\nRunning TetGen..." << std::endl;
    try {
        // TetGen switches:
        // p - Tetrahedralize a piecewise linear complex (PLC)
        // q - Quality mesh generation
        // a - Maximum tetrahedron volume constraint
        // A - Assign attributes to regions
        // V - Verbose output
        tetrahedralize("pq1.2a0.1V", &in, &out);
        
        std::cout << "\n=== TetGen Results ===" << std::endl;
        std::cout << "Generated " << out.numberofpoints << " nodes" << std::endl;
        std::cout << "Generated " << out.numberoftetrahedra << " tetrahedra" << std::endl;
        std::cout << "Generated " << out.numberoftrifaces << " boundary faces" << std::endl;
        
        // Write output
        std::cout << "\nWriting output files..." << std::endl;
        writeTetgenOutput(out, "output_mesh");
        
        std::cout << "\nMesh generation completed successfully!" << std::endl;
        std::cout << "Output files: output_mesh.node, output_mesh.ele, output_mesh.face" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during tetrahedralization: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
