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

// Compute bounding box and return center and maximum length
void computeBoundingBox(const std::vector<Point3D>& points, 
                        double& xmin, double& xmax,
                        double& ymin, double& ymax,
                        double& zmin, double& zmax,
                        double& maxLength) {
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
    maxLength = std::max({dx, dy, dz});
    
    std::cout << "Bounding box: [" << xmin << ", " << xmax << "] x [" 
              << ymin << ", " << ymax << "] x [" << zmin << ", " << zmax << "]" << std::endl;
    std::cout << "Maximum length: " << maxLength << std::endl;
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
    // Note: Face normals should point outward from the box
    boxFaces.clear();
    
    // Bottom face (z = zmin): 0-1-2, 0-2-3 (normal pointing down)
    boxFaces.push_back({0, 1, 2});
    boxFaces.push_back({0, 2, 3});
    
    // Top face (z = zmax): 4-7-6, 4-6-5 (normal pointing up)
    boxFaces.push_back({4, 7, 6});
    boxFaces.push_back({4, 6, 5});
    
    // Front face (y = ymin): 0-4-5, 0-5-1 (normal pointing forward)
    boxFaces.push_back({0, 4, 5});
    boxFaces.push_back({0, 5, 1});
    
    // Back face (y = ymax): 3-2-6, 3-6-7 (normal pointing backward)
    boxFaces.push_back({3, 2, 6});
    boxFaces.push_back({3, 6, 7});
    
    // Left face (x = xmin): 0-3-7, 0-7-4 (normal pointing left)
    boxFaces.push_back({0, 3, 7});
    boxFaces.push_back({0, 7, 4});
    
    // Right face (x = xmax): 1-5-6, 1-6-2 (normal pointing right)
    boxFaces.push_back({1, 5, 6});
    boxFaces.push_back({1, 6, 2});
}

// Compute center point of boundary geometry (for hole point)
Point3D computeBoundaryCenter(const std::vector<Point3D>& points) {
    Point3D center(0, 0, 0);
    if (points.empty()) return center;
    
    for (const auto& p : points) {
        center.x += p.x;
        center.y += p.y;
        center.z += p.z;
    }
    
    center.x /= points.size();
    center.y /= points.size();
    center.z /= points.size();
    
    return center;
}

// Create TetGen input from boundary faces and box
// The boundary surface defines a hole inside the box
void createTetGenInput(const std::vector<Point3D>& boundaryPoints,
                       const std::vector<Triangle>& boundaryTriangles,
                       const std::vector<Point3D>& boxPoints,
                       const std::vector<std::vector<int> >& boxFaces,
                       const Point3D& holePoint,
                       tetgenio& in) {
    
    int numBoundaryPoints = boundaryPoints.size();
    int numBoxPoints = boxPoints.size();
    int totalPoints = numBoundaryPoints + numBoxPoints;
    
    // Allocate point list
    in.numberofpoints = totalPoints;
    in.pointlist = new REAL[totalPoints * 3];
    
    // Copy boundary points first (these will be preserved with -Y flag)
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
    
    // Add boundary triangles (marker = 1) - these define the hole boundary
    // Note: The boundary surface should be oriented with normals pointing outward from the hole
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
    
    // Add box faces (marker = 2) - outer boundary
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
    
    // Define regions: 
    // Region 1: Outside the box (should not exist, but defined for safety)
    // Region 2: Inside the box but outside the hole (the domain we want to mesh)
    // Region 3: Inside the hole (should be empty)
    in.numberofregions = 1;
    in.regionlist = new REAL[in.numberofregions * 5];  // x, y, z, attribute, max volume
    
    // Region 2: Point inside the box but outside the hole
    // Use a point near a box corner (guaranteed to be outside the boundary)
    // Use box corner 6 (max x, max y, max z) but move slightly inward
    Point3D boxCorner = boxPoints[6];  // (xmax, ymax, zmax)
    Point3D boxCorner0 = boxPoints[0]; // (xmin, ymin, zmin)
    
    // Place region point near the box corner but slightly inward (99% of the way)
    // This ensures it's definitely inside the box but far from the boundary
    Point3D regionPoint(
        boxCorner0.x + (boxCorner.x - boxCorner0.x) * 0.99,
        boxCorner0.y + (boxCorner.y - boxCorner0.y) * 0.99,
        boxCorner0.z + (boxCorner.z - boxCorner0.z) * 0.99
    );
    
    in.regionlist[0] = regionPoint.x;
    in.regionlist[1] = regionPoint.y;
    in.regionlist[2] = regionPoint.z;
    in.regionlist[3] = 2.0;  // Region attribute (marker)
    in.regionlist[4] = 0.0;  // Max volume constraint (0 = no constraint)
    
    // Define hole point (inside the boundary surface)
    in.numberofholes = 1;
    in.holelist = new REAL[in.numberofholes * 3];
    in.holelist[0] = holePoint.x;
    in.holelist[1] = holePoint.y;
    in.holelist[2] = holePoint.z;
    
    std::cout << "Created TetGen input: " << totalPoints << " points, " 
              << totalFacets << " facets (" << numBoundaryFacets 
              << " boundary + " << numBoxFacets << " box)" << std::endl;
    std::cout << "Hole point: (" << holePoint.x << ", " << holePoint.y << ", " << holePoint.z << ")" << std::endl;
    std::cout << "Region point: (" << regionPoint.x << ", " << regionPoint.y << ", " << regionPoint.z << ")" << std::endl;
}

// Check if a point is inside a closed surface using ray casting
// This is a simple implementation that counts intersections with triangles
bool isPointInsideBoundary(const Point3D& point, 
                           const std::vector<Point3D>& boundaryPoints,
                           const std::vector<Triangle>& boundaryTriangles) {
    // Cast a ray in +X direction and count intersections
    int intersections = 0;
    const double epsilon = 1e-12;
    
    for (const auto& tri : boundaryTriangles) {
        const Point3D& v0 = boundaryPoints[tri.v1];
        const Point3D& v1 = boundaryPoints[tri.v2];
        const Point3D& v2 = boundaryPoints[tri.v3];
        
        // Ray-triangle intersection test (Möller-Trumbore algorithm)
        // Ray: point + t * (1, 0, 0)
        Point3D edge1(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
        Point3D edge2(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);
        Point3D rayDir(1.0, 0.0, 0.0);  // Ray direction (positive X)
        Point3D s(point.x - v0.x, point.y - v0.y, point.z - v0.z);
        
        // Cross product: edge1 x rayDir
        Point3D h(
            edge1.y * rayDir.z - edge1.z * rayDir.y,
            edge1.z * rayDir.x - edge1.x * rayDir.z,
            edge1.x * rayDir.y - edge1.y * rayDir.x
        );
        
        double a = edge2.x * h.x + edge2.y * h.y + edge2.z * h.z;
        
        if (std::abs(a) < epsilon) continue;  // Ray is parallel to triangle
        
        double f = 1.0 / a;
        double u = f * (s.x * h.x + s.y * h.y + s.z * h.z);
        if (u < 0.0 || u > 1.0) continue;
        
        // Cross product: s x edge2
        Point3D q(
            s.y * edge2.z - s.z * edge2.y,
            s.z * edge2.x - s.x * edge2.z,
            s.x * edge2.y - s.y * edge2.x
        );
        
        double v = f * (rayDir.x * q.x + rayDir.y * q.y + rayDir.z * q.z);
        if (v < 0.0 || u + v > 1.0) continue;
        
        double t = f * (edge1.x * q.x + edge1.y * q.y + edge1.z * q.z);
        if (t > epsilon) {  // Intersection found (ray goes forward)
            intersections++;
        }
    }
    
    // Odd number of intersections means point is inside
    return (intersections % 2) == 1;
}

// Compute centroid of a tetrahedron
Point3D computeTetrahedronCentroid(const Point3D& p0, const Point3D& p1, 
                                    const Point3D& p2, const Point3D& p3) {
    return Point3D(
        (p0.x + p1.x + p2.x + p3.x) / 4.0,
        (p0.y + p1.y + p2.y + p3.y) / 4.0,
        (p0.z + p1.z + p2.z + p3.z) / 4.0
    );
}

// Write TetGen output to VTK format, filtering out tetrahedra inside the boundary
void writeVTKOutput(tetgenio& out, const std::string& filename,
                    const std::vector<Point3D>& boundaryPoints,
                    const std::vector<Triangle>& boundaryTriangles) {
    std::ofstream vtkFile(filename);
    if (!vtkFile.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << " for writing" << std::endl;
        return;
    }
    
    // Filter tetrahedra: only keep those in the external domain
    std::vector<int> validTetrahedra;
    const double EXTERNAL_REGION_ATTRIBUTE = 2.0;
    
    std::cout << "Filtering tetrahedra to keep only external domain..." << std::endl;
    
    // Pre-compute tetrahedron centroids for filtering
    std::vector<Point3D> centroids;
    centroids.reserve(out.numberoftetrahedra);
    
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        int v0 = out.tetrahedronlist[i*4];
        int v1 = out.tetrahedronlist[i*4+1];
        int v2 = out.tetrahedronlist[i*4+2];
        int v3 = out.tetrahedronlist[i*4+3];
        
        Point3D p0(out.pointlist[v0*3], out.pointlist[v0*3+1], out.pointlist[v0*3+2]);
        Point3D p1(out.pointlist[v1*3], out.pointlist[v1*3+1], out.pointlist[v1*3+2]);
        Point3D p2(out.pointlist[v2*3], out.pointlist[v2*3+1], out.pointlist[v2*3+2]);
        Point3D p3(out.pointlist[v3*3], out.pointlist[v3*3+1], out.pointlist[v3*3+2]);
        
        centroids.push_back(computeTetrahedronCentroid(p0, p1, p2, p3));
    }
    
    // Use geometric test as primary method - more reliable than region attributes
    int keptByGeometry = 0;
    int keptByAttribute = 0;
    
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        bool keep = false;
        
        // Primary method: geometric test - check if centroid is outside the boundary
        bool isOutside = !isPointInsideBoundary(centroids[i], boundaryPoints, boundaryTriangles);
        
        if (isOutside) {
            keep = true;
            keptByGeometry++;
        }
        
        // Optional: verify with region attributes for diagnostics
        if (out.tetrahedronattributelist != NULL && out.numberoftetrahedronattributes > 0) {
            double regionAttr = out.tetrahedronattributelist[i * out.numberoftetrahedronattributes];
            if (std::abs(regionAttr - EXTERNAL_REGION_ATTRIBUTE) < 1e-6) {
                keptByAttribute++;
            }
        }
        
        if (keep) {
            validTetrahedra.push_back(i);
        }
    }
    
    std::cout << "  Kept by geometry (outside boundary): " << keptByGeometry << std::endl;
    if (out.tetrahedronattributelist != NULL && out.numberoftetrahedronattributes > 0) {
        std::cout << "  Kept by region attribute (2.0): " << keptByAttribute << std::endl;
    }
    
    std::cout << "Kept " << validTetrahedra.size() << " out of " 
              << out.numberoftetrahedra << " tetrahedra (external domain only)" << std::endl;
    
    // Write VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Volume mesh from TetGen (external domain only)\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write points (all points, even if some aren't used)
    vtkFile << "POINTS " << out.numberofpoints << " float\n";
    for (int i = 0; i < out.numberofpoints; i++) {
        vtkFile << out.pointlist[i*3] << " "
                << out.pointlist[i*3+1] << " "
                << out.pointlist[i*3+2] << "\n";
    }
    
    // Write cells (only valid tetrahedra)
    int totalSize = validTetrahedra.size() * 5;  // 4 vertices + 1 count
    vtkFile << "\nCELLS " << validTetrahedra.size() << " " << totalSize << "\n";
    for (int idx : validTetrahedra) {
        vtkFile << "4 "  // 4 vertices for tetrahedron
                << out.tetrahedronlist[idx*4] << " "
                << out.tetrahedronlist[idx*4+1] << " "
                << out.tetrahedronlist[idx*4+2] << " "
                << out.tetrahedronlist[idx*4+3] << "\n";
    }
    
    // Write cell types (10 = VTK_TETRA)
    vtkFile << "\nCELL_TYPES " << validTetrahedra.size() << "\n";
    for (size_t i = 0; i < validTetrahedra.size(); i++) {
        vtkFile << "10\n";  // VTK_TETRA = 10
    }
    
    // Write region attributes if available
    if (out.tetrahedronattributelist != NULL && out.numberoftetrahedronattributes > 0) {
        vtkFile << "\nCELL_DATA " << validTetrahedra.size() << "\n";
        vtkFile << "SCALARS RegionAttribute float 1\n";
        vtkFile << "LOOKUP_TABLE default\n";
        for (int idx : validTetrahedra) {
            vtkFile << out.tetrahedronattributelist[idx * out.numberoftetrahedronattributes] << "\n";
        }
    }
    
    vtkFile.close();
    std::cout << "Written VTK output to " << filename << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <boundary.vtk> <output.vtk> [quality]" << std::endl;
        std::cerr << "  boundary.vtk: Input VTK file with watertight boundary surface" << std::endl;
        std::cerr << "  output.vtk: Output VTK file with volume mesh" << std::endl;
        std::cerr << "  quality: TetGen quality parameter (default: 1.2)" << std::endl;
        std::cerr << std::endl;
        std::cerr << "This program generates a volume mesh between a bounding box" << std::endl;
        std::cerr << "and the boundary surface, creating a hole inside the domain." << std::endl;
        return 1;
    }
    
    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    double quality = 1.2;
    if (argc > 3) {
        quality = std::atof(argv[3]);
    }
    
    std::cout << "=== Volume Mesh Generator using TetGen ===" << std::endl;
    std::cout << "Input boundary surface: " << inputFile << std::endl;
    std::cout << "Output mesh: " << outputFile << std::endl;
    std::cout << "Quality parameter: " << quality << std::endl;
    std::cout << std::endl;
    
    // Read boundary faces from VTK file
    std::vector<Point3D> boundaryPoints;
    std::vector<Triangle> boundaryTriangles;
    
    std::cout << "Reading boundary surface..." << std::endl;
    if (!readVTKFile(inputFile, boundaryPoints, boundaryTriangles)) {
        return 1;
    }
    
    if (boundaryPoints.empty() || boundaryTriangles.empty()) {
        std::cerr << "Error: No points or triangles found in input file" << std::endl;
        return 1;
    }
    
    // Compute bounding box
    std::cout << "\nComputing bounding box..." << std::endl;
    double xmin, xmax, ymin, ymax, zmin, zmax, maxLength;
    computeBoundingBox(boundaryPoints, xmin, xmax, ymin, ymax, zmin, zmax, maxLength);
    
    // Calculate box dimensions (5x the maximum length)
    double boxSize = 5.0 * maxLength;
    double centerX = (xmin + xmax) / 2.0;
    double centerY = (ymin + ymax) / 2.0;
    double centerZ = (zmin + zmax) / 2.0;
    
    double boxXmin = centerX - boxSize / 2.0;
    double boxXmax = centerX + boxSize / 2.0;
    double boxYmin = centerY - boxSize / 2.0;
    double boxYmax = centerY + boxSize / 2.0;
    double boxZmin = centerZ - boxSize / 2.0;
    double boxZmax = centerZ + boxSize / 2.0;
    
    std::cout << "\nCreating bounding box (5x maximum length)..." << std::endl;
    std::cout << "Box size: " << boxSize << std::endl;
    std::cout << "Box bounds: [" << boxXmin << ", " << boxXmax << "] x [" 
              << boxYmin << ", " << boxYmax << "] x [" 
              << boxZmin << ", " << boxZmax << "]" << std::endl;
    
    // Compute hole point (inside the boundary surface)
    Point3D holePoint = computeBoundaryCenter(boundaryPoints);
    std::cout << "\nHole point (inside boundary): (" 
              << holePoint.x << ", " << holePoint.y << ", " << holePoint.z << ")" << std::endl;
    
    // Verify hole point is inside boundary (diagnostic)
    bool holePointInside = isPointInsideBoundary(holePoint, boundaryPoints, boundaryTriangles);
    std::cout << "  Verification: hole point is " << (holePointInside ? "INSIDE" : "OUTSIDE") 
              << " the boundary (should be INSIDE)" << std::endl;
    
    // Create box faces
    std::vector<Point3D> boxPoints;
    std::vector<std::vector<int> > boxFaces;
    createBoxFaces(boxXmin, boxXmax, boxYmin, boxYmax, boxZmin, boxZmax, 
                   boxPoints, boxFaces);
    
    // Verify region point will be outside boundary (diagnostic)
    // Region point will be near box corner, so compute it here for verification
    Point3D boxCorner = boxPoints[6];
    Point3D boxCorner0 = boxPoints[0];
    Point3D testRegionPoint(
        boxCorner0.x + (boxCorner.x - boxCorner0.x) * 0.99,
        boxCorner0.y + (boxCorner.y - boxCorner0.y) * 0.99,
        boxCorner0.z + (boxCorner.z - boxCorner0.z) * 0.99
    );
    bool regionPointInside = isPointInsideBoundary(testRegionPoint, boundaryPoints, boundaryTriangles);
    std::cout << "  Region point will be: (" << testRegionPoint.x << ", " 
              << testRegionPoint.y << ", " << testRegionPoint.z << ")" << std::endl;
    std::cout << "  Verification: region point is " << (regionPointInside ? "INSIDE" : "OUTSIDE") 
              << " the boundary (should be OUTSIDE)" << std::endl;
    
    // Create TetGen input
    std::cout << "\nPreparing TetGen input..." << std::endl;
    tetgenio in, out;
    createTetGenInput(boundaryPoints, boundaryTriangles, boxPoints, boxFaces, holePoint, in);
    
    // Run TetGen
    std::cout << "\nRunning TetGen..." << std::endl;
    try {
        // TetGen switches:
        // p - Tetrahedralize a piecewise linear complex (PLC)
        // q - Quality mesh generation with specified quality
        // Y - Preserve input surface mesh (fixed surface tessellation)
        // A - Assign attributes to regions
        // V - Verbose output
        std::ostringstream switches;
        switches << "pq" << quality << "YAV";
        
        std::cout << "TetGen switches: " << switches.str() << std::endl;
        std::cout << "  p: Tetrahedralize PLC" << std::endl;
        std::cout << "  q" << quality << ": Quality mesh generation" << std::endl;
        std::cout << "  Y: Preserve input surface mesh (fixed boundary tessellation)" << std::endl;
        std::cout << "  A: Assign region attributes" << std::endl;
        std::cout << "  V: Verbose output" << std::endl;
        
        tetrahedralize(const_cast<char*>(switches.str().c_str()), &in, &out);
        
        std::cout << "\n=== TetGen Results ===" << std::endl;
        std::cout << "Generated " << out.numberofpoints << " nodes" << std::endl;
        std::cout << "Generated " << out.numberoftetrahedra << " tetrahedra" << std::endl;
        std::cout << "Generated " << out.numberoftrifaces << " boundary faces" << std::endl;
        
        // Write output (filtering out internal tetrahedra)
        std::cout << "\nWriting output..." << std::endl;
        writeVTKOutput(out, outputFile, boundaryPoints, boundaryTriangles);
        
        std::cout << "\nMesh generation completed successfully!" << std::endl;
        std::cout << "The volume mesh has a hole defined by the boundary surface." << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during tetrahedralization: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

