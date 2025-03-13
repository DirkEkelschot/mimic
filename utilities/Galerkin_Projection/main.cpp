#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <CGAL/Nef_polyhedron_3.h>
// #include <CGAL/Nef_3>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
// Define the kernel type
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Define the Delaunay triangulation type
typedef CGAL::Delaunay_triangulation_3<Kernel> Delaunay;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

// #include <CGAL/Point_location.h>
// Define a simple 2D point structure
struct Point {
    double x, y, z;
};

// Define a triangular element with vertices and cell-centered value

Nef_polyhedron create_tetrahedron(Kernel::Point_3 p1, Kernel::Point_3 p2, Kernel::Point_3 p3, Kernel::Point_3 p4) {
    Nef_polyhedron tetrahedron;
    tetrahedron = create_tetrahedron(p1, p2, p3, p4);
    return tetrahedron;
}



// Function to compute the intersection volume of two tetrahedra
void compute_intersection_volume(Kernel::Point_3 p1[], Kernel::Point_3 p2[]) 
{
    // Create Nef Polyhedra for the tetrahedra
    // Nef_polyhedron tetrahedron1 = create_tetrahedron(p1[0], p1[1], p1[2], p1[3]);
    // Nef_polyhedron tetrahedron2 = create_tetrahedron(p2[0], p2[1], p2[2], p2[3]);
    // Nef_polyhedron tetrahedron1;
    // Nef_polyhedron tetrahedron2;
    Polyhedron pm1;
    pm1.make_tetrahedron(p1[0], p1[1], p1[2], p1[3]);
    // Convert the Polyhedron to a Nef_polyhedron
    Nef_polyhedron nef1(pm1);
    Polyhedron pm2;
    pm2.make_tetrahedron(p2[0], p2[1], p2[2], p2[3]);
    // Convert the Polyhedron to a Nef_polyhedron
    Nef_polyhedron nef2(pm2);

    //Nef_polyhedron intersection;
    //Compute the intersection
    // Nef_polyhedron intersection = pm1.intersection(pm2);
    Nef_polyhedron intersection = nef1 * nef2;
    // Calculate the volume of the intersection
    // Note: This step is simplified and might require additional handling
    // depending on the actual implementation of volume calculation for Nef Polyhedra
    if (!intersection.is_empty()) {
        // Assuming a method to calculate volume exists
        // This is a placeholder; actual implementation may vary
        // double volume = intersection.volume();
        Polyhedron pm_intersection;
        CGAL::convert_nef_polyhedron_to_polygon_mesh(intersection, pm_intersection);
        Polyhedron pm1s;
        CGAL::convert_nef_polyhedron_to_polygon_mesh(nef1, pm1s);
        Polyhedron pm2s;
        CGAL::convert_nef_polyhedron_to_polygon_mesh(nef1, pm2s);
        // Calculate volume of intersection
        double volume  = CGAL::Polygon_mesh_processing::volume(pm_intersection);
        double volume1 = CGAL::Polygon_mesh_processing::volume(pm1s);
        double volume2 = CGAL::Polygon_mesh_processing::volume(pm2s);
        // double volume = 0.0;
        // Nef_polyhedron pm;
        // CGAL::convert_nef_polyhedron_to_polygon_mesh(intersection, pm);
        // double volume = CGAL::Polygon_mesh_processing::volume(pm);
        
        //std::cout << "Intersection Volume: " << volume << " " << volume1 << " " << volume2 << std::endl;
    } else {
        std::cout << "No intersection." << std::endl;
    }
}


struct Tetrahedron {
    Point vertices[4];
    double value; // Cell-centered value

    // Constructor to initialize Triangle objects
    Tetrahedron(const Point& v1, const Point& v2, const Point& v3, const Point& v4, double val)
        : value(val) {
        vertices[0] = v1;
        vertices[1] = v2;
        vertices[2] = v3;
        vertices[3] = v4;
    }
};

// Function to compute volume of a tetrahedron
double computeVolume(const Tetrahedron& tetra) {
    // Simplified formula for volume computation
    // In practice, use a more robust method like determinant calculation
    double a = tetra.vertices[0].x - tetra.vertices[3].x;
    double b = tetra.vertices[0].y - tetra.vertices[3].y;
    double c = tetra.vertices[0].z - tetra.vertices[3].z;
    double d = tetra.vertices[1].x - tetra.vertices[3].x;
    double e = tetra.vertices[1].y - tetra.vertices[3].y;
    double f = tetra.vertices[1].z - tetra.vertices[3].z;
    double g = tetra.vertices[2].x - tetra.vertices[3].x;
    double h = tetra.vertices[2].y - tetra.vertices[3].y;
    double i = tetra.vertices[2].z - tetra.vertices[3].z;
    // std::cout << "coords " << tetra.vertices[0].x << " " << tetra.vertices[3].x << std::endl;
    // std::cout << "abc " << a << " " << b << " " << c << std::endl;
    double det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    // std::cout << "det " << det << std::endl; 

    return std::abs(det) / 6.0;
}

// Function to compute overlap volume between two tetrahedra (simplified for this example)
double computeOverlapVolume(const Tetrahedron& tetra1, const Tetrahedron& tetra2) {
    // For simplicity, assume full overlap (in practice, use geometric intersection algorithms)
    // std::cout << "vol T1 " << computeVolume(tetra1) << " T2 " << computeVolume(tetra2) << std::endl;
    return std::min(computeVolume(tetra1), computeVolume(tetra2));
}


void galerkinProjectionVolume(const std::vector<Tetrahedron>& donorMesh,
    const std::vector<Tetrahedron>& targetMesh,
    std::vector<double>& interpolatedValues) {
size_t numDonor = donorMesh.size();
size_t numTarget = targetMesh.size();

// Assemble mass matrix and RHS vector
Eigen::MatrixXd M(numTarget, numTarget);
Eigen::VectorXd b(numTarget);
M.setZero();
b.setZero();
std::cout << numTarget << " " << numTarget << std::endl;
for (size_t i = 0; i < numTarget; ++i) 
{
    // for(int zz = 0;zz<4;zz++)
    // {
    //     std::cout << targetMesh[i].vertices[zz].x << " " << targetMesh[i].vertices[zz].y << " " << targetMesh[i].vertices[zz].z << std::endl;
    // }
    // std::cout << "VOLUMES "<< " "  << computeVolume(targetMesh[i]) << std::endl;
    for (size_t j = 0; j < numDonor; ++j) 
    {
        // std::cout << "VOLUMES "<< " "  << computeVolume(targetMesh[i]) << ", " << computeVolume(donorMesh[j]) << std::endl;
        //std::cout << "test " << targetMesh[i].vertices[0].x << " " << targetMesh[i].vertices[0].y << std::endl;
        double overlapVolume = computeOverlapVolume(targetMesh[i], donorMesh[j]);
        // std::cout << "overlapVolume " << overlapVolume << std::endl;
        M(i, i) += overlapVolume; // Diagonal entries of mass matrix
        b(i) += overlapVolume * donorMesh[j].value; // RHS vector
        std::cout << M(i, j) << " ";
    }

    std::cout << std::endl;
}

for (size_t i = 0; i < numTarget; ++i) 
{
    std::cout << b(i) << std::endl;
}

// Solve linear system M * u = b
Eigen::VectorXd u = M.colPivHouseholderQr().solve(b);

// Store interpolated values
interpolatedValues.resize(numTarget);
for (size_t i = 0; i < numTarget; ++i) 
{
    interpolatedValues[i] = u(i);
    std::cout <<  u(i) << std::endl;
}
}



// struct Triangle {
//     Point vertices[3];
//     double value;

//     // Constructor to initialize Triangle objects
//     Triangle(const Point& v1, const Point& v2, const Point& v3, double val)
//         : value(val) {
//         vertices[0] = v1;
//         vertices[1] = v2;
//         vertices[2] = v3;
//     }
// };


// // Function to compute area of a triangle
// double computeArea(const Triangle& tri) {
//     const auto& [x1, y1, z1] = tri.vertices[0];
//     const auto& [x2, y2, z1] = tri.vertices[1];
//     const auto& [x3, y3, z1] = tri.vertices[2];
//     return 0.5 * std::abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
// }

// Function to compute overlap area between two triangles (simplified for this example)
// double computeOverlapArea(const Triangle& tri1, const Triangle& tri2) {
//     // For simplicity, assume full overlap (in practice, use geometric intersection algorithms)
//     return std::min(computeArea(tri1), computeArea(tri2));
// }

// // Galerkin projection function
// void galerkinProjection(const std::vector<Triangle>& donorMesh,
//                         const std::vector<Triangle>& targetMesh,
//                         std::vector<double>& interpolatedValues) {
//     size_t numDonor = donorMesh.size();
//     size_t numTarget = targetMesh.size();

//     // Assemble mass matrix and RHS vector
//     Eigen::MatrixXd M(numTarget, numTarget);
//     Eigen::VectorXd b(numTarget);
//     M.setZero();
//     b.setZero();

//     for (size_t i = 0; i < numTarget; ++i) {
//         for (size_t j = 0; j < numDonor; ++j) {
//             double overlapArea = computeOverlapArea(targetMesh[i], donorMesh[j]);
//             M(i, i) += overlapArea; // Diagonal entries of mass matrix
//             b(i) += overlapArea * donorMesh[j].value; // RHS vector
//         }
//     }

//     // Solve linear system M * u = b
//     Eigen::VectorXd u = M.colPivHouseholderQr().solve(b);

//     // Store interpolated values
//     interpolatedValues.resize(numTarget);
//     for (size_t i = 0; i < numTarget; ++i) {
//         interpolatedValues[i] = u(i);
//     }
// }

int main() {
    // Define donor mesh: Two triangles forming a square (top-left to bottom-right split)
    // std::vector<Triangle> donorMesh = {
    //     {Point{0.0, 0.0}, Point{1.0, 0.0}, Point{0.0, 1.0}, 1.0},
    //     {Point{1.0, 1.0}, Point{1.0, 0.0}, Point{0.0, 1.0}, 2.0}
    // };

    // std::vector<Triangle> targetMesh = {
    //     {Point{0.0, 0.0}, Point{1.0, 0.0}, Point{1.0, 1.0}, 0.0},
    //     {Point{0.0, 0.0}, Point{1.0, 1.0}, Point{0.0, 1.0}, 0.0}
    // };

    // // Perform Galerkin projection
    // std::vector<double> interpolatedValues;
    // galerkinProjection(donorMesh, targetMesh, interpolatedValues);

    // // Output results
    // std::cout << "Interpolated values on target mesh:\n";
    // for (size_t i = 0; i < interpolatedValues.size(); ++i) {
    //     std::cout << "Triangle " << i + 1 << ": " << interpolatedValues[i] << "\n";
    // }


    std::vector<Kernel::Point_3> points = {
        Kernel::Point_3(0, 0, 0),
        Kernel::Point_3(1, 0, 0),
        Kernel::Point_3(0, 1, 0),
        Kernel::Point_3(0, 0, 1)
    };

    // Define donor mesh: A box consisting of 6 tetrahedra
    std::vector<Tetrahedron> donorMesh = {
        {Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0}, Point{1.0, 1.0, 0.0}, Point{1.0, 1.0, 1.0}, 1.0},
        {Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0}, Point{1.0, 0.0, 1.0}, Point{1.0, 1.0, 1.0}, 2.0},
        {Point{0.0, 0.0, 0.0}, Point{0.0, 1.0, 0.0}, Point{1.0, 1.0, 0.0}, Point{1.0, 1.0, 1.0}, 3.0},
        {Point{0.0, 0.0, 0.0}, Point{0.0, 1.0, 0.0}, Point{0.0, 1.0, 1.0}, Point{1.0, 1.0, 1.0}, 4.0},
        {Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0}, Point{0.0, 0.0, 1.0}, Point{1.0, 1.0, 1.0}, 5.0},
        {Point{0.0, 0.0, 0.0}, Point{0.0, 1.0, 0.0}, Point{0.0, 0.0, 1.0}, Point{1.0, 1.0, 1.0}, 6.0}
    };

    // Define target mesh: A box with different orientation, also 6 tetrahedra
    std::vector<Tetrahedron> targetMesh = {
        {Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0}, Point{0.0, 1.0, 0.0}, Point{0.0, 0.0, 1.0}, 0.0},
        {Point{1.0, 1.0, 1.0}, Point{1.0, 1.0, 0.0}, Point{1.0, 0.0, 1.0}, Point{0.0, 1.0, 1.0}, 0.0},
        {Point{0.0, 0.0, 0.0}, Point{1.0, 1.0, 0.0}, Point{0.0, 1.0, 1.0}, Point{1.0, 0.0, 1.0}, 0.0},
        {Point{0.0, 0.0, 0.0}, Point{1.0, 1.0, 1.0}, Point{1.0, 0.0, 0.0}, Point{0.0, 1.0, 0.0}, 0.0},
        {Point{0.0, 0.0, 0.0}, Point{1.0, 1.0, 1.0}, Point{0.0, 0.0, 1.0}, Point{1.0, 0.0, 0.0}, 0.0},
        {Point{0.0, 0.0, 0.0}, Point{1.0, 1.0, 1.0}, Point{0.0, 1.0, 0.0}, Point{0.0, 0.0, 1.0}, 0.0}
    };

    // Perform Galerkin projection
    std::vector<double> interpolatedValues;
    galerkinProjectionVolume(donorMesh, targetMesh, interpolatedValues);

    // Output results
    std::cout << "Interpolated values on target mesh:\n";
    for (size_t i = 0; i < interpolatedValues.size(); ++i) {
        std::cout << "Tetrahedron " << i + 1 << ": " << interpolatedValues[i] << "\n";
    }


    Kernel::Point_3 p1[] = {Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0), Kernel::Point_3(0, 1, 0), Kernel::Point_3(0, 0, 1)};
    Kernel::Point_3 p2[] = {Kernel::Point_3(0.5, 0, 0), Kernel::Point_3(1.5, 0, 0), Kernel::Point_3(0.5, 1, 0), Kernel::Point_3(0.5, 0, 1)};

    clock_t start = clock();

    for(int i=0;i<10000;i++)
    {
        compute_intersection_volume(p1, p2);
    }

    clock_t end = clock();
    double time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
    std::cout << time_taken << std::endl;

    return 0;
}
