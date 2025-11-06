#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include "egads.h"

/*
 * Enhanced tessellation with vertex merging to eliminate gaps
 *
 * Compile with:
 *   g++ -o tessellate_nogaps tessellate_nogaps.cpp -I${ESP_ROOT}/include \
 *       -L${ESP_ROOT}/lib -legads -lm -lstdc++
 *
 * Usage:
 *   ./tessellate_nogaps <input.iges> <output.vtk> [merge_tolerance]
 *   Default merge_tolerance = 1e-6
 */

struct Point3D {
    double x, y, z;
};

// Merge duplicate vertices
void mergeVertices(std::vector<Point3D>& points,
                   std::vector<int>& triangles,
                   double tolerance) {
    
    int n_orig = points.size();
    std::vector<int> vertex_map(n_orig);
    std::vector<Point3D> unique_points;
    
    printf("\nMerging duplicate vertices...\n");
    printf("  Original vertices: %d\n", n_orig);
    printf("  Tolerance: %.10f\n", tolerance);
    
    for (int i = 0; i < n_orig; i++) {
        bool found = false;
        
        for (size_t j = 0; j < unique_points.size(); j++) {
            double dx = points[i].x - unique_points[j].x;
            double dy = points[i].y - unique_points[j].y;
            double dz = points[i].z - unique_points[j].z;
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            
            if (dist < tolerance) {
                vertex_map[i] = j;
                found = true;
                break;
            }
        }
        
        if (!found) {
            vertex_map[i] = unique_points.size();
            unique_points.push_back(points[i]);
        }
    }
    
    printf("  Unique vertices: %d\n", (int)unique_points.size());
    printf("  Removed %d duplicates\n", n_orig - (int)unique_points.size());
    
    // Remap triangles
    for (size_t i = 0; i < triangles.size(); i++) {
        triangles[i] = vertex_map[triangles[i]];
    }
    
    points = unique_points;
}

int main(int argc, char *argv[])
{
    ego    context, model, geom, *bodies, *dum;
    int    status, oclass, mtype, nbody, *senses;
    int    major, minor;
    double size, box[6], params[3];
    const char *OCCrev;
    double merge_tol = 1e-6;
    
    if (argc < 3) {
        printf("Usage: %s <input.iges> <output.vtk> [merge_tolerance]\n", argv[0]);
        printf("Default merge_tolerance = 1e-6\n");
        return 1;
    }
    
    if (argc >= 4) {
        merge_tol = atof(argv[3]);
    }
    
    // Initialize EGADS
    status = EG_open(&context);
    if (status != EGADS_SUCCESS) {
        printf("Error: EG_open returned %d\n", status);
        return 1;
    }
    
    EG_revision(&major, &minor, &OCCrev);
    printf("EGADS Version: %d.%d\n", major, minor);
    
    // Load IGES
    printf("Loading IGES file: %s\n", argv[1]);
    status = EG_loadModel(context, 0, argv[1], &model);
    if (status != EGADS_SUCCESS) {
        printf("Error: EG_loadModel returned %d\n", status);
        EG_close(context);
        return 1;
    }
    
    // Get bodies
    status = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
                            &bodies, &senses);
    if (status != EGADS_SUCCESS) {
        printf("Error: EG_getTopology returned %d\n", status);
        EG_close(context);
        return 1;
    }
    
    printf("Number of bodies: %d\n", nbody);
    
    if (nbody == 0) {
        printf("Error: No bodies found\n");
        EG_close(context);
        return 1;
    }
    
    // Get bounding box
    status = EG_getBoundingBox(model, box);
    if (status != EGADS_SUCCESS) {
        printf("Error: EG_getBoundingBox returned %d\n", status);
        EG_close(context);
        return 1;
    }
    
    printf("Bounding box: [%.3f, %.3f, %.3f] to [%.3f, %.3f, %.3f]\n",
           box[0], box[1], box[2], box[3], box[4], box[5]);
    
    size = sqrt((box[3]-box[0])*(box[3]-box[0]) +
                (box[4]-box[1])*(box[4]-box[1]) +
                (box[5]-box[2])*(box[5]-box[2]));
    
    // Fine tessellation parameters
    params[0] = size * 0.02;
    params[1] = size * 0.0005;
    params[2] = 10.0;
    
    printf("\nTessellation parameters:\n");
    printf("  Max edge length: %.6f\n", params[0]);
    printf("  Max sag:         %.6f\n", params[1]);
    printf("  Max angle:       %.1f degrees\n", params[2]);
    
    // Tessellate all bodies
    std::vector<ego> tessellations(nbody);
    
    for (int ibody = 0; ibody < nbody; ibody++) {
        printf("\nTessellating body %d...\n", ibody + 1);
        status = EG_makeTessBody(bodies[ibody], params, &tessellations[ibody]);
        if (status != EGADS_SUCCESS) {
            printf("Error: EG_makeTessBody returned %d\n", status);
            tessellations[ibody] = NULL;
        }
    }
    
    // Collect all points and triangles
    std::vector<Point3D> all_points;
    std::vector<int> all_triangles;
    
    for (int ibody = 0; ibody < nbody; ibody++) {
        if (tessellations[ibody] == NULL) continue;
        
        int nface, nedge;
        status = EG_statusTessBody(tessellations[ibody], &geom, &nface, &nedge);
        if (status != EGADS_SUCCESS) continue;
        
        for (int iface = 1; iface <= nface; iface++) {
            int len, ntri;
            const int    *ptype, *pindex, *tris, *tric;
            const double *xyz, *uv;
            
            status = EG_getTessFace(tessellations[ibody], iface, &len, &xyz, &uv,
                                    &ptype, &pindex, &ntri, &tris, &tric);
            if (status != EGADS_SUCCESS) continue;
            
            int base_idx = all_points.size();
            
            // Add points
            for (int i = 0; i < len; i++) {
                Point3D p;
                p.x = xyz[3*i];
                p.y = xyz[3*i+1];
                p.z = xyz[3*i+2];
                all_points.push_back(p);
            }
            
            // Add triangles
            for (int i = 0; i < ntri; i++) {
                all_triangles.push_back(base_idx + tris[3*i] - 1);
                all_triangles.push_back(base_idx + tris[3*i+1] - 1);
                all_triangles.push_back(base_idx + tris[3*i+2] - 1);
            }
        }
    }
    
    printf("\nBefore merging:\n");
    printf("  Vertices: %d\n", (int)all_points.size());
    printf("  Triangles: %d\n", (int)all_triangles.size()/3);
    
    // Merge duplicate vertices
    mergeVertices(all_points, all_triangles, merge_tol);
    
    printf("\nAfter merging:\n");
    printf("  Vertices: %d\n", (int)all_points.size());
    printf("  Triangles: %d\n", (int)all_triangles.size()/3);
    
    // Write VTK
    printf("\nWriting VTK: %s\n", argv[2]);
    FILE *fp = fopen(argv[2], "w");
    if (!fp) {
        printf("Error: Cannot open output file\n");
        return 1;
    }
    
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Merged Surface Tessellation\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n");
    
    fprintf(fp, "POINTS %d float\n", (int)all_points.size());
    for (size_t i = 0; i < all_points.size(); i++) {
        fprintf(fp, "%.9f %.9f %.9f\n", all_points[i].x, all_points[i].y, all_points[i].z);
    }
    
    fprintf(fp, "\nPOLYGONS %d %d\n", (int)all_triangles.size()/3, (int)all_triangles.size()*4/3);
    for (size_t i = 0; i < all_triangles.size(); i += 3) {
        fprintf(fp, "3 %d %d %d\n", all_triangles[i], all_triangles[i+1], all_triangles[i+2]);
    }
    
    fclose(fp);
    
    // Cleanup
    for (int ibody = 0; ibody < nbody; ibody++) {
        if (tessellations[ibody] != NULL) {
            EG_deleteObject(tessellations[ibody]);
        }
    }
    EG_deleteObject(model);
    EG_close(context);
    
    printf("\n✅ Success! Gaps should be eliminated.\n");
    return 0;
}
