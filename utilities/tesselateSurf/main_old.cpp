#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "egads.h"

/*
 * Program to tessellate a hemisphere from an IGES file using ESP/EGADS
 * and output directly to VTK format
 *
 * Compile with:
 *   g++ -o tessellate_hemisphere tessellate_hemisphere.cpp -I${ESP_ROOT}/include \
 *       -L${ESP_ROOT}/lib -legads -lm -ldl -lstdc++
 *
 * Usage:
 *   ./tessellate_hemisphere <input.iges> <output.vtk>
 */

int main(int argc, char *argv[])
{
    ego    context, model, geom, *bodies, *dum, tess;
    int    status, oclass, mtype, nbody, *senses;
    int    nface, nedge, nnode;
    int    major, minor;
    double size, box[6], params[3];
    const char *OCCrev;
    
    if (argc != 3) {
        printf("Usage: %s <input.iges> <output.vtk>\n", argv[0]);
        return 1;
    }
    
    // Initialize EGADS
    status = EG_open(&context);
    if (status != EGADS_SUCCESS) {
        printf("Error: EG_open returned %d\n", status);
        return 1;
    }
    
    // Get EGADS and OCC version
    EG_revision(&major, &minor, &OCCrev);
    printf("EGADS Version: %d.%d\n", major, minor);
    printf("OpenCASCADE Version: %s\n", OCCrev);
    
    // Load the IGES file
    printf("Loading IGES file: %s\n", argv[1]);
    status = EG_loadModel(context, 0, argv[1], &model);
    if (status != EGADS_SUCCESS) {
        printf("Error: EG_loadModel returned %d\n", status);
        EG_close(context);
        return 1;
    }
    
    // Get the bodies from the model
    status = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
                            &bodies, &senses);
    if (status != EGADS_SUCCESS) {
        printf("Error: EG_getTopology returned %d\n", status);
        EG_close(context);
        return 1;
    }
    
    printf("Number of bodies: %d\n", nbody);
    
    if (nbody == 0) {
        printf("Error: No bodies found in model\n");
        EG_close(context);
        return 1;
    }
    
    // Arrays to store tessellations for each body
    ego *tessellations = (ego *)malloc(nbody * sizeof(ego));
    if (tessellations == NULL) {
        printf("Error: Memory allocation failed\n");
        EG_close(context);
        return 1;
    }
    
    // Tessellation parameters
    // params[0] = maximum edge length (chord length)
    // params[1] = maximum surface curvature (sag)
    // params[2] = maximum angle between surface normals (degrees)
    
    // Get overall bounding box to determine tessellation size
    status = EG_getBoundingBox(model, box);
    if (status != EGADS_SUCCESS) {
        printf("Error: EG_getBoundingBox returned %d\n", status);
        free(tessellations);
        EG_close(context);
        return 1;
    }
    
    printf("Overall bounding box: [%.3f, %.3f, %.3f] to [%.3f, %.3f, %.3f]\n",
           box[0], box[1], box[2], box[3], box[4], box[5]);
    
    // Calculate a reasonable tessellation size based on bounding box
    size = sqrt((box[3]-box[0])*(box[3]-box[0]) +
                (box[4]-box[1])*(box[4]-box[1]) +
                (box[5]-box[2])*(box[5]-box[2]));
    
    params[0] = size * 0.05;  // 5% of diagonal for max edge length
    params[1] = size * 0.001; // 0.1% of diagonal for sag
    params[2] = 15.0;          // 15 degrees max angle
    
    printf("Tessellation parameters:\n");
    printf("  Max edge length: %.6f\n", params[0]);
    printf("  Max sag:         %.6f\n", params[1]);
    printf("  Max angle:       %.1f degrees\n", params[2]);
    
    // Tessellate all bodies
    int total_faces = 0;
    for (int ibody = 0; ibody < nbody; ibody++) {
        printf("\nTessellating body %d...\n", ibody + 1);
        status = EG_makeTessBody(bodies[ibody], params, &tessellations[ibody]);
        if (status != EGADS_SUCCESS) {
            printf("Error: EG_makeTessBody returned %d for body %d\n", status, ibody + 1);
            tessellations[ibody] = NULL;
            continue;
        }
        
        // Get tessellation statistics for this body
        int nface_body, nedge_body;
        status = EG_statusTessBody(tessellations[ibody], &geom, &nface_body, &nedge_body);
        if (status != EGADS_SUCCESS) {
            printf("Error: EG_statusTessBody returned %d for body %d\n", status, ibody + 1);
            continue;
        }
        
        printf("  Body %d: %d faces, %d edges\n", ibody + 1, nface_body, nedge_body);
        total_faces += nface_body;
    }
    
    printf("\nTotal faces across all bodies: %d\n", total_faces);
    
    // First pass: count total points and triangles across all bodies
    int total_points = 0;
    int total_triangles = 0;
    
    for (int ibody = 0; ibody < nbody; ibody++) {
        if (tessellations[ibody] == NULL) continue;
        
        int nface_body, nedge_body;
        status = EG_statusTessBody(tessellations[ibody], &geom, &nface_body, &nedge_body);
        if (status != EGADS_SUCCESS) continue;
        
        for (int iface = 1; iface <= nface_body; iface++) {
            int len, ntri;
            const int    *ptype, *pindex, *tris, *tric;
            const double *xyz, *uv;
            
            status = EG_getTessFace(tessellations[ibody], iface, &len, &xyz, &uv,
                                    &ptype, &pindex, &ntri, &tris, &tric);
            if (status != EGADS_SUCCESS) {
                printf("Error: EG_getTessFace returned %d for body %d, face %d\n",
                       status, ibody + 1, iface);
                continue;
            }
            
            printf("  Body %d, Face %d: %d vertices, %d triangles\n",
                   ibody + 1, iface, len, ntri);
            total_points += len;
            total_triangles += ntri;
        }
    }
    
    printf("\nTotal vertices: %d\n", total_points);
    printf("Total triangles: %d\n", total_triangles);
    
    // Write VTK file
    printf("\nWriting VTK file: %s\n", argv[2]);
    FILE *fp = fopen(argv[2], "w");
    if (fp == NULL) {
        printf("Error: Cannot open output file %s\n", argv[2]);
        for (int ibody = 0; ibody < nbody; ibody++) {
            if (tessellations[ibody] != NULL) {
                EG_deleteObject(tessellations[ibody]);
            }
        }
        free(tessellations);
        EG_close(context);
        return 1;
    }
    
    // VTK header
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Hemisphere Surface Tessellation from %s\n", argv[1]);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n");
    
    // Write all vertices from all bodies
    fprintf(fp, "POINTS %d float\n", total_points);
    
    for (int ibody = 0; ibody < nbody; ibody++) {
        if (tessellations[ibody] == NULL) continue;
        
        int nface_body, nedge_body;
        status = EG_statusTessBody(tessellations[ibody], &geom, &nface_body, &nedge_body);
        if (status != EGADS_SUCCESS) continue;
        
        for (int iface = 1; iface <= nface_body; iface++) {
            int len, ntri;
            const int    *ptype, *pindex, *tris, *tric;
            const double *xyz, *uv;
            
            status = EG_getTessFace(tessellations[ibody], iface, &len, &xyz, &uv,
                                    &ptype, &pindex, &ntri, &tris, &tric);
            if (status != EGADS_SUCCESS) continue;
            
            for (int i = 0; i < len; i++) {
                fprintf(fp, "%.9f %.9f %.9f\n", xyz[3*i], xyz[3*i+1], xyz[3*i+2]);
            }
        }
    }
    
    // Write triangles with correct vertex indexing across all bodies
    fprintf(fp, "\nPOLYGONS %d %d\n", total_triangles, 4*total_triangles);
    
    int vertex_offset = 0;
    for (int ibody = 0; ibody < nbody; ibody++) {
        if (tessellations[ibody] == NULL) continue;
        
        int nface_body, nedge_body;
        status = EG_statusTessBody(tessellations[ibody], &geom, &nface_body, &nedge_body);
        if (status != EGADS_SUCCESS) continue;
        
        for (int iface = 1; iface <= nface_body; iface++) {
            int len, ntri;
            const int    *ptype, *pindex, *tris, *tric;
            const double *xyz, *uv;
            
            status = EG_getTessFace(tessellations[ibody], iface, &len, &xyz, &uv,
                                    &ptype, &pindex, &ntri, &tris, &tric);
            if (status != EGADS_SUCCESS) continue;
            
            // Write triangles for this face with global vertex indices
            for (int i = 0; i < ntri; i++) {
                // tris indices are 1-based, convert to 0-based and add offset
                fprintf(fp, "3 %d %d %d\n",
                        vertex_offset + tris[3*i] - 1,
                        vertex_offset + tris[3*i+1] - 1,
                        vertex_offset + tris[3*i+2] - 1);
            }
            
            vertex_offset += len;
        }
    }
    
    fclose(fp);
    printf("VTK file successfully written!\n");
    printf("Visualization tip: Open in ParaView and set representation to 'Surface With Edges'\n");
    
    // Clean up
    for (int ibody = 0; ibody < nbody; ibody++) {
        if (tessellations[ibody] != NULL) {
            EG_deleteObject(tessellations[ibody]);
        }
    }
    free(tessellations);
    EG_deleteObject(model);
    EG_close(context);
    
    printf("\nSuccess!\n");
    return 0;
}

