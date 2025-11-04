#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include "egads.h"

/*
 * Program to tessellate a hemisphere from an IGES file using ESP/EGADS
 * and output to VTK format with prismatic layer extrusion with growth rate
 *
 * Compile with:
 *   g++ -o tessellate_hemisphere tessellate_hemisphere.cpp -I${ESP_ROOT}/include \
 *       -L${ESP_ROOT}/lib -legads -lm -ldl -lstdc++
 *
 * Usage:
 *   ./tessellate_hemisphere <input.iges> <output.vtk> [num_layers] [first_layer_thickness] [growth_rate]
 *   Default: 5 layers, first layer = 0.01, growth rate = 1.2
 */

// Structure to store vertex information
struct Vertex {
    double x, y, z;
    double nx, ny, nz;  // Normal vector
    int global_id;
};

// Function to calculate triangle normal
void calculateTriangleNormal(const double* v0, const double* v1, const double* v2,
                            double* normal) {
    double edge1[3], edge2[3];
    
    // Calculate edge vectors
    edge1[0] = v1[0] - v0[0];
    edge1[1] = v1[1] - v0[1];
    edge1[2] = v1[2] - v0[2];
    
    edge2[0] = v2[0] - v0[0];
    edge2[1] = v2[1] - v0[1];
    edge2[2] = v2[2] - v0[2];
    
    // Cross product
    normal[0] = edge1[1] * edge2[2] - edge1[2] * edge2[1];
    normal[1] = edge1[2] * edge2[0] - edge1[0] * edge2[2];
    normal[2] = edge1[0] * edge2[1] - edge1[1] * edge2[0];
    
    // Normalize
    double len = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    if (len > 1e-10) {
        normal[0] /= len;
        normal[1] /= len;
        normal[2] /= len;
    }
}

int main(int argc, char *argv[])
{
    ego    context, model, geom, *bodies, *dum;
    int    status, oclass, mtype, nbody, *senses;
    int    nface, nedge, nnode;
    int    major, minor;
    double size, box[6], params[3];
    const char *OCCrev;
    
    // Parse command line arguments
    int num_layers = 5;
    double first_layer_thickness = 0.01;
    double growth_rate = 1.2;  // Each layer is 1.2x thicker than previous
    
    if (argc < 3) {
        printf("Usage: %s <input.iges> <output.vtk> [num_layers] [first_layer_thickness] [growth_rate]\n", argv[0]);
        printf("Default: 5 layers, first layer = 0.01, growth rate = 1.2\n");
        printf("\nGrowth rate examples:\n");
        printf("  1.0 = uniform thickness\n");
        printf("  1.2 = each layer 20%% thicker (recommended)\n");
        printf("  1.5 = each layer 50%% thicker (aggressive)\n");
        return 1;
    }
    
    if (argc >= 4) {
        num_layers = atoi(argv[3]);
        if (num_layers < 1 || num_layers > 100) {
            printf("Error: num_layers must be between 1 and 100\n");
            return 1;
        }
    }
    
    if (argc >= 5) {
        first_layer_thickness = atof(argv[4]);
        if (first_layer_thickness <= 0) {
            printf("Error: first_layer_thickness must be positive\n");
            return 1;
        }
    }
    
    if (argc >= 6) {
        growth_rate = atof(argv[5]);
        if (growth_rate < 1.0 || growth_rate > 2.0) {
            printf("Error: growth_rate must be between 1.0 and 2.0\n");
            return 1;
        }
    }
    
    printf("Prismatic layer settings:\n");
    printf("  Number of layers: %d\n", num_layers);
    printf("  First layer thickness: %.6f\n", first_layer_thickness);
    printf("  Growth rate: %.3f\n", growth_rate);
    
    // Calculate total thickness
    double total_thickness = 0.0;
    double current_thickness = first_layer_thickness;
    for (int i = 0; i < num_layers; i++) {
        total_thickness += current_thickness;
        current_thickness *= growth_rate;
    }
    printf("  Total thickness: %.6f\n\n", total_thickness);
    
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
    
    // Tessellation parameters
    params[0] = size * 0.05;  // 5% of diagonal for max edge length
    params[1] = size * 0.0005; // 0.1% of diagonal for sag
    params[2] = 7.5;          // 15 degrees max angle
    
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
    
    // Collect all surface vertices and triangles
    std::vector<Vertex> surface_vertices;
    std::vector<int> surface_triangles;  // Indices into surface_vertices
    
    // First pass: collect vertices and calculate vertex normals
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
            
            int vertex_base = surface_vertices.size();
            
            // Add vertices for this face
            for (int i = 0; i < len; i++) {
                Vertex v;
                v.x = xyz[3*i];
                v.y = xyz[3*i+1];
                v.z = xyz[3*i+2];
                v.nx = 0.0;
                v.ny = 0.0;
                v.nz = 0.0;
                v.global_id = surface_vertices.size();
                surface_vertices.push_back(v);
            }
            
            // Add triangles and accumulate normals
            for (int i = 0; i < ntri; i++) {
                int i0 = vertex_base + tris[3*i] - 1;
                int i1 = vertex_base + tris[3*i+1] - 1;
                int i2 = vertex_base + tris[3*i+2] - 1;
                
                surface_triangles.push_back(i0);
                surface_triangles.push_back(i1);
                surface_triangles.push_back(i2);
                
                // Calculate triangle normal
                double v0[3] = {surface_vertices[i0].x, surface_vertices[i0].y, surface_vertices[i0].z};
                double v1[3] = {surface_vertices[i1].x, surface_vertices[i1].y, surface_vertices[i1].z};
                double v2[3] = {surface_vertices[i2].x, surface_vertices[i2].y, surface_vertices[i2].z};
                double normal[3];
                calculateTriangleNormal(v0, v1, v2, normal);
                
                // Accumulate normal at each vertex
                surface_vertices[i0].nx += normal[0];
                surface_vertices[i0].ny += normal[1];
                surface_vertices[i0].nz += normal[2];
                
                surface_vertices[i1].nx += normal[0];
                surface_vertices[i1].ny += normal[1];
                surface_vertices[i1].nz += normal[2];
                
                surface_vertices[i2].nx += normal[0];
                surface_vertices[i2].ny += normal[1];
                surface_vertices[i2].nz += normal[2];
            }
        }
    }
    
    // Calculate geometric center of all vertices
    double center[3] = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < surface_vertices.size(); i++) {
        center[0] += surface_vertices[i].x;
        center[1] += surface_vertices[i].y;
        center[2] += surface_vertices[i].z;
    }
    center[0] /= surface_vertices.size();
    center[1] /= surface_vertices.size();
    center[2] /= surface_vertices.size();
    
    printf("Geometric center: [%.6f, %.6f, %.6f]\n", center[0], center[1], center[2]);
    
    // Normalize vertex normals and ensure they point outward from center
    int flipped_count = 0;
    for (size_t i = 0; i < surface_vertices.size(); i++) {
        // Normalize the accumulated normal
        double len = sqrt(surface_vertices[i].nx * surface_vertices[i].nx +
                         surface_vertices[i].ny * surface_vertices[i].ny +
                         surface_vertices[i].nz * surface_vertices[i].nz);
        if (len > 1e-10) {
            surface_vertices[i].nx /= len;
            surface_vertices[i].ny /= len;
            surface_vertices[i].nz /= len;
        }
        
        // Vector from center to vertex
        double to_vertex[3];
        to_vertex[0] = surface_vertices[i].x - center[0];
        to_vertex[1] = surface_vertices[i].y - center[1];
        to_vertex[2] = surface_vertices[i].z - center[2];
        
        // Dot product to check if normal points outward
        double dot = surface_vertices[i].nx * to_vertex[0] +
                     surface_vertices[i].ny * to_vertex[1] +
                     surface_vertices[i].nz * to_vertex[2];
        
        // If dot product is negative, normal points inward - flip it
        if (dot < 0) {
            surface_vertices[i].nx = -surface_vertices[i].nx;
            surface_vertices[i].ny = -surface_vertices[i].ny;
            surface_vertices[i].nz = -surface_vertices[i].nz;
            flipped_count++;
        }
    }
    
    printf("Flipped %d normals to point outward\n", flipped_count);
    
    int num_surface_vertices = surface_vertices.size();
    int num_surface_triangles = surface_triangles.size() / 3;
    
    printf("\nSurface mesh statistics (before merging):\n");
    printf("  Surface vertices: %d\n", num_surface_vertices);
    printf("  Surface triangles: %d\n", num_surface_triangles);
    
    // Merge duplicate vertices to eliminate gaps at surface boundaries
    printf("\nMerging duplicate vertices...\n");
    double merge_tolerance = size * 1e-6;  // Very small tolerance based on model size
    printf("  Merge tolerance: %.10f\n", merge_tolerance);
    
    std::vector<int> vertex_map(surface_vertices.size());  // Maps old index to new index
    std::vector<Vertex> unique_vertices;
    
    // Simple O(n^2) merging - could be optimized with spatial hashing for large meshes
    for (size_t i = 0; i < surface_vertices.size(); i++) {
        bool found_duplicate = false;
        
        // Check if this vertex is close to any existing unique vertex
        for (size_t j = 0; j < unique_vertices.size(); j++) {
            double dx = surface_vertices[i].x - unique_vertices[j].x;
            double dy = surface_vertices[i].y - unique_vertices[j].y;
            double dz = surface_vertices[i].z - unique_vertices[j].z;
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            
            if (dist < merge_tolerance) {
                // Found a duplicate - map to existing vertex
                vertex_map[i] = j;
                
                // Average the normals (will be renormalized later)
                unique_vertices[j].nx += surface_vertices[i].nx;
                unique_vertices[j].ny += surface_vertices[i].ny;
                unique_vertices[j].nz += surface_vertices[i].nz;
                
                found_duplicate = true;
                break;
            }
        }
        
        if (!found_duplicate) {
            // This is a new unique vertex
            vertex_map[i] = unique_vertices.size();
            unique_vertices.push_back(surface_vertices[i]);
        }
    }
    
    printf("  Merged %d vertices into %d unique vertices\n",
           (int)surface_vertices.size(), (int)unique_vertices.size());
    printf("  Eliminated %d duplicate vertices\n",
           (int)(surface_vertices.size() - unique_vertices.size()));
    
    // Renormalize the averaged normals in unique vertices
    for (size_t i = 0; i < unique_vertices.size(); i++) {
        double len = sqrt(unique_vertices[i].nx * unique_vertices[i].nx +
                         unique_vertices[i].ny * unique_vertices[i].ny +
                         unique_vertices[i].nz * unique_vertices[i].nz);
        if (len > 1e-10) {
            unique_vertices[i].nx /= len;
            unique_vertices[i].ny /= len;
            unique_vertices[i].nz /= len;
        }
    }
    
    // Replace surface_vertices with unique_vertices
    surface_vertices = unique_vertices;
    
    // Remap triangle indices to use merged vertices
    for (size_t i = 0; i < surface_triangles.size(); i++) {
        surface_triangles[i] = vertex_map[surface_triangles[i]];
    }
    
    // Update counts
    num_surface_vertices = surface_vertices.size();
    
    printf("\nSurface mesh statistics (after merging):\n");
    printf("  Surface vertices: %d\n", num_surface_vertices);
    printf("  Surface triangles: %d\n", num_surface_triangles);
    
    // Create extruded layers with growth rate
    std::vector<Vertex> all_vertices;
    std::vector<int> prism_cells;  // 6 indices per prism (wedge)
    
    // Create all vertex layers with growth rate
    double cumulative_offset = 0.0;
    std::vector<double> layer_offsets;
    layer_offsets.push_back(0.0);  // Surface layer
    
    current_thickness = first_layer_thickness;
    for (int layer = 0; layer < num_layers; layer++) {
        cumulative_offset += current_thickness;
        layer_offsets.push_back(cumulative_offset);
        current_thickness *= growth_rate;  // Grow for next layer
    }
    
    printf("\nLayer offset details:\n");
    for (int layer = 0; layer <= num_layers; layer++) {
        if (layer == 0) {
            printf("  Layer %d (surface): offset = %.6f\n", layer, layer_offsets[layer]);
        } else {
            double thickness = layer_offsets[layer] - layer_offsets[layer-1];
            printf("  Layer %d: offset = %.6f, thickness = %.6f\n",
                   layer, layer_offsets[layer], thickness);
        }
    }
    
    // Create vertices for each layer
    for (int layer = 0; layer <= num_layers; layer++) {
        double offset = layer_offsets[layer];
        for (int i = 0; i < num_surface_vertices; i++) {
            Vertex v;
            v.x = surface_vertices[i].x + offset * surface_vertices[i].nx;
            v.y = surface_vertices[i].y + offset * surface_vertices[i].ny;
            v.z = surface_vertices[i].z + offset * surface_vertices[i].nz;
            v.nx = surface_vertices[i].nx;
            v.ny = surface_vertices[i].ny;
            v.nz = surface_vertices[i].nz;
            all_vertices.push_back(v);
        }
    }
    
    // Create prism cells
    for (int layer = 0; layer < num_layers; layer++) {
        int base_layer = layer * num_surface_vertices;
        int next_layer = (layer + 1) * num_surface_vertices;
        
        for (int i = 0; i < num_surface_triangles; i++) {
            int i0 = surface_triangles[3*i];
            int i1 = surface_triangles[3*i+1];
            int i2 = surface_triangles[3*i+2];
            
            // Prism (wedge) connectivity: bottom triangle + top triangle
            prism_cells.push_back(base_layer + i0);
            prism_cells.push_back(base_layer + i1);
            prism_cells.push_back(base_layer + i2);
            prism_cells.push_back(next_layer + i0);
            prism_cells.push_back(next_layer + i1);
            prism_cells.push_back(next_layer + i2);
        }
    }
    
    int total_vertices = all_vertices.size();
    int total_prisms = prism_cells.size() / 6;
    
    printf("Extruded mesh statistics:\n");
    printf("  Total vertices: %d\n", total_vertices);
    printf("  Total prisms: %d\n", total_prisms);
    
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
    fprintf(fp, "Hemisphere with %d Prismatic Layers from %s\n", num_layers, argv[1]);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    
    // Write all vertices
    fprintf(fp, "POINTS %d float\n", total_vertices);
    for (size_t i = 0; i < all_vertices.size(); i++) {
        fprintf(fp, "%.9f %.9f %.9f\n", all_vertices[i].x, all_vertices[i].y, all_vertices[i].z);
    }
    
    // Write prism cells
    // VTK cell type 13 = WEDGE (prism)
    fprintf(fp, "\nCELLS %d %d\n", total_prisms, total_prisms * 7);  // 7 = 1 count + 6 indices
    for (size_t i = 0; i < prism_cells.size(); i += 6) {
        fprintf(fp, "6 %d %d %d %d %d %d\n",
                prism_cells[i], prism_cells[i+1], prism_cells[i+2],
                prism_cells[i+3], prism_cells[i+4], prism_cells[i+5]);
    }
    
    fprintf(fp, "\nCELL_TYPES %d\n", total_prisms);
    for (int i = 0; i < total_prisms; i++) {
        fprintf(fp, "13\n");  // VTK_WEDGE
    }
    
    fclose(fp);
    printf("VTK file successfully written!\n");
    printf("\nVisualization tips:\n");
    printf("  - Open in ParaView\n");
    printf("  - Use 'Extract Surface' filter to see outer surface\n");
    printf("  - Use 'Slice' filter to see internal structure\n");
    
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
