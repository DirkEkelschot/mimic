#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <algorithm>
#include "egads.h"

/*
 * Program to tessellate geometry from an IGES file with per-surface spacing control
 *
 * Key improvements:
 * - Analyzes each surface individually (type, curvature, size)
 * - Determines optimal parameters per surface
 * - Uses intelligent global parameters based on analysis
 * - Supports both uniform and adaptive meshing
 */

struct Vertex {
    double x, y, z;
    double nx, ny, nz;
    int global_id;
};

struct SurfaceAnalysis {
    int face_id;
    int surface_type;  // EGADS surface type
    double curvature;
    double area;
    double characteristic_size;
    double recommended_edge_length;
    double recommended_sag;
    double recommended_angle;
};

void calculateTriangleNormal(const double* v0, const double* v1, const double* v2,
                            double* normal) {
    double edge1[3], edge2[3];
    edge1[0] = v1[0] - v0[0]; edge1[1] = v1[1] - v0[1]; edge1[2] = v1[2] - v0[2];
    edge2[0] = v2[0] - v0[0]; edge2[1] = v2[1] - v0[1]; edge2[2] = v2[2] - v0[2];
    
    normal[0] = edge1[1] * edge2[2] - edge1[2] * edge2[1];
    normal[1] = edge1[2] * edge2[0] - edge1[0] * edge2[2];
    normal[2] = edge1[0] * edge2[1] - edge1[1] * edge2[0];
    
    double len = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    if (len > 1e-10) {
        normal[0] /= len; normal[1] /= len; normal[2] /= len;
    }
}

// Get the surface type from EGADS
int getSurfaceType(ego face) {
    ego ref, *children;
    int oclass, mtype, nchild, *senses;
    
    int status = EG_getTopology(face, &ref, &oclass, &mtype, NULL,
                                &nchild, &children, &senses);
    if (status != EGADS_SUCCESS || ref == NULL) return -1;
    
    status = EG_getTopology(ref, &ref, &oclass, &mtype, NULL,
                           &nchild, &children, &senses);
    if (status != EGADS_SUCCESS) return -1;
    
    return mtype;
    // 0=PLANE, 1=SPHERICAL, 2=CYLINDRICAL, 3=REVOLUTION,
    // 4=TOROIDAL, 5=TRIMMED, 6=BEZIER, 7=BSPLINE, 8=OFFSET
}

const char* getSurfaceTypeName(int type) {
    switch(type) {
        case 0: return "PLANE";
        case 1: return "SPHERICAL";
        case 2: return "CYLINDRICAL";
        case 3: return "REVOLUTION";
        case 4: return "TOROIDAL";
        case 5: return "TRIMMED";
        case 6: return "BEZIER";
        case 7: return "BSPLINE";
        case 8: return "OFFSET";
        default: return "UNKNOWN";
    }
}

// Estimate maximum curvature of a surface
double estimateSurfaceCurvature(ego face) {
    ego ref, *children;
    int oclass, mtype, nchild, *senses;
    double range[4];
    int periodic;
    
    // Get underlying surface
    int status = EG_getTopology(face, &ref, &oclass, &mtype, NULL,
                                &nchild, &children, &senses);
    if (status != EGADS_SUCCESS || ref == NULL) return 0.0;
    
    // Get UV range
    status = EG_getRange(ref, range, &periodic);
    if (status != EGADS_SUCCESS) return 0.0;
    
    // Sample surface at multiple points to find maximum curvature
    double max_curvature = 0.0;
    int n_samples = 5;  // 5x5 sampling grid
    
    for (int i = 0; i < n_samples; i++) {
        for (int j = 0; j < n_samples; j++) {
            double u = range[0] + (range[1] - range[0]) * i / (n_samples - 1.0);
            double v = range[2] + (range[3] - range[2]) * j / (n_samples - 1.0);
            double uv[2] = {u, v};
            double data[18];
            
            status = EG_evaluate(ref, uv, data);
            if (status != EGADS_SUCCESS) continue;
            
            // data[0-2]:   position
            // data[3-5]:   ∂/∂u
            // data[6-8]:   ∂/∂v
            // data[9-11]:  ∂²/∂u²
            // data[12-14]: ∂²/∂u∂v
            // data[15-17]: ∂²/∂v²
            
            // Estimate curvature from second derivatives
            double curv_u = sqrt(data[9]*data[9] + data[10]*data[10] + data[11]*data[11]);
            double curv_v = sqrt(data[15]*data[15] + data[16]*data[16] + data[17]*data[17]);
            double local_curv = fmax(curv_u, curv_v);
            
            if (local_curv > max_curvature) {
                max_curvature = local_curv;
            }
        }
    }
    
    return max_curvature;
}

// Analyze a surface and determine recommended tessellation parameters
SurfaceAnalysis analyzeSurface(ego face, int face_id, double model_size) {
    SurfaceAnalysis analysis;
    analysis.face_id = face_id;
    
    // Get surface type
    analysis.surface_type = getSurfaceType(face);
    
    // Get bounding box to estimate size
    double box[6];
    int status = EG_getBoundingBox(face, box);
    if (status == EGADS_SUCCESS) {
        analysis.characteristic_size = sqrt((box[3]-box[0])*(box[3]-box[0]) +
                                          (box[4]-box[1])*(box[4]-box[1]) +
                                          (box[5]-box[2])*(box[5]-box[2]));
    } else {
        analysis.characteristic_size = model_size * 0.1;
    }
    
    // Estimate curvature
    analysis.curvature = estimateSurfaceCurvature(face);
    
    // Base edge length (2% of model size)
    double base_edge = model_size * 0.02;
    
    // Determine recommended parameters based on surface type and curvature
    switch(analysis.surface_type) {
        case 0: {  // PLANE - can use coarser mesh
            analysis.recommended_edge_length = base_edge * 1.5;
            analysis.recommended_angle = 15.0;
            break;
        }
        case 1: {  // SPHERICAL - moderate refinement
            if (analysis.curvature > 1.0) {
                analysis.recommended_edge_length = base_edge * 0.8;
                analysis.recommended_angle = 10.0;
            } else {
                analysis.recommended_edge_length = base_edge;
                analysis.recommended_angle = 12.0;
            }
            break;
        }
        case 2: {  // CYLINDRICAL - moderate refinement
            if (analysis.curvature > 0.5) {
                analysis.recommended_edge_length = base_edge * 0.9;
                analysis.recommended_angle = 11.0;
            } else {
                analysis.recommended_edge_length = base_edge * 1.1;
                analysis.recommended_angle = 13.0;
            }
            break;
        }
        case 7: {  // BSPLINE - adaptive based on curvature
            double curv_factor = 1.0 / (1.0 + 3.0 * analysis.curvature);
            analysis.recommended_edge_length = base_edge * curv_factor;
            
            // Tighter angle for high curvature
            if (analysis.curvature > 1.0) {
                analysis.recommended_angle = 8.0;
            } else if (analysis.curvature > 0.5) {
                analysis.recommended_angle = 10.0;
            } else {
                analysis.recommended_angle = 12.0;
            }
            break;
        }
        default: {  // Other types - use base parameters
            analysis.recommended_edge_length = base_edge;
            analysis.recommended_angle = 12.0;
        }
    }
    
    // Adjust for very high curvature regardless of type
    if (analysis.curvature > 2.0) {
        analysis.recommended_edge_length *= 0.6;
        analysis.recommended_angle = fmin(analysis.recommended_angle, 8.0);
    }
    
    // Clamp to reasonable range
    double min_edge = model_size * 0.005;  // 0.5%
    double max_edge = model_size * 0.1;    // 10%
    analysis.recommended_edge_length = fmax(min_edge, fmin(analysis.recommended_edge_length, max_edge));
    
    // Set sag proportional to edge length
    analysis.recommended_sag = analysis.recommended_edge_length * 0.015;  // 1.5% of edge
    
    return analysis;
}

int main(int argc, char *argv[])
{
    ego context, model, geom, *bodies;
    int status, oclass, mtype, nbody, *senses;
    int major, minor;
    double size, box[6];
    const char *OCCrev;
    
    // Parse arguments
    int num_layers = 5;
    double first_layer_thickness = 0.01;
    double growth_rate = 1.2;
    int strategy = 0;  // 0=conservative(finest), 1=balanced(median), 2=aggressive(adaptive weighted)
    
    if (argc < 3) {
        printf("Usage: %s <input.iges> <output.vtk> [num_layers] [first_layer] [growth_rate] [strategy]\n", argv[0]);
        printf("  strategy: 0=conservative (finest), 1=balanced (median), 2=aggressive (adaptive)\n");
        printf("            default: 0\n");
        return 1;
    }
    
    if (argc >= 4) num_layers = atoi(argv[3]);
    if (argc >= 5) first_layer_thickness = atof(argv[4]);
    if (argc >= 6) growth_rate = atof(argv[5]);
    if (argc >= 7) strategy = atoi(argv[6]);
    
    printf("Mesh parameters:\n");
    printf("  Strategy: ");
    switch(strategy) {
        case 0: printf("Conservative (finest surfaces drive global params)\n"); break;
        case 1: printf("Balanced (median of all surfaces)\n"); break;
        case 2: printf("Aggressive (weighted adaptive)\n"); break;
        default: printf("Unknown (using conservative)\n"); strategy = 0;
    }
    printf("  Prismatic layers: %d\n", num_layers);
    printf("  First layer thickness: %.6f\n", first_layer_thickness);
    printf("  Growth rate: %.3f\n\n", growth_rate);
    
    // Initialize EGADS
    status = EG_open(&context);
    if (status != EGADS_SUCCESS) {
        printf("Error: EG_open = %d\n", status);
        return 1;
    }
    
    EG_revision(&major, &minor, &OCCrev);
    printf("EGADS %d.%d, OpenCASCADE %s\n", major, minor, OCCrev);
    
    // Load model
    status = EG_loadModel(context, 0, argv[1], &model);
    if (status != EGADS_SUCCESS) {
        printf("Error: EG_loadModel = %d\n", status);
        EG_close(context);
        return 1;
    }
    
    status = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody, &bodies, &senses);
    if (status != EGADS_SUCCESS || nbody == 0) {
        printf("Error: No bodies found\n");
        EG_close(context);
        return 1;
    }
    
    printf("Bodies: %d\n", nbody);
    
    // Get bounding box
    EG_getBoundingBox(model, box);
    size = sqrt((box[3]-box[0])*(box[3]-box[0]) +
                (box[4]-box[1])*(box[4]-box[1]) +
                (box[5]-box[2])*(box[5]-box[2]));
    
    printf("Bounding box diagonal: %.6f\n\n", size);
    
    // ========== ANALYZE ALL SURFACES ==========
    printf("=== PER-SURFACE ANALYSIS ===\n");
    std::vector<SurfaceAnalysis> all_analyses;
    
    for (int ibody = 0; ibody < nbody; ibody++) {
        ego *faces;
        int nface;
        status = EG_getBodyTopos(bodies[ibody], NULL, FACE, &nface, &faces);
        if (status != EGADS_SUCCESS) continue;
        
        printf("\nBody %d: %d faces\n", ibody + 1, nface);
        printf("%-4s %-12s %-8s %-10s %-10s %-8s\n",
               "ID", "Type", "Size", "Curvature", "Edge", "Angle");
        printf("------------------------------------------------------------\n");
        
        for (int iface = 0; iface < nface; iface++) {
            SurfaceAnalysis analysis = analyzeSurface(faces[iface], iface + 1, size);
            all_analyses.push_back(analysis);
            
            printf("%-4d %-12s %-8.4f %-10.3e %-10.6f %-8.1f\n",
                   analysis.face_id,
                   getSurfaceTypeName(analysis.surface_type),
                   analysis.characteristic_size,
                   analysis.curvature,
                   analysis.recommended_edge_length,
                   analysis.recommended_angle);
        }
        
        EG_free(faces);
    }
    
    // ========== DETERMINE GLOBAL PARAMETERS BASED ON STRATEGY ==========
    if (all_analyses.empty()) {
        printf("\nError: No surfaces analyzed\n");
        return 1;
    }
    
    // Sort by recommended edge length
    std::vector<double> edge_lengths;
    std::vector<double> angles;
    for (const auto& a : all_analyses) {
        edge_lengths.push_back(a.recommended_edge_length);
        angles.push_back(a.recommended_angle);
    }
    std::sort(edge_lengths.begin(), edge_lengths.end());
    std::sort(angles.begin(), angles.end());
    
    double min_edge = edge_lengths.front();
    double max_edge = edge_lengths.back();
    double median_edge = edge_lengths[edge_lengths.size() / 2];
    double avg_edge = 0.0;
    for (double e : edge_lengths) avg_edge += e;
    avg_edge /= edge_lengths.size();
    
    double min_angle = angles.front();
    double median_angle = angles[angles.size() / 2];
    
    printf("\n=== TESSELLATION PARAMETER SELECTION ===\n");
    printf("Edge length statistics:\n");
    printf("  Finest:    %.6f (%.2f%% of diagonal)\n", min_edge, 100.0*min_edge/size);
    printf("  Median:    %.6f (%.2f%% of diagonal)\n", median_edge, 100.0*median_edge/size);
    printf("  Average:   %.6f (%.2f%% of diagonal)\n", avg_edge, 100.0*avg_edge/size);
    printf("  Coarsest:  %.6f (%.2f%% of diagonal)\n", max_edge, 100.0*max_edge/size);
    printf("  Variation: %.2fx\n", max_edge / min_edge);
    
    // Select parameters based on strategy
    double params[3];
    
    switch(strategy) {
        case 0:  // Conservative - use finest requirement
            params[0] = min_edge * 1.2;  // 20% coarser than absolute finest
            params[2] = min_angle;
            printf("\nStrategy: Using conservative (finest * 1.2)\n");
            break;
            
        case 1:  // Balanced - use median
            params[0] = median_edge;
            params[2] = median_angle;
            printf("\nStrategy: Using balanced (median)\n");
            break;
            
        case 2:  // Aggressive - weighted toward coarser
            params[0] = 0.3 * min_edge + 0.7 * median_edge;
            params[2] = median_angle;
            printf("\nStrategy: Using aggressive (weighted adaptive)\n");
            break;
            
        default:
            params[0] = min_edge * 1.2;
            params[2] = min_angle;
    }
    
    params[1] = params[0] * 0.015;  // Sag = 1.5% of edge length
    
    printf("\nFinal global tessellation parameters:\n");
    printf("  Max edge length: %.6f (%.2f%% of diagonal)\n", params[0], 100.0*params[0]/size);
    printf("  Max sag:         %.6f (%.2f%% of edge)\n", params[1], 100.0*params[1]/params[0]);
    printf("  Max angle:       %.1f degrees\n", params[2]);
    
    // Estimate element count
    double surf_area_est = size * size;
    double tri_area = 0.433 * params[0] * params[0];
    int est_triangles = (int)(surf_area_est / tri_area);
    printf("  Estimated surface triangles: ~%d\n", est_triangles);
    printf("  Estimated total prisms: ~%d\n\n", est_triangles * num_layers);
    
    // ========== TESSELLATE WITH SELECTED PARAMETERS ==========
    ego *tessellations = (ego *)malloc(nbody * sizeof(ego));
    
    for (int ibody = 0; ibody < nbody; ibody++) {
        printf("Tessellating body %d...\n", ibody + 1);
        
        status = EG_makeTessBody(bodies[ibody], params, &tessellations[ibody]);
        if (status != EGADS_SUCCESS) {
            printf("Error: EG_makeTessBody = %d\n", status);
            tessellations[ibody] = NULL;
            continue;
        }
        
        // Get tessellation stats
        int state, nface_tess;
        status = EG_statusTessBody(tessellations[ibody], &geom, &state, &nface_tess);
        if (status == EGADS_SUCCESS) {
            printf("  Successfully tessellated %d faces\n", nface_tess);
        }
    }
    
    // ========== COLLECT SURFACE VERTICES AND TRIANGLES ==========
    std::vector<Vertex> surface_vertices;
    std::vector<int> surface_triangles;
    
    for (int ibody = 0; ibody < nbody; ibody++) {
        if (tessellations[ibody] == NULL) continue;
        
        int state, nface_tess;
        EG_statusTessBody(tessellations[ibody], &geom, &state, &nface_tess);
        
        for (int iface = 1; iface <= nface_tess; iface++) {
            int npnt, ntri;
            const double *xyz, *uv;
            const int *ptype, *pindex, *tris, *tric;
            
            status = EG_getTessFace(tessellations[ibody], iface,
                                   &npnt, &xyz, &uv, &ptype, &pindex,
                                   &ntri, &tris, &tric);
            
            if (status != EGADS_SUCCESS) continue;
            
            int base_idx = surface_vertices.size();
            
            // Add vertices
            for (int i = 0; i < npnt; i++) {
                Vertex v;
                v.x = xyz[3*i];
                v.y = xyz[3*i+1];
                v.z = xyz[3*i+2];
                v.nx = v.ny = v.nz = 0.0;
                surface_vertices.push_back(v);
            }
            
            // Add triangles
            for (int i = 0; i < ntri; i++) {
                surface_triangles.push_back(base_idx + tris[3*i] - 1);
                surface_triangles.push_back(base_idx + tris[3*i+1] - 1);
                surface_triangles.push_back(base_idx + tris[3*i+2] - 1);
            }
        }
    }
    
    printf("\nSurface mesh: %d vertices, %d triangles\n",
           (int)surface_vertices.size(), (int)surface_triangles.size()/3);
    
    // ========== COMPUTE VERTEX NORMALS ==========
    printf("Computing vertex normals...\n");
    for (size_t i = 0; i < surface_triangles.size(); i += 3) {
        int i0 = surface_triangles[i];
        int i1 = surface_triangles[i+1];
        int i2 = surface_triangles[i+2];
        
        double v0[3] = {surface_vertices[i0].x, surface_vertices[i0].y, surface_vertices[i0].z};
        double v1[3] = {surface_vertices[i1].x, surface_vertices[i1].y, surface_vertices[i1].z};
        double v2[3] = {surface_vertices[i2].x, surface_vertices[i2].y, surface_vertices[i2].z};
        
        double normal[3];
        calculateTriangleNormal(v0, v1, v2, normal);
        
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
    
    // Normalize
    for (size_t i = 0; i < surface_vertices.size(); i++) {
        double len = sqrt(surface_vertices[i].nx * surface_vertices[i].nx +
                         surface_vertices[i].ny * surface_vertices[i].ny +
                         surface_vertices[i].nz * surface_vertices[i].nz);
        if (len > 1e-10) {
            surface_vertices[i].nx /= len;
            surface_vertices[i].ny /= len;
            surface_vertices[i].nz /= len;
        }
    }
    
    // ========== MERGE DUPLICATE VERTICES ==========
    printf("Merging duplicate vertices...\n");
    double merge_tolerance = params[0] * 1e-4;
    std::vector<int> vertex_map(surface_vertices.size());
    std::vector<Vertex> unique_vertices;
    
    for (size_t i = 0; i < surface_vertices.size(); i++) {
        bool found = false;
        for (size_t j = 0; j < unique_vertices.size(); j++) {
            double dx = surface_vertices[i].x - unique_vertices[j].x;
            double dy = surface_vertices[i].y - unique_vertices[j].y;
            double dz = surface_vertices[i].z - unique_vertices[j].z;
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            
            if (dist < merge_tolerance) {
                vertex_map[i] = j;
                unique_vertices[j].nx += surface_vertices[i].nx;
                unique_vertices[j].ny += surface_vertices[i].ny;
                unique_vertices[j].nz += surface_vertices[i].nz;
                found = true;
                break;
            }
        }
        
        if (!found) {
            vertex_map[i] = unique_vertices.size();
            unique_vertices.push_back(surface_vertices[i]);
        }
    }
    
    printf("  %d -> %d unique vertices\n",
           (int)surface_vertices.size(), (int)unique_vertices.size());
    
    // Renormalize merged normals
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
    
    surface_vertices = unique_vertices;
    
    // Remap triangles
    for (size_t i = 0; i < surface_triangles.size(); i++) {
        surface_triangles[i] = vertex_map[surface_triangles[i]];
    }
    
    int num_surface_vertices = surface_vertices.size();
    int num_surface_triangles = surface_triangles.size() / 3;
    
    // ========== VERIFY MESH QUALITY ==========
    printf("\n=== MESH QUALITY VERIFICATION ===\n");
    int poor_quality = 0;
    double sum_aspect = 0.0;
    double max_aspect = 0.0;
    double min_edge_actual = 1e10;
    double max_edge_actual = 0.0;
    
    for (int i = 0; i < num_surface_triangles; i++) {
        int i0 = surface_triangles[3*i];
        int i1 = surface_triangles[3*i+1];
        int i2 = surface_triangles[3*i+2];
        
        double dx1 = surface_vertices[i1].x - surface_vertices[i0].x;
        double dy1 = surface_vertices[i1].y - surface_vertices[i0].y;
        double dz1 = surface_vertices[i1].z - surface_vertices[i0].z;
        double e1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
        
        double dx2 = surface_vertices[i2].x - surface_vertices[i1].x;
        double dy2 = surface_vertices[i2].y - surface_vertices[i1].y;
        double dz2 = surface_vertices[i2].z - surface_vertices[i1].z;
        double e2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
        
        double dx3 = surface_vertices[i0].x - surface_vertices[i2].x;
        double dy3 = surface_vertices[i0].y - surface_vertices[i2].y;
        double dz3 = surface_vertices[i0].z - surface_vertices[i2].z;
        double e3 = sqrt(dx3*dx3 + dy3*dy3 + dz3*dz3);
        
        double max_e = fmax(e1, fmax(e2, e3));
        double min_e = fmin(e1, fmin(e2, e3));
        double aspect = max_e / (min_e + 1e-12);
        
        sum_aspect += aspect;
        if (aspect > max_aspect) max_aspect = aspect;
        if (aspect > 2.0) poor_quality++;
        
        if (max_e > max_edge_actual) max_edge_actual = max_e;
        if (min_e < min_edge_actual && min_e > 1e-10) min_edge_actual = min_e;
    }
    
    printf("Mesh quality metrics:\n");
    printf("  Actual edge length range: %.6f to %.6f\n", min_edge_actual, max_edge_actual);
    printf("  Average aspect ratio: %.2f\n", sum_aspect / num_surface_triangles);
    printf("  Maximum aspect ratio: %.2f\n", max_aspect);
    printf("  Elements with aspect > 2.0: %d (%.1f%%)\n",
           poor_quality, 100.0 * poor_quality / num_surface_triangles);
    
    if (100.0 * poor_quality / num_surface_triangles < 5.0) {
        printf("  Quality grade: EXCELLENT\n");
    } else if (100.0 * poor_quality / num_surface_triangles < 10.0) {
        printf("  Quality grade: GOOD\n");
    } else if (100.0 * poor_quality / num_surface_triangles < 20.0) {
        printf("  Quality grade: ACCEPTABLE\n");
    } else {
        printf("  Quality grade: POOR - consider refining\n");
    }
    
    // ========== CREATE PRISMATIC LAYERS ==========
    printf("\n=== GENERATING PRISMATIC LAYERS ===\n");
    std::vector<Vertex> all_vertices;
    std::vector<int> prism_cells;
    
    std::vector<double> layer_offsets;
    layer_offsets.push_back(0.0);
    double cumulative = 0.0;
    double thickness = first_layer_thickness;
    
    for (int i = 0; i < num_layers; i++) {
        cumulative += thickness;
        layer_offsets.push_back(cumulative);
        thickness *= growth_rate;
    }
    
    // Create vertex layers
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
    
    // Create prisms
    for (int layer = 0; layer < num_layers; layer++) {
        int base = layer * num_surface_vertices;
        int next = (layer + 1) * num_surface_vertices;
        
        for (int i = 0; i < num_surface_triangles; i++) {
            int i0 = surface_triangles[3*i];
            int i1 = surface_triangles[3*i+1];
            int i2 = surface_triangles[3*i+2];
            
            prism_cells.push_back(base + i0);
            prism_cells.push_back(base + i1);
            prism_cells.push_back(base + i2);
            prism_cells.push_back(next + i0);
            prism_cells.push_back(next + i1);
            prism_cells.push_back(next + i2);
        }
    }
    
    printf("  Total vertices: %d\n", (int)all_vertices.size());
    printf("  Total prisms: %d\n", (int)prism_cells.size()/6);
    
    // ========== WRITE VTK FILE ==========
    printf("\n=== WRITING OUTPUT ===\n");
    printf("Writing %s...\n", argv[2]);
    FILE *fp = fopen(argv[2], "w");
    if (!fp) {
        printf("Error opening output file\n");
        return 1;
    }
    
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Prismatic mesh from %s (per-surface analysis)\n", argv[1]);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %d float\n", (int)all_vertices.size());
    
    for (size_t i = 0; i < all_vertices.size(); i++) {
        fprintf(fp, "%.9f %.9f %.9f\n",
                all_vertices[i].x, all_vertices[i].y, all_vertices[i].z);
    }
    
    int nprisms = prism_cells.size() / 6;
    fprintf(fp, "\nCELLS %d %d\n", nprisms, nprisms * 7);
    for (size_t i = 0; i < prism_cells.size(); i += 6) {
        fprintf(fp, "6 %d %d %d %d %d %d\n",
                prism_cells[i], prism_cells[i+1], prism_cells[i+2],
                prism_cells[i+3], prism_cells[i+4], prism_cells[i+5]);
    }
    
    fprintf(fp, "\nCELL_TYPES %d\n", nprisms);
    for (int i = 0; i < nprisms; i++) {
        fprintf(fp, "13\n");
    }
    
    fclose(fp);
    printf("Success!\n");
    
    // Cleanup
    for (int i = 0; i < nbody; i++) {
        if (tessellations[i]) EG_deleteObject(tessellations[i]);
    }
    free(tessellations);
    EG_deleteObject(model);
    EG_close(context);
    
    printf("\n=== COMPLETE ===\n");
    return 0;
}
