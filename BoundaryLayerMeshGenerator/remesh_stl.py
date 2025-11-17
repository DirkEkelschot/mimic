#!/usr/bin/env python3
"""
Script to generate an isotropic surface mesh from an STL file using Gmsh.

Usage:
    python remesh_stl.py input.stl output.msh [target_edge_length]
    
Arguments:
    input.stl          - Input STL file
    output.msh         - Output mesh file (can be .msh, .vtk, .stl, etc.)
    target_edge_length - Optional target edge length (default: auto-detect)
"""

import sys
import gmsh
import math

def compute_bounding_box(stl_file):
    """Compute bounding box of STL file to estimate appropriate mesh size."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    
    # Load STL temporarily to get bounding box
    gmsh.open(stl_file)
    
    # Get all entities and compute overall bounding box
    entities = gmsh.model.getEntities()
    if not entities:
        gmsh.finalize()
        return 0.0, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    
    # Get bounding box of first entity
    dim, tag = entities[0]
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    
    # Expand to include all entities
    for dim, tag in entities[1:]:
        x1, y1, z1, x2, y2, z2 = gmsh.model.getBoundingBox(dim, tag)
        xmin = min(xmin, x1)
        ymin = min(ymin, y1)
        zmin = min(zmin, z1)
        xmax = max(xmax, x2)
        ymax = max(ymax, y2)
        zmax = max(zmax, z2)
    
    gmsh.finalize()
    
    diagonal = math.sqrt((xmax-xmin)**2 + (ymax-ymin)**2 + (zmax-zmin)**2)
    return diagonal, (xmin, ymin, zmin, xmax, ymax, zmax)

def remesh_stl(input_stl, output_file, target_edge_length=None):
    """
    Generate isotropic surface mesh from STL file.
    
    Args:
        input_stl: Path to input STL file
        output_file: Path to output mesh file
        target_edge_length: Target edge length (None for auto-detect)
    """
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    
    # Load STL file
    print(f"Loading STL file: {input_stl}")
    gmsh.open(input_stl)
    
    # Get bounding box
    entities = gmsh.model.getEntities()
    if not entities:
        print("Error: No entities found in STL file")
        gmsh.finalize()
        return
    
    # Get bounding box of first entity
    dim, tag = entities[0]
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    
    # Expand to include all entities
    for dim, tag in entities[1:]:
        x1, y1, z1, x2, y2, z2 = gmsh.model.getBoundingBox(dim, tag)
        xmin = min(xmin, x1)
        ymin = min(ymin, y1)
        zmin = min(zmin, z1)
        xmax = max(xmax, x2)
        ymax = max(ymax, y2)
        zmax = max(zmax, z2)
    
    diagonal = math.sqrt((xmax-xmin)**2 + (ymax-ymin)**2 + (zmax-zmin)**2)
    
    print(f"Bounding box: [{xmin:.6f}, {xmax:.6f}] x [{ymin:.6f}, {ymax:.6f}] x [{zmin:.6f}, {zmax:.6f}]")
    print(f"Diagonal: {diagonal:.6f}")
    
    # Auto-detect target edge length if not provided
    if target_edge_length is None:
        # Use 2% of diagonal as default edge length
        target_edge_length = diagonal * 0.02
        print(f"Auto-detected target edge length: {target_edge_length:.6f} (2% of diagonal)")
    else:
        print(f"Using specified target edge length: {target_edge_length:.6f}")
    
    # Get all surfaces
    surfaces = gmsh.model.getEntities(2)
    print(f"Found {len(surfaces)} surface(s)")
    
    # Set mesh size field for isotropic meshing
    # Create a simple size field based on target edge length
    gmsh.model.mesh.field.add("MathEval", 1)
    gmsh.model.mesh.field.setString(1, "F", str(target_edge_length))
    gmsh.model.mesh.field.setAsBackgroundMesh(1)
    
    # Set mesh options for isotropic surface meshing
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay for 2D
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # HXT for 3D (if needed)
    gmsh.option.setNumber("Mesh.Smoothing", 1)  # Enable smoothing
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 0)  # Linear elements
    
    # Generate 2D mesh (surface mesh)
    print("Generating surface mesh...")
    gmsh.model.mesh.generate(2)
    
    # Get mesh statistics
    num_nodes = gmsh.model.mesh.getNodes()[0].shape[0]
    num_elements = len(gmsh.model.mesh.getElements(2)[1][0]) if len(gmsh.model.mesh.getElements(2)[1]) > 0 else 0
    
    print(f"Mesh generated: {num_nodes} nodes, {num_elements} triangles")
    
    # Save mesh
    print(f"Saving mesh to: {output_file}")
    gmsh.write(output_file)
    
    # Finalize
    gmsh.finalize()
    print("Done!")

def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    
    input_stl = sys.argv[1]
    output_file = sys.argv[2]
    target_edge_length = None
    
    if len(sys.argv) >= 4:
        try:
            target_edge_length = float(sys.argv[3])
        except ValueError:
            print(f"Error: Invalid target edge length: {sys.argv[3]}")
            sys.exit(1)
    
    try:
        remesh_stl(input_stl, output_file, target_edge_length)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

