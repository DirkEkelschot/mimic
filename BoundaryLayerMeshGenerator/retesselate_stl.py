#!/usr/bin/env python3
"""
Script to retesselate (remesh) surface geometry from an STL file using Gmsh.

This script properly remeshes the surface with a new triangulation, not just
copying the original mesh.

Usage:
    python retesselate_stl.py input.stl output.vtk [target_edge_length]
    
Arguments:
    input.stl          - Input STL file
    output.vtk         - Output mesh file (VTK format)
    target_edge_length - Optional target edge length (default: auto-detect)
"""

import sys
import gmsh
import math

def retesselate_stl(input_stl, output_file, target_edge_length=None):
    """
    Retesselate surface geometry from STL file with new triangulation.
    
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
    
    # Get all entities and compute bounding box
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
    
    # IMPORTANT: Remove the original mesh to force remeshing
    print("Removing original mesh...")
    gmsh.model.mesh.clear()
    
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
    
    # Force remeshing by setting mesh size constraints
    # Set minimum and maximum mesh sizes
    gmsh.option.setNumber("Mesh.MeshSizeMin", target_edge_length * 0.5)
    gmsh.option.setNumber("Mesh.MeshSizeMax", target_edge_length * 2.0)
    
    # Set mesh size for all surfaces
    for dim, tag in surfaces:
        gmsh.model.mesh.setSize([(dim, tag)], target_edge_length)
    
    # Generate 2D mesh (surface mesh) - this will create a new triangulation
    print("Generating new surface mesh (retesselating)...")
    gmsh.model.mesh.generate(2)
    
    # Optional: Optimize the mesh
    print("Optimizing mesh...")
    gmsh.model.mesh.optimize("Netgen")
    
    # Get mesh statistics
    nodeTags, coords, _ = gmsh.model.mesh.getNodes()
    num_nodes = len(nodeTags)
    
    # Get 2D elements (triangles)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2)
    num_elements = 0
    if len(elemTags) > 0 and len(elemTags[0]) > 0:
        num_elements = len(elemTags[0])
    
    print(f"New mesh generated: {num_nodes} nodes, {num_elements} triangles")
    
    # Save mesh
    print(f"Saving remeshed surface to: {output_file}")
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
        retesselate_stl(input_stl, output_file, target_edge_length)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

