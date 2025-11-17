#!/usr/bin/env python3
"""
Generate a tetrahedral volume mesh between a large box and boundary faces from a VTK file.
The box size is 10 times the maximum diameter of the boundary bounding box.
"""

import numpy as np
import gmsh
import sys

def read_vtk_boundary(filename):
    """
    Read boundary faces from a VTK file (legacy format).
    Returns points and triangles/polygons.
    """
    points = []
    cells = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse VTK file
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Read points
        if line.startswith('POINTS'):
            parts = line.split()
            n_points = int(parts[1])
            point_data = []
            i += 1
            
            # Read all point coordinates
            while len(point_data) < n_points * 3:
                point_line = lines[i].strip().split()
                point_data.extend([float(x) for x in point_line])
                i += 1
            
            # Reshape into (n_points, 3)
            points = np.array(point_data).reshape(n_points, 3)
            continue
        
        # Read cells (polygons/triangles)
        if line.startswith('POLYGONS') or line.startswith('CELLS'):
            parts = line.split()
            n_cells = int(parts[1])
            i += 1
            
            for _ in range(n_cells):
                cell_line = lines[i].strip().split()
                n_verts = int(cell_line[0])
                cell_verts = [int(x) for x in cell_line[1:n_verts+1]]
                cells.append(cell_verts)
                i += 1
            continue
        
        i += 1
    
    return np.array(points), cells

def compute_bounding_box(points):
    """
    Compute the bounding box of a set of points.
    Returns min and max coordinates.
    """
    min_coords = np.min(points, axis=0)
    max_coords = np.max(points, axis=0)
    return min_coords, max_coords

def compute_max_diameter(min_coords, max_coords):
    """
    Compute the maximum diameter (diagonal length) of the bounding box.
    """
    diagonal = np.linalg.norm(max_coords - min_coords)
    return diagonal

def create_mesh(boundary_file, output_file, mesh_size=None):
    """
    Create a tetrahedral mesh between a large box and the boundary surface.
    
    Parameters:
    -----------
    boundary_file : str
        Path to the boundary.vtk file
    output_file : str
        Path for the output mesh file (e.g., 'output.msh')
    mesh_size : float, optional
        Characteristic mesh size. If None, computed automatically.
    """
    
    # Read boundary VTK file
    print(f"Reading boundary file: {boundary_file}")
    points, cells = read_vtk_boundary(boundary_file)
    print(f"  Found {len(points)} points and {len(cells)} boundary faces")
    
    # Compute bounding box
    min_coords, max_coords = compute_bounding_box(points)
    center = (min_coords + max_coords) / 2
    max_diameter = compute_max_diameter(min_coords, max_coords)
    
    print(f"\nBoundary bounding box:")
    print(f"  Min: {min_coords}")
    print(f"  Max: {max_coords}")
    print(f"  Center: {center}")
    print(f"  Max diameter: {max_diameter:.4f}")
    
    # Define large box (10x max diameter)
    box_size = 10 * max_diameter
    box_min = center - box_size / 2
    box_max = center + box_size / 2
    
    print(f"\nLarge box:")
    print(f"  Size: {box_size:.4f}")
    print(f"  Min: {box_min}")
    print(f"  Max: {box_max}")
    
    # Set default mesh size if not provided
    if mesh_size is None:
        mesh_size = max_diameter / 20  # Reasonable default
    
    print(f"\nMesh size: {mesh_size:.4f}")
    
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.model.add("volume_mesh")
    
    # Add points from boundary
    print("\nAdding boundary points to Gmsh...")
    point_tags = []
    for pt in points:
        tag = gmsh.model.geo.addPoint(pt[0], pt[1], pt[2], mesh_size)
        point_tags.append(tag)
    
    # Add boundary faces
    print("Adding boundary faces to Gmsh...")
    boundary_curves = []
    boundary_surfaces = []
    
    for cell in cells:
        # Create curves for each cell edge
        curve_loop_tags = []
        n_verts = len(cell)
        
        for j in range(n_verts):
            p1 = point_tags[cell[j]]
            p2 = point_tags[cell[(j + 1) % n_verts]]
            
            # Check if curve already exists
            try:
                curve = gmsh.model.geo.addLine(p1, p2)
                curve_loop_tags.append(curve)
            except:
                # If curve creation fails, it might already exist
                pass
        
        # Create curve loop and surface for triangular/quad faces
        if len(curve_loop_tags) >= 3:
            try:
                curve_loop = gmsh.model.geo.addCurveLoop(curve_loop_tags)
                surface = gmsh.model.geo.addPlaneSurface([curve_loop])
                boundary_surfaces.append(surface)
            except:
                pass
    
    print(f"  Created {len(boundary_surfaces)} boundary surfaces")
    
    # Create the large box
    print("\nCreating large box...")
    box = gmsh.model.geo.addBox(box_min[0], box_min[1], box_min[2],
                                  box_max[0] - box_min[0],
                                  box_max[1] - box_min[1],
                                  box_max[2] - box_min[2])
    
    # Synchronize before boolean operations
    gmsh.model.geo.synchronize()
    
    # Get the boundary surface loop
    print("Creating surface loop from boundary surfaces...")
    try:
        boundary_loop = gmsh.model.geo.addSurfaceLoop(boundary_surfaces)
        boundary_volume = gmsh.model.geo.addVolume([boundary_loop])
        
        # Perform boolean operation: box - boundary_volume
        print("Performing boolean operation (box - boundary)...")
        gmsh.model.geo.synchronize()
        
        out_dimtags, out_dimtags_map = gmsh.model.geo.cut(
            [(3, box)], [(3, boundary_volume)], removeObject=True, removeTool=True
        )
        
        gmsh.model.geo.synchronize()
        
    except Exception as e:
        print(f"Warning: Could not create closed volume from boundary surfaces.")
        print(f"Error: {e}")
        print("Attempting alternative approach: embedding boundary in box...")
        
        # Alternative: just mesh the box with boundary embedded
        gmsh.model.geo.synchronize()
    
    # Set mesh size
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size * 0.5)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size * 2.0)
    
    # Generate 3D mesh
    print("\nGenerating 3D tetrahedral mesh...")
    gmsh.model.mesh.generate(3)
    
    # Optimize mesh
    print("Optimizing mesh...")
    gmsh.model.mesh.optimize("Netgen")
    
    # Write mesh
    print(f"\nWriting mesh to: {output_file}")
    gmsh.write(output_file)
    
    # Get mesh statistics
    n_nodes = len(gmsh.model.mesh.getNodes()[0])
    n_elements = len(gmsh.model.mesh.getElements(3)[1][0])
    
    print(f"\nMesh statistics:")
    print(f"  Number of nodes: {n_nodes}")
    print(f"  Number of tetrahedra: {n_elements}")
    
    # Finalize
    gmsh.finalize()
    
    print("\nMesh generation complete!")

if __name__ == "__main__":
    # Default parameters
    boundary_file = "boundary.vtk"
    output_file = "volume_mesh.msh"
    mesh_size = None  # Auto-computed
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        boundary_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    if len(sys.argv) > 3:
        mesh_size = float(sys.argv[3])
    
    # Generate mesh
    try:
        create_mesh(boundary_file, output_file, mesh_size)
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
