import pyvista as pv
import numpy as np

# ============================================
# CONFIGURATION PARAMETERS
# ============================================
FARFIELD_DISTANCE = 50.0        # Distance from boundary to farfield
FARFIELD_REFINEMENT = 2         # Subdivision level (0-5)
                                # 0 = coarse (12 faces)
                                # 1 = fine (48 faces)  
                                # 2 = finer (192 faces)
                                # 3 = very fine (768 faces)
                                # 4 = extremely fine (3072 faces)

# Alternative: Use target edge length instead of subdivision
USE_TARGET_EDGE_LENGTH = False  # Set to True to use edge length instead
TARGET_EDGE_LENGTH = 30.0       # Target edge length for farfield faces

print("="*60)
print("VOLUME MESH GENERATION WITH REFINED FARFIELD")
print("="*60)

# Read wall surface
wall_surf = pv.read("boundary.vtk")

# If it's an UnstructuredGrid, extract the surface first
if isinstance(wall_surf, pv.UnstructuredGrid):
    print("Converting UnstructuredGrid to PolyData surface...")
    wall_surf = wall_surf.extract_surface()

# Create farfield box
bounds = wall_surf.bounds
farfield_box = pv.Box(bounds=[
    bounds[0] - FARFIELD_DISTANCE, bounds[1] + FARFIELD_DISTANCE,
    bounds[2] - FARFIELD_DISTANCE, bounds[3] + FARFIELD_DISTANCE,
    bounds[4] - FARFIELD_DISTANCE, bounds[5] + FARFIELD_DISTANCE
])

max_dist = np.sqrt((bounds[0]-bounds[1])*(bounds[0]-bounds[1])
                  +(bounds[2]-bounds[3])*(bounds[2]-bounds[3])
                  +(bounds[4]-bounds[5])*(bounds[4]-bounds[5]))

print(f"Max distance: {max_dist}")




# Check if wall surface is already triangulated
def is_triangulated(mesh):
    """Check if mesh contains only triangular faces."""
    if mesh.n_cells == 0:
        return False
    for i in range(mesh.n_cells):
        cell = mesh.get_cell(i)
        if cell.n_points != 3:
            return False
    return True

# Only triangulate if necessary
if not is_triangulated(wall_surf):
    print("Warning: wall_surf contains non-triangular faces. Triangulating...")
    wall_surf = wall_surf.triangulate()
else:
    print("Wall surface is already triangulated. Preserving original tessellation.")

# ============================================
# REFINE FARFIELD BOX TESSELLATION
# ============================================
print("\n=== Refining Farfield Box ===")
print(f"Farfield distance: {FARFIELD_DISTANCE}")

# Triangulate the box first
farfield_box = farfield_box.triangulate()
print(f"Initial farfield: {farfield_box.n_cells} faces")

if USE_TARGET_EDGE_LENGTH:
    # Method A: Adaptive refinement based on edge length
    print(f"Using target edge length: {TARGET_EDGE_LENGTH}")
    
    # Subdivide until all edges are smaller than target
    farfield_refined = farfield_box
    max_iterations = 5
    for iteration in range(max_iterations):
        edges = farfield_refined.extract_all_edges()
        edge_lengths = edges.compute_cell_sizes()['Length']
        max_edge = edge_lengths.max()
        
        if max_edge <= TARGET_EDGE_LENGTH:
            print(f"  Converged after {iteration} iterations")
            break
        
        farfield_refined = farfield_refined.subdivide(nsub=1, subfilter='butterfly')
        print(f"  Iteration {iteration+1}: max edge = {max_edge:.1f}, faces = {farfield_refined.n_cells}")
    
    farfield_box_refined = farfield_refined
else:
    # Method B: Fixed subdivision level
    if FARFIELD_REFINEMENT > 0:
        print(f"Subdivision level: {FARFIELD_REFINEMENT}")
        farfield_box_refined = farfield_box.subdivide(
            nsub=FARFIELD_REFINEMENT, 
            subfilter='butterfly'
        )
    else:
        farfield_box_refined = farfield_box
    
print(f"Final farfield: {farfield_box_refined.n_points} points, {farfield_box_refined.n_cells} faces")

# Calculate average edge length
edges = farfield_box_refined.extract_all_edges()
avg_edge = edges.compute_cell_sizes()['Length'].mean()
print(f"Average farfield edge length: {avg_edge:.1f}")

# ============================================
# FIND A POINT TRULY INSIDE THE SURFACE
# ============================================
print("\n=== Finding hole point ===")

# Create a grid of test points inside the bounding box
n_samples = 30
test_grid = pv.ImageData(
    dimensions=(n_samples, n_samples, n_samples),
    spacing=(
        (bounds[1]-bounds[0])/(n_samples-1),
        (bounds[3]-bounds[2])/(n_samples-1),
        (bounds[5]-bounds[4])/(n_samples-1)
    ),
    origin=(bounds[0], bounds[2], bounds[4])
)

# Use PyVista's select_enclosed_points to find points inside the surface
selection = test_grid.select_enclosed_points(wall_surf, tolerance=0.0, check_surface=False)
inside_mask = selection['SelectedPoints'].view(bool)

if inside_mask.sum() > 0:
    inside_points = test_grid.points[inside_mask]
    hole_point = inside_points[len(inside_points)//2]
    print(f"Hole point: {hole_point}")
else:
    print("WARNING: Could not find points inside surface, using centroid")
    hole_point = wall_surf.center

# ============================================
# EXTRACT POINTS AND FACES
# ============================================
def extract_faces(mesh):
    """Extract triangular faces from a mesh."""
    faces = mesh.faces
    if faces is None or len(faces) == 0:
        face_list = []
        for i in range(mesh.n_cells):
            cell = mesh.get_cell(i)
            if cell.n_points == 3:
                face_list.append([cell.point_ids[0], cell.point_ids[1], cell.point_ids[2]])
        return np.array(face_list) if face_list else np.array([])
    else:
        return faces.reshape(-1, 4)[:, 1:]

farfield_points = farfield_box_refined.points
wall_points = wall_surf.points

farfield_faces = extract_faces(farfield_box_refined)
wall_faces = extract_faces(wall_surf)

# Combine points
points = np.vstack([farfield_points, wall_points])

# Adjust wall face indices
n_farfield_points = len(farfield_points)
wall_faces_adjusted = wall_faces + n_farfield_points

all_faces = np.vstack([farfield_faces, wall_faces_adjusted])
n_farfield_faces = len(farfield_faces)
n_wall_faces = len(wall_faces)

print(f"\nFarfield: {len(farfield_points)} points, {n_farfield_faces} faces")
print(f"Wall: {len(wall_points)} points, {n_wall_faces} faces")
print(f"Combined: {len(points)} points, {len(all_faces)} faces")

# ============================================
# WRITE .POLY FILE
# ============================================
poly_filename = 'flow_domain.poly'
print(f"\n=== Writing {poly_filename} ===")

with open(poly_filename, 'w') as f:
    f.write(f"{len(points)} 3 0 1\n")
    
    for i, pt in enumerate(points, 1):
        f.write(f"{i} {pt[0]:.10f} {pt[1]:.10f} {pt[2]:.10f}\n")
    
    f.write(f"{len(all_faces)} 1\n")
    
    for i in range(n_farfield_faces):
        face = farfield_faces[i]
        f.write("1 0 1\n")
        f.write(f"3 {face[0]+1} {face[1]+1} {face[2]+1}\n")
    
    for i in range(n_wall_faces):
        face = wall_faces_adjusted[i]
        f.write("1 0 2\n")
        f.write(f"3 {face[0]+1} {face[1]+1} {face[2]+1}\n")
    
    f.write("1\n")
    f.write(f"1 {hole_point[0]:.10f} {hole_point[1]:.10f} {hole_point[2]:.10f}\n")

print(f"✓ Wrote {poly_filename}")

# Save refined farfield for inspection
farfield_box_refined.save('farfield_refined.vtk')
print("✓ Saved farfield_refined.vtk for inspection")

print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Configuration:")
print(f"  Farfield distance: {FARFIELD_DISTANCE}")
if USE_TARGET_EDGE_LENGTH:
    print(f"  Refinement method: Target edge length = {TARGET_EDGE_LENGTH}")
else:
    print(f"  Refinement method: Subdivision level = {FARFIELD_REFINEMENT}")
print(f"  Farfield faces: {n_farfield_faces}")
print(f"  Wall faces: {n_wall_faces} (preserved from original)")
print(f"\nNext steps:")
print(f"  tetgen -pq1.4/15YAa5.0V {poly_filename}")
print("="*60)