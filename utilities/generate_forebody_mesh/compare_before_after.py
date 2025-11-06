"""
Compare mesh deviation BEFORE and AFTER projection to CAD surface.
This version writes the mesh, modifies coordinates in memory, and re-imports.
"""
import gmsh
import sys
import numpy as np
import matplotlib.pyplot as plt



def project_to_nearest_surface(point, surfaces_list):
    best_distance = float('inf')
    best_point = None
    
    for surf_dim, surf_tag in surfaces_list:
        try:
            closest_pt, uv = gmsh.model.getClosestPoint(surf_dim, surf_tag, point)
            distance = np.linalg.norm(np.array(closest_pt) - np.array(point))
            
            if distance < best_distance:
                best_distance = distance
                best_point = closest_pt
        except:
            pass
    
    return best_point


def project_line(start_pt,end_pt,all_surfaces,n_steps_line):
    line_pts = []
    p0_tag = gmsh.model.occ.addPoint(start_pt[0], start_pt[1], start_pt[2])
    p1_tag = gmsh.model.occ.addPoint(end_pt[0]  , end_pt[1]  , end_pt[2])
    line_pts.append(p0_tag)
    current_pt = start_pt.copy()

    for step in range(1, n_steps_line):
        direction = end_pt - current_pt
        direction_norm = np.linalg.norm(direction)
        
        if direction_norm < 1e-6:
            break
        
        direction = direction / direction_norm
        step_size = direction_norm / (n_steps_line - step + 1)
        next_guess = current_pt + direction * step_size
        
        closest_pt = project_to_nearest_surface(next_guess.tolist(), all_surfaces)
        
        if closest_pt is not None:
            dist = np.linalg.norm(np.array(closest_pt) - current_pt)
            if dist > 1e-6:
                current_pt = np.array(closest_pt)
                pt = gmsh.model.occ.addPoint(closest_pt[0], closest_pt[1], closest_pt[2])
                line_pts.append(pt)

    line_pts.append(p1_tag)
    return line_pts


def circle_line(n_points, x_plane, R, start_alpha, alpha, startc, endc):
    circle_points = []
    delta = alpha/(n_points-1)
    theta = start_alpha
    
    for i in range(n_points):
        x = x_plane
        y = R * np.cos(theta)
        z = R * np.sin(theta)
        
        if i == 0:
            pt_tag = gmsh.model.occ.addPoint(startc[0], startc[1], startc[2])
        elif i == n_points-1:
            pt_tag = gmsh.model.occ.addPoint(endc[0], endc[1], endc[2])
        else:
            pt_tag = gmsh.model.occ.addPoint(x, y, z)
        circle_points.append(pt_tag)
        theta = theta + delta

    return circle_points


def generateTransfiniteSurfaceMesh(curves,nmesh0,nmesh1):
    wire = gmsh.model.occ.addWire(curves, checkClosed=True)
    stag = gmsh.model.occ.addSurfaceFilling(wire)
    gmsh.model.occ.synchronize()
    boundary = gmsh.model.getBoundary([(2, stag)], oriented=True, recursive=False)
    boundary_curve_tags = [abs(c[1]) for c in boundary]
    
    if len(boundary_curve_tags) == 4:
        gmsh.model.mesh.setTransfiniteCurve(boundary_curve_tags[0], nmesh0)
        gmsh.model.mesh.setTransfiniteCurve(boundary_curve_tags[1], nmesh1)
        gmsh.model.mesh.setTransfiniteCurve(boundary_curve_tags[2], nmesh0)
        gmsh.model.mesh.setTransfiniteCurve(boundary_curve_tags[3], nmesh1)
    
    gmsh.model.mesh.setTransfiniteSurface(stag)
    gmsh.model.mesh.setRecombine(2, stag)

    return stag


def measure_deviation(surface_tags, cad_surfaces):
    """Measure deviation from CAD"""
    all_distances = []
    interior_distances = []
    
    for surf_tag in surface_tags:
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes(2, surf_tag, includeBoundary=True)
        
        # Get boundary nodes
        boundary_entities = gmsh.model.getBoundary([(2, surf_tag)], oriented=False, recursive=False)
        boundary_node_set = set()
        for dim, tag in boundary_entities:
            bnodes, _, _ = gmsh.model.mesh.getNodes(abs(dim), abs(tag))
            boundary_node_set.update(bnodes)
        
        for i, node_tag in enumerate(node_tags):
            x = node_coords[3*i]
            y = node_coords[3*i + 1]
            z = node_coords[3*i + 2]
            
            projected = project_to_nearest_surface([x, y, z], cad_surfaces)
            
            if projected is not None:
                distance = np.linalg.norm(np.array([x,y,z]) - np.array(projected))
                all_distances.append(distance)
                
                if node_tag not in boundary_node_set:
                    interior_distances.append(distance)
    
    return np.array(all_distances), np.array(interior_distances)


def project_mesh_nodes_to_cad(surface_tags, cad_surfaces):
    """
    Project mesh nodes to CAD by modifying them in-place.
    Uses the lower-level node access API.
    """
    print("\n🔧 Projecting mesh nodes to CAD surface...")
    
    total_nodes = 0
    total_moved = 0
    
    # Collect all nodes that need to be moved
    all_node_updates = {}  # node_tag -> (x, y, z)
    
    for surf_tag in surface_tags:
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes(2, surf_tag, includeBoundary=True)
        total_nodes += len(node_tags)
        
        for i, node_tag in enumerate(node_tags):
            x = node_coords[3*i]
            y = node_coords[3*i + 1]
            z = node_coords[3*i + 2]
            
            # Find closest point on CAD
            projected = project_to_nearest_surface([x, y, z], cad_surfaces)
            
            if projected is not None:
                distance = np.linalg.norm(np.array([x, y, z]) - np.array(projected))
                if distance > 1e-10:
                    all_node_updates[node_tag] = projected
                    total_moved += 1
    
    # Now update all nodes
    # We need to use setNode (singular) for each node
    print(f"  Updating {total_moved}/{total_nodes} nodes...")
    
    for node_tag, (px, py, pz) in all_node_updates.items():
        try:
            # setNode(tag, x, y, z) - singular, updates one node
            gmsh.model.mesh.setNode(int(node_tag), [px, py, pz], [])
        except Exception as e:
            print(f"    Warning: Could not update node {node_tag}: {e}")
            pass
    
    print(f"✓ Projected {total_moved} nodes to CAD surface")


if __name__ == "__main__":
    
    gmsh.initialize()
    gmsh.model.add("comparison")

    iges_file = "your_model.iges"
    gmsh.merge(iges_file)
    gmsh.model.occ.synchronize()
    surfaces = gmsh.model.getEntities(2)
    all_surfaces = surfaces.copy()
    
    # Get bounds
    x_min = y_min = z_min = float('inf')
    x_max = y_max = z_max = float('-inf')

    for surf_dim, surf_tag in surfaces:
        bbox = gmsh.model.occ.getBoundingBox(surf_dim, surf_tag)
        x_min = min(x_min, bbox[0])
        y_min = min(y_min, bbox[1])
        z_min = min(z_min, bbox[2])
        x_max = max(x_max, bbox[3])
        y_max = max(y_max, bbox[4])
        z_max = max(z_max, bbox[5])

    gmsh.model.occ.synchronize()

    square_coords = [
        np.array([0, 0, 0]),
        np.array([0, 0, 1]),
        np.array([0, 1, 1]),
        np.array([0, 1, 0])
    ]

    R = y_max
    theta  = 45*np.pi/180
    twotheta = 2*theta
    
    top = [x_max,y_max,0]
    p45 = [x_max,R*np.cos(theta), R*np.sin(theta)]
    p90 = [x_max,R*np.cos(twotheta), R*np.sin(twotheta)]
    
    n_steps_line   = 50
    n_mesh0        = 50
    n_mesh1        = 50
    n_mesh2        = 50
    
    p0_phys = project_to_nearest_surface(square_coords[0], all_surfaces)
    p1_phys = project_to_nearest_surface(square_coords[1], all_surfaces)
    p2_phys = project_to_nearest_surface(square_coords[2], all_surfaces)
    p3_phys = project_to_nearest_surface(square_coords[3], all_surfaces)
    p4_phys = project_to_nearest_surface(top, all_surfaces)
    p5_phys = project_to_nearest_surface(p45, all_surfaces)
    p6_phys = project_to_nearest_surface(p90, all_surfaces)
    
    # Project lines
    line0 = project_line(np.array(p0_phys), np.array(p1_phys), all_surfaces, n_steps_line)
    line1 = project_line(np.array(p1_phys), np.array(p2_phys), all_surfaces, n_steps_line)
    line2 = project_line(np.array(p2_phys), np.array(p3_phys), all_surfaces, n_steps_line)
    line3 = project_line(np.array(p3_phys), np.array(p0_phys), all_surfaces, n_steps_line)
    line4 = project_line(np.array(p3_phys), np.array(p4_phys), all_surfaces, n_steps_line)
    line5 = circle_line(n_steps_line, x_max, R, 0.0, theta, np.array(p4_phys), np.array(p5_phys))
    line6 = project_line(np.array(p5_phys), np.array(p2_phys), all_surfaces, n_steps_line)
    line7 = circle_line(n_steps_line, x_max, R, theta, theta, np.array(p5_phys), np.array(p6_phys))
    line8 = project_line(np.array(p6_phys), np.array(p1_phys), all_surfaces, n_steps_line)
    
    line0_spline = gmsh.model.occ.addBSpline(line0)
    line1_spline = gmsh.model.occ.addBSpline(line1)
    line2_spline = gmsh.model.occ.addBSpline(line2)
    line3_spline = gmsh.model.occ.addBSpline(line3)
    line4_spline = gmsh.model.occ.addBSpline(line4)
    line5_spline = gmsh.model.occ.addBSpline(line5)
    line6_spline = gmsh.model.occ.addBSpline(line6)
    line7_spline = gmsh.model.occ.addBSpline(line7)
    line8_spline = gmsh.model.occ.addBSpline(line8)
    
    gmsh.model.occ.synchronize()
    
    # Set transfinite
    gmsh.model.mesh.setTransfiniteCurve(line0_spline, n_mesh0)
    gmsh.model.mesh.setTransfiniteCurve(line1_spline, n_mesh1)
    gmsh.model.mesh.setTransfiniteCurve(line2_spline, n_mesh0)
    gmsh.model.mesh.setTransfiniteCurve(line3_spline, n_mesh1)
    gmsh.model.mesh.setTransfiniteCurve(line4_spline, n_mesh2)
    gmsh.model.mesh.setTransfiniteCurve(line5_spline, n_mesh1)
    gmsh.model.mesh.setTransfiniteCurve(line6_spline, n_mesh2)
    gmsh.model.mesh.setTransfiniteCurve(line7_spline, n_mesh1)
    gmsh.model.mesh.setTransfiniteCurve(line8_spline, n_mesh2)
    
    # Create surfaces
    curves  = [line0_spline, line1_spline, line2_spline, line3_spline]
    surface_tag = generateTransfiniteSurfaceMesh(curves, n_mesh0, n_mesh1)

    curves2 = [line4_spline, line5_spline, line6_spline, line2_spline]
    surface_tag2 = generateTransfiniteSurfaceMesh(curves2, n_mesh0, n_mesh1)
    
    curves3 = [line7_spline, line8_spline, line1_spline, line6_spline]
    surface_tag3 = generateTransfiniteSurfaceMesh(curves3, n_mesh1, n_mesh2)
    
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setCompound(2, [surface_tag, surface_tag2, surface_tag3])
    
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Smoothing", 100)

    gmsh.model.addPhysicalGroup(2, [surface_tag, surface_tag2, surface_tag3], 1)
    gmsh.model.setPhysicalName(2, 1, "StructuredPatch")

    # Generate mesh
    print("\n🔧 Generating mesh...")
    gmsh.model.mesh.generate(2)
    
    print("\n" + "="*70)
    print("BEFORE PROJECTION")
    print("="*70)
    
    # Measure BEFORE projection
    all_before, interior_before = measure_deviation([surface_tag, surface_tag2, surface_tag3], all_surfaces)
    
    print(f"Interior nodes:")
    print(f"  Mean deviation:   {np.mean(interior_before):.6e}")
    print(f"  Max deviation:    {np.max(interior_before):.6e}")
    print(f"  Median deviation: {np.median(interior_before):.6e}")
    
    gmsh.write("mesh_BEFORE_projection.vtk")
    print("✓ Saved: mesh_BEFORE_projection.vtk")
    
    # NOW PROJECT
    project_mesh_nodes_to_cad([surface_tag, surface_tag2, surface_tag3], all_surfaces)
    
    print("\n" + "="*70)
    print("AFTER PROJECTION")
    print("="*70)
    
    # Measure AFTER projection
    all_after, interior_after = measure_deviation([surface_tag, surface_tag2, surface_tag3], all_surfaces)
    
    print(f"Interior nodes:")
    print(f"  Mean deviation:   {np.mean(interior_after):.6e}")
    print(f"  Max deviation:    {np.max(interior_after):.6e}")
    print(f"  Median deviation: {np.median(interior_after):.6e}")
    
    gmsh.write("mesh_AFTER_projection.vtk")
    print("✓ Saved: mesh_AFTER_projection.vtk")
    
    # Calculate improvement
    print("\n" + "="*70)
    print("IMPROVEMENT")
    print("="*70)
    
    improvement_mean = (1 - np.mean(interior_after)/np.mean(interior_before)) * 100
    improvement_max = (1 - np.max(interior_after)/np.max(interior_before)) * 100
    
    print(f"Mean deviation reduced by:  {improvement_mean:.1f}%")
    print(f"Max deviation reduced by:   {improvement_max:.1f}%")
    print(f"\nBefore max: {np.max(interior_before):.6e}")
    print(f"After max:  {np.max(interior_after):.6e}")
    print(f"Ratio:      {np.max(interior_after)/np.max(interior_before):.2e}x")
    
    # Create comparison plot
    plt.figure(figsize=(14, 5))
    
    plt.subplot(1, 3, 1)
    plt.hist(interior_before, bins=50, alpha=0.7, color='red', edgecolor='black', label='Before')
    plt.xlabel('Distance from CAD')
    plt.ylabel('Number of nodes')
    plt.title('BEFORE Projection')
    plt.yscale('log')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.subplot(1, 3, 2)
    plt.hist(interior_after, bins=50, alpha=0.7, color='green', edgecolor='black', label='After')
    plt.xlabel('Distance from CAD')
    plt.ylabel('Number of nodes')
    plt.title('AFTER Projection')
    plt.yscale('log')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.subplot(1, 3, 3)
    plt.hist(interior_before, bins=50, alpha=0.5, color='red', edgecolor='black', label='Before')
    plt.hist(interior_after, bins=50, alpha=0.5, color='green', edgecolor='black', label='After')
    plt.xlabel('Distance from CAD')
    plt.ylabel('Number of nodes')
    plt.title('Comparison')
    plt.yscale('log')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('before_after_comparison.png', dpi=150)
    print(f"\n✓ Saved: before_after_comparison.png")
    
    gmsh.finalize()
    
    print("\n" + "="*70)
    print("FILES CREATED")
    print("="*70)
    print("1. mesh_BEFORE_projection.vtk")
    print("2. mesh_AFTER_projection.vtk")
    print("3. before_after_comparison.png")
    print("\n💡 Load both VTK files in ParaView side-by-side!")
    print("   Color by 'Quality' or use 'Curvature' filter to see differences.")
