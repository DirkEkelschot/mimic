import gmsh
import sys
import numpy as np



# Function to project a point to the nearest surface
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

        # theta = i / n_points * 0.25 * np.pi   # Full circle (don't include last point to avoid duplicate)
        
        x = x_plane
        y = R * np.cos(theta)
        z = R * np.sin(theta)
        print(theta*180/np.pi)
        
        # This point should already be on the sphere surface (or very close)
        
        print(x,y,z)
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
    
    print(f"\nSurface {stag} boundary curves: {boundary_curve_tags}")
    print(f"Original input curves: {curves}")
    
    # Set transfinite on the actual boundary curves
    # The boundary should have 4 curves corresponding to our wire
    if len(boundary_curve_tags) == 4:
        # Match boundary curves to original curves based on order
        # The wire was created with order: [line0, line1, line2, line3]
        # So boundary should be in same order
        gmsh.model.mesh.setTransfiniteCurve(boundary_curve_tags[0], n_mesh0)  # corresponds to line0
        gmsh.model.mesh.setTransfiniteCurve(boundary_curve_tags[1], n_mesh1)  # corresponds to line1
        gmsh.model.mesh.setTransfiniteCurve(boundary_curve_tags[2], n_mesh0)  # corresponds to line2
        gmsh.model.mesh.setTransfiniteCurve(boundary_curve_tags[3], n_mesh1)  # corresponds to line3
        
        print(f"Set transfinite: curve {boundary_curve_tags[0]} -> {n_mesh0} nodes")
        print(f"Set transfinite: curve {boundary_curve_tags[1]} -> {n_mesh1} nodes")
        print(f"Set transfinite: curve {boundary_curve_tags[2]} -> {n_mesh0} nodes")
        print(f"Set transfinite: curve {boundary_curve_tags[3]} -> {n_mesh1} nodes")
    
    # Set transfinite surface and recombine
    gmsh.model.mesh.setTransfiniteSurface(stag)
    gmsh.model.mesh.setRecombine(2, stag)

    return stag


if __name__ == "__main__":
    
    gmsh.initialize()
    gmsh.model.add("projection")

    iges_file = "your_model.iges"
    gmsh.merge(iges_file)
    gmsh.model.occ.synchronize()
    surfaces = gmsh.model.getEntities(2)
    all_surfaces = surfaces.copy()
    # Compute overall bounds
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

    print(f"X: [{x_min}, {x_max}]")
    print(f"Y: [{y_min}, {y_max}]")
    print(f"Z: [{z_min}, {z_max}]")

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
    
    n_steps_line  = 10
    n_mesh0        = 10
    n_mesh1        = 10
    n_mesh2        = 10
    
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
    
    # Set transfinite curves on input splines
    
    # Surface 1
    gmsh.model.mesh.setTransfiniteCurve(line0_spline, n_mesh0)
    gmsh.model.mesh.setTransfiniteCurve(line1_spline, n_mesh1)
    gmsh.model.mesh.setTransfiniteCurve(line2_spline, n_mesh0)
    gmsh.model.mesh.setTransfiniteCurve(line3_spline, n_mesh1)
    
    # Surface 2
    gmsh.model.mesh.setTransfiniteCurve(line4_spline, n_mesh2)
    gmsh.model.mesh.setTransfiniteCurve(line5_spline, n_mesh1)
    gmsh.model.mesh.setTransfiniteCurve(line6_spline, n_mesh2)
    
    # Surface 3
    gmsh.model.mesh.setTransfiniteCurve(line7_spline, n_mesh1)
    gmsh.model.mesh.setTransfiniteCurve(line8_spline, n_mesh2)
    
    # Create wire and surface
    curves  = [line0_spline,
               line1_spline,
               line2_spline,
               line3_spline]
               
    surface_tag = generateTransfiniteSurfaceMesh(curves,n_mesh0,n_mesh1)

    curves2 = [line4_spline,
               line5_spline,
               line6_spline,
               line2_spline]
               
    surface_tag2 = generateTransfiniteSurfaceMesh(curves2,n_mesh0,n_mesh1)
    
    curves3 = [line7_spline,
               line8_spline,
               line1_spline,
               line6_spline]
               
    surface_tag3 = generateTransfiniteSurfaceMesh(curves3,n_mesh1,n_mesh2)
    
    # Optional: Set mesh algorithm
    gmsh.option.setNumber("Mesh.Algorithm", 8)  # Frontal-Delaunay for quads
    gmsh.option.setNumber("Mesh.RecombineAll", 1)

    all_curve_tags = [line0_spline, line1_spline, line2_spline, line3_spline,
                  line4_spline, line5_spline, line6_spline, line7_spline, line8_spline]

    gmsh.model.addPhysicalGroup(1, all_curve_tags, 2)
    # Create a physical group for export so only this surface is written
    gmsh.model.addPhysicalGroup(2, [surface_tag,surface_tag2,surface_tag3], 1)
    gmsh.model.setPhysicalName(2, 1, "StructuredPatch")

    gmsh.model.mesh.generate(2)

    gmsh.write("projected_square_patch.msh")
    gmsh.write("projected_square_patch.vtk")
    gmsh.finalize()
    
    print("\nMesh generated successfully!")
