import numpy as np
import meshio
from collections import defaultdict

# Extended element type mapping
SU2_TO_MESHIO_TYPE = {
    3: "line",
    5: "triangle",
    9: "quad",
    10: "tetra",
    12: "hexahedron",
    13: "wedge",
    14: "pyramid",
    21: "triangle6",   # Quadratic triangle
    22: "quad8",       # Quadratic quad
    23: "tetra10",     # Quadratic tetra
    24: "hexa20",      # Quadratic hexa
    25: "wedge15",     # Quadratic wedge
    26: "pyramid13",   # Quadratic pyramid
    7: "line3"         # Quadratic line
}

def parse_su2_boundaries(file_path):
    boundaries = {}
    current_boundary = None
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith("NMARK="):
            nmark = int(line.split('=')[1].strip())
        
        elif line.startswith("MARKER_TAG="):
            marker_name = line.split('=')[1].strip()
            current_boundary = {
                'name': marker_name,
                'elements': []
            }
            i += 1
            marker_elems_line = lines[i].strip()
            if not marker_elems_line.startswith("MARKER_ELEMS="):
                raise ValueError("Expected MARKER_ELEMS after MARKER_TAG")
            num_elems = int(marker_elems_line.split('=')[1].strip())
            i += 1
            
            # Read elements
            for _ in range(num_elems):
                element_line = lines[i].strip()
                parts = element_line.split()
                element_type = int(parts[0])
                vertices = list(map(int, parts[1:]))
                
                cell_type = SU2_TO_MESHIO_TYPE.get(element_type)
                if not cell_type:
                    raise ValueError(f"Unsupported element type: {element_type}")
                
                current_boundary['elements'].append((cell_type, vertices))
                i += 1
            
            boundaries[marker_name] = current_boundary
            continue
        
        i += 1
    
    return boundaries

def write_boundary_vtu(boundary, points, output_path):
    # Extract all unique points
    all_vertices = set()
    for cell_type, vertices in boundary['elements']:
        all_vertices.update(vertices)
    
    glob2loc={}
    loc_vid = 0
    local_points = []
    for cell_type, vertices in boundary['elements']:
        # print(vertices)
        vertices_loc = [0] * len(vertices)
        for i in range(0,len(vertices)):
            vid = vertices[i]
    
            if vid not in glob2loc:

                glob2loc[vid] = loc_vid
                local_points.append(points[vid])
                vertices_loc[i] = loc_vid
                loc_vid=loc_vid+1
            else:
                loc_vid_tmp = glob2loc[vid]
                vertices_loc[i] = loc_vid_tmp
        all_vertices.update(vertices_loc)

    # Create global to local mapping
    #  global_to_local = {global_idx: local_idx for local_idx, global_idx in enumerate(sorted(all_vertices))}
    # local_points    = points[list(all_vertices)]
    global_to_local = glob2loc
    # Group cells by type
    cell_groups = defaultdict(list)
    for cell_type, vertices in boundary['elements']:
        local_vertices = [global_to_local[v] for v in vertices]
        
        cell_groups[cell_type].append(local_vertices)
    
    # Create cell blocks with uniform arrays
    cell_blocks = []
    for cell_type, cells_list in cell_groups.items():
        # Create array with consistent shape
        max_vertices = max(len(cell) for cell in cells_list)
        padded_cells = [cell + [-1] * (max_vertices - len(cell)) for cell in cells_list]
        cell_blocks.append((cell_type, np.array(padded_cells, dtype=int)))
    
    # Create mesh
    boundary_mesh = meshio.Mesh(
        points=local_points,
        cells=cell_blocks
    )
    
    # Write to VTU
    boundary_mesh.write(output_path)

# Main processing
input_file = "mimic_hemi.su2"

# 1. Read points from the mesh
mesh = meshio.read(input_file)
points = mesh.points

# 2. Parse boundaries directly from SU2 file
boundaries = parse_su2_boundaries(input_file)

# 3. Write each boundary to separate VTU file
for name, boundary in boundaries.items():
    output_file = f"{name}_boundary.vtu"
    write_boundary_vtu(boundary, points, output_file)
    print(f"Created {output_file} with {len(boundary['elements'])} elements")

# Verification
print("\nBoundary Verification:")
for name, boundary in boundaries.items():
    print(f"{name}: {len(boundary['elements'])} elements")
