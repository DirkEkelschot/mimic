from OCC.Core.STEPControl import STEPControl_Writer
from OCC.Core.StlAPI import StlAPI_Reader
from OCC.Core.TopoDS import TopoDS_Shape
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Sewing
from OCC.Core.ShapeFix import ShapeFix_Shape
from OCC.Core.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.GProp import GProp_GProps
from OCC.Core import BRepGProp
import sys

# Get input/output files from command line or use defaults
if len(sys.argv) >= 2:
    input_stl = sys.argv[1]
else:
    input_stl = "100031.stl"

if len(sys.argv) >= 3:
    output_step = sys.argv[2]
else:
    output_step = "output.step"

print(f"Reading STL file: {input_stl}")

# Read STL
stl_reader = StlAPI_Reader()
shape = TopoDS_Shape()
status = stl_reader.Read(shape, input_stl)

if not status:
    print(f"Error: Failed to read STL file: {input_stl}")
    sys.exit(1)

# Count initial faces
def count_faces(shape):
    count = 0
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        count += 1
        explorer.Next()
    return count

initial_face_count = count_faces(shape)
print(f"Initial number of faces: {initial_face_count}")

print("Sewing faces together with increased tolerance...")

# Sew faces together with larger tolerance to merge more surfaces
# Increased tolerance helps merge surfaces that are close but not exactly coincident
# For complex STL files, use even larger tolerance
sewing_tolerance = 1.0e-2  # Very large tolerance for complex geometries (was 1.0e-3)
sewing = BRepBuilderAPI_Sewing(sewing_tolerance)
sewing.Add(shape)
sewing.Perform()
sewn_shape = sewing.SewedShape()

sewn_face_count = count_faces(sewn_shape)
print(f"Faces after sewing: {sewn_face_count}")

print("Fixing and cleaning shape...")

# Fix the shape (removes small edges, fixes topology)
# Use more aggressive fixing for complex geometries
shape_fix = ShapeFix_Shape(sewn_shape)
# Set tolerance for fixing (if available)
try:
    if hasattr(shape_fix, 'SetPrecision'):
        shape_fix.SetPrecision(1.0e-3)  # Larger precision for more aggressive fixing
    if hasattr(shape_fix, 'SetMaxTolerance'):
        shape_fix.SetMaxTolerance(1.0e-2)  # Allow larger tolerance
except:
    pass
shape_fix.Perform()
fixed_shape = shape_fix.Shape()

fixed_face_count = count_faces(fixed_shape)
print(f"Faces after fixing: {fixed_face_count}")

print("Unifying adjacent surfaces (merging coplanar and connected faces)...")

# Unify same domain (merge adjacent coplanar faces) - this simplifies the geometry
# This merges adjacent faces that are on the same plane, reducing complexity
# Iterate multiple times to merge as much as possible
# Also try with different tolerances and settings
# For complex geometries, do more iterations
current_shape = fixed_shape
max_iterations = 10  # Increased from 5 to 10 for complex geometries
prev_face_count = fixed_face_count

for iteration in range(max_iterations):
    # Use more aggressive unification settings
    # Parameters: (shape, UnifyFaces, UnifyEdges, UnifyVertices)
    unifier = ShapeUpgrade_UnifySameDomain(current_shape, True, True, True)
    
    # Try to set tolerance for merging (if available)
    # Some versions of OpenCASCADE allow setting tolerance
    # Use very large tolerances for complex geometries
    try:
        # Set tolerance for geometric comparison (if method exists)
        if hasattr(unifier, 'SetLinearTolerance'):
            unifier.SetLinearTolerance(1.0e-2)  # Very large tolerance for complex STL files
        if hasattr(unifier, 'SetAngularTolerance'):
            unifier.SetAngularTolerance(5.0e-2)  # Large angular tolerance (about 3 degrees)
        if hasattr(unifier, 'SetTolerance'):
            unifier.SetTolerance(1.0e-2)  # General tolerance
    except:
        pass
    
    unifier.Build()
    unified_shape = unifier.Shape()
    
    current_face_count = count_faces(unified_shape)
    print(f"  Iteration {iteration + 1}: {current_face_count} faces")
    
    # If no more faces were merged, stop iterating
    if current_face_count >= prev_face_count:
        break
    
    prev_face_count = current_face_count
    current_shape = unified_shape

unified_shape = current_shape
final_unified_count = count_faces(unified_shape)
print(f"Final unified face count: {final_unified_count} (reduced from {initial_face_count})")

print("Removing small surfaces and edges...")

# Calculate total surface area to determine threshold
props = GProp_GProps()
total_area = None

try:
    # Try different possible function names for surface area calculation
    if hasattr(BRepGProp, 'SurfaceProperties_s'):
        BRepGProp.SurfaceProperties_s(unified_shape, props)
        total_area = props.Mass()
    elif hasattr(BRepGProp, 'surfaceProperties'):
        BRepGProp.surfaceProperties(unified_shape, props)
        total_area = props.Mass()
    else:
        # Try to get available methods
        methods = [m for m in dir(BRepGProp) if 'Surface' in m and not m.startswith('_')]
        if methods:
            print(f"Available BRepGProp methods: {methods[:5]}")
        # Fallback: use a reasonable default
        total_area = 1000.0
        print("Warning: Could not find surface area calculation method, using default threshold")
except Exception as e:
    # If all else fails, use a fixed threshold
    print(f"Warning: Could not calculate surface area ({e}), using fixed threshold")
    total_area = 1000.0  # Default threshold

if total_area is None or total_area <= 0:
    total_area = 1000.0
    print("Warning: Invalid surface area, using default threshold")

# Set minimum area threshold (adjust this percentage as needed)
# Keep only surfaces that are at least this percentage of total area
# For complex STL files, use a higher threshold to remove more small faces
min_area_percentage = 0.05  # 5% of total area - more aggressive for complex geometries (was 0.02)
min_area_threshold = total_area * min_area_percentage
print(f"Total surface area: {total_area:.6f}")
print(f"Minimum area threshold: {min_area_threshold:.6f} ({min_area_percentage*100}% of total)")
print("Note: Using aggressive filtering for complex geometries")

# Remove small faces by building a new shape with only large faces
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeShell, BRepBuilderAPI_MakeSolid
from OCC.Core.TopTools import TopTools_ListOfShape
from OCC.Core.TopAbs import TopAbs_SHELL, TopAbs_SOLID

# Try to build a shell from large faces
large_faces = TopTools_ListOfShape()
face_count = 0
kept_faces = 0

explorer = TopExp_Explorer(unified_shape, TopAbs_FACE)
while explorer.More():
    face = explorer.Current()
    face_count += 1
    
    face_props = GProp_GProps()
    try:
        if hasattr(BRepGProp, 'SurfaceProperties_s'):
            BRepGProp.SurfaceProperties_s(face, face_props)
        elif hasattr(BRepGProp, 'surfaceProperties'):
            BRepGProp.surfaceProperties(face, face_props)
        else:
            # Skip this face if we can't calculate area
            print(f"  Skipping face {face_count}: cannot calculate area")
            explorer.Next()
            continue
        face_area = face_props.Mass()
    except Exception as e:
        # Skip faces where area calculation fails
        print(f"  Skipping face {face_count}: error calculating area - {e}")
        explorer.Next()
        continue
    
    if face_area >= min_area_threshold:
        large_faces.Append(face)
        kept_faces += 1
        print(f"  Keeping face {face_count}: area = {face_area:.6f}")
    else:
        print(f"  Removing small face {face_count}: area = {face_area:.6f} (threshold: {min_area_threshold:.6f})")
    
    explorer.Next()

print(f"Kept {kept_faces} out of {face_count} faces")

# After filtering, try to unify again to merge any newly adjacent faces
# Do multiple passes for complex geometries
if large_faces.Size() > 0:
    print("Unifying remaining faces after filtering (multiple passes)...")
    # Build a temporary shell to unify
    temp_shell_builder = BRepBuilderAPI_MakeShell()
    for i in range(1, large_faces.Size() + 1):
        temp_shell_builder.Add(large_faces.Value(i))
    
    if temp_shell_builder.IsDone():
        temp_shell = temp_shell_builder.Shell()
        current_shell = temp_shell
        prev_face_count = large_faces.Size()
        
        # Multiple unification passes - more for complex geometries
        for unify_pass in range(5):  # Increased from 3 to 5
            # Unify again with aggressive settings
            unifier2 = ShapeUpgrade_UnifySameDomain(current_shell, True, True, True)
            try:
                if hasattr(unifier2, 'SetLinearTolerance'):
                    unifier2.SetLinearTolerance(1.0e-2)
                if hasattr(unifier2, 'SetAngularTolerance'):
                    unifier2.SetAngularTolerance(5.0e-2)
            except:
                pass
            unifier2.Build()
            unified_after_filter = unifier2.Shape()
            
            # Re-extract faces from unified shape
            large_faces = TopTools_ListOfShape()
            explorer2 = TopExp_Explorer(unified_after_filter, TopAbs_FACE)
            while explorer2.More():
                large_faces.Append(explorer2.Current())
                explorer2.Next()
            
            current_count = large_faces.Size()
            print(f"  Pass {unify_pass + 1}: {current_count} faces")
            
            if current_count >= prev_face_count:
                break
            
            prev_face_count = current_count
            
            # Rebuild shell for next iteration
            if current_count > 0:
                temp_shell_builder2 = BRepBuilderAPI_MakeShell()
                for i in range(1, large_faces.Size() + 1):
                    temp_shell_builder2.Add(large_faces.Value(i))
                if temp_shell_builder2.IsDone():
                    current_shell = temp_shell_builder2.Shell()
                else:
                    break
            else:
                break
        
        print(f"After post-filter unification: {large_faces.Size()} faces")

# Build shell from large faces
if large_faces.Size() > 0:
    shell_builder = BRepBuilderAPI_MakeShell()
    for i in range(1, large_faces.Size() + 1):
        shell_builder.Add(large_faces.Value(i))
    
    if shell_builder.IsDone():
        shell_shape = shell_builder.Shell()
        
        # Try to build a solid if the shell is closed - this might merge more surfaces
        from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeSolid
        solid_builder = BRepBuilderAPI_MakeSolid()
        solid_builder.Add(shell_shape)
        
        if solid_builder.IsDone():
            # Check if the solid is valid
            from OCC.Core.BRepCheck import BRepCheck_Analyzer
            analyzer = BRepCheck_Analyzer(solid_builder.Solid())
            if analyzer.IsValid():
                final_shape = solid_builder.Solid()
                print(f"Built solid from large faces: {count_faces(final_shape)} faces")
            else:
                final_shape = shell_shape
                print(f"Built shell from large faces: {count_faces(final_shape)} faces (solid was invalid)")
        else:
            final_shape = shell_shape
            print(f"Built shell from large faces: {count_faces(final_shape)} faces")
    else:
        # If shell building fails, use the unified shape
        print("Warning: Could not build shell, using unified shape")
        final_shape = unified_shape
else:
    print("Warning: No large faces found, using unified shape")
    final_shape = unified_shape

# Final unification pass to merge any remaining mergeable surfaces
# This is critical - do multiple passes with different approaches
# For complex geometries, do even more iterations
print("\nPerforming final unification pass...")
final_unified = final_shape
for final_iteration in range(10):  # Increased from 5 to 10 for complex geometries
    unifier_final = ShapeUpgrade_UnifySameDomain(final_unified, True, True, True)
    
    # Try to set more aggressive tolerances for complex geometries
    try:
        if hasattr(unifier_final, 'SetLinearTolerance'):
            unifier_final.SetLinearTolerance(1.0e-2)  # Very large tolerance
        if hasattr(unifier_final, 'SetAngularTolerance'):
            unifier_final.SetAngularTolerance(5.0e-2)  # Large angular tolerance
        if hasattr(unifier_final, 'SetTolerance'):
            unifier_final.SetTolerance(1.0e-2)
    except:
        pass
    
    unifier_final.Build()
    temp_shape = unifier_final.Shape()
    temp_count = count_faces(temp_shape)
    current_count = count_faces(final_unified)
    
    if temp_count < current_count:
        print(f"  Final iteration {final_iteration + 1}: {temp_count} faces (reduced from {current_count})")
        final_unified = temp_shape
    else:
        # Try one more time with shape fixing before unification
        if final_iteration < 9:  # Increased from 4 to 9
            print(f"  No reduction, trying shape fix before unification...")
            shape_fix_final = ShapeFix_Shape(final_unified)
            shape_fix_final.Perform()
            fixed_final = shape_fix_final.Shape()
            
            unifier_final2 = ShapeUpgrade_UnifySameDomain(fixed_final, True, True, True)
            try:
                if hasattr(unifier_final2, 'SetLinearTolerance'):
                    unifier_final2.SetLinearTolerance(1.0e-2)  # Very large tolerance
                if hasattr(unifier_final2, 'SetAngularTolerance'):
                    unifier_final2.SetAngularTolerance(5.0e-2)  # Large angular tolerance
                if hasattr(unifier_final2, 'SetTolerance'):
                    unifier_final2.SetTolerance(1.0e-2)
            except:
                pass
            unifier_final2.Build()
            temp_shape2 = unifier_final2.Shape()
            temp_count2 = count_faces(temp_shape2)
            
            if temp_count2 < current_count:
                print(f"  After shape fix: {temp_count2} faces (reduced from {current_count})")
                final_unified = temp_shape2
                continue
        break

final_shape = final_unified

# Final face count
final_face_count = count_faces(final_shape)
print(f"\nFinal number of faces in STEP file: {final_face_count}")
print(f"Reduction: {initial_face_count} -> {final_face_count} ({100.0 * (1.0 - final_face_count / initial_face_count):.1f}% reduction)")

# Write STEP
print(f"Writing STEP file: {output_step}")
step_writer = STEPControl_Writer()
status = step_writer.Transfer(final_shape, 1)  # 1 = STEPControl_AsIs

if status != 1:
    print(f"Warning: Transfer status: {status}")

status = step_writer.Write(output_step)

if status != 1:
    print(f"Error: Failed to write STEP file: {output_step}")
    sys.exit(1)

print("Conversion successful!")