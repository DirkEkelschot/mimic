#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <ShapeFix_Shape.hxx>
#include <gp_Ax2.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <TopoDS_Shape.hxx>
#include <STEPControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include <Interface_Static.hxx>

#include <iostream>
#include <string>

void createSimpleAirplane(const std::string& filename) {
    std::cout << "Creating simple airplane geometry (unified solid)..." << std::endl;
    
    // Strategy: Build everything as one compound, then use BRepBuilderAPI_Sewing
    // to merge into a single shell, then make it a solid
    
    // Main fuselage body (cylinder)
    gp_Ax2 axis(gp_Pnt(-3, 0, 0), gp_Dir(1, 0, 0));
    TopoDS_Shape fuselageBody = BRepPrimAPI_MakeCylinder(axis, 0.5, 6.0).Shape();
    
    // Nose cone
    gp_Ax2 noseAxis(gp_Pnt(3, 0, 0), gp_Dir(1, 0, 0));
    TopoDS_Shape nose = BRepPrimAPI_MakeCone(noseAxis, 0.5, 0.1, 1.5).Shape();
    
    // Tail cone
    gp_Ax2 tailAxis(gp_Pnt(-3, 0, 0), gp_Dir(-1, 0, 0));
    TopoDS_Shape tail = BRepPrimAPI_MakeCone(tailAxis, 0.5, 0.2, 1.0).Shape();
    
    // Combine fuselage parts first
    TopoDS_Shape fuselageComplete = BRepAlgoAPI_Fuse(fuselageBody, nose).Shape();
    fuselageComplete = BRepAlgoAPI_Fuse(fuselageComplete, tail).Shape();
    
    // Make wings thicker so they properly intersect with fuselage
    // This ensures clean Boolean operations
    TopoDS_Shape leftWing = BRepPrimAPI_MakeBox(gp_Pnt(-1, 0.3, -0.15),
                                                  gp_Pnt(1, 4.5, 0.15)).Shape();
    TopoDS_Shape rightWing = BRepPrimAPI_MakeBox(gp_Pnt(-1, -4.5, -0.15),
                                                   gp_Pnt(1, -0.3, 0.15)).Shape();
    
    // Combine wings into fuselage - this should create a unified solid
    TopoDS_Shape airplane = BRepAlgoAPI_Fuse(fuselageComplete, leftWing).Shape();
    airplane = BRepAlgoAPI_Fuse(airplane, rightWing).Shape();
    
    // Tail wings (horizontal stabilizer) - also made thicker
    TopoDS_Shape leftTailWing = BRepPrimAPI_MakeBox(gp_Pnt(-3.5, 0.2, -0.08),
                                                      gp_Pnt(-2.5, 1.5, 0.08)).Shape();
    TopoDS_Shape rightTailWing = BRepPrimAPI_MakeBox(gp_Pnt(-3.5, -1.5, -0.08),
                                                       gp_Pnt(-2.5, -0.2, 0.08)).Shape();
    
    airplane = BRepAlgoAPI_Fuse(airplane, leftTailWing).Shape();
    airplane = BRepAlgoAPI_Fuse(airplane, rightTailWing).Shape();
    
    // Vertical stabilizer
    TopoDS_Shape vertStabilizer = BRepPrimAPI_MakeBox(gp_Pnt(-3.5, -0.08, 0),
                                                        gp_Pnt(-2.5, 0.08, 1.2)).Shape();
    
    airplane = BRepAlgoAPI_Fuse(airplane, vertStabilizer).Shape();
    
    // CRITICAL: Unify same-domain faces to remove internal boundaries
    // This merges coplanar/coincident faces and removes internal edges
    std::cout << "  Unifying solid (removing internal faces)..." << std::endl;
    ShapeUpgrade_UnifySameDomain unifier(airplane);
    unifier.AllowInternalEdges(Standard_False);  // Don't allow internal edges
    unifier.Build();
    airplane = unifier.Shape();
    
    // Fix the shape if needed
    Handle(ShapeFix_Shape) shapeFixer = new ShapeFix_Shape(airplane);
    shapeFixer->Perform();
    airplane = shapeFixer->Shape();
    
    std::cout << "  Airplane built as unified solid (internal faces removed)" << std::endl;
    
    // Export to STEP
    STEPControl_Writer writer;
    Interface_Static::SetCVal("write.step.schema", "AP214");
    writer.Transfer(airplane, STEPControl_AsIs);
    IFSelect_ReturnStatus status = writer.Write(filename.c_str());
    
    if (status == IFSelect_RetDone) {
        std::cout << "✓ Airplane saved to: " << filename << std::endl;
    } else {
        std::cout << "✗ Error writing file" << std::endl;
    }
}
    if (status == IFSelect_RetDone) {
        std::cout << "✓ Airplane saved to: " << filename << std::endl;
    } else {
        std::cout << "✗ Error writing file" << std::endl;
    }
}

void createSimpleCar(const std::string& filename) {
    std::cout << "Creating simple car geometry..." << std::endl;
    
    // Car body (main box)
    TopoDS_Shape body = BRepPrimAPI_MakeBox(gp_Pnt(-2, -1, 0),
                                             gp_Pnt(2, 1, 0.8)).Shape();
    
    // Cabin (smaller box on top)
    TopoDS_Shape cabin = BRepPrimAPI_MakeBox(gp_Pnt(-0.5, -0.8, 0.8),
                                              gp_Pnt(1.0, 0.8, 1.5)).Shape();
    
    // Combine body and cabin
    TopoDS_Shape carBody = BRepAlgoAPI_Fuse(body, cabin).Shape();
    
    // Wheels (cylinders)
    gp_Ax2 wheelAxis(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0));
    
    TopoDS_Shape wheel1 = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(-1.2, -1.2, 0.3), gp_Dir(0, 1, 0)), 0.4, 0.2).Shape();
    TopoDS_Shape wheel2 = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(-1.2, 1.0, 0.3), gp_Dir(0, 1, 0)), 0.4, 0.2).Shape();
    TopoDS_Shape wheel3 = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(1.2, -1.2, 0.3), gp_Dir(0, 1, 0)), 0.4, 0.2).Shape();
    TopoDS_Shape wheel4 = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(1.2, 1.0, 0.3), gp_Dir(0, 1, 0)), 0.4, 0.2).Shape();
    
    // Combine all
    TopoDS_Shape car = BRepAlgoAPI_Fuse(carBody, wheel1).Shape();
    car = BRepAlgoAPI_Fuse(car, wheel2).Shape();
    car = BRepAlgoAPI_Fuse(car, wheel3).Shape();
    car = BRepAlgoAPI_Fuse(car, wheel4).Shape();
    
    // Export to STEP
    STEPControl_Writer writer;
    Interface_Static::SetCVal("write.step.schema", "AP214");
    writer.Transfer(car, STEPControl_AsIs);
    IFSelect_ReturnStatus status = writer.Write(filename.c_str());
    
    if (status == IFSelect_RetDone) {
        std::cout << "✓ Car saved to: " << filename << std::endl;
    } else {
        std::cout << "✗ Error writing file" << std::endl;
    }
}

void createSimpleSubmarine(const std::string& filename) {
    std::cout << "Creating simple submarine geometry..." << std::endl;
    
    // Main hull (cylinder with rounded ends)
    gp_Ax2 axis(gp_Pnt(-3, 0, 0), gp_Dir(1, 0, 0));
    TopoDS_Shape hull = BRepPrimAPI_MakeCylinder(axis, 0.8, 6.0).Shape();
    
    // Nose hemisphere
    TopoDS_Shape nose = BRepPrimAPI_MakeSphere(gp_Pnt(3, 0, 0), 0.8).Shape();
    TopoDS_Shape noseKeep = BRepPrimAPI_MakeBox(gp_Pnt(3, -1, -1),
                                                  gp_Pnt(4, 1, 1)).Shape();
    nose = BRepAlgoAPI_Cut(nose, noseKeep).Shape();
    
    // Tail hemisphere
    TopoDS_Shape tailSphere = BRepPrimAPI_MakeSphere(gp_Pnt(-3, 0, 0), 0.8).Shape();
    TopoDS_Shape tailCut = BRepPrimAPI_MakeBox(gp_Pnt(-4, -1, -1),
                                                 gp_Pnt(-3, 1, 1)).Shape();
    TopoDS_Shape tailEnd = BRepAlgoAPI_Cut(tailSphere, tailCut).Shape();
    
    // Conning tower (vertical cylinder)
    gp_Ax2 towerAxis(gp_Pnt(0, 0, 0.8), gp_Dir(0, 0, 1));
    TopoDS_Shape tower = BRepPrimAPI_MakeCylinder(towerAxis, 0.3, 0.6).Shape();
    
    // Combine all parts
    TopoDS_Shape submarine = BRepAlgoAPI_Fuse(hull, nose).Shape();
    submarine = BRepAlgoAPI_Fuse(submarine, tailEnd).Shape();
    submarine = BRepAlgoAPI_Fuse(submarine, tower).Shape();
    
    // Export to STEP
    STEPControl_Writer writer;
    Interface_Static::SetCVal("write.step.schema", "AP214");
    writer.Transfer(submarine, STEPControl_AsIs);
    IFSelect_ReturnStatus status = writer.Write(filename.c_str());
    
    if (status == IFSelect_RetDone) {
        std::cout << "✓ Submarine saved to: " << filename << std::endl;
    } else {
        std::cout << "✗ Error writing file" << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::cout << "=== Test STEP File Generator ===" << std::endl;
    std::cout << "Creates simple geometries for testing mesh generation" << std::endl << std::endl;
    
    std::string choice = "all";
    if (argc >= 2) {
        choice = argv[1];
    }
    
    if (choice == "airplane" || choice == "all") {
        createSimpleAirplane("test_airplane.step");
    }
    
    if (choice == "car" || choice == "all") {
        createSimpleCar("test_car.step");
    }
    
    if (choice == "submarine" || choice == "all") {
        createSimpleSubmarine("test_submarine.step");
    }
    
    std::cout << "\n=== Done ===" << std::endl;
    std::cout << "You can now test with:" << std::endl;
    std::cout << "  ./mesh_around_cad test_airplane.step" << std::endl;
    std::cout << "  ./mesh_around_cad test_car.step" << std::endl;
    std::cout << "  ./mesh_around_cad test_submarine.step" << std::endl;
    
    return 0;
}
