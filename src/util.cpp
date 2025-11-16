#include "util.h"
#include "types.h"

#include <STEPCAFControl_Reader.hxx>
#include <TDocStd_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <XCAFApp_Application.hxx>

#include <V3d_View.hxx>
#include <V3d_Viewer.hxx>
#include <AIS_InteractiveContext.hxx>
#include <AIS_Shape.hxx>
#include <STEPControl_Reader.hxx>
#include <Aspect_DisplayConnection.hxx>
#include <OpenGl_GraphicDriver.hxx>

#include <Aspect_NeutralWindow.hxx>
#include <Aspect_Handle.hxx>
#include <Xw_Window.hxx>

#include <X11/Xlib.h>

#include <TopoDS_Shape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Pnt.hxx>
#include <TopoDS_Face.hxx>
#include <Geom_Line.hxx>
#include <AIS_Line.hxx>
#include <Geom_CartesianPoint.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <Geom_Axis2Placement.hxx>
#include <AIS_Axis.hxx>
#include <AIS_Trihedron.hxx>
#include <BRepBndLib.hxx>

using namespace types;

TopoDS_Shape util::loadStepShape(const std::string& path)
{
    Handle(XCAFApp_Application) app = XCAFApp_Application::GetApplication();
    Handle(TDocStd_Document) doc;
    app->NewDocument("MDTV-XCAF", doc);

    STEPCAFControl_Reader reader;
    IFSelect_ReturnStatus status = reader.ReadFile(path.c_str());

    if (status != IFSelect_RetDone)
        throw std::runtime_error("Failed STEP read");

    reader.Transfer(doc);

    Handle(XCAFDoc_ShapeTool) shapeTool =
        XCAFDoc_DocumentTool::ShapeTool(doc->Main());

    TDF_LabelSequence labels;
    shapeTool->GetFreeShapes(labels);

    if (labels.Length() < 1)
        throw std::runtime_error("No shapes in STEP");

    TDF_Label L = labels.Value(1);

    // Extract geometry + its full location
    TopoDS_Shape shape = shapeTool->GetShape(L);
    TopLoc_Location loc = shapeTool->GetLocation(L);

    // Apply instance transform
    shape.Move(loc);

    return shape;
}


int util::displayShape(const TopoDS_Shape shape, const PointSet samples)
{
    // ----------------------------
    // Create X11 window
    // ----------------------------
    Display* xDisplay = XOpenDisplay(nullptr);
    if (!xDisplay) {
        std::cerr << "Cannot open X11 display." << std::endl;
        return 1;
    }

    int screen = DefaultScreen(xDisplay);
    Window root = RootWindow(xDisplay, screen);

    Window win = XCreateSimpleWindow(
        xDisplay, root, 0, 0, 800, 600, 0,
        BlackPixel(xDisplay, screen),
        WhitePixel(xDisplay, screen)
    );
    XStoreName(xDisplay, win, "OpenCascade Viewer");
    XMapWindow(xDisplay, win);
    XFlush(xDisplay);

    // ----------------------------
    // Wrap X11 window into OCCT
    // ----------------------------
    Handle(Aspect_DisplayConnection) displayConnection = new Aspect_DisplayConnection();
    Handle(Xw_Window) window = new Xw_Window(displayConnection, (Aspect_Drawable)win);

    if (!window->IsMapped())
        window->Map();

    Handle(OpenGl_GraphicDriver) driver = new OpenGl_GraphicDriver(displayConnection);
    Handle(V3d_Viewer) viewer = new V3d_Viewer(driver);
    Handle(V3d_View) view = viewer->CreateView();
    view->SetWindow(window);

    // --------------------------
    // AIS context + shape
    // --------------------------
    Handle(AIS_InteractiveContext) context = new AIS_InteractiveContext(viewer);

    Handle(AIS_Shape) aisShape = new AIS_Shape(shape);
    context->Display(aisShape, Standard_True);
    
    // ----------------------------
    // ADD AXES
    // ----------------------------

    // Display coordinate trihedron (X/Y/Z axes with arrows)
    Handle(Geom_Axis2Placement) originA2P = 
        new Geom_Axis2Placement(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));

    Handle(AIS_Trihedron) trihedron = new AIS_Trihedron(originA2P);
    trihedron->SetDatumDisplayMode(Prs3d_DM_Shaded);
    trihedron->Attributes()->SetDatumAspect(new Prs3d_DatumAspect);

    // Colors
    trihedron->SetXAxisColor(Quantity_NOC_GREEN);
    trihedron->SetYAxisColor(Quantity_NOC_RED);
    trihedron->SetAxisColor(Quantity_NOC_YELLOW);

    trihedron->SetSize(10);

    context->Display(trihedron, Standard_False);

    // ----------------------------
    // Draw sample points only
    // ----------------------------
    for (const auto& s : samples.getSamples())
    {
        arraydouble3 p = s.getPosition();

        gp_Pnt gpP(p.x, p.y, p.z);

        // Draw small sphere at the point
        TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(gpP, 0.5).Shape();
        Handle(AIS_Shape) aisPoint = new AIS_Shape(sphere);
        aisPoint->SetColor(Quantity_NOC_RED);

        context->Display(aisPoint, Standard_False);
    }

    view->FitAll();
    view->Redraw();

    Bnd_Box bb;
    BRepBndLib::Add(shape, bb);
    Standard_Real xmin, ymin, zmin, xmax, ymax, zmax;
    bb.Get(xmin, ymin, zmin, xmax, ymax, zmax);

    gp_Pnt center(
        0.5 * (xmin + xmax),
        0.5 * (ymin + ymax),
        0.5 * (zmin + zmax)
    );

    double angle = 0.0;
    double radius = 300;

    while (true)
    {
        angle += 0.01;
        if (angle > 2*M_PI) angle = 0.0;

        // Update camera
        view->SetEye(center.X() + radius * cos(angle), center.Y() + radius * sin(angle), center.Z() + radius * 0.15);
        view->SetAt(center.X(), center.Y(), center.Z());
        view->SetUp(0,0,1);

        view->Redraw();
    }

    return 0;
}



void util::extractMeshData(const TopoDS_Shape& shape, std::vector<arraydouble3>& outPoints, std::vector<std::array<int, 3>>& outTriangles)
{
    outPoints.clear();
    outTriangles.clear();

    int vertexOffset = 0;

    for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next())
    {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        
        TopLoc_Location loc;
        Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, loc);

        if (triangulation.IsNull())
            continue;

        gp_Trsf triTrsf  = loc.Transformation();
        gp_Trsf faceTrsf = face.Location().Transformation();

        for (int i = 1; i <= triangulation->NbNodes(); i++)
        {
            gp_Pnt p = triangulation->Node(i);

            p.Transform(triTrsf);

            p.Transform(faceTrsf);

            outPoints.push_back({ 
                p.X(), 
                p.Y(), 
                p.Z() 
            });
        }

        for (int i = 1; i <= triangulation->NbTriangles(); i++)
        {
            int n1, n2, n3;
            triangulation->Triangle(i).Get(n1, n2, n3);

            outTriangles.push_back({
                vertexOffset + n1 - 1,
                vertexOffset + n2 - 1,
                vertexOffset + n3 - 1
            });
        }

        vertexOffset += triangulation->NbNodes();
    }
}

