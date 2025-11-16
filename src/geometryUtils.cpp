#include "geometryUtils.h"
#include "types.h"
#include "debugVars.h"

#include <BRep_Tool.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepLProp_SLProps.hxx>
#include <ShapeAnalysis_Surface.hxx>

#include <algorithm>
#include <vector>
#include <array>
#include <numbers>

#include <omp.h>

using namespace types;

namespace geometryUtils
{
    std::array<arraydouble3, 3> dot3x3(const std::array<arraydouble3, 3>& m1, const std::array<arraydouble3, 3>& m2)
    {
        std::array<arraydouble3, 3> m3{};

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    m3[i][j] += m1[i][k] * m2[k][j];
        
        return m3;
    }

    arraydouble3 dot3x1(const std::array<arraydouble3, 3>& m1, const arraydouble3& m2)
    {
        arraydouble3 m3{};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m3[i] += m1[i][j] * m2[j];
        return m3;
    }

    std::array<arraydouble3, 3> rotationMatrixToZ(arraydouble3 normal)
    {   
        normal = normal / sqrtf(normal.dot(normal));

        arraydouble3 zAxis {0.0f, 0.0f, 1.0f};
        double epsilon {0.001f};
        double theta {0.0f};
        arraydouble3 k {};
        arraydouble3 arb {};

        double normalDot = normal.dot(zAxis);

        if ( normalDot >= 1 - epsilon )
            return {{ {1.0f,0.0f,0.0f},
                    {0.0f,1.0f,0.0f},
                    {0.0f,0.0f,1.0f} }};
        if ( normalDot <= -1 + epsilon )
        {
            arb = {1.0f,0.0f,0.0f};
            if (fabs(normal.dot(arb)) > 1.0f - epsilon)
                arb = {0.0f, 1.0f, 0.0f};
            k = normal.cross(arb);
            k = k / sqrtf(k.dot(k));
            theta = std::numbers::pi_v<double>;
        }
        else
        {
            k = normal.cross(zAxis);
            double kLen = sqrtf(k.dot(k));
            k = {k[0] / kLen, k[1] / kLen, k[2] / kLen};
            double abDot = normal.dot(zAxis);
            if (abDot < -1)
                abDot = -1;
            else if (abDot > 1)
                abDot = 1;
            theta = acos(abDot);
        }

        std::array<arraydouble3, 3> skewK {{ 
            { 0.0f, -k[2], k[1] },
            { k[2], 0.0f, -k[0] },
            { -k[1], k[0], 0.0f } 
        }};

        std::array<arraydouble3, 3> R {{ 
            { 1.0f, 0.0f, 0.0f },
            { 0.0f, 1.0f, 0.0f },
            { 0.0f, 0.0f, 1.0f } 
        }};

        std::array<arraydouble3, 3> skewK2{dot3x3(skewK, skewK)};

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                R[i][j] += sinf(theta) * skewK[i][j] + (1 - cosf(theta)) * skewK2[i][j];

        return R;
    }

    std::array<std::array<double, 4>, 4> transformMatrixToZAxis(arraydouble3 normal, arraydouble3 point)
    {
        std::array<arraydouble3, 3> R = rotationMatrixToZ(normal);

        arraydouble3 t{};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                t[i] -= R[i][j] * point[j];

        std::array<std::array<double, 4>, 4> T {{
            {{ R[0][0], R[0][1], R[0][2], t[0] }},
            {{ R[1][0], R[1][1], R[1][2], t[1] }},
            {{ R[2][0], R[2][1], R[2][2], t[2] }},
            {{ 0.0f,    0.0f,    0.0f,    1.0f }}
        }};

        return T;
    }

    bool isPointInXYTriangle(std::array<double, 2> point, std::array<arraydouble3, 3> triangle, double epsilon)
    {
        double px = point[0];
        double py = point[1];
        arraydouble3 v1 = triangle[0];
        arraydouble3 v2 = triangle[1];
        arraydouble3 v3 = triangle[2];

        double areaV1V2V3 = 0.5 * fabs((v2[0] - v1[0]) * (v3[1] - v1[1]) - (v3[0] - v1[0]) * (v2[1] - v1[1]));

        double pv1[2] { px - v1[0], py - v1[1] };
        double pv2[2] { px - v2[0], py - v2[1] };
        double pv3[2] { px - v3[0], py - v3[1] };

        double areaPV2V3 = 0.5 * fabs(pv2[0] * pv3[1] - pv3[0] * pv2[1]);
        double areaPV1V3 = 0.5 * fabs(pv1[0] * pv3[1] - pv3[0] * pv1[1]);
        double areaPV1V2 = 0.5 * fabs(pv1[0] * pv2[1] - pv2[0] * pv1[1]);

        if ((fabs(areaPV2V3 + areaPV1V3 + areaPV1V2) - areaV1V2V3) < epsilon)
            return true;
        return false;
    }

    arraydouble3 triangleNormal(std::array<arraydouble3, 3> triangle)
    {
        arraydouble3 v1 = triangle[0];
        arraydouble3 v2 = triangle[1];
        arraydouble3 v3 = triangle[2];
        arraydouble3 v1v2 = v2-v1;
        arraydouble3 v1v3 = v3-v1;
        return v1v2.cross(v1v3);
    }

    double triangleD(std::array<arraydouble3, 3> triangle, arraydouble3 normal)
    {
        return -normal.dot(triangle[0]);
    }

    double triangleD(std::array<arraydouble3, 3> triangle)
    {
        arraydouble3 normal = triangleNormal(triangle);
        return normal.dot(triangle[0]);
    }

    double zFromXY(std::array<double, 2> P, arraydouble3 normal, double D)
    {
        return -(normal.x * P[0] + normal.y * P[1] + D) / normal.z;
    }

    std::vector<std::array<double, 2>> generateUniformPointsOnFace(TopoDS_Face face, double resolution)
    {
        BRepAdaptor_Surface adaptorSurface = BRepAdaptor_Surface(face);

        double uMin = adaptorSurface.FirstUParameter();
        double uMax = adaptorSurface.LastUParameter();
        double vMin = adaptorSurface.FirstVParameter();
        double vMax = adaptorSurface.LastVParameter();
        double uRange = uMax - uMin;
        double vRange = vMax - vMin;

        GeomAbs_SurfaceType surfaceType = adaptorSurface.GetType();
        double radius{};
        if (surfaceType == GeomAbs_Cylinder)
        {
            radius = adaptorSurface.Cylinder().Radius();
            uRange = (radius * uRange);
        }
        else if (surfaceType == GeomAbs_Sphere)
        {
            radius = adaptorSurface.Sphere().Radius();
            uRange = radius * uRange;
            vRange = radius * vRange;
        }
        else if (surfaceType == GeomAbs_Torus)
        {
            gp_Torus torus = adaptorSurface.Torus();
            uRange = torus.MajorRadius() * uRange;
            vRange = torus.MinorRadius() * vRange;
        }

        if (surfaceType == GeomAbs_Cone || surfaceType == GeomAbs_SurfaceOfRevolution)
        {
            return generateNonUniformRevolutionUV(adaptorSurface, uMin, vMin, vMax, uRange, vRange, resolution); 
        }

        return generateUniformUV(uMin, uMax, vMin, vMax, uRange, vRange, resolution); 
    }   

    std::vector<std::array<double, 2>> generateNonUniformRevolutionUV(
        BRepAdaptor_Surface adaptorSurface, 
        double uMin, double vMin, double vMax, double uRange, 
        double vRange, double resolution, double epsilon)
    {
        const int minNumberOfPoints = 3;
        const int maxNumberOfPoints = 200;

        int nV = std::clamp(
            static_cast<int>(std::round(vRange / resolution)),
            minNumberOfPoints,
            maxNumberOfPoints);
        nV += (nV % 2 == 0);
        
        double vStepSize = (vMax - vMin - 2 * epsilon) / (nV - 1);
        gp_Ax1 axis{};
        
        if (adaptorSurface.GetType() == GeomAbs_Cone)
            axis = adaptorSurface.Cone().Axis();
        else
            axis = adaptorSurface.AxeOfRevolution();

        gp_Pnt origin = axis.Location();
        gp_Vec direction = static_cast<gp_Vec>(axis.Direction());

        std::vector<int> nUs;
        nUs.reserve(nV);
        
        gp_Pnt p{};
        gp_Vec vec{};
        double radius{};
        double circumference{};
        int nU{};
        for (int i = 0; i < nV; i++)
        {
            p = adaptorSurface.Value(0, vMin + epsilon + i * vStepSize);
            vec = {origin, p};
            radius = vec.Crossed(direction).Magnitude() / direction.Magnitude();
            circumference = 2 * std::numbers::pi_v<double> * radius;
            nU = std::clamp(static_cast<int>(round(circumference / resolution)), minNumberOfPoints, maxNumberOfPoints);
            if (nU % 2 == 0) nU++;
            nUs.push_back(nU);
        }

        std::vector<std::array<double, 2>> uvPoints;
        double uStepSize{};
        for (int i = 0; i < nV; i++)
        {    
            uStepSize = (uRange - 2 * epsilon) / (nUs[i] - 1);
            // This loop can be moved into the previous loop
            for (int j = 0; j < nUs[i]; j++)
                uvPoints.push_back({uMin + epsilon + j * uStepSize, vMin + epsilon + i * vStepSize});
        }

        return uvPoints;
    }

    std::vector<std::array<double, 2>> filterPointsOutsideFace(std::vector<std::array<double, 2>> uvPoints, TopoDS_Face face)
    {
        Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
        std::vector<std::array<double, 2>> insidePoints;

        BRepClass_FaceClassifier classifier;

        for (auto& uv : uvPoints)
        {
            gp_Pnt P = surface->Value(uv[0], uv[1]);

            classifier.Perform(face, P, 1e-9);

            if (classifier.State() == TopAbs_IN)
                insidePoints.push_back(uv);
        }
        
        return insidePoints;
    }

    void setPositionAndNormals(
        const std::vector<std::array<double, 2>>& uvSamples,
        const TopoDS_Face& face,
        std::vector<arraydouble3>& outPoints,
        std::vector<arraydouble3>& outNormals)
    {
        outPoints.clear();
        outNormals.clear();

        for (auto& uv : uvSamples)
        {
            outPoints.push_back(uvPositionToGlobal(face, uv[0], uv[1]));
            outNormals.push_back(uvNormal(face, uv[0], uv[1]));
        }
    }

    arraydouble3 uvPositionToGlobal(TopoDS_Face face, double u, double v)
    {
        Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
        ShapeAnalysis_Surface sas = ShapeAnalysis_Surface(surface);

        gp_Pnt p = sas.Value(u, v);

        return {p.X(), p.Y(), p.Z()};
    }

    arraydouble3 uvNormal(TopoDS_Face face, double u, double v)
    {
        BRepAdaptor_Surface adaptorSurface = BRepAdaptor_Surface(face);
        BRepLProp_SLProps props = BRepLProp_SLProps(adaptorSurface, u, v, 1, 1e-6);
        if (!props.IsNormalDefined())
            throw std::invalid_argument("Normal is not defined.");

        gp_Dir normal = props.Normal();

        if (face.Orientation() == TopAbs_REVERSED)
            normal.Reverse();

        return {normal.X(), normal.Y(), normal.Z()};
    }


    std::vector<std::array<double, 2>> generateUniformUV(
    double uMin, double uMax, double vMin, 
    double vMax, double uRange, double vRange, double resolution, double epsilon)
    {
        std::vector<std::array<double, 2>> points;

        const int minNumberOfPoints = 3;
        const int maxNumberOfPoints = 200;

        int uCount = std::clamp(static_cast<int>(std::round(uRange / resolution)), minNumberOfPoints, maxNumberOfPoints);
        int vCount = std::clamp(static_cast<int>(std::round(vRange / resolution)), minNumberOfPoints, maxNumberOfPoints);
        
        uCount += (uCount % 2 == 0);
        vCount += (vCount % 2 == 0);

        for (int i = 0; i < uCount; ++i)
        {
            for (int j = 0; j < vCount; ++j)
            {
                double u = uMin + epsilon + (static_cast<double>(i) / (uCount - 1)) * ((uMax - epsilon) - (uMin + epsilon));
                double v = vMin + epsilon + (static_cast<double>(j) / (vCount - 1)) * ((vMax - epsilon) - (vMin + epsilon));
                points.push_back({u, v});
            }
        }

        return points;
    }

}

types::ZMap geometryUtils::meshZMap(arraydouble3 point, arraydouble3 normal, 
                                    double size, Mesh mesh, int samplingRate, 
                                    double zMaxThreshold)
{
    std::array<std::array<double, 4>, 4> tMatrix = transformMatrixToZAxis(normal, point);

    std::vector<double> samples(samplingRate);
    double stepSize{size / samplingRate};

    for (int i = 0; i < samplingRate; i++)
        samples[i] = (i - samplingRate/2) * stepSize;

    mesh.removeTinyTriangles(0.5f);

    std::vector<std::vector<double>> zMaxMap(samplingRate, std::vector<double>(samplingRate));
    
    std::vector<std::array<arraydouble3, 3>> tTriangles = mesh.transformedTriangles(tMatrix);

    mesh.shiftTriangles(tTriangles, {0, 0, -zMaxThreshold});
    mesh.insideTriangles(tTriangles, 0, size/2, size/2);
    double epsilon{1e-7};
    mesh.nonOrthogonalTriangles(tTriangles, epsilon);

    std::vector<arraydouble3> normals;
    std::vector<double> Ds;
    int N = samples.size();
    int tN = tTriangles.size();
    for (int i = 0; i < tN; i++)
    {
        arraydouble3 triNormal = triangleNormal(tTriangles[i]);
        normals.push_back(triNormal);
        Ds.push_back(triangleD(tTriangles[i], triNormal));
    }
    double z = 0.0f;
    std::vector<double> zMapData(N*N, 0.0f);
    for (int i = 0; i < tN; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
            {
                z = zFromXY({samples[j], samples[k]}, normals[i], Ds[i]);
                z += zMaxThreshold;
                if (z > zMapData[j*N + k])
                    zMapData[j*N + k] = z;
            }
    return ZMap(std::move(zMapData), tMatrix, size, samplingRate);
}