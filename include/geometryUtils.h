#pragma once
#include <array>
#include <vector>
#include <TopoDS_Face.hxx>
#include "types.h"
#include <BRepAdaptor_Surface.hxx>

using namespace types;

namespace geometryUtils
{
types::ZMap meshZMap( arraydouble3 point, arraydouble3 normal, 
                      double size, Mesh mesh, int samplingRate, 
                      double zMaxThreshold );

std::array<arraydouble3, 3> dot3x3(const std::array<arraydouble3, 3>& m1, const std::array<arraydouble3, 3>& m2);

arraydouble3 dot3x1(const std::array<arraydouble3, 3>& m1, const arraydouble3& m2);

std::array<arraydouble3, 3> rotationMatrixToZ(arraydouble3 normal);

std::array<std::array<double, 4>, 4> transformMatrixToZAxis(arraydouble3 normal, arraydouble3 point);

arraydouble3 triangleNormal(std::array<arraydouble3, 3> triangle);

double triangleD(std::array<arraydouble3, 3> triangle);

bool isPointInXYTriangle(std::array<double, 2> point, std::array<arraydouble3, 3> triangle, double epsilon = 1e-7);

std::vector<std::array<double, 2>> generateUniformPointsOnFace(TopoDS_Face face, double resolution);

double zFromXY(std::array<double, 2> P, arraydouble3 normal, double D);

std::vector<std::array<double, 2>> generateNonUniformRevolutionUV(
        BRepAdaptor_Surface adaptorSurface, double uMin,
        double vMin, double vMax, 
        double uRange, double vRange, 
        double resolution, double epsilon = 1e-6);

std::vector<std::array<double, 2>> generateUniformUV(
    double uMin, double uMax, double vMin, 
    double vMax, double uRange, double vRange, 
    double resolution, double epsilon = 1e-6);

void setPositionAndNormals(
        const std::vector<std::array<double, 2>>& uvSamples,
        const TopoDS_Face& face,
        std::vector<arraydouble3>& outPoints,
        std::vector<arraydouble3>& outNormals);

arraydouble3 uvPositionToGlobal(TopoDS_Face face, double u, double v);

arraydouble3 uvNormal(TopoDS_Face face, double u, double v);

std::vector<std::array<double, 2>> filterPointsOutsideFace(std::vector<std::array<double, 2>> uvPoints, TopoDS_Face face);

}

