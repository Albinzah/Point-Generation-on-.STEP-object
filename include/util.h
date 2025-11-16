#pragma once
#include "types.h"
#include <string>
#include <TopoDS_Shape.hxx>

using namespace types;

namespace util
{
    TopoDS_Shape loadStepShape(const std::string& path);

    int displayShape(const TopoDS_Shape shape, const PointSet samples);

    void extractMeshData(const TopoDS_Shape& shape, std::vector<arraydouble3>& outPoints, std::vector<std::array<int, 3>>& outTriangles);
}