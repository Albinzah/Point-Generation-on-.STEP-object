#include <iostream>
#include <string>
#include <BRepMesh_IncrementalMesh.hxx>

#include "include/util.h"
#include "include/geometryUtils.h"
#include "include/types.h"
#include "include/debugVars.h"

#include <omp.h>


int main()
{
    double fullstart = omp_get_wtime();
    std::cout << "---" << std::endl;

    std::string filePath = "/mnt/c/Users/Azahn/Desktop/Random/Version2.0/assets/cads/silver_box.stp";
    TopoDS_Shape rootShape = util::loadStepShape(filePath);

    BRepMesh_IncrementalMesh mesher(rootShape, 0.5);

    std::vector<arraydouble3> points;
    std::vector<std::array<int, 3>> triangles;

    util::extractMeshData(rootShape, points, triangles);

    types::Mesh mesh(points, triangles);

    const double pointDistance = 7;
    const double mass = 0;
    const arraydouble3 cog{0.0f, 0.0f, 0.0f};
    
    double shapeStart = omp_get_wtime();
    Shape shape(rootShape, mesh, pointDistance, mass, cog);
    double shapeEnd = omp_get_wtime();

    PointSet sampleSet = shape.getPointSet();
    std::vector<types::PointSample> samples = sampleSet.getSamples();


    std::cout << "========= Main Summary =========" << '\n';
    std::cout << "Full time:  " << omp_get_wtime() - fullstart << "s\n";
    std::cout << "Shape time: " << shapeEnd - shapeStart << "s\n";
    std::cout << "================================" << '\n';

    // util::displayShape(rootShape, sampleSet);

    return 0;
}