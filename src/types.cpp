#include "types.h"
#include "geometryUtils.h"

#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

#include <stdexcept>
#include <vector>
#include <array>
#include <numbers>
#include <iostream>

#include <omp.h>

// TODO: Validate that 'data' has samplingRate x samplingRate shape.

using namespace geometryUtils;

namespace types {
using Matrix4X4 = std::array<std::array<double, 4>, 4>;

ZMap::ZMap(std::vector<double> data, const Matrix4X4 matrix, double size, int samplingRate)
    : data(std::move(data)), transformationMatrix(matrix), size(size), samplingRate(samplingRate) 
    {
        stepSize = size / samplingRate;
    }

Matrix4X4 ZMap::rotate(double theta)
{
    double Rz[4][4] = {
            {cosf(theta),-sinf(theta), 0, 0}, 
            {sinf(theta), cosf(theta), 0, 0},
            {0,             0,             1, 0},
            {0,             0,             0, 1}
    };

    Matrix4X4 newTMatrix{};
    
    // This should use loop reduction + take a look at loop rearangement.
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            for (int k=0; k<4; k++)
                newTMatrix[i][j] += transformationMatrix[i][k] * Rz[k][j];

    return newTMatrix;
}

std::vector<double> ZMap::extractPatch(double radius)
{
    if (data.empty())
        throw std::invalid_argument("extractPatch: 'data' is empty.");
    
    if (radius > (size / 2))
        throw std::invalid_argument("extractPatch: 'radius' exceeds half the z-map size.");

    double center   {(samplingRate - 1) / 2.0f};

    std::vector<double> maskedData;
    maskedData.reserve(samplingRate * samplingRate);

    double distance {0.0f};

    for (int i=0; i<samplingRate; i++)
        for (int j=0; j<samplingRate; j++)
        {
            distance = sqrtf((i - center)*(i - center) + (j - center) * (j - center)) * stepSize;
            if (distance <= radius)
                maskedData.emplace_back(data[i * samplingRate + j]);
        }

    return maskedData;
}

std::vector<std::vector<double>> ZMap::extractRectanglePatches(const double width, const double height, const int rotations)
{
    int patchWidthPx  = static_cast<int>(std::ceil(width  / stepSize));
    int patchHeightPx = static_cast<int>(std::ceil(height / stepSize));

    if (patchWidthPx > samplingRate)
        throw std::invalid_argument("extractRectanglePatches: Patch 'width' larger than sampling size.");
    if (patchHeightPx > samplingRate)
        throw std::invalid_argument("extractRectanglePatches: Patch 'height' larger than sampling size.");

    // This should probably round up instead?
    int widthOffset  = (samplingRate - patchWidthPx) / 2;
    int heightOffset = (samplingRate - patchHeightPx) / 2;

    std::vector<std::vector<bool>> mask(samplingRate, std::vector<bool>(samplingRate, false));

    for (int i=widthOffset; i<samplingRate-widthOffset; i++)
        for (int j=heightOffset; j<samplingRate-heightOffset; j++)
            mask[j][i] = true;

    double rotationsStep {2 * std::numbers::pi_v<double> / rotations};

    std::vector<std::vector<double>> patches(rotations);
    
    for (auto& v : patches)
        v.reserve(samplingRate*samplingRate);
    
    for (int idx=0; idx < rotations; idx++)
    {    
        double angle {rotationsStep * idx};
        auto rotatedMask = rotateMask(mask, angle);
        for (int i=0; i<samplingRate; i++)
            for (int j=0; j<samplingRate; j++)
                if (rotatedMask[i][j])
                    patches[idx].emplace_back(data[i * samplingRate + j]);
    }
    
    return patches;
}

bool ZMap::isFlat(double penetrationThreshold, std::vector<double> patch)
{
    for (auto& v : patch)
        if (std::fabs(v) > penetrationThreshold)
            return false;
    return true;
}

bool ZMap::isCovered(double coverageThreshold, double penetrationThreshold, std::vector<double> patch)
{
    double peakValue{patch[0]};
    for (auto& v : patch)
        if (v > peakValue)
            peakValue = v;

    if (peakValue > penetrationThreshold)
        return false;

    std::vector<double> difference(patch);

    const int N = patch.size();

    for (int i=0; i < N; i++)
        difference[i] = std::fabs(patch[i] - peakValue);

    int validSamples = 0;

    for (auto& d : difference)
        if (d <= penetrationThreshold)
            validSamples++;
    
    double fraction = static_cast<double>(validSamples) / patch.size();

    return (fraction >= coverageThreshold);
}

std::vector<std::vector<bool>> ZMap::rotateMask(const std::vector<std::vector<bool>> matrix, double theta)
{
    if (matrix.size() != matrix[0].size())
        throw std::invalid_argument("rotateMatrix: 'matrix' is not square.");
    if (matrix.empty())
        throw std::invalid_argument("rotateMatrix: 'matrix' is empty.");
    
    int N = static_cast<int>(matrix.size());
    double center = (N - 1) / 2.0f;

    std::vector<std::vector<bool>> rotatedMatrix(N, std::vector<bool>(N, false));

    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
        {
            int x = j - center;
            int y = i - center;

            double newX = x * cos(theta) - y * sin(theta);
            double newY = x * sin(theta) + y * cos(theta);

            int newI = static_cast<int>(std::round(newX + center));
            int newJ = static_cast<int>(std::round(newY + center));

            if (0 <= newI && newI < N && 0 <= newJ && newJ < N)
                rotatedMatrix[newI][newJ] = matrix[i][j];
        }
    
    return rotatedMatrix;
}

// PointSample
PointSample::PointSample(arraydouble3 position, arraydouble3 normal, ZMap zMap) 
    : position(position), normal(normal), zMap(zMap) {};

double PointSample::radialDistanceToCog(arraydouble3 cog)
{
    arraydouble3 displacement = cog - position;
    double d = displacement.dot(normal);
    arraydouble3 projected = displacement - normal * d;
    return std::sqrt(projected.dot(projected));
}

const arraydouble3& PointSample::getPosition() const { return position; };
const arraydouble3& PointSample::getNormal() const { return normal; };


// PointSet

PointSet::PointSet(std::vector<arraydouble3> points,
                   std::vector<arraydouble3> normals,
                   double zMapSize, int samplingRate,
                   const Mesh& mesh, arraydouble3 cog, double mass)
    : samples(), cog(cog), mass(mass), zMapSize(zMapSize), samplingRate(samplingRate)
{
    const int N = points.size();
    samples.reserve(N);
    for (int i=0; i<N; i++)
    {
        ZMap zmap = meshZMap(points[i], normals[i], zMapSize, mesh, samplingRate, zThreshold);
        samples.emplace_back(points[i], normals[i], zmap);
    };
}

const std::vector<PointSample>& PointSet::getSamples() const { return samples; };


// Shape
Shape::Shape(TopoDS_Shape shape, Mesh mesh, double pointDistance, double mass, arraydouble3 cog)
    : shape(std::move(shape)),
      mesh(std::move(mesh)),
      pointDistance(pointDistance),
      mass(mass),
      cog(cog),
      samples()
{
    samples = generateSamples();
}

PointSet Shape::generateSamples()
{
    std::vector<arraydouble3> points;
    std::vector<arraydouble3> normals;

    int loopCount = 0;
    double startLoopT;
    double loopTime = 0;
    double genStart;
    double genTime = 0;
    double filterStart;
    double filterTime = 0;

    for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {

        const TopoDS_Face& face = TopoDS::Face(exp.Current());

        std::vector<std::array<double, 2>> uvSamples = geometryUtils::generateUniformPointsOnFace(face, pointDistance);
        
        std::vector<std::array<double, 2>> filteredSamples = geometryUtils::filterPointsOutsideFace(uvSamples, face);

        std::vector<arraydouble3> samplePoints;
        std::vector<arraydouble3> sampleNormals;

        geometryUtils::setPositionAndNormals(filteredSamples, face, samplePoints, sampleNormals);

        points.insert(points.end(), samplePoints.begin(), samplePoints.end());
        normals.insert(normals.end(), sampleNormals.begin(), sampleNormals.end());
    }

    const double zMapSize{10.0f};
    const int samplingRate{20};

    return PointSet(points, normals, zMapSize, samplingRate, mesh, cog, mass);
}

// Mesh

Mesh::Mesh(std::vector<arraydouble3> points, std::vector<std::array<int, 3>> triangles)
    : _points(points), _triangles(triangles) {}

[[nodiscard]] int Mesh::nbTriangles() const noexcept { return static_cast<int>(_triangles.size()); } 

double Mesh::triangleArea(std::array<int, 3> triangle) 
{
    arraydouble3 v1 = _points[triangle[0]];
    arraydouble3 v2 = _points[triangle[1]];
    arraydouble3 v3 = _points[triangle[2]];

    arraydouble3 v1v2 = v2 - v1;
    arraydouble3 v1v3 = v3 - v1;

    arraydouble3 c = v1v2.cross(v1v3);
    
    double area = 0.5f * sqrtf(c.dot(c));

    return area;
}

void Mesh::removeTinyTriangles(double minArea)
{
    std::vector<std::array<int, 3>> largeTriangles{};
    for (auto& t : _triangles)
        if (triangleArea(t) > minArea)
            largeTriangles.emplace_back(t);

    _triangles = largeTriangles;
}

std::vector<std::array<arraydouble3, 3>> Mesh::transformedTriangles(const std::array<std::array<double, 4>, 4>& transformMatrix)
{   
    const int N = _triangles.size();
    std::array<arraydouble3, 3> rotationMatrix = {{ {transformMatrix[0][0], transformMatrix[0][1], transformMatrix[0][2]},
                                                   {transformMatrix[1][0], transformMatrix[1][1], transformMatrix[1][2]},
                                                   {transformMatrix[2][0], transformMatrix[2][1], transformMatrix[2][2]}
                                                }};

    std::vector<std::array<arraydouble3, 3>> transformedTriangles(N);
    arraydouble3 p1, p2, p3;

    for (int i = 0; i < N; i++)
    {
        p1 = dot3x1(rotationMatrix, _points[_triangles[i][0]]);
        p2 = dot3x1(rotationMatrix, _points[_triangles[i][1]]);
        p3 = dot3x1(rotationMatrix, _points[_triangles[i][2]]);

        transformedTriangles[i] = {p1, p2, p3};
    }

    return transformedTriangles;
}

void Mesh::shiftTriangles(std::vector<std::array<arraydouble3, 3>>& triangles, arraydouble3 direction)
{
    for (auto& tri : triangles)
    {
        tri[0] = tri[0] + direction;
        tri[1] = tri[1] + direction;
        tri[2] = tri[2] + direction;
    }
}

void Mesh::insideTriangles(std::vector<std::array<arraydouble3, 3>>& triangles, double xRange, double yRange, double zRange)
{
    std::erase_if(triangles, [&](const auto& t) {
        for (auto& v : t)
            if (fabs(v.x) > xRange || fabs(v.y) > yRange || fabs(v.z) > zRange)
                return true;
        return false;
    });
}
void Mesh::nonOrthogonalTriangles(std::vector<std::array<arraydouble3, 3>>& triangles, double epsilon)
{
    std::erase_if(triangles, [&](const std::array<arraydouble3, 3>& t)
    {
        arraydouble3 normal = triangleNormal(t);
        return std::fabs(normal.z) <= epsilon;
    });
}
}