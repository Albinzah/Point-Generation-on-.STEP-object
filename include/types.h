#pragma once

#include <stdexcept>
#include <vector>
#include <array>
#include <cmath>
#include <ios>
#include <iomanip>
#include <TopoDS_Shape.hxx>

namespace types
{

struct arraydouble3
{
    double x{};
    double y{};
    double z{};

    arraydouble3() = default;
    arraydouble3(double x, double y, double z) : x(x), y(y), z(z) {}

    operator std::array<double, 3>() const {
        return {x, y, z};
    }

    explicit arraydouble3(const std::array<double, 3>& arr)
        : x(arr[0]), y(arr[1]), z(arr[2]) {}

    double dot(const arraydouble3& other)
    {
        return (x * other.x + y * other.y + z * other.z);
    }

    arraydouble3 cross(const arraydouble3& other) const
    {
        return { y*other.z - z*other.y,
                 z*other.x - x*other.z,
                 x*other.y - y*other.x };
    }
    
    double& operator[](int index)
    {
        switch (index)
        {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default:
                throw std::out_of_range("arraydouble3 index out of range");
        }
    }
    const double& operator[](int index) const
    {
        switch (index)
        {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default:
                throw std::out_of_range("arraydouble3 index out of range");
        }
    }
    arraydouble3 operator+(const arraydouble3& other) const
    {
        return {x + other.x, y + other.y, z + other.z};
    }
    arraydouble3 operator-(const arraydouble3& other) const
    {
        return {x - other.x, y - other.y, z - other.z};
    }
    arraydouble3 operator*(double s) const
    {
        return {x * s, y * s, z * s};
    }
    arraydouble3 operator/(double s) const
    {
        return {x / s, y / s, z / s};
    }
    arraydouble3& operator+=(const arraydouble3& other)
    {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }
    arraydouble3& operator-=(const arraydouble3& other)
    {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }
    arraydouble3& operator*=(double s)
    {
        x *= s; y *= s; z *= s;
        return *this;
    }
    arraydouble3& operator/=(double s)
    {
        x /= s; y /= s; z /= s;
        return *this;
    }
    arraydouble3 operator-() const
    {
        return {-x, -y, -z};
    }

    bool equals(const arraydouble3& other, double epsilon=0.0f) const
    {
        return std::fabs(x - other.x) < epsilon &&
               std::fabs(y - other.y) < epsilon &&
               std::fabs(z - other.z) < epsilon;
    }
};

inline std::ostream& operator<<(std::ostream& os, const arraydouble3& v) {
    os << std::fixed;
    os.precision(3);
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

class Mesh
{
private:
    std::vector<arraydouble3> _points;
    std::vector<std::array<int, 3>> _triangles;
public:
    Mesh(std::vector<arraydouble3> points, std::vector<std::array<int, 3>> triangles);

    [[nodiscard]] int nbTriangles() const noexcept;
    void removeTinyTriangles(double minArea);
    void shiftTriangles(std::vector<std::array<arraydouble3, 3>>& triangles, arraydouble3 direction);
    void insideTriangles(std::vector<std::array<arraydouble3, 3>>& triangles, double xRange, double yRange, double zRange);
    void nonOrthogonalTriangles(std::vector<std::array<arraydouble3, 3>>& triangles, double epsilon);
    double triangleArea(std::array<int, 3> triangle);
    std::vector<std::array<arraydouble3, 3>> transformedTriangles(const std::array<std::array<double, 4>, 4>& transformMatrix);
};

class ZMap
{
using Matrix4X4 = std::array<std::array<double, 4>, 4>;
private:
    std::vector<double> data;
    Matrix4X4 transformationMatrix;
    double size;
    int samplingRate;
    double stepSize;

    std::vector<std::vector<bool>> rotateMask(const std::vector<std::vector<bool>> matrix, double theta);

public:
    // Constructor
    ZMap(std::vector<double> data, const Matrix4X4 matrix, double size, int samplingRate);

    // Getters
    [[nodiscard]] const std::vector<double>& getData() const noexcept { return data; }
    [[nodiscard]] std::vector<double>& getData() noexcept { return data; }

    // Operations
    Matrix4X4 rotate(double theta);
    std::vector<double> extractPatch(double radius);
    std::vector<std::vector<double>> extractRectanglePatches(const double width, const double height, const int rotations = 4);
    bool isFlat(double penetrationThreshold, std::vector<double> patch);
    bool isCovered(double coverageThreshold, double penetrationThreshold, std::vector<double> patch);
};

class PointSample
{
private:
    arraydouble3 position;
    arraydouble3 normal;
    ZMap zMap;
public:
    PointSample(arraydouble3 position, arraydouble3 normal, ZMap zMap);
    const arraydouble3& getPosition() const;
    const arraydouble3& getNormal() const;

    double radialDistanceToCog(arraydouble3 cog);
};

class PointSet
{
private:
    std::vector<PointSample> samples;
    arraydouble3 cog;
    double mass;
    double zMapSize;
    int samplingRate;
    double zThreshold = 1e5f;
public: 
    PointSet() = default;
    PointSet( std::vector<arraydouble3> points, std::vector<arraydouble3> normals, 
              double zMapSize, int samplingRate, const Mesh& mesh, 
              arraydouble3 cog, double mass);
    const std::vector<PointSample>& getSamples() const;
};

class Shape
{
private:
    TopoDS_Shape shape;
    Mesh mesh;
    double pointDistance;
    double mass;
    arraydouble3 cog;
    PointSet samples;

    PointSet generateSamples();
public:
    Shape(TopoDS_Shape shape, Mesh mesh, double pointDistance, double mass, arraydouble3 cog);
    const PointSet& getPointSet() const { return samples; }
    
};

}