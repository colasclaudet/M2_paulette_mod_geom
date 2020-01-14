#ifndef DELAUNAY_H
#define DELAUNAY_H
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QDesktopServices>
#include <QProcess>
#include <QDebug>
#include "triangle.h"
#include "mainwindow.h"

struct Edge
{
    MyMesh::Point p1;
    MyMesh::Point p2;
    bool isBad = false;
};
class Delaunay
{
public:
    Delaunay(MyMesh * _mesh);
    std::vector<Triangle> triangulate();
    std::vector<Triangle> triangulate2();
    std::vector<Triangle> getTriangles() const;


    MyMesh draw_mesh();

private:
    MyMesh * mesh;
    std::vector<Triangle> triangles;

    std::vector<Edge> edges;
    std::vector<MyMesh::Point> _vertices;
};


#endif // DELAUNAY_H
