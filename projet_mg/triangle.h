#ifndef TRIANGLE_H
#define TRIANGLE_H
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QDesktopServices>
#include <QProcess>
#include <QDebug>

#include "mainwindow.h"
using namespace OpenMesh;
using namespace OpenMesh::Attributes;
class Triangle
{

public:
    Triangle(MyMesh::Point v1, MyMesh::Point v2, MyMesh::Point v3);
    bool containsVertex(MyMesh::Point v);
    MyMesh::Point getNorm();
    MyMesh::Point getVect(MyMesh::Point p1, MyMesh::Point p2);
    float * equation_plane(MyMesh::Point P1
                        , MyMesh::Point P2
                        , MyMesh::Point P3);
    bool circumCircleContains3D(MyMesh::Point v);
    bool circumCircleContains2D(MyMesh::Point v);
    bool operator ==(const Triangle tr);
    bool containVertex(MyMesh::Point p);
    bool almost_equal(const Triangle t1, const Triangle t2);
    bool almost_equal(const double x, const double y, int ulp=2);
    bool almost_equal(MyMesh::Point v1, MyMesh::Point v2, int ulp=2);

    MyMesh::Point getPoint1() const;
    void setPoint1(const MyMesh::Point &value);

    MyMesh::Point getPoint2() const;
    void setPoint2(const MyMesh::Point &value);

    MyMesh::Point getPoint3() const;
    void setPoint3(const MyMesh::Point &value);
    double dist2(MyMesh::Point v1,MyMesh::Point v2);
    MyMesh::Point return_gravity_center();

    bool isBad = false;

private:
    MyMesh::Point point1;
    MyMesh::Point point2;
    MyMesh::Point point3;

};

#endif // TRIANGLE_H
