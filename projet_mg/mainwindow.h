#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QDesktopServices>
#include <QProcess>
#include <QMessageBox>
#include <cstdlib>
#include "delaunay.h"

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

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

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void displayMeshes(QVector<MyMesh *> vector_meshes, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);
    MyMesh * cut(MyMesh * _mesh, MyMesh::Point P1, MyMesh::Point P2, MyMesh::Point P3);
    float * equation_plane(MyMesh::Point P1
                        , MyMesh::Point P2
                        , MyMesh::Point P3);
    void showBorder(MyMesh* _mesh,float * eq_plane) ;
    //void changeValuePointCut(OpenMesh::Vec3f first, OpenMesh::Vec3f last);

    void boite_englobante();
    bool pointIsInMesh(MyMesh::Point point);
    bool createSeeds();
    int intersect(FaceHandle triangle, MyMesh::Point point_extern, MyMesh::Point point);
public slots:
    void changeValuePointCut(OpenMesh::Vec3f first, OpenMesh::Vec3f last);
private slots:
    void on_pushButton_chargement_clicked();
    void on_pushButton_generer_clicked();
    void on_pushButton_couper_clicked();

    void on_pushButton_clicked();

    void on_showSelections_clicked();

    void on_vertexSelection_valueChanged(int arg1);

    void on_pushButton_save_clicked();

private:

    MyMesh mesh;

    Ui::MainWindow *ui;
    float be_min_x = 0.0f;
    float be_max_x = 0.0f;
    float be_min_y = 0.0f;
    float be_max_y = 0.0f;
    float be_min_z = 0.0f;
    float be_max_z = 0.0f;
    int vertexSelection = -1;

    std::vector<MyMesh::Point> seeds_impact;
    std::vector<MyMesh> list_meshes_to_save;
    unsigned nb_seeds = 100;
};

#endif // MAINWINDOW_H
