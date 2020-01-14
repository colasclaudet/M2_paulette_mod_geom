#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QDebug>
#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QDesktopServices>
#include <QProcess>

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

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    /* ***Affichage et objets*** */
    void resetAllColorsAndThickness(MyMesh* _mesh);
    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void displayMeshes(QVector<MyMesh *> vector_meshes, bool isTemperatureMap = false, float mapRange = -1);
    void showBorder(MyMesh* _mesh/*,float * eq_plane*/) ;
    void chooseFragToDisplay();
    void chooseFragToSave();

    /* ***Calcul et manipulation*** */
    double *equation_droite(double *eq1, double *eq2);
    double *equation_plane(MyMesh::Point P1, MyMesh::Point P2, MyMesh::Point P3);
    MyMesh::Point resol_3eq_3inc(double *eq1, double *eq2, double *eq3);
    //MyMesh * cut(MyMesh * _mesh, MyMesh::Point P1, MyMesh::Point P2, MyMesh::Point P3); //obsolete
    MyMesh *cutV2(MyMesh *_mesh, MyMesh::Point P1, MyMesh::Point P2, MyMesh::Point P3);
    void boite_englobante();
    int intersect(FaceHandle triangle, MyMesh::Point point_extern, MyMesh::Point point);
    bool pointIsInMesh(MyMesh::Point point);
    void createSeeds();
    std::vector<double> crossProduct(MyMesh::Point p1, MyMesh::Point p2);
    void divide();
    void regular_divide(int divisor);
private slots:
    /* * MENU  * */
    void importFile();
    void generer_cube();
    /* *FIN MENU * */
    void changeValuePointCut(OpenMesh::Vec3f first, OpenMesh::Vec3f last);
    void on_pushButton_couper_clicked();
    void on_pushButton_clicked(); //création seed
    void on_checkBox_display_frag1_clicked();
    void on_checkBox_display_frag2_clicked();
    void on_checkBox_save_frag1_clicked();
    void on_checkBox_save_frag2_clicked();
    void on_pushButton_save_clicked();
    //inutile
    void on_showSelections_clicked();
    void on_vertexSelection_valueChanged(int arg1);
    void on_pushButton_chargement_clicked();
    void on_pushButton_generer_clicked();
private:
    MyMesh mesh; //maillage de base
    QVector<MyMesh> cubes; //liste des cubes apres la division
    MyMesh cube;
    MyMesh seed; //maillage comprennant seulement
    Ui::MainWindow *ui;
    //valeurs de la boite englobante
    float be_min_x = 0.0f;
    float be_max_x = 0.0f;
    float be_min_y = 0.0f;
    float be_max_y = 0.0f;
    float be_min_z = 0.0f;
    float be_max_z = 0.0f;

    int vertexSelection = -1;
    std::vector<MyMesh> list_meshes_to_save;
    QVector<MyMesh*> list_meshes_to_display;
    MyMesh frag1;
    MyMesh frag2;
    MyMesh frag3; //mesh comprenant seulement les points d'intersections

    std::vector<MyMesh::Point> seeds_impact;
    unsigned nb_seeds = 1000;
};

#endif // MAINWINDOW_H
