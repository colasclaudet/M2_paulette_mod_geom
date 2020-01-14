#ifndef MESHVIEWERWIDGET_H
#define MESHVIEWERWIDGET_H

#define GL_GLEXT_PROTOTYPES
#include <QDebug>
#include <QGLWidget>
#include "QMouseEvent"
#include <QMainWindow>
#include <OpenMesh/Core/Geometry/VectorT.hh>


#ifdef __APPLE__
    #include <gl.h>
    #include <glu.h>
    #include <glext.h>
#else
    #include <GL/gl.h>
    #include <GL/glu.h>
    #include <GL/glx.h>
    #include <GL/glext.h>
    #include <GL/glut.h>
#endif



using namespace OpenMesh;

const double TRACKBALL_RADIUS = 0.6;

class MeshViewerWidget : public QGLWidget
{
    Q_OBJECT
public:
    MeshViewerWidget(QWidget* _parent=0);
    MeshViewerWidget(QGLFormat& _fmt, QWidget* _parent );

    // events GL
    void initializeGL();
    void resizeGL(int _w, int _h);
    void paintGL();

    // gestion de la vue et de la trackball
    void update_projection_matrix();
    void view_all();
    void set_scene_pos( const OpenMesh::Vec3f& _cog, float _radius );
    bool map_to_sphere(const QPoint& _point, OpenMesh::Vec3f& _result);
    void translate( const OpenMesh::Vec3f& _trans );
    void rotate( const OpenMesh::Vec3f& _axis, float _angle );

    void reloadPOV();
    void loadMesh(GLfloat* verts, GLfloat* colors, int nVerts, GLuint* triangles, int nTriangles);
    void loadLines(GLfloat* verts, GLfloat* colors, int nVerts, GLuint* lines, int nLines, QList<QPair<float, int> > es);
    void loadPoints(GLfloat* verts, GLfloat* colors, int nVerts, GLuint* points, int nPoints, QList<QPair<float, int> > vs);


    // buffer pour les VBO
    GLuint TriDataBuffers[2];
    // Nombre de triangles du mesh (pour le call du draw)
    int triToDraw;

    GLuint LinesDataBuffers[2];
    int linesToDraw;
    QList<QPair<float, int> > edgeSizes;

    GLuint PointsDataBuffers[2];
    int pointsToDraw;
    QList<QPair<float, int> > vertsSizes;

    // pour savoir si les buffer sont init
    bool init;

    // variables de gestion de la vue et de la trackball
    OpenMesh::Vec3f  center_;
    OpenMesh::Vec3f  user_position;
    float            radius_;

    GLdouble    projection_matrix_[16], modelview_matrix_[16];

    QPoint           last_point_2D_;
    OpenMesh::Vec3f  last_point_3D_;
    bool             last_point_ok_;

    //coordonnées des points pour la decoupe
    OpenMesh::Vec3f  first_point_to_cut3D;
    OpenMesh::Vec3f  last_point_to_cut3D;
    OpenMesh::Vec3f center_cut3D;
    float Rx[3][3];
    float Ry[3][3];
    float Rz[3][3];
    float angle_h = 0.0;
    float angle_v = 0.0;
    float angle_p = 0.0;
    Vec3f axis_h;
    Vec3f axis_v;
    Vec3f axis_p;
    QWidget* parent;
protected:

    // Qt mouse events
    virtual void mousePressEvent( QMouseEvent* );
    virtual void mouseReleaseEvent( QMouseEvent* );
    virtual void mouseMoveEvent( QMouseEvent* );
    virtual void wheelEvent( QWheelEvent* );
signals:
    void sig_changePoint(OpenMesh::Vec3f, OpenMesh::Vec3f);

};

#endif // MESHVIEWERWIDGET_H
