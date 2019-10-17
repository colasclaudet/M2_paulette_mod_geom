#include "mainwindow.h"
#include "ui_mainwindow.h"


/* **** début de la partie boutons et IHM **** */


// exemple pour charger un fichier .obj
void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

// exemple pour construire un mesh face par face
void MainWindow::on_pushButton_generer_clicked()
{
    MyMesh mesh;

    // on construit une liste de sommets
    MyMesh::VertexHandle sommets[8];
    sommets[0] = mesh.add_vertex(MyMesh::Point(-1, -1,  1));
    sommets[1] = mesh.add_vertex(MyMesh::Point( 1, -1,  1));
    sommets[2] = mesh.add_vertex(MyMesh::Point( 1,  1,  1));
    sommets[3] = mesh.add_vertex(MyMesh::Point(-1,  1,  1));
    sommets[4] = mesh.add_vertex(MyMesh::Point(-1, -1, -1));
    sommets[5] = mesh.add_vertex(MyMesh::Point( 1, -1, -1));
    sommets[6] = mesh.add_vertex(MyMesh::Point( 1,  1, -1));
    sommets[7] = mesh.add_vertex(MyMesh::Point(-1,  1, -1));


    // on construit des faces à partir des sommets

    std::vector<MyMesh::VertexHandle> uneNouvelleFace;

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[3]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[7]);
    uneNouvelleFace.push_back(sommets[6]);
    uneNouvelleFace.push_back(sommets[5]);
    uneNouvelleFace.push_back(sommets[4]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[4]);
    uneNouvelleFace.push_back(sommets[5]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[5]);
    uneNouvelleFace.push_back(sommets[6]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[3]);
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[6]);
    uneNouvelleFace.push_back(sommets[7]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[3]);
    uneNouvelleFace.push_back(sommets[7]);
    uneNouvelleFace.push_back(sommets[4]);
    mesh.add_face(uneNouvelleFace);

    mesh.update_normals();
    this->mesh = mesh;
    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);

}

/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

float * MainWindow::equation_plane(MyMesh::Point P1
                    , MyMesh::Point P2
                    , MyMesh::Point P3)
{
    float * eq_plane = new float[4];
    float a1 = P2[0] - P1[0];
    float b1 = P2[1] - P1[1];
    float c1 = P2[2] - P1[2];
    float a2 = P3[0] - P1[0];
    float b2 = P3[1] - P1[1];
    float c2 = P3[2] - P1[2];
    float a = b1 * c2 - b2 * c1;
    float b = a2 * c1 - a1 * c2;
    float c = a1 * b2 - b1 * a2;
    float d = (- a * P1[0] - b * P1[1] - c * P1[2]);

    qDebug() << "equation of plane is " << a << " x + " << b
        << " y + " << c << " z + " << d << " = 0.";
    eq_plane[0]=a;
    eq_plane[1]=b;
    eq_plane[2]=c;
    eq_plane[3]=d;
    return eq_plane;
}

MyMesh * MainWindow::cut(MyMesh * _mesh, MyMesh::Point P1, MyMesh::Point P2, MyMesh::Point P3)
{
    MyMesh frag1;
    MyMesh frag2;
    /*MyMesh::VertexHandle vh1 = _mesh->add_vertex(P1);
    MyMesh::VertexHandle vh2 = _mesh->add_vertex(P2);
    MyMesh::VertexHandle vh3 = _mesh->add_vertex(P3);*/

    float * eq_plane = equation_plane(P1,P2,P3);
    /*for(MyMesh::VertexIter vi;vi!=_mesh->vertices_end(); vi++)
    {

    }*/
    /*frag1.add_vertex(P1);
    frag1.add_vertex(P2);
    frag1.add_vertex(P3);

    frag2.add_vertex(P1);
    frag2.add_vertex(P2);
    frag2.add_vertex(P3);
    */

    for(int i=0;i < _mesh->n_vertices();i++)
    {
        MyMesh::Point PM = _mesh->point(_mesh->vertex_handle(i));

        if((PM[0]*eq_plane[0]+PM[1]*eq_plane[1]+PM[2]*eq_plane[2]+eq_plane[3])>=0)
        {
            frag1.add_vertex(PM);
        }
        else
        {
            frag2.add_vertex(PM);
        }
    }
    MyMesh * tabmesh = new MyMesh[2];
    tabmesh[0] = frag1;
    tabmesh[1] = frag2;
    return tabmesh;
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->displayWidget, SIGNAL(sig_changePoint(OpenMesh::Vec3f, OpenMesh::Vec3f)), this, SLOT(changeValuePointCut(OpenMesh::Vec3f, OpenMesh::Vec3f)));
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::showBorder(MyMesh* _mesh, float * eq_plane) //<affiche bordure> met plusieurs seconde à s'executer
{
    resetAllColorsAndThickness(_mesh);
    std::vector<MyMesh::VertexHandle> vector_vertex;
    int count = 0;
    for(MyMesh::EdgeIter currentEdge = _mesh->edges_begin(); currentEdge != _mesh->edges_end(); currentEdge ++)
    {
        EdgeHandle e = *currentEdge;
        //_mesh->set_color(e, MyMesh::Color(0,255,0));*/
        //On parcourt toutes les faces du mesh
        for(MyMesh::FaceIter f = _mesh->faces_begin(); f != _mesh->faces_end(); f++)
        {
            //On parcourt toutes les aretes sur chacune des faces
            for(MyMesh::FaceEdgeIter current_edgeFace = _mesh->fe_iter(f); current_edgeFace.is_valid(); current_edgeFace ++)
            {
                EdgeHandle ce = *current_edgeFace;
                if(e == ce) //si notre arrête selectionné égale à l'arrete courante de l'une des faces
                    count ++;
            }

        }
        if(count == 1)
        {
            MyMesh::HalfedgeHandle heh1 = _mesh->halfedge_handle(e,0);

            MyMesh::VertexHandle V1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(e,0));
            MyMesh::VertexHandle V2 = _mesh->to_vertex_handle(_mesh->halfedge_handle(e,1));
            vector_vertex.push_back(V1);
            vector_vertex.push_back(V2);
            //_mesh->set_color(e, MyMesh::Color(200,0,170));
            //_mesh->data(e).thickness = 6;
        }
        count = 0;
    }
    // on affiche le nouveau maillage
    /*
     std::vector<MyMesh::VertexHandle> vector_vertex_temp;
    for(int i = 0; i<vector_vertex.size()/3-3; i++)
    {
        for(int j = 0;j<3;j++)
        {
            vector_vertex_temp.push_back(vector_vertex.at(i*3+j));
        }
        _mesh->add_face(vector_vertex_temp);
        vector_vertex_temp.clear();
    }*/
    _mesh->add_face(vector_vertex);
    displayMesh(_mesh);
}

void MainWindow::changeValuePointCut(OpenMesh::Vec3f first, OpenMesh::Vec3f last)
{
    ui->lb_firstX->setNum(first[0]);
    ui->lb_firstY->setNum(first[1]);
    ui->lb_firstZ->setNum(first[2]);

    ui->lb_lastX->setNum(last[0]);
    ui->lb_lastY->setNum(last[1]);
    ui->lb_lastZ->setNum(last[2]);
}
void MainWindow::on_pushButton_couper_clicked()
{
    //MyMesh::Point P1 = mesh.point(mesh.vertex_handle(1));
    //MyMesh::Point P2 = mesh.point(mesh.vertex_handle(3));
    //MyMesh::Point P3 = mesh.point(mesh.vertex_handle(mesh.n_vertices()));

    /*MyMesh::Point P1 = MyMesh::Point(0,0,0);
    MyMesh::Point P2 = MyMesh::Point(3,0,0);
    MyMesh::Point P3 = MyMesh::Point(0,3,0);*/

    MyMesh::Point P1 = MyMesh::Point(ui->displayWidget->center_[0],ui->displayWidget->center_[1],ui->displayWidget->center_[2]);
    MyMesh::Point P2 = MyMesh::Point(ui->displayWidget->first_point_to_cut3D[0],ui->displayWidget->first_point_to_cut3D[1],ui->displayWidget->first_point_to_cut3D[2]);
    MyMesh::Point P3 = MyMesh::Point(ui->displayWidget->last_point_to_cut3D[0],ui->displayWidget->last_point_to_cut3D[1],ui->displayWidget->last_point_to_cut3D[2]);
    //MyMesh frag1 = cut(&mesh,P1,P2,P3)[0];
    //MyMesh frag2 = cut(&mesh,P1,P2,P3)[1];
    //MyMesh cpy = mesh;
    MyMesh cpy = mesh;



    //MyMesh cpy = mesh;
    float * eq_plane = equation_plane(P1,P2,P3);


    for(MyMesh::FaceIter f = mesh.faces_begin(); f != mesh.faces_end(); f++)
    {

        QVector<MyMesh::Point> facePoints;
        for(MyMesh::FaceVertexIter fv = mesh.fv_begin(*f);fv.is_valid();fv++)
        {
            facePoints.push_back(mesh.point(*fv));
        }
        float x = (facePoints.at(0)[0]+facePoints.at(1)[0]+facePoints.at(2)[0])/3;
        float y = (facePoints.at(0)[1]+facePoints.at(1)[1]+facePoints.at(2)[1])/3;
        float z = (facePoints.at(0)[2]+facePoints.at(1)[2]+facePoints.at(2)[2])/3;

        //MyMesh::Point PM = mesh.point(mesh.vertex_handle(vem->idx()));
        if((x*eq_plane[0]+y*eq_plane[1]+z*eq_plane[2]+eq_plane[3]) < 0.0)
        {
            /*for(MyMesh::VertexFaceIter vf = mesh.vf_begin(*vem);vf.is_valid();vf ++)
            {
                mesh.delete_face(mesh.face_handle(vf->idx()),false);
            }*/

            /*for(MyMesh::FaceEdgeIter fe = mesh.fe_begin(*f);fe->is_valid();fe++)
            {

                mesh.delete_edge(*fe,false);
            }*/

            mesh.delete_face(*f,true);

            /*for(MyMesh::FaceVertexIter fv = mesh.fv_begin(*f);fv->is_valid();fv++)
            {
                mesh.delete_vertex(*fv,false);
            }*/
        }



    }
    std::vector<MyMesh::VertexHandle> vectvertices;

    mesh.garbage_collection();
    for(MyMesh::FaceIter f = cpy.faces_begin(); f != cpy.faces_end(); f++)
    {

        QVector<MyMesh::Point> facePoints;
        for(MyMesh::FaceVertexIter fv = cpy.fv_begin(*f);fv.is_valid();fv++)
        {
            facePoints.push_back(cpy.point(*fv));
        }
        float x = (facePoints.at(0)[0]+facePoints.at(1)[0]+facePoints.at(2)[0])/3;
        float y = (facePoints.at(0)[1]+facePoints.at(1)[1]+facePoints.at(2)[1])/3;
        float z = (facePoints.at(0)[2]+facePoints.at(1)[2]+facePoints.at(2)[2])/3;

        //MyMesh::Point PM = mesh.point(mesh.vertex_handle(vem->idx()));
        if((x*eq_plane[0]+y*eq_plane[1]+z*eq_plane[2]+eq_plane[3]) >= 0.0)
        {


            cpy.delete_face(*f,true);
            MyMesh::VertexHandle verth = cpy.add_vertex(MyMesh::Point(x,y,z));
            //if((x*eq_plane[0]+y*eq_plane[1]+z*eq_plane[2]+eq_plane[3]) >= -0.0025 && (x*eq_plane[0]+y*eq_plane[1]+z*eq_plane[2]+eq_plane[3]) <= 0.0025)
            if((x*eq_plane[0]+y*eq_plane[1]+z*eq_plane[2]+eq_plane[3]) <= 0.0025)
            {
                cpy.set_color(verth, MyMesh::Color(200,0,170));
                cpy.data(verth).thickness = 6;
                vectvertices.push_back(verth);
            }


        }
        /*int nb_face = 0;
        for(MyMesh::FaceFaceIter ff = cpy.ff_begin(*f);ff.is_valid();ff++)
        {
            nb_face ++;
            if(nb_face>2)
            {
                break;
            }

        }
        if(nb_face<2)
        {
            mesh.delete_face(*f,true);
        }
        cpy.garbage_collection();*/

    }

    cpy.garbage_collection();
    for(MyMesh::VertexIter vi = cpy.vertices_begin(); vi != cpy.vertices_end(); vi++)
    {
        float x = cpy.point(*vi)[0];
        float y = cpy.point(*vi)[1];
        float z = cpy.point(*vi)[2];
        if((x*eq_plane[0]+y*eq_plane[1]+z*eq_plane[2]+eq_plane[3]) >= 0.0)
        {
            cpy.delete_vertex(*vi);
        }
    }

    /*if(vectvertices.size()>=3)
    {
        cpy.add_face(vectvertices);
    }*/
        //cpy.add_face(vectvertices);
    QString fileName = "../mesh1";

            qDebug()<<fileName;
    if (true)
    {
         fileName+=".obj";
         QFile file(fileName);
         file.open(QIODevice::WriteOnly|QIODevice::Text);
         QTextStream out(&file);
         out.setCodec("UTF-8");
         //out.setVersion(QDataStream::Qt_4_5);
         out<<"g"<<endl;
         for(MyMesh::VertexIter ve = mesh.vertices_begin(); ve != mesh.vertices_end();ve++)
         {
            out<<"v ";
            out<<cpy.point(*ve)[0]<<" ";
            out<<cpy.point(*ve)[1]<<" ";
            out<<cpy.point(*ve)[2];
            out<<endl;
         }
         if (!file.open(QIODevice::WriteOnly))
         {
            //QMessageBox::information(this, tr("Unable to open file"),
            file.errorString();

         }
         QProcess process1;
         qDebug()<<qApp->applicationDirPath();
         process1.start("./../starterLight/remesh/build/scanLIDAR ../mesh1.obj");
         //process1.waitForFinished();

         //QDesktopServices::openUrl(QUrl("file:./remesh/build/scanLIDAR"));
    }
    vectvertices.clear();
    //cpy.garbage_collection();

    //showBorder(&cpy,eq_plane);
    //resetAllColorsAndThickness(&mesh);
    MyMesh result;
    QString fileNameobj = qApp->applicationDirPath() + "mesh.obj";

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(result, fileName.toUtf8().constData());

    result.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&result);

    // on affiche le maillage
    displayMesh(&cpy);
    //displayMesh(&cpy);
}
