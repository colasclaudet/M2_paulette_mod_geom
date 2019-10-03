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
}

MainWindow::~MainWindow()
{
    delete ui;
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

    MyMesh frag1 = cut(&mesh,P1,P2,P3)[0];
    MyMesh frag2 = cut(&mesh,P1,P2,P3)[1];
    MyMesh cpy = mesh;
    mesh.clean();
    for(MyMesh::VertexIter ve = frag1.vertices_begin(); ve != frag1.vertices_end(); ve++)
    {
        MyMesh::VertexHandle vh = mesh.add_vertex(frag1.point(frag1.vertex_handle(ve->idx())));
        mesh.set_color(mesh.vertex_handle(vh.idx()), MyMesh::Color(0, 0, 255));
        mesh.data(mesh.vertex_handle(vh.idx())).thickness = 2;
    }
    for(MyMesh::VertexIter ve = frag2.vertices_begin(); ve != frag2.vertices_end(); ve++)
    {
        MyMesh::VertexHandle vh = mesh.add_vertex(frag2.point(frag2.vertex_handle(ve->idx())));
        mesh.set_color(mesh.vertex_handle(vh.idx()), MyMesh::Color(255, 0, 0));
        mesh.data(mesh.vertex_handle(vh.idx())).thickness = 2;
    }
    for(MyMesh::FaceIter f = cpy.faces_begin();f!=cpy.faces_end();f++)
    {

        if(true)
        {
        std::vector<MyMesh::VertexHandle> uneNouvelleFace;
        std::vector<MyMesh::Point> face_points;
        for(MyMesh::FaceVertexIter fv = cpy.fv_begin(*f); fv.is_valid(); fv++)
        {
            //uneNouvelleFace.push_back(cpy.vertex_handle(fv->idx()));
            face_points.push_back(cpy.point(cpy.vertex_handle(fv->idx())));
        }
        //qDebug()<<face_points.size();

        uneNouvelleFace.push_back(cpy.vertex_handle(0));
        uneNouvelleFace.push_back(cpy.vertex_handle(0));
        uneNouvelleFace.push_back(cpy.vertex_handle(0));
        for (MyMesh::VertexIter v = mesh.vertices_begin(); v!=mesh.vertices_end();v++)
        {
            MyMesh::Point P =  mesh.point(mesh.vertex_handle(v->idx()));
            if(face_points.at(0) == P)
            {
                //uneNouvelleFace.push_back(mesh.vertex_handle(v->idx()));
                uneNouvelleFace.at(0) = mesh.vertex_handle(v->idx());
            }//verifier point de cpy est dans mesh
            if(face_points.at(1) == P)
            {
                //uneNouvelleFace.push_back(mesh.vertex_handle(v->idx()));
                uneNouvelleFace.at(1) = mesh.vertex_handle(v->idx());
            }//verifier point de cpy est dans mesh
            if(face_points.at(2) == P)
            {
                //uneNouvelleFace.push_back(mesh.vertex_handle(v->idx()));
                uneNouvelleFace.at(2) = mesh.vertex_handle(v->idx());
            }//verifier point de cpy est dans mesh
        }
        //qDebug()<<uneNouvelleFace.size();
        //mesh.add_face(uneNouvelleFace);
        if(mesh.is_valid_handle(mesh.add_face(uneNouvelleFace))) //123
        {

        }
        else
        {
            MyMesh::VertexHandle vh = uneNouvelleFace.at(0);
            uneNouvelleFace.at(0) = uneNouvelleFace.at(2);
            uneNouvelleFace.at(2) = vh;
            //mesh.add_face(uneNouvelleFace);
            if(mesh.is_valid_handle(mesh.add_face(uneNouvelleFace)))//321
            {

            }
            else
            {
                vh = uneNouvelleFace.at(0);
                uneNouvelleFace.at(0) = uneNouvelleFace.at(1);
                uneNouvelleFace.at(1) = vh;
                //mesh.add_face(uneNouvelleFace);
                if(mesh.is_valid_handle(mesh.add_face(uneNouvelleFace)))//231
                {

                }
                else
                {
                    vh = uneNouvelleFace.at(1);
                    uneNouvelleFace.at(1) = uneNouvelleFace.at(2);
                    uneNouvelleFace.at(2) = vh;
                    if(mesh.is_valid_handle(mesh.add_face(uneNouvelleFace)))//213
                    {

                    }
                    else
                    {
                        vh = uneNouvelleFace.at(0);
                        uneNouvelleFace.at(0) = uneNouvelleFace.at(2);
                        uneNouvelleFace.at(2) = vh;
                        if(mesh.is_valid_handle(mesh.add_face(uneNouvelleFace)))//312
                        {

                        }
                        else
                        {
                            vh = uneNouvelleFace.at(1);
                            uneNouvelleFace.at(1) = uneNouvelleFace.at(0);
                            uneNouvelleFace.at(0) = vh;

                            if(mesh.is_valid_handle(mesh.add_face(uneNouvelleFace)))//132
                            {

                            }
                            else
                            {
                                vh = uneNouvelleFace.at(2);
                                uneNouvelleFace.at(2) = uneNouvelleFace.at(1);
                                uneNouvelleFace.at(1) = vh;
                                mesh.add_face(uneNouvelleFace);
                            }


                        }

                    }
                }
            }
        }

        face_points.clear();
        uneNouvelleFace.clear();
        }
    }

    /*for(MyMesh::VertexIter v = frag1.vertices_begin(); v!=frag1.vertices_end();v++)
    {
        MyMesh::Point pf1 = frag1.point(frag1.vertex_handle(v->idx()));
        for(MyMesh::VertexIter v2 = frag2.vertices_begin(); v2!=frag2.vertices_end();v2++)
        {
            MyMesh::Point pf2 = frag2.point(frag2.vertex_handle(v2->idx()));
            if(pf1 == pf2)
            {
                for (MyMesh::VertexIter vm = cpy.vertices_begin(); vm!=cpy.vertices_end();vm++)
                {
                    MyMesh::Point pfm = cpy.point(cpy.vertex_handle(vm->idx()));
                    if(pfm == pf1)
                    {
                        cpy.delete_vertex(cpy.vertex_handle(vm->idx()));
                    }
                }
            }
        }
    }*/

    /*std::vector<MyMesh::VertexHandle> uneNouvelleFace;
    int i=1;
    for(MyMesh::VertexIter ve = frag1.vertices_begin(); ve != frag1.vertices_end(); ve++)
    {
        uneNouvelleFace.push_back(mesh.vertex_handle(ve->idx()));
        i++;
        if(i%4 == 0)
        {
            i=1;
            qDebug()<<uneNouvelleFace.size();
            mesh.add_face(uneNouvelleFace);
            uneNouvelleFace.clear();
        }
    }*/


    displayMesh(&mesh);
}
