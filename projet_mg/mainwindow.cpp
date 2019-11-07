#include "mainwindow.h"
#include "ui_mainwindow.h"

#define SMALL_NUM   0.00000001 // anything that avoids division overflow
#define dot(u,v)   (u[0] * v[0] + u[1] * v[1] + u[2] * v[2])


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

void MainWindow::displayMeshes(QVector<MyMesh *> vector_meshes, bool isTemperatureMap, float mapRange)
{
    unsigned int nb_faces_,nb_vertices_,nb_edges_ = 0;
    for (int j = 0;j<vector_meshes.size();j++)
    {
        nb_faces_+=vector_meshes.at(j)->n_faces();
        nb_edges_+=vector_meshes.at(j)->n_edges();
        nb_vertices_+=vector_meshes.at(j)->n_vertices();
    }
    GLuint* triIndiceArray = new GLuint[nb_faces_ * 3];
    GLfloat* triCols = new GLfloat[nb_faces_ * 3 * 3];
    GLfloat* triVerts = new GLfloat[nb_faces_ * 3 * 3];
    int i = 0;
    for (int j = 0;j<vector_meshes.size();j++)
    {

        if(isTemperatureMap)
        {
            QVector<float> values;

            if(mapRange == -1)
            {
                for (MyMesh::VertexIter curVert = vector_meshes.at(j)->vertices_begin(); curVert != vector_meshes.at(j)->vertices_end(); curVert++)
                    values.append(fabs(vector_meshes.at(j)->data(*curVert).value));
                qSort(values);
                mapRange = values.at(values.size()*0.8);
                qDebug() << "mapRange" << mapRange;
            }

            float range = mapRange;
            MyMesh::ConstFaceIter fIt(vector_meshes.at(j)->faces_begin()), fEnd(vector_meshes.at(j)->faces_end());
            MyMesh::ConstFaceVertexIter fvIt;

            for (; fIt!=fEnd; ++fIt)
            {
                fvIt = vector_meshes.at(j)->cfv_iter(*fIt);
                if(vector_meshes.at(j)->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0);}
                else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0);}
                triVerts[3*i+0] = vector_meshes.at(j)->point(*fvIt)[0]; triVerts[3*i+1] = vector_meshes.at(j)->point(*fvIt)[1]; triVerts[3*i+2] = vector_meshes.at(j)->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++; ++fvIt;
                if(vector_meshes.at(j)->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0);}
                else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0);}
                triVerts[3*i+0] = vector_meshes.at(j)->point(*fvIt)[0]; triVerts[3*i+1] = vector_meshes.at(j)->point(*fvIt)[1]; triVerts[3*i+2] = vector_meshes.at(j)->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++; ++fvIt;
                if(vector_meshes.at(j)->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0);}
                else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-vector_meshes.at(j)->data(*fvIt).value/range) * 255.0, 255.0);}
                triVerts[3*i+0] = vector_meshes.at(j)->point(*fvIt)[0]; triVerts[3*i+1] = vector_meshes.at(j)->point(*fvIt)[1]; triVerts[3*i+2] = vector_meshes.at(j)->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++;
            }
        }
        else
        {
            MyMesh::ConstFaceIter fIt(vector_meshes.at(j)->faces_begin()), fEnd(vector_meshes.at(j)->faces_end());
            MyMesh::ConstFaceVertexIter fvIt;
            for (; fIt!=fEnd; ++fIt)
            {
                fvIt = vector_meshes.at(j)->cfv_iter(*fIt);
                triCols[3*i+0] = vector_meshes.at(j)->color(*fIt)[0]; triCols[3*i+1] = vector_meshes.at(j)->color(*fIt)[1]; triCols[3*i+2] = vector_meshes.at(j)->color(*fIt)[2];
                triVerts[3*i+0] = vector_meshes.at(j)->point(*fvIt)[0]; triVerts[3*i+1] = vector_meshes.at(j)->point(*fvIt)[1]; triVerts[3*i+2] = vector_meshes.at(j)->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++; ++fvIt;
                triCols[3*i+0] = vector_meshes.at(j)->color(*fIt)[0]; triCols[3*i+1] = vector_meshes.at(j)->color(*fIt)[1]; triCols[3*i+2] = vector_meshes.at(j)->color(*fIt)[2];
                triVerts[3*i+0] = vector_meshes.at(j)->point(*fvIt)[0]; triVerts[3*i+1] = vector_meshes.at(j)->point(*fvIt)[1]; triVerts[3*i+2] = vector_meshes.at(j)->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++; ++fvIt;
                triCols[3*i+0] = vector_meshes.at(j)->color(*fIt)[0]; triCols[3*i+1] = vector_meshes.at(j)->color(*fIt)[1]; triCols[3*i+2] = vector_meshes.at(j)->color(*fIt)[2];
                triVerts[3*i+0] = vector_meshes.at(j)->point(*fvIt)[0]; triVerts[3*i+1] = vector_meshes.at(j)->point(*fvIt)[1]; triVerts[3*i+2] = vector_meshes.at(j)->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++;
            }
        }
     }


    ui->displayWidget->loadMesh(triVerts, triCols, nb_faces_ * 3 * 3, triIndiceArray, nb_faces_ * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;
    GLuint* linesIndiceArray = new GLuint[nb_edges_ * 2];
    GLfloat* linesCols = new GLfloat[nb_edges_ * 2 * 3];
    GLfloat* linesVerts = new GLfloat[nb_edges_ * 2 * 3];

    i = 0;

    QList<QPair<float, int> > edgeSizes;
    for (int j = 0;j<vector_meshes.size();j++)
    {

        QHash<float, QList<int> > edgesIDbyThickness;
        for (MyMesh::EdgeIter eit = vector_meshes.at(j)->edges_begin(); eit != vector_meshes.at(j)->edges_end(); ++eit)
        {
            float t = vector_meshes.at(j)->data(*eit).thickness;
            if(t > 0)
            {
                if(!edgesIDbyThickness.contains(t))
                    edgesIDbyThickness[t] = QList<int>();
                edgesIDbyThickness[t].append((*eit).idx());
            }
        }
        QHashIterator<float, QList<int> > it(edgesIDbyThickness); //?????

        while (it.hasNext())
        {
            it.next();

            for(int e = 0; e < it.value().size(); e++)
            {
                int eidx = it.value().at(e);

                MyMesh::VertexHandle vh1 = vector_meshes.at(j)->to_vertex_handle(vector_meshes.at(j)->halfedge_handle(vector_meshes.at(j)->edge_handle(eidx), 0));
                linesVerts[3*i+0] = vector_meshes.at(j)->point(vh1)[0];
                linesVerts[3*i+1] = vector_meshes.at(j)->point(vh1)[1];
                linesVerts[3*i+2] = vector_meshes.at(j)->point(vh1)[2];
                linesCols[3*i+0] = vector_meshes.at(j)->color(vector_meshes.at(j)->edge_handle(eidx))[0];
                linesCols[3*i+1] = vector_meshes.at(j)->color(vector_meshes.at(j)->edge_handle(eidx))[1];
                linesCols[3*i+2] = vector_meshes.at(j)->color(vector_meshes.at(j)->edge_handle(eidx))[2];
                linesIndiceArray[i] = i;
                i++;

                MyMesh::VertexHandle vh2 = vector_meshes.at(j)->from_vertex_handle(vector_meshes.at(j)->halfedge_handle(vector_meshes.at(j)->edge_handle(eidx), 0));
                linesVerts[3*i+0] = vector_meshes.at(j)->point(vh2)[0];
                linesVerts[3*i+1] = vector_meshes.at(j)->point(vh2)[1];
                linesVerts[3*i+2] = vector_meshes.at(j)->point(vh2)[2];
                linesCols[3*i+0] = vector_meshes.at(j)->color(vector_meshes.at(j)->edge_handle(eidx))[0];
                linesCols[3*i+1] = vector_meshes.at(j)->color(vector_meshes.at(j)->edge_handle(eidx))[1];
                linesCols[3*i+2] = vector_meshes.at(j)->color(vector_meshes.at(j)->edge_handle(eidx))[2];
                linesIndiceArray[i] = i;
                i++;
            }
            edgeSizes.append(qMakePair(it.key(), it.value().size()));
        }
    }
    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[nb_vertices_];
    GLfloat* pointsCols = new GLfloat[nb_vertices_ * 3];
    GLfloat* pointsVerts = new GLfloat[nb_vertices_ * 3];
    i = 0;
    QList<QPair<float, int> > vertsSizes;
    for (int j = 0;j<vector_meshes.size();j++)
    {

        QHash<float, QList<int> > vertsIDbyThickness;
        for (MyMesh::VertexIter vit = vector_meshes.at(j)->vertices_begin(); vit != vector_meshes.at(j)->vertices_end(); ++vit)
        {
            float t = vector_meshes.at(j)->data(*vit).thickness;
            if(t > 0)
            {
                if(!vertsIDbyThickness.contains(t))
                    vertsIDbyThickness[t] = QList<int>();
                vertsIDbyThickness[t].append((*vit).idx());
            }
        }
        QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);


        while (vitt.hasNext())
        {
            vitt.next();

            for(int v = 0; v < vitt.value().size(); v++)
            {
                int vidx = vitt.value().at(v);

                pointsVerts[3*i+0] = vector_meshes.at(j)->point(vector_meshes.at(j)->vertex_handle(vidx))[0];
                pointsVerts[3*i+1] = vector_meshes.at(j)->point(vector_meshes.at(j)->vertex_handle(vidx))[1];
                pointsVerts[3*i+2] = vector_meshes.at(j)->point(vector_meshes.at(j)->vertex_handle(vidx))[2];
                pointsCols[3*i+0] = vector_meshes.at(j)->color(vector_meshes.at(j)->vertex_handle(vidx))[0];
                pointsCols[3*i+1] = vector_meshes.at(j)->color(vector_meshes.at(j)->vertex_handle(vidx))[1];
                pointsCols[3*i+2] = vector_meshes.at(j)->color(vector_meshes.at(j)->vertex_handle(vidx))[2];
                pointsIndiceArray[i] = i;
                i++;
            }
            vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
        }
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

    /*for(int i=0;i < _mesh->n_vertices();i++)
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
    }*/
    for(MyMesh::FaceIter f = _mesh->faces_begin(); f != _mesh->faces_end(); f++)
    {

        std::vector<MyMesh::Point> facePoints;

        for(MyMesh::FaceVertexIter fv = _mesh->fv_begin(*f);fv.is_valid();fv++)
        {
            facePoints.push_back(_mesh->point(*fv));
            //vertices.push_back(*fv);
        }
        float x = (facePoints.at(0)[0]+facePoints.at(1)[0]+facePoints.at(2)[0])/3;
        float y = (facePoints.at(0)[1]+facePoints.at(1)[1]+facePoints.at(2)[1])/3;
        float z = (facePoints.at(0)[2]+facePoints.at(1)[2]+facePoints.at(2)[2])/3;

        //MyMesh::Point PM = mesh.point(mesh.vertex_handle(vem->idx()));
        if((x*eq_plane[0]+y*eq_plane[1]+z*eq_plane[2]+eq_plane[3]) >= 0.0)
        {

            //cpy.delete_face(*f,true);
            std::vector<VertexHandle> vertices;
            for(unsigned int i = 0;i<facePoints.size();i++)
            {
                vertices.push_back(frag1.add_vertex(facePoints.at(i)));
            }
            MyMesh::FaceHandle fh = frag1.add_face(vertices);
            //MyMesh::VertexHandle verth = cpy.add_vertex(MyMesh::Point(x,y,z));
            frag1.set_color(fh,MyMesh::Color(0,0,255));
            vertices.clear();
        }
        else
        {
            std::vector<VertexHandle> vertices;
            for(unsigned int i = 0;i<facePoints.size();i++)
            {
                vertices.push_back(frag2.add_vertex(facePoints.at(i)));
            }
            MyMesh::FaceHandle fh = frag2.add_face(vertices);
            //MyMesh::VertexHandle verth = cpy.add_vertex(MyMesh::Point(x,y,z));
            frag2.set_color(fh,MyMesh::Color(255,0,0));
            vertices.clear();
        }

    }
    for (MyMesh::EdgeIter ei = frag1.edges_begin();ei != frag1.edges_end();ei++)
    {
        frag1.set_color(*ei,MyMesh::Color(150,150,0));
        frag1.data(*ei).thickness = 1;
    }
    for (MyMesh::EdgeIter ei = frag2.edges_begin();ei != frag2.edges_end();ei++)
    {
        frag2.set_color(*ei,MyMesh::Color(0,150,150));
        frag2.data(*ei).thickness = 1;
    }
    MyMesh * tabmesh = new MyMesh[2];
    tabmesh[0] = frag1;
    tabmesh[1] = frag2;
    resetAllColorsAndThickness(&frag1);
    resetAllColorsAndThickness(&frag2);
    return tabmesh;
}
void MainWindow::on_pushButton_couper_clicked()
{
    //MyMesh::Point P1 = mesh.point(mesh.vertex_handle(1));
    //MyMesh::Point P2 = mesh.point(mesh.vertex_handle(3));
    //MyMesh::Point P3 = mesh.point(mesh.vertex_handle(mesh.n_vertices()));

    /*MyMesh::Point P1 = MyMesh::Point(0,0,0);
    MyMesh::Point P2 = MyMesh::Point(3,0,0);
    MyMesh::Point P3 = MyMesh::Point(0,3,0);*/

    //MyMesh::Point P1 = MyMesh::Point(ui->displayWidget->center_[0],ui->displayWidget->center_[1],ui->displayWidget->center_[2]);
    MyMesh::Point P1 = MyMesh::Point(ui->displayWidget->center_cut3D[0],ui->displayWidget->center_cut3D[1],ui->displayWidget->center_cut3D[2]);
    MyMesh::Point P2 = MyMesh::Point(ui->displayWidget->first_point_to_cut3D[0],ui->displayWidget->first_point_to_cut3D[1],ui->displayWidget->first_point_to_cut3D[2]);
    MyMesh::Point P3 = MyMesh::Point(ui->displayWidget->last_point_to_cut3D[0],ui->displayWidget->last_point_to_cut3D[1],ui->displayWidget->last_point_to_cut3D[2]);
    MyMesh frag1 = cut(&mesh,P1,P2,P3)[0];
    MyMesh frag2 = cut(&mesh,P1,P2,P3)[1];
    //MyMesh cpy = mesh;
    //________________________________________________________________________________________________________



    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    //resetAllColorsAndThickness(&frag1);

    // on affiche le maillage
    QVector<MyMesh * > listmeshes;
    listmeshes.push_back(&frag1);
    listmeshes.push_back(&frag2);
    displayMeshes(listmeshes);
    //displayMesh(&frag1);
    //displayMesh(&cpy);
}

void MainWindow::boite_englobante(){
    for (MyMesh::VertexIter curVert = mesh.vertices_begin(); curVert != mesh.vertices_end(); curVert++)
    {
        if(mesh.point(curVert)[0] < be_min_x) be_min_x = mesh.point(curVert)[0];
        if(mesh.point(curVert)[0] > be_max_x) be_max_x = mesh.point(curVert)[0];
        if(mesh.point(curVert)[1] < be_min_y) be_min_y = mesh.point(curVert)[1];
        if(mesh.point(curVert)[1] > be_max_y) be_max_y = mesh.point(curVert)[1];
        if(mesh.point(curVert)[2] < be_min_z) be_min_z = mesh.point(curVert)[2];
        if(mesh.point(curVert)[2] > be_max_z) be_max_z = mesh.point(curVert)[2];
    }
    qDebug() << "Boîte englobante :";
    qDebug() << "    min x : " << be_min_x << "; min y : " << be_min_y << "; min z : " << be_min_z;
    qDebug() << "    max x : " << be_max_x << "; max y : " << be_max_y << "; max z : " << be_max_z;
}

int MainWindow::intersect(FaceHandle triangle, MyMesh::Point point_extern, MyMesh::Point point){
    //aidé par http://geomalgorithms.com/a06-_intersect-2.html
    MyMesh::Normal n;
    MyMesh::Point u, v;
    //Vector    u, v, n;              // triangle vectors : normal + 2 edge
    MyMesh::Point dir, w0, w;           // ray vectors :
    float     r, a, b;              // params to calc ray-plane intersect

    MyMesh::Point tv0;
    MyMesh::Point tv1;
    MyMesh::Point tv2;
    // get triangle edge vectors and plane normal

    MyMesh::FaceVertexIter curVertex = mesh.fv_iter(triangle);
    if (curVertex.is_valid()) tv2 = mesh.point(*curVertex); //verifier que c'est bien une edge
    else qDebug() << "-----------------------------curVertex 1 not valid";
    curVertex ++;
    if (curVertex.is_valid()) tv1 = mesh.point(*curVertex);
    else qDebug() << "-----------------------------curEdge 2 not valid";
    curVertex ++;
    if (curVertex.is_valid()) tv0 = mesh.point(*curVertex);
    else qDebug() << "-----------------------------curEdge 3 not valid";

   /* VectorT<float, 3> u;
    VectorT<float, 3> v;
    */
    //u = T.V1 - T.V0;
    //v = T.V2 - T.V0;
    u = tv1 - tv0;
    v = tv2 - tv0;

    //n = u * v;              // cross product

    n = mesh.calc_face_normal(triangle);
    if (n == MyMesh::Normal(0,0,0))             // triangle is degenerate
        return -1;                  // do not deal with this case
    /*VectorT<float, 3> n;
    n[0] = n_tmp[0];
    n[1] = n_tmp[1];
    n[2] = n_tmp[2];*/

    //dir = R.P1 - R.P0;              // ray direction vector
    dir = point_extern - point;

    //w0 = R.P0 - T.V0;
    w0 = point - tv0;

    a = -dot(n,w0);
//    qDebug() << "a" << a;
//    float ab = n[0] * w0[0] + n[1] * w0[1] + n[2] * w0[2];
//    qDebug() << "ab" << ab;
    b = dot(n,dir);
//    qDebug() << "b" << b;
//    float bb = n[0] * dir[0] + n[1] * dir[1] + n[2] * dir[2];
//    qDebug() << "bb" << bb;
    if (fabs(b) < SMALL_NUM) {     // ray is  parallel to triangle plane
        if (a == 0.0f)                 // ray lies in triangle plane
            return 2;
        else return 0;              // ray disjoint from plane
    }

    // get intersect point of ray with triangle plane
    r = a / b;
    if (r < 0.0f && r > 1.0f)                    // ray goes away from triangle
        return 0;                   // => no intersect
    // for a segment, also test if (r > 1.0) => no intersect

    //MyMesh::Point I = R.P0 + r * dir;            // intersect point of ray and plane
    MyMesh::Point I = point + r * dir;

    // is I inside T?
    float    uu, uv, vv, wu, wv, D;
    uu = dot(u,u);
    uv = dot(u,v);
    vv = dot(v,v);
    w = I - point;
    wu = dot(w,u);
    wv = dot(w,v);
    D = uv * uv - uu * vv;

    // get and test parametric coords
    float s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0f || s > 1.0f)         // I is outside T
        return 0;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0f || (s + t) > 1.0f)  // I is outside T
        return 0;

    return 1;                       // I is in T
}
bool MainWindow::pointIsInMesh(MyMesh::Point point){

    MyMesh::Point point_extern = MyMesh::Point(be_max_x+42, be_max_y+42, be_max_z+42);
    int intersections = 0;
    for (MyMesh::FaceIter curFace = mesh.faces_begin(); curFace != mesh.faces_end(); curFace++)
    {
        FaceHandle fh = *curFace;
        if(intersect(fh, point_extern, point) == 1)
            intersections++;
    }

    if(intersections % 2 == 1) {
        return true;
    } else {
        return false;
    }
    
}

float randomFloat(float min, float max)
{
    float random = ((float) rand()) / (float) RAND_MAX;
    float range = max - min;
    return (random*range) + min;
}

bool MainWindow::createSeeds(){
    boite_englobante();
    for (unsigned i = 0; i < nb_seeds;) {
        //crée un point dans la boite englobante
        qDebug() <<"x :" << be_min_x << be_max_x;
        qDebug() <<"y :" << be_min_y << be_max_y;
        qDebug() <<"z :" << be_min_z << be_max_z;
        //float tmp_x = rand()% ((int)be_max_x+1) + be_min_x;
        //float tmp_y = rand()% ((int)be_max_y+1) + be_min_y;
        //float tmp_z = rand()% ((int)be_max_z+1) + be_min_z;
        float tmp_x = randomFloat(be_min_x, be_max_x);
        float tmp_y = randomFloat(be_min_y, be_max_y);
        float tmp_z = randomFloat(be_min_z, be_max_z);
        MyMesh::Point point = MyMesh::Point(tmp_x, tmp_y, tmp_z);
        if(pointIsInMesh(point)) {
            seeds_impact.push_back(point);
            VertexHandle v = mesh.add_vertex(point);
            mesh.data(v).thickness = 4;
            mesh.set_color(v, MyMesh::Color(255, 0, 0));
            i++;
            qDebug() <<"point " << v.idx() << " is in :" << point[0] << point[1] << point[2];
        }
        else qDebug() << " is out :" << point[0] << point[1] << point[2];
    }
    displayMesh(&mesh);
}
void MainWindow::on_pushButton_clicked()
{
    createSeeds();
}

void MainWindow::on_showSelections_clicked()
{
    resetAllColorsAndThickness(&mesh);

    /* **** à compléter ! ****
     * cette fonction utilise les vatiables de sélection vertexSelection, edgeSelection et faceSelection
     * qui sont les ID des élements sélectionnés et qui sont égales à -1 si la sélection est vide
     */

    if(vertexSelection >= 0 && vertexSelection < mesh.n_vertices()){
        VertexHandle vh = mesh.vertex_handle(vertexSelection);
        mesh.set_color(vh, MyMesh::Color(255,0,0));
        mesh.data(vh).thickness = 8;
    }

    displayMesh(&mesh);
}

void MainWindow::on_vertexSelection_valueChanged(int arg1)
{
    vertexSelection = arg1;
}
