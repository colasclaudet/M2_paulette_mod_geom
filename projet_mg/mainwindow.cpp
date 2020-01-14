#include "mainwindow.h"
#include "ui_mainwindow.h"

#define SMALL_NUM   0.00000001 // anything that avoids division overflow
#define dot(u,v)   (u[0] * v[0] + u[1] * v[1] + u[2] * v[2])

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->displayWidget, SIGNAL(sig_changePoint(OpenMesh::Vec3f, OpenMesh::Vec3f)), this, SLOT(changeValuePointCut(OpenMesh::Vec3f, OpenMesh::Vec3f)));
    //menu
    connect(ui->actionCharger_OBJ, SIGNAL(triggered()), this, SLOT(importFile()));
    connect(ui->actionGenerer_cube, SIGNAL(triggered()), this, SLOT(generer_cube()));
    list_meshes_to_display.push_back(&frag1);
    list_meshes_to_display.push_back(&frag2);

    list_meshes_to_save.push_back(frag1);
    list_meshes_to_save.push_back(frag2);
}

MainWindow::~MainWindow()
{
    delete ui;
}
/*
 *
 * A SUPPRIMER - pb avec le fichier .moc lorsqu'on supprime ces fonctions (dû à l'ui)
 *
 */
void MainWindow::on_pushButton_chargement_clicked(){

}
void MainWindow::on_pushButton_generer_clicked(){

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

/*
 *
 * FIN SUPPRIMER
 *
 */
/* **** début de la partie boutons et IHM **** */

/* * MENU  * */
// exemple pour charger un fichier .obj
void MainWindow::importFile(){
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
void MainWindow::generer_cube()
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

/* * FIN MENU  * */
//modifie les valeurs l'affichage des points pour le découpage en 2
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
    //MyMesh::Point P1 = MyMesh::Point(ui->displayWidget->center_[0],ui->displayWidget->center_[1],ui->displayWidget->center_[2]);
    MyMesh::Point P1 = MyMesh::Point(ui->displayWidget->center_cut3D[0],ui->displayWidget->center_cut3D[1],ui->displayWidget->center_cut3D[2]);
    MyMesh::Point P2 = MyMesh::Point(ui->displayWidget->first_point_to_cut3D[0],ui->displayWidget->first_point_to_cut3D[1],ui->displayWidget->first_point_to_cut3D[2]);
    MyMesh::Point P3 = MyMesh::Point(ui->displayWidget->last_point_to_cut3D[0],ui->displayWidget->last_point_to_cut3D[1],ui->displayWidget->last_point_to_cut3D[2]);

    cutV2(&mesh,P1,P2,P3);
    chooseFragToSave();
    displayMeshes(list_meshes_to_display);
}

//création des seeds
void MainWindow::on_pushButton_clicked()
{
    createSeeds();
}

//modifie la liste des fichiers affichés
void MainWindow::on_checkBox_display_frag1_clicked()
{
    chooseFragToDisplay();
}
void MainWindow::on_checkBox_display_frag2_clicked()
{
    chooseFragToDisplay();
}

//modifie la liste des fichiers à sauvegarder
void MainWindow::on_checkBox_save_frag1_clicked()
{
    chooseFragToSave();
}
void MainWindow::on_checkBox_save_frag2_clicked()
{
    chooseFragToSave();
}
void MainWindow::on_pushButton_save_clicked()
{
    if(list_meshes_to_save.size() == 0){
        ui->statusBar->showMessage("Aucun élément à sauvegarder.");
        return;
    }
    ui->statusBar->showMessage("Sauvegarde en cours...");
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save file to obj"), "",tr("File"));
    //qDebug()<<fileName;
    if (true)
    {
        QDir dir(fileName);
        if (!dir.exists())
            dir.mkpath(".");

         //QFile file(fileName);
         //file.open(QIODevice::WriteOnly|QIODevice::Text);
        QString cp = fileName;
        for (unsigned int i=0;i<this->list_meshes_to_save.size();i++)
        {

            std::stringstream st;
            st<< i;

            std::string s = st.str();
            QString qs = QString::fromStdString(s);
            fileName = cp + "/" +"frag" + qs +".obj";
            //qDebug()<<fileName;
            OpenMesh::IO::write_mesh(this->list_meshes_to_save.at(i), fileName.toStdString());
            //qDebug()<<"save : frag"<<i;
        }
        this->list_meshes_to_save.clear();
    }
    ui->statusBar->showMessage("Sauvegarde terminée.");
}

/* **** fin de la partie boutons et IHM **** */

/* **** fonctions supplémentaires **** */
/* ***Affichage et objets*** */
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

//affichage d'un maillage (ou d'un vecteur de mesh)
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
            //qDebug() << "mapRange" << mapRange;
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

    if(vector_meshes.size() == 0){
        qDebug() << "aucun affichage a faire";
        //TODO enlever tous les affichages de la fenetre
        return;
    }
    ui->statusBar->showMessage("Affichage en cours...");
    unsigned int nb_faces_ =0,nb_vertices_=0,nb_edges_ = 0;
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
                //qDebug() << "mapRange" << mapRange;
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

    ui->statusBar->showMessage("Affichage réalisé.");
}

//<affiche bordure> met plusieurs seconde à s'executer (n'est plus utilisé)
void MainWindow::showBorder(MyMesh* _mesh/*, float * eq_plane*/)
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
            //MyMesh::HalfedgeHandle heh1 = _mesh->halfedge_handle(e,0);
            _mesh->halfedge_handle(e,0);
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

void MainWindow::chooseFragToDisplay(){
    list_meshes_to_display.clear();
    if(ui->checkBox_display_frag1->checkState()){
        list_meshes_to_display.push_back(&frag1);
    }
    if(ui->checkBox_display_frag2->checkState()){
        list_meshes_to_display.push_back(&frag2);
    }
    list_meshes_to_display.push_back(&frag3);
    displayMeshes(list_meshes_to_display);
}
void MainWindow::chooseFragToSave(){
    list_meshes_to_save.clear();
    if(ui->checkBox_save_frag1->checkState()){
        //qDebug() << "frag1 choose";
        list_meshes_to_save.push_back(frag1);
    }
    if(ui->checkBox_save_frag2->checkState()){
        //qDebug() << "frag2 choose";
        list_meshes_to_save.push_back(frag2);
    }
}

/* ***Calcul et manipulation*** */
//Calcul l'équation d'une droite à partir de deux équations de plans (n'est plus utilisé)
double * MainWindow::equation_droite(double * eq1, double * eq2){
    double * eq_droite = new double[4];
    eq_droite[0] = eq2[0] - eq1[0];
    eq_droite[1] = eq2[1] - eq1[1];
    eq_droite[2] = eq2[2] - eq1[2];
    eq_droite[3] = eq2[3] - eq1[3];
    qDebug() << "-------------------------------------------";
    qDebug() << "equation of eq1 is " << eq1[0] << " x + " << eq1[1]
        << " y + " << eq1[2] << " z + " << eq1[3] << " = 0.";
    qDebug() << "equation of eq2 is " << eq2[0] << " x + " << eq2[1]
        << " y + " << eq2[2] << " z + " << eq2[3] << " = 0.";

    qDebug() << "equation of droite is " << eq_droite[0] << " x + " << eq_droite[1]
        << " y + " << eq_droite[2] << " z + " << eq_droite[3] << " = 0.";
    qDebug() << "-------------------------------------------";
    return eq_droite;

}

//calcul l'équation d'un plan à partir de 3 points
double * MainWindow::equation_plane(MyMesh::Point P1, MyMesh::Point P2, MyMesh::Point P3)
{
    double * eq_plane = new double[4];
    /*float a1 = P2[0] - P1[0];
    float b1 = P2[1] - P1[1];
    float c1 = P2[2] - P1[2];
    float a2 = P3[0] - P1[0];
    float b2 = P3[1] - P1[1];
    float c2 = P3[2] - P1[2];
    float a = b1 * c2 - b2 * c1;
    float b = a2 * c1 - a1 * c2;
    float c = a1 * b2 - b1 * a2;
    float d = (- a * P1[0] - b * P1[1] - c * P1[2]);*/
    double a = static_cast<double>((P2[1]-P1[1])*(P3[2] - P1[2]) - (P2[2]-P1[2])*(P3[1]-P1[2]));
    double b = static_cast<double>(-((P2[0]-P1[0])*(P3[2]-P1[2])-(P2[2]-P1[2])*(P3[0]-P1[0])));
    double c = static_cast<double>((P2[0]-P1[0])*(P3[1]-P1[1]) - (P2[1]-P1[1])*(P3[0]-P1[0]));
    double d = static_cast<double>(-( a * P1[0] + b * P1[1] + c * P1[2]));
    //qDebug() << "-( a * P1[0] + b * P1[1] + c * P1[2]) :" << a << "*" << P1[0] << "+" << b << "*" << P1[1] << "+" << c << "*" << P1[2];
    //qDebug() << "equation of plane is " << a << " x + " << b
        //<< " y + " << c << " z + " << d << " = 0.";
    eq_plane[0]=a;
    eq_plane[1]=b;
    eq_plane[2]=c;
    eq_plane[3]=d;
    return eq_plane;
}

//Résolution d'un système de 3 équations à 3 inconnues
//Réutilisation d'un code du TER de M1 (Fonction ThreePlanIntersection) et modification
MyMesh::Point MainWindow::resol_3eq_3inc(double * eq1, double * eq2, double * eq3){

    double matrice[3][4];
    double coefficient,x,y,z;
    int i=0;
    //qDebug() << "Resolution d'un systeme de 3 equations a trois inconnues\n";
    //qDebug() << "Premiere equation, entrez en ordre respectif x,y,z et la constante\n";
    //cin >> matrice[0][0]>>matrice[0][1]>>matrice[0][2]>>matrice[0][3];
    matrice[0][0] = eq1[0];
    matrice[0][1] = eq1[1];
    matrice[0][2] = eq1[2];
    matrice[0][3] = eq1[3];
    //cout << "Seconde equation, entrez en ordre respectif x,y,z et la constante\n";
    //cin >> matrice[1][0]>>matrice[1][1]>>matrice[1][2]>>matrice[1][3];
    matrice[1][0] = eq2[0];
    matrice[1][1] = eq2[1];
    matrice[1][2] = eq2[2];
    matrice[1][3] = eq2[3];
    //cout << "Premiere equation, entrez en ordre respectif x,y,z et la constante\n";
    //cin >> matrice[2][0]>>matrice[2][1]>>matrice[2][2]>>matrice[2][3];
    matrice[2][0] = eq3[0];
    matrice[2][1] = eq3[1];
    matrice[2][2] = eq3[2];
    matrice[2][3] = eq3[3];
    coefficient=(-1.0*matrice[1][0]/matrice[0][0]);
    for(int p = 0; p<3;p++)
        for(int l = 0; l<4; l++)
            //qDebug() <<"matrice : "<<matrice[p][l]<<endl;
    for(;i<=3;i++)
    {
        matrice[1][i]=(coefficient*matrice[0][i])+matrice[1][i];
    }
    coefficient=(-1.0*matrice[2][0]/matrice[0][0]);
    i=0;
    for(;i<=3;i++)
    {
        matrice[2][i]=(coefficient*matrice[0][i])+matrice[2][i];
    }
    coefficient=(-1.0*matrice[2][1]/matrice[1][1]);
    i=1;
    for(;i<=3;i++)
    {
        matrice[2][i]=(coefficient*matrice[1][i])+matrice[2][i];
    }
    z=matrice[2][3]/matrice[2][2];
    y=(matrice[1][3]-(matrice[1][2]*z))/matrice[1][1];
    x=(matrice[0][3]-((matrice[0][1]*y)+(matrice[0][2]*z)))/matrice[0][0];
    //qDebug()  << "X est egal a " << x << "\n";
    //qDebug()  << "Y est egal a " << y << "\n";
    //qDebug()  << "Z est egal a " << z << "\n";
    MyMesh::Point point(x,y,z);
    //system("PAUSE");
    return point;

}

//---Découpage en 2 - fonctionne avec des maillages conséquents (difficile avec cube)
/*MyMesh * MainWindow::cut(MyMesh * _mesh, MyMesh::Point P1, MyMesh::Point P2, MyMesh::Point P3)
{
    MyMesh frag1;
    MyMesh frag2;
//    MyMesh::VertexHandle vh1 = _mesh->add_vertex(P1);
//    MyMesh::VertexHandle vh2 = _mesh->add_vertex(P2);
//    MyMesh::VertexHandle vh3 = _mesh->add_vertex(P3);

    float * eq_plane = equation_plane(P1,P2,P3);
//    for(MyMesh::VertexIter vi;vi!=_mesh->vertices_end(); vi++)
//    {

//    }
//    frag1.add_vertex(P1);
//    frag1.add_vertex(P2);
//    frag1.add_vertex(P3);

//    frag2.add_vertex(P1);
//    frag2.add_vertex(P2);
//    frag2.add_vertex(P3);


//    for(int i=0;i < _mesh->n_vertices();i++)
//    {
//        MyMesh::Point PM = _mesh->point(_mesh->vertex_handle(i));

//        if((PM[0]*eq_plane[0]+PM[1]*eq_plane[1]+PM[2]*eq_plane[2]+eq_plane[3])>=0)
//        {
//            frag1.add_vertex(PM);
//        }
//        else
//        {
//            frag2.add_vertex(PM);
//        }
//    }
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
}*/
MyMesh * MainWindow::cutV2(MyMesh * _mesh, MyMesh::Point P1, MyMesh::Point P2, MyMesh::Point P3)
{
    double * eq_plane = equation_plane(P2,P3,P1);

    //parours des sommets
    //tri les sommets dans 2 mesh différents en fonction du plan de coupe
    for (MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end(); v++){
        MyMesh::Point sommet = mesh.point(*v);
        double x = static_cast<double>(sommet[0]);
        double y = static_cast<double>(sommet[1]);
        double z = static_cast<double>(sommet[2]);

        if((x*eq_plane[0]+y*eq_plane[1]+z*eq_plane[2]+eq_plane[3]) >= 0.0)
        {
            MyMesh::VertexHandle vh = frag1.add_vertex(sommet);
            frag1.set_color(vh,MyMesh::Color(0,0,255));
            frag1.data(vh).thickness = 3;
        }
        else
        {
            MyMesh::VertexHandle vh = frag2.add_vertex(sommet);
            frag2.set_color(vh,MyMesh::Color(255,0,0));
            frag2.data(vh).thickness = 3;
        }
    }

    /*
     * Ajout des points d'intersections dans les deux frag
     */
    for(MyMesh::EdgeIter e = _mesh->edges_begin(); e != _mesh->edges_end(); e++){
        EdgeHandle eh = *e;
        HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
        VertexHandle vh1 = mesh.from_vertex_handle(heh);
        VertexHandle vh2 =  mesh.to_vertex_handle(heh);
        MyMesh::Point sommet1 = mesh.point(vh1);
        MyMesh::Point sommet2 = mesh.point(vh2);

        double result1 = eq_plane[0] * sommet1[0] + eq_plane[1] * sommet1[1] + eq_plane[2] * sommet1[2] + eq_plane[3];
        double result2 = eq_plane[0] * sommet2[0] + eq_plane[1] * sommet2[1] + eq_plane[2] * sommet2[2] + eq_plane[3];

        //TODO gérer le = 0
        if (!((result1 < 0 && result2 > 0) || (result1 > 0 && result2 < 0))){ //pas intersection
            //qDebug() << "pas intersecion";
            continue;
        }

        //prendre un 3eme point c
        MyMesh::Point Pc = MyMesh::Point(sommet1[0] + 0.2, sommet1[1] + 0.2, sommet1[2] + 0.2); //valeur arbitraire
        //plan sommet2 sommet2 c
        double * eq_plane_pc = equation_plane(sommet1,sommet2,Pc);
        //produit vectorielle  -> 4eme point tmp
        std::vector<double> d = crossProduct(sommet2-sommet1, Pc-sommet1); //produit vectoriel permettant de calculer un plan perpendiculaire au plan de Pc, s1, s2
        MyMesh::Point Pd(sommet1[0]+d[0], sommet1[1]+d[1], sommet1[2]+d[2]);
        //plan sommet1 sommet2 d
        double * eq_plane_pd = equation_plane(sommet1,sommet2,Pd);

        /*  INUTILE -> prefere garder 3 équations pour permettre la resolution
        //intersection plan 1 et plan 2 --> equation droite
        float * eq_droite = equation_droite(eq_plane_pc, eq_plane_pd);
        */

        //eq_plane_pc VS eq_plane_pd VS eq_plane -> systeme 3 equation 3 inconnu
        MyMesh::Point I = resol_3eq_3inc(eq_plane, eq_plane_pc, eq_plane_pd);
        //les points d'intersections sont inversés par rapport au maillage
        I[0] = -I[0];
        I[1] = -I[1];
        I[2] = -I[2];

        //ajout du sommet dans les 2 fragments
        frag1.add_vertex(I);
        frag2.add_vertex(I);
        VertexHandle vh = frag3.add_vertex(I);
        frag3.set_color(vh, MyMesh::Color(255,255,0));
        frag3.data(vh).thickness = 5;
    }

    //  TODO REMAILLAGE

    MyMesh * tabmesh = new MyMesh[3];
    tabmesh[0] = frag1;
    tabmesh[1] = frag2;
    tabmesh[2] = frag3;

    return tabmesh;
}

//calcul de la boite englobante du maillage
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
//    qDebug() << "Boîte englobante :";
//    qDebug() << "    min x : " << be_min_x << "; min y : " << be_min_y << "; min z : " << be_min_z;
//    qDebug() << "    max x : " << be_max_x << "; max y : " << be_max_y << "; max z : " << be_max_z;
}

//calcul l'intersection entre une face triangulaire et une droite passant par 2 points - fonction invalide
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
    if (r < 0.0f || r > 1.0f)                    // ray goes away from triangle
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
//vérifie si un point est à l'intérieur du maillage
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

//retourne un float aléatoire compris entre min et max
float randomFloat(float min, float max)
{
    float random = (static_cast<float>(rand())) / (static_cast<float>(RAND_MAX));
    float range = max - min;
    return (random*range) + min;
}

//créer des points disposés aléatoirement à l'intérieur du maillage
//version 1 basique
/*void MainWindow::createSeeds(){
    boite_englobante();
    for (unsigned i = 0; i < nb_seeds;) {
        //crée un point dans la boite englobante
//        qDebug() <<"x :" << be_min_x << be_max_x;
//        qDebug() <<"y :" << be_min_y << be_max_y;
//        qDebug() <<"z :" << be_min_z << be_max_z;
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
//            qDebug() <<"point " << v.idx() << " is in :" << point[0] << point[1] << point[2];
        }
//        else qDebug() << " is out :" << point[0] << point[1] << point[2];
    }
    displayMesh(&mesh);
}*/
//version 2 - ajout dela division par cube
void MainWindow::createSeeds(){

    boite_englobante();

    for (unsigned i = 0; i < nb_seeds;) {
        //crée un point dans la boite englobante
        float tmp_x = randomFloat(be_min_x, be_max_x);
        float tmp_y = randomFloat(be_min_y, be_max_y);
        float tmp_z = randomFloat(be_min_z, be_max_z);
        MyMesh::Point point = MyMesh::Point(tmp_x, tmp_y, tmp_z);
        //if(pointIsInMesh(point)) { //fonction invalide
            seeds_impact.push_back(point);
            VertexHandle v = seed.add_vertex(point);
            seed.data(v).thickness = 4;
            seed.set_color(v, MyMesh::Color(255, 0, 0));
            i++;
            qDebug() <<"point " << v.idx() << " is in :" << point[0] << point[1] << point[2];
        //}
        //else qDebug() << " is out :" << point[0] << point[1] << point[2];
    }

    divide();
    qDebug()<<"Draw frag : ";

    //probleme de conversion reference
    QVector<MyMesh*> cubes_tmp;
    for(int i = 0; i < cubes.size();i++){
        resetAllColorsAndThickness(&cubes[i]);
        cubes_tmp.push_back(&cubes[i]);
    }
    //cubes_tmp.push_back(&mesh);
    displayMeshes(cubes_tmp);
}
//cross product
std::vector<double> MainWindow::crossProduct(MyMesh::Point p1, MyMesh::Point p2)
{
    std::vector<double> Res;
    Res.push_back(static_cast<double>(p1[1]*p2[2] - p1[2]*p2[1]));
    Res.push_back(static_cast<double>(p1[2]*p2[0] - p1[0]*p2[2]));
    Res.push_back(static_cast<double>(p1[0]*p2[1] - p1[1]*p2[0]));
    return Res;
}

//division en cube a partir des graines
void MainWindow::divide()
{
    //parcours des seeds
    int nbiter = 0;
    //MyMesh cube;

    for (MyMesh::VertexIter v = this->seed.vertices_begin(); v != this->seed.vertices_end();v++)
    {
        cube.clean();
        nbiter++;
        //MyMesh::VertexHandle vh;
        double dist = 10000000.0;
        MyMesh::Point p1 = this->seed.point(*v);
        MyMesh::VertexHandle vertexMesh;
        for (MyMesh::VertexIter vm = this->mesh.vertices_begin(); vm != this->mesh.vertices_end();vm++)
        {
            MyMesh::Point p2 = this->mesh.point(*vm);
            double dist_temp = sqrt(pow(p2[0] - p1[0],2)+pow(p2[1] - p1[1],2)+pow(p2[2] - p1[2],2));
            //double dist_temp = 1;
            if(dist>dist_temp)
            {
                dist = dist_temp;
                vertexMesh = *vm;
                qDebug()<<"distance :"<<dist_temp;
            }
        }
        qDebug()<<"end_dist";

        MyMesh::Point p2 = mesh.point(vertexMesh);

        MyMesh::VertexHandle sommets[8];
        if(p2[0]<p1[0])
        {
            sommets[0] = cube.add_vertex(MyMesh::Point(p2[0], p2[1],  p1[2]));
            sommets[1] = cube.add_vertex(MyMesh::Point( p1[0], p2[1],  p1[2]));
            sommets[2] = cube.add_vertex(MyMesh::Point( p1[0],  p1[1],  p1[2]));
            sommets[3] = cube.add_vertex(MyMesh::Point(p2[0],  p1[1],  p1[2]));
            sommets[4] = cube.add_vertex(MyMesh::Point(p2[0], p2[1], p2[2]));
            sommets[5] = cube.add_vertex(MyMesh::Point( p1[0], p2[1], p2[2]));
            sommets[6] = cube.add_vertex(MyMesh::Point( p1[0],  p1[1], p2[2]));
            sommets[7] = cube.add_vertex(MyMesh::Point(p2[0],  p1[1], p2[2]));

            this->mesh.add_vertex(MyMesh::Point(p2[0], p2[1],  p1[2]));
            this->mesh.add_vertex(MyMesh::Point( p1[0], p2[1],  p1[2]));
            this->mesh.add_vertex(MyMesh::Point(p2[0],  p1[1],  p1[2]));
            this->mesh.add_vertex(MyMesh::Point(p1[0], p2[1], p2[2]));
            this->mesh.add_vertex(MyMesh::Point(p1[0],  p1[1], p2[2]));
            this->mesh.add_vertex(MyMesh::Point(p2[0],  p1[1], p2[2]));
        }
        else
        {
            MyMesh::Point tmp = p1;
            p1 = p2;
            p2 = tmp;
            sommets[0] = cube.add_vertex(MyMesh::Point(p2[0], p2[1],  p1[2]));
            sommets[1] = cube.add_vertex(MyMesh::Point( p1[0], p2[1],  p1[2]));
            sommets[2] = cube.add_vertex(MyMesh::Point( p1[0],  p1[1],  p1[2]));
            sommets[3] = cube.add_vertex(MyMesh::Point(p2[0],  p1[1],  p1[2]));
            sommets[4] = cube.add_vertex(MyMesh::Point(p2[0], p2[1], p2[2]));
            sommets[5] = cube.add_vertex(MyMesh::Point( p1[0], p2[1], p2[2]));
            sommets[6] = cube.add_vertex(MyMesh::Point( p1[0],  p1[1], p2[2]));
            sommets[7] = cube.add_vertex(MyMesh::Point(p2[0],  p1[1], p2[2]));

            this->mesh.add_vertex(MyMesh::Point(p2[0], p2[1],  p1[2]));
            this->mesh.add_vertex(MyMesh::Point( p1[0], p2[1],  p1[2]));
            this->mesh.add_vertex(MyMesh::Point(p2[0],  p1[1],  p1[2]));
            this->mesh.add_vertex(MyMesh::Point(p1[0], p2[1], p2[2]));
            this->mesh.add_vertex(MyMesh::Point(p1[0],  p1[1], p2[2]));
            this->mesh.add_vertex(MyMesh::Point(p2[0],  p1[1], p2[2]));
        }

        // on construit des faces à partir des sommets

        std::vector<MyMesh::VertexHandle> uneNouvelleFace;

        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(sommets[0]);
        uneNouvelleFace.push_back(sommets[1]);
        uneNouvelleFace.push_back(sommets[2]);
        uneNouvelleFace.push_back(sommets[3]);
        cube.add_face(uneNouvelleFace);

        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(sommets[7]);
        uneNouvelleFace.push_back(sommets[6]);
        uneNouvelleFace.push_back(sommets[5]);
        uneNouvelleFace.push_back(sommets[4]);
        cube.add_face(uneNouvelleFace);

        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(sommets[1]);
        uneNouvelleFace.push_back(sommets[0]);
        uneNouvelleFace.push_back(sommets[4]);
        uneNouvelleFace.push_back(sommets[5]);
        cube.add_face(uneNouvelleFace);

        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(sommets[2]);
        uneNouvelleFace.push_back(sommets[1]);
        uneNouvelleFace.push_back(sommets[5]);
        uneNouvelleFace.push_back(sommets[6]);
        cube.add_face(uneNouvelleFace);

        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(sommets[3]);
        uneNouvelleFace.push_back(sommets[2]);
        uneNouvelleFace.push_back(sommets[6]);
        uneNouvelleFace.push_back(sommets[7]);
        cube.add_face(uneNouvelleFace);

        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(sommets[0]);
        uneNouvelleFace.push_back(sommets[3]);
        uneNouvelleFace.push_back(sommets[7]);
        uneNouvelleFace.push_back(sommets[4]);
        cube.add_face(uneNouvelleFace);

        //cube.update_normals();
        //this->mesh = mesh;
        //cube.clear();
        cubes.push_back(cube);
        this->mesh.add_vertex(this->seed.point(*v));

        this->seed.delete_vertex(*v);

        // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
        //resetAllColorsAndThickness(&mesh);
        //delete seed
        qDebug()<<"NB_ITER"<<nbiter;
    }
    QVector<MyMesh> cpy = cubes;
    for(int i = 0; i<cpy.size();i++)
    {
        double distance = 10000000000.0;
        int cubeselect = 0;
        for (int j = 0;j<cpy.size();j++)
        {
            if(i!=j)
            {
                for(MyMesh::VertexIter vcube1 = cpy.at(i).vertices_begin();vcube1 != cpy.at(i).vertices_end();vcube1++)
                {
                    MyMesh::Point p1 = cpy.at(i).point(*vcube1);
                    for(MyMesh::VertexIter vcube2 = cpy.at(j).vertices_begin();vcube2 != cpy.at(j).vertices_end();vcube2++)
                    {
                        MyMesh::Point p2 = cpy.at(j).point(*vcube2);
                        double dist_temp = sqrt(pow(p2[0] - p1[0],2)+pow(p2[1] - p1[1],2)+pow(p2[2] - p1[2],2));
                        if(distance>dist_temp)
                        {
                            distance = dist_temp;
                            cubeselect = j;
                        }
                    }
                }
            }
        }
        qDebug()<<i<<" with "<<cubeselect;

    }
    qDebug() << "after for vertices cubes" << cubes.last().n_vertices();
    qDebug() << "after for vertices cube" << cube.n_vertices();

}
