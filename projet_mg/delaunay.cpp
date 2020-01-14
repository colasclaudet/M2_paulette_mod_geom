#include "delaunay.h"

//contructeur
Delaunay::Delaunay(MyMesh * _mesh)
{
    this->mesh = _mesh;
}
/*
 * ancienne fonction de triangulation (que pour 2 dimentions)
*/
std::vector<Triangle> Delaunay::triangulate2()
{

    float minx = this->mesh->point(this->mesh->vertex_handle(0))[0];
            float maxx = this->mesh->point(this->mesh->vertex_handle(0))[0];
    float miny = this->mesh->point(this->mesh->vertex_handle(0))[1];
            float maxy = this->mesh->point(this->mesh->vertex_handle(0))[1];

    float minz = this->mesh->point(this->mesh->vertex_handle(0))[2];
            float maxz = this->mesh->point(this->mesh->vertex_handle(0))[2];

    //calcul des mins et max du maillages
    for(MyMesh::VertexIter v = this->mesh->vertices_begin(); v != this->mesh->vertices_end(); v++)
    {
        if(this->mesh->point(*v)[0]<minx)
            minx = this->mesh->point(*v)[0];
        if(this->mesh->point(*v)[1]<miny)
            miny = this->mesh->point(*v)[1];
        if(this->mesh->point(*v)[2]<minz)
            minz = this->mesh->point(*v)[2];

        if(this->mesh->point(*v)[0]>maxx)
            maxx = this->mesh->point(*v)[0];
        if(this->mesh->point(*v)[1]>maxy)
            maxy = this->mesh->point(*v)[1];  
        if(this->mesh->point(*v)[2]>maxz)
            maxz = this->mesh->point(*v)[2];
    }
    float dx = maxx - minx;
    float dy = maxy - miny;
    float dz = maxz - minz;

    double dmax = std::max(dx,dy);
    dmax = std::max((float)dmax,dz);
    //for Z
    //dmax = std::max(dmax,dz);
    //
    float midx = 0.5*(minx + maxx);
    float midy = 0.5*(miny + maxy);
    float midz = 0.5*(minz + maxz);

    /*MyMesh::Point p1(midx - 20 * dmax, midy - dmax, midz - dmax);
    MyMesh::Point p2(midx,midy+20*dmax,midz);
    MyMesh::Point p3(midx+20*dmax,midy-dmax,midz-dmax);*/
    MyMesh::Point p1(midx -  dmax, midy - dmax, midz - dmax);
    MyMesh::Point p2(midx,midy+dmax,midz);
    MyMesh::Point p3(midx+dmax,midy-dmax,midz-dmax);
    //Triangle englobant tout le maillage
    Triangle tr(p1,p2,p3);
    triangles.push_back(tr);
    //on parcourt tout les vertex et on regarde ceux qui appartiennent au triangle global
    for(MyMesh::VertexIter v = this->mesh->vertices_begin(); v != this->mesh->vertices_end(); v++)
    {
        std::vector<Edge> geom;
        for(unsigned int i = 0; i<triangles.size(); i++)
        {
            if(triangles.at(i).circumCircleContains3D(this->mesh->point(*v)))
            {
                triangles.at(i).isBad = true;
                Edge e;
                e.p1 = triangles.at(i).getPoint1();
                e.p2 = triangles.at(i).getPoint2();
                geom.push_back(e);
                e.p1 = triangles.at(i).getPoint2();
                e.p2 = triangles.at(i).getPoint3();
                geom.push_back(e);
                e.p1 = triangles.at(i).getPoint3();
                e.p2 = triangles.at(i).getPoint1();
                geom.push_back(e);
                //edges.push_back(e);
                //MyMesh::Edge e()
                //geom.push_back()
            }
        }
        /*triangles.erase(std::remove_if(begin(triangles),end(triangles),[](Triangle &t)
        {
            return t.isBad;
        }),end(triangles));*/
        //test d'égalité entre les edges
        for(unsigned int i = 0; i<geom.size();i++)
        {
            for(unsigned int j = 1; j<geom.size();j++)
            {

                if(tr.almost_equal(geom.at(i).p1,geom.at(j).p1) ||
                        tr.almost_equal(geom.at(i).p2,geom.at(j).p1)||
                        tr.almost_equal(geom.at(i).p2,geom.at(j).p2))
                {
                    geom.at(i).isBad = true;
                    geom.at(j).isBad = true;
                }
            }
        }
        geom.erase(std::remove_if(begin(geom), end(geom), [](Edge &e){
                    return e.isBad;
                }), end(geom));
        //création des triangles
        for(unsigned int i = 0; i<geom.size();i++)
        {
            Triangle tri(geom.at(i).p1,geom.at(i).p2,mesh->point(*v));
            triangles.push_back(tri);
        }
    }
    /*triangles.erase(std::remove_if(begin(triangles), end(triangles), [p1, p2, p3](Triangle &t){
            return t.containsVertex(p1) || t.containsVertex(p2) || t.containsVertex(p3);
        }), end(triangles));*/
    //création des edges
    for(unsigned int i = 0; i<triangles.size();i++)
    {
        Edge e;
        e.p1 = triangles.at(i).getPoint1();
        e.p2 = triangles.at(i).getPoint2();
        edges.push_back(e);
        e.p1 = triangles.at(i).getPoint2();
        e.p2 = triangles.at(i).getPoint3();
        edges.push_back(e);
        e.p1 = triangles.at(i).getPoint3();
        e.p2 = triangles.at(i).getPoint1();
        edges.push_back(e);
    }
    return triangles;
}

/*
 * Fonction de triangulation de Delaunay
 * on commence par faire une triangulation globale des seeds
 * Ensuite on va vérifier pour chaque triangle obtenu si les vertex appartiennent à l'un des triangles de la triangulation
 * si il appartient au triangle on forme alors 3 nouveaux triangles qui ont pour sommet commun le vextex en question
*/

std::vector<Triangle> Delaunay::triangulate()
{

    float minx = this->mesh->point(this->mesh->vertex_handle(0))[0];
            float maxx = this->mesh->point(this->mesh->vertex_handle(0))[0];
    float miny = this->mesh->point(this->mesh->vertex_handle(0))[1];
            float maxy = this->mesh->point(this->mesh->vertex_handle(0))[1];

    float minz = this->mesh->point(this->mesh->vertex_handle(0))[2];
            float maxz = this->mesh->point(this->mesh->vertex_handle(0))[2];
    //calcul des mins et max du maillages
    for(MyMesh::VertexIter v = this->mesh->vertices_begin(); v != this->mesh->vertices_end(); v++)
    {
        if(this->mesh->point(*v)[0]<minx)
            minx = this->mesh->point(*v)[0];
        if(this->mesh->point(*v)[1]<miny)
            miny = this->mesh->point(*v)[1];
        if(this->mesh->point(*v)[2]<minz)
            minz = this->mesh->point(*v)[2];

        if(this->mesh->point(*v)[0]>maxx)
            maxx = this->mesh->point(*v)[0];
        if(this->mesh->point(*v)[1]>maxy)
            maxy = this->mesh->point(*v)[1];
        if(this->mesh->point(*v)[2]>maxz)
            maxz = this->mesh->point(*v)[2];
    }
    float dx = maxx - minx;
    float dy = maxy - miny;
    float dz = maxz - minz;

    double dmax = std::max(dx,dy);
    dmax = std::max((float)dmax,dz);
    //for Z
    //dmax = std::max(dmax,dz);
    //
    float midx = 0.5*(minx + maxx);
    float midy = 0.5*(miny + maxy);
    float midz = 0.5*(minz + maxz);

    MyMesh::Point p1(midx -  dmax, midy - dmax, 0);
    MyMesh::Point p2(midx,midy+dmax,0);
    MyMesh::Point p3(midx+dmax,midy-dmax,0);

    //Triangle englobant tout le maillage
    Triangle triangle_main(p1,p2,p3);
    triangles.push_back(triangle_main);
    MyMesh::Point gravity_center_main = triangle_main.return_gravity_center();

    unsigned int nb_iteration = 0;
    //on parcourt tout les vertex et on regarde ceux qui appartiennent au triangle global
    for(MyMesh::VertexIter v = this->mesh->vertices_begin(); v != this->mesh->vertices_end(); v++)
    {
        nb_iteration++;
        qDebug()<<nb_iteration;
        std::vector<Edge> geom;
        //problèmede boucle infinie
        if(nb_iteration>600)
        {
            break;
        }
        for(unsigned int i = 0; i < triangles.size(); i++)
        {
            if(triangles.at(i).circumCircleContains3D(mesh->point(*v))) //on vérifie si le point V peut être projeté sur le triangle
            {
                qDebug()<<"new triangle";
                MyMesh::Point P = mesh->point(*v);
                MyMesh::Point P1 = triangles.at(i).getPoint1();
                MyMesh::Point P2 = triangles.at(i).getPoint2();
                MyMesh::Point P3 = triangles.at(i).getPoint3();

                Triangle T1(P1,P2,P);
                Triangle T2(P1,P3,P);
                Triangle T3(P2,P3,P);
                triangles.push_back(T1);
                triangles.push_back(T2);
                triangles.push_back(T3);

                //triangles.erase(triangles.begin()+i);
                //triangles.at(i).isBad = true; //si le point est bon , on supprime le triangle initial
                nb_iteration++;
                qDebug()<<nb_iteration;

                Edge e;
                e.p1 = triangles.at(i).getPoint1();
                e.p2 = triangles.at(i).getPoint2();
                geom.push_back(e);
                e.p1 = triangles.at(i).getPoint2();
                e.p2 = triangles.at(i).getPoint3();
                geom.push_back(e);
                e.p1 = triangles.at(i).getPoint3();
                e.p2 = triangles.at(i).getPoint1();
                geom.push_back(e);

                if(nb_iteration>300)
                {
                    break;
                }
            }

        }
        for(unsigned int i = 0; i< triangles.size(); i++)
        {
            if(triangles.at(i).isBad)
            {
                triangles.erase(triangles.begin()+i); //suppression des triangles non valides
                qDebug()<<"erase "<< i;
            }
        }

        for(unsigned int i = 0; i<geom.size();i++)
        {
            for(unsigned int j = 1; j<geom.size();j++)
            {

                if(triangle_main.almost_equal(geom.at(i).p1,geom.at(j).p1) ||
                        triangle_main.almost_equal(geom.at(i).p2,geom.at(j).p1)||
                        triangle_main.almost_equal(geom.at(i).p2,geom.at(j).p2))
                {
                    geom.at(i).isBad = true;
                    geom.at(j).isBad = true;
                }
            }
        }
        for(unsigned int i = 0; i < geom.size(); i++)
        {
            if(geom.at(i).isBad)
            {
                geom.erase(geom.begin()+i);
            }
        }

        for(unsigned int i = 0; i<geom.size();i++)
        {
            Triangle tri(geom.at(i).p1,geom.at(i).p2,mesh->point(*v));
            triangles.push_back(tri);
        }
    }

    for(unsigned int i = 0; i<triangles.size();i++)
    {
        Edge e;
        e.p1 = triangles.at(i).getPoint1();
        e.p2 = triangles.at(i).getPoint2();
        edges.push_back(e);
        e.p1 = triangles.at(i).getPoint2();
        e.p2 = triangles.at(i).getPoint3();
        edges.push_back(e);
        e.p1 = triangles.at(i).getPoint3();
        e.p2 = triangles.at(i).getPoint1();
        edges.push_back(e);
    }
    return triangles;
}

std::vector<Triangle> Delaunay::getTriangles() const
{
    return triangles;
}
/*
 * fonction de génération du maillage à afficher après la triangulation
*/
MyMesh Delaunay::draw_mesh()
{
    MyMesh _mesh;

    int b = 10;
    int r = 255;
    int g = 125;
    for(unsigned int i = 0; i < this->edges.size()-1;i++)
        {
            MyMesh::VertexHandle v1 = _mesh.add_vertex(edges.at(i).p1);
            MyMesh::VertexHandle v2 = _mesh.add_vertex(edges.at(i).p2);
            MyMesh::VertexHandle v3 = _mesh.add_vertex(edges.at(i+1).p1);
            _mesh.data(v1).thickness = 10;
            _mesh.set_color(v1, MyMesh::Color(r%255, g%255, b%255));
            _mesh.data(v2).thickness = 10;
            _mesh.set_color(v2, MyMesh::Color(r%255, g%255, b%255));
            b+=10;
            g+=20;
            r+=5;
            std::vector <MyMesh::VertexHandle> pf;
            pf.push_back(v1);
            pf.push_back(v2);
            pf.push_back(v3);
            //_mesh.add_face(pf);


            if(pf.size()>=3)
                //_mesh.add_face(pf);
            qDebug()<<"FACE";
        }
    return _mesh;
    //on retourne un mesh composé de triangles
}
