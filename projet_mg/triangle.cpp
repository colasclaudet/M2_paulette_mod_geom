#include "triangle.h"
//contructeur
Triangle::Triangle(MyMesh::Point v1, MyMesh::Point v2, MyMesh::Point v3)
{
    this->point1 = v1;
    this->point2 = v2;
    this->point3 = v3;

}
//On vérifie si un point appartient au triangle
bool Triangle::containsVertex(MyMesh::Point v)
{
    // return p1 == v || p2 == v || p3 == v;
    return almost_equal(this->point1, v) || almost_equal(this->point2, v) || almost_equal(this->point3, v);
}

//retourne la normale au triangle
MyMesh::Point Triangle::getNorm()
{
    MyMesh::Point vectAB = point2 - point1;
    MyMesh::Point vectBC = point3 - point2;

    float vx = vectAB[1] * vectBC[2] - vectAB[2] * vectBC[1];

    float vy = vectAB[2] * vectBC[0] - vectAB[0] * vectBC[2];
    float vz = vectAB[0] * vectBC[1] - vectAB[1] * vectBC[0];
    MyMesh::Point vn(vx,vy,vz);
    MyMesh::Point n = vn/(sqrt(pow(vx,2)+pow(vy,2)+pow(vz,2)));
    return n;
}

//retourne le veteur directeur entre 2 points
MyMesh::Point Triangle::getVect(MyMesh::Point p1, MyMesh::Point p2)
{
    return p2-p1;
}

//retourne la distance entre 2 points
double Triangle::dist2(MyMesh::Point v1, MyMesh::Point v2)
{
    const double dx = v1[0] - v2[0];
    const double dy = v1[1] - v2[1];
    const double dz = v1[2] - v2[2];
    return dx * dx + dy * dy + dz*dz;
}

//retourne le centre de gravité du triangle
MyMesh::Point Triangle::return_gravity_center()
{
    float x = (this->point1[0]+this->point2[0]+this->point3[0])/3;
    float y = (this->point1[1]+this->point2[1]+this->point3[1])/3;
    float z = (this->point1[2]+this->point2[2]+this->point3[2])/3;
    return MyMesh::Point(x,y,z);
}

//retourne l'équation de plan du triangle
float * Triangle::equation_plane(MyMesh::Point P1
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

/*
 * on cherche à vérifier si un point appartient à la sphère englobante du triangle
 * Pour cela on simplifie les calculs en se servant simplement du centre de gravité
 * du triangle pour former une pyramide dans la direction de la normale au triangle
 * principal et ensuite vérifier si le point est contenu dans cette pyramide
*/
bool Triangle::circumCircleContains3D(MyMesh::Point v)
{
    const double ab = this->getPoint1().norm();
    const double cd = this->getPoint2().norm();
    const double ef = this->getPoint3().norm();

    const double ax = this->getPoint1()[0];
    const double ay = this->getPoint1()[1];
    const double az = this->getPoint1()[2];

    const double bx = this->getPoint2()[0];
    const double by = this->getPoint2()[1];
    const double bz = this->getPoint2()[2];

    const double cx = this->getPoint3()[0];
    const double cy = this->getPoint3()[1];
    const double cz = this->getPoint3()[2];
    MyMesh::Point pcenter(0.0,0.0,-10);

    MyMesh::Point triangle_norm = this->getNorm();
    const double ox = triangle_norm[0] * - (sqrt(pow(this->point2[0] - this->point1[0],2) + pow(this->point2[1] - this->point1[1],2) + pow(this->point2[2] - this->point1[2],2))); //10
    const double oy = triangle_norm[1] * - (sqrt(pow(this->point3[0] - this->point2[0],2) + pow(this->point3[1] - this->point2[1],2) + pow(this->point3[2] - this->point2[2],2))); //10
    const double oz = triangle_norm[2] * - (sqrt(pow(this->point3[0] - this->point1[0],2) + pow(this->point3[1] - this->point1[1],2) + pow(this->point3[2] - this->point1[2],2))); //10

    /*const double ox = 0.0;
    const double oy = 0.0;
    const double oz = -10.0;*/
    MyMesh::Point o(ox,oy,oz);
    const double vpx = v[0];
    const double vpy = v[1];
    const double vpz = v[2];

    MyMesh::Point OP = getVect(o,v);
    MyMesh::Point g = this->return_gravity_center();
    MyMesh::Point OG = getVect(o,g);
    Triangle OAC(o,point1,point3);
    Triangle OBC(o,point2,point3);
    Triangle OAB(o,point1,point2);

    MyMesh::Point calc_for_OAB = OAB.getNorm() * OP * OAB.getNorm() * OG;
    MyMesh::Point calc_for_OAC = OAC.getNorm() * OP * OAC.getNorm() * OG;
    MyMesh::Point calc_for_OBC = OBC.getNorm() * OP * OBC.getNorm() * OG;

    if(calc_for_OAB[0]>=0 && calc_for_OAB[1]>=0 && calc_for_OAB[2]>=0
            && calc_for_OBC[0]>=0 && calc_for_OBC[1]>=0 && calc_for_OBC[2]>=0
            && calc_for_OAC[0]>=0 && calc_for_OAC[1]>=0 && calc_for_OAC[2]>=0)
    {
        return true;
    }
    else
    {
        return false;
    }
    //création des vecteurs OP et OG maintenant il faut trouver la normale de tout les triangles et les push dans un vecteur
    /*const double circum_x = (ab * (cy - by) + cd * (ay - cy) + ef * (by - ay)) / (ax * (cy - by) + bx * (ay - cy) + cx * (by - ay));
    const double circum_y = (ab * (cx - bx) + cd * (ax - cx) + ef * (bx - ax)) / (ay * (cx - bx) + by * (ax - cx) + cy * (bx - ax));
    const double circum_z = (ab * (cx - bx) + cd * (ax - cx) + ef * (bx - ax)) / (az * (cx - bx) + bz * (ax - cx) + cz * (bx - ax));

    const MyMesh::Point circum(0.5*circum_x, 0.5*circum_y,0.5*circum_z);
    const double circum_radius = dist2(this->point1,circum);
    const double dist = dist2(v,circum);
    return dist <= circum_radius;*/
}

bool Triangle::circumCircleContains2D(MyMesh::Point v)
{
    const double ab = this->getPoint1().norm();
    const double cd = this->getPoint2().norm();
    const double ef = this->getPoint3().norm();


    const double ax = this->getPoint1()[0];
    const double ay = this->getPoint1()[1];
    //const double az = this->getPoint1()[2];

    const double bx = this->getPoint2()[0];
    const double by = this->getPoint2()[1];
    //const double bz = this->getPoint2()[2];

    const double cx = this->getPoint3()[0];
    const double cy = this->getPoint3()[1];
    //const double cz = this->getPoint3()[2];

    const double circum_x = (ab * (cy - by) + cd * (ay - cy) + ef * (by - ay)) / (ax * (cy - by) + bx * (ay - cy) + cx * (by - ay));
    const double circum_y = (ab * (cx - bx) + cd * (ax - cx) + ef * (bx - ax)) / (ay * (cx - bx) + by * (ax - cx) + cy * (bx - ax));
    //const double circum_z = (ab * (cx - bx) + cd * (ax - cx) + ef * (bx - ax)) / (az * (cx - bx) + bz * (ax - cx) + cz * (bx - ax));

    const MyMesh::Point circum(0.5*circum_x, 0.5*circum_y,0);//,0.5*circum_z);
    const double circum_radius = dist2(this->point1,circum);
    const double dist = dist2(v,circum);
    return dist <= circum_radius;
}

bool Triangle::operator ==(const Triangle tr)
{
    return	(this->point1 == tr.point1 || this->point1 == tr.point2 || this->point1 == tr.point3) &&
                (this->point2 == tr.point1 || this->point2 == tr.point2 || this->point2 == tr.point3) &&
            (this->point3 == tr.point1|| this->point3 == tr.point2 || this->point3 == tr.point3);
}

MyMesh::Point Triangle::getPoint1() const
{
    return point1;
}

void Triangle::setPoint1(const MyMesh::Point &value)
{
    point1 = value;
}

MyMesh::Point Triangle::getPoint2() const
{
    return point2;
}

void Triangle::setPoint2(const MyMesh::Point &value)
{
    point2 = value;
}

MyMesh::Point Triangle::getPoint3() const
{
    return point3;
}

void Triangle::setPoint3(const MyMesh::Point &value)
{
    point3 = value;
}


bool Triangle::almost_equal(const double x, const double y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<double>::epsilon() * std::abs(x+y) * static_cast<double>(ulp)
    // unless the result is subnormal
        || std::abs(x-y) < std::numeric_limits<double>::min();
}

bool Triangle::almost_equal(MyMesh::Point v1, MyMesh::Point v2, int ulp)
{
    return almost_equal(v1[0], v2[0], ulp) && almost_equal(v1[1], v2[1], ulp) && almost_equal(v1[2], v2[2], ulp); //Z
}

bool Triangle::almost_equal(const Triangle t1, const Triangle t2)
{
    return	(almost_equal(t1.getPoint1() , t2.getPoint1()) || almost_equal(t1.getPoint1() , t2.getPoint2()) || almost_equal(t1.getPoint1() , t2.getPoint3())) &&
            (almost_equal(t1.getPoint2(), t2.getPoint1()) || almost_equal(t1.getPoint2() , t2.getPoint2()) || almost_equal(t1.getPoint2() , t2.getPoint3())) &&
            (almost_equal(t1.getPoint3() , t2.getPoint1()) || almost_equal(t1.getPoint3() , t2.getPoint2()) || almost_equal(t1.getPoint3() , t2.getPoint3()));
}
bool Triangle::containVertex(MyMesh::Point p)
{
    return almost_equal(getPoint1(),p)||almost_equal(getPoint2(),p)||almost_equal(getPoint3(),p);
}
