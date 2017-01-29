///////////////////////            CIRCLES             //////////////////////////////////////////
///////////////////////            CIRCLES             //////////////////////////////////////////


// for all the unknown functions & variables refer to other files of this repository.
 


//structure of a circle- integer version
struct circle{ point_i c; int r;};


//to check whther a point (p) is inside and on the circumference of a circle with a radius(r) and center at point (c) 
int inside&onCircle ( point p, point c, double r )
{
    double dx=p.x - c.x, dy = p.y - c.y;
    double Euc = dx*dx + dy*dy, rSq = r*r;
    return Euc <= rSq ? 1 : 0;
}


//to check whther a point (p) is inside a circle with a radius(r) and center at point (c) 
int insideCircle ( point p, point c, double r )
{
    double dx=p.x - c.x, dy = p.y - c.y;
    double Euc = dx*dx + dy*dy, rSq = r*r;
    return Euc < rSq ? 1 : 0;
}


// to get the circle made from the two given points (p1) & (p2) and given radius (r)
// return false if no such circle exist else true
// to get the other possible circle just pass the two points (p1) & (p2) in opposite order
bool circle2PtsRad(point p1, point p2, double r, point &c)
{
    double d2 = (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y) *(p1.y-p2.y);
    double det=r*r / d2 -0.25;

    if(det < 0.0) return false;
    double h=sqrt(det);
    c.x=(p1.x +p2.x) * 0.5 + (p1.y - p2.y) *h;
    c.y = (p1.y + p2.y) * 0.5 +(p2.x - p1.x)* h;
    return true;
}
