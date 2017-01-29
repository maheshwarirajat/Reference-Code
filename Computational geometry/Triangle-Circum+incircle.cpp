//////////////////////                TRIANGLES           ////////////////////////////////////////

double heron( double ab, double bc, double ca)
{
    double s= (ab+bc+ca)/2;
    return sqrt(s*(s-ab)*(s-bc)*(s-ca));
}

double rInCircle(double ab, double bc, double ca) { return heron(ab,bc,ca)/ (0.5 * (ab+bc+ca)); }
double rInCircle(point a, point b, point c) { return rInCircle(dist(a,b),dist(b,c), dist(c,a) ); }

int inCircle ( point p1, point p2, point p3, point &ctr, double &r )
{
    r = rInCircle(p1,p2,p3);
    if (fabs(r) < EPS ) return 0;

    line l1,l2;
    double ratio = dist(p1, p2) / dist(p2,p3);
    point p=translate(p2, scale(toVec(p2,p3), ratio / (1+ratio)));

    ratio = dist(p2,p1) / dist(p2,p3);
    p= translate(p1, scale(toVec(p1,p3), ratio / (1 + ratio)));
    pointsToLine(p2, p, l2);

    areIntersect(l1, l2, ctr);
    return 1;
}

double rCircumCircle(double ab, double bc, double ca) { return ab* bc * ca / (4.0 * heron(ab, bc, ca)); }
double rCircumCircle( point a, point b, point c) { return rCircumCircle(dist(a, b), dist(b, c), dist(c, a)); }

double perimeterTriangle(point a, point b, point c) { return dist(a,b)+dist(b,c)+dist(c,a); }

bool check_triangle (double a, double b, double c) { if((a+b > c) && (b+c > a) && (a+c > b) ) return true; else return false; }

int circumCircle(point p1, point p2, point p3, point &ctr, double &r )
{
    r = rCircumCircle(p1,p2,p3);
    if(fabs(r) < EPS) return 0;

    line l1,l2;
    point mid;
    mid.x=(p1.x+p2.x)/2; mid.y = (p1.y + p2.y)/2;
    p2.x-=mid.x; p2.y-=mid.y;
    double temp;
    temp=p2.x;
    p2.x=-p2.y; p2.y=temp;
    p2.x+=mid.x; p2.y+=mid.y;
    pointsToLine(mid, p2, l1);

    mid.x=(p1.x+p3.x)/2; mid.y = (p1.y + p3.y)/2;
    p3.x-=mid.x; p3.y-=mid.y;
    temp=p3.x;
    p3.x=-p3.y; p3.y=temp;
    p3.x+=mid.x; p3.y+=mid.y;
    pointsToLine(mid, p3, l2);

    areIntersect(l1, l2, ctr);
    return 1;
}
