//////////////////////                VECTORS              ///////////////////////
struct vec
{
    double x,y;
    vec () { x=y=0.0; }
    vec(double _x, double _y) : x(_x), y(_y) {}
};

vec toVec(point a, point b) { return vec(b.x - a.x, b.y - a.y); }

vec scale ( vec v, double s ) { return vec(v.x * s, v.y * s); }

point translate (point p, vec v) { return point(p.x +v.x, p.y+v.y); }

double dot(vec a, vec b) { return (a.x *b.x + a.y*b.y); }

double norm_sq(vec v) { return v.x*v.x + v.y*v.y ;}

double distToLine( point p,point a, point b, point &c)
{
    vec ap=toVec(a,p), ab=toVec(a,b);
    double u= dot(ap, ab)/ norm_sq(ab);
    c= translate(a,scale(ab,u));

    return dist(p,c);
}

double distToLineSegment (point p, point a, point b, point &c)
{
    vec ap=toVec(a,p), ab=toVec(a,b);
    double u=dot(ap,ab)/norm_sq(ab);

    if(u<0.0)
    {
        c=point(a.x,a.y);
        return dist(p,a);
    }

    if(u>1.0)
    {
        c=point(b.x,b.y);
        return dist(p,b);
    }

    return distToLine(p,a,b,c);
}

//remember about order of points
double angle (point a, point o, point b)
{
    vec oa=toVec(o,a), ob=toVec(o,b);
    return acos(dot(oa,ob)/ sqrt(norm_sq(oa)* norm_sq(ob)));
}

double cross(vec a, vec b) { return a.x * b.y - a.y * b.x; }

int ccw (point p, point q, point r)
{
    double t=cross(toVec(p,q), toVec(p,r) );
    if(t>0) return true;
    else return false;
}

bool onSegment(point p, point q, point r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
       return true;

    return false;
}

bool collinear( point p, point q, point r) { return fabs(cross(toVec(p,q), toVec(p,r))) < EPS; }
