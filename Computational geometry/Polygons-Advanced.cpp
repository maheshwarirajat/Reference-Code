///////////////////////////          POLYGONS               //////////////////////////////

double perimeter (const vector<point> &P)
{
    double result  =0.0;
    for(int i=0 ;i< (int)P.size() -1 ; i++)
        result +=dist(P[i], P[i+1]);
    return result;
}

double area(const vector<point> &P)
{
    double result =0.0, x1,y1,x2,y2;
    for(int i=0; i< (int)P.size()-1; i++)
    {
        x1=P[i].x; x2=P[i+1].x;
        y1=P[i].y; y2=P[i+1].y;
        result+= (x1*y2 - x2*y1);
    }

    return fabs(result)/2.0;
}

bool isConvex( const vector<point> &P)
{
    int sz =(int)P.size();
    if(sz <=3) return false;
    bool isLeft =ccw(P[0], P[1],P[2]);
    for(int i=1; i<sz-1; i++)
        if(ccw(P[i], P[i+1], P[(i+2) == sz ? 1 : i+2])!= isLeft) return false;
    return true;
}



bool inPolygon( point pt, const vector<point> &P )
{
    if((int)P.size()==0) return false;
    double sum=0;
    for(int i=0; i< (int)P.size()-1; i++)
    {
        if(ccw(pt, P[i], P[i+1]))
            sum+= angle(P[i], pt, P[i+1]);
        else sum-= angle(P[i], pt, P[i+1]);
    }

    return fabs(fabs(sum)- 2*PI) < EPS;
}


point lineIntersectSeg(point p, point q, point A, point B)
{
    double a = B.y - A.y;
    double b = A.x - B.x;
    double c =  B.x * A.y - A.x * B.y;
    double u = fabs(a * p.x + b * p.y + c);
    double v = fabs(a * q.x + b * q.y + c);
    return point( (p.x * v + q.x *u) / (u+v), (p.y * v + q.y * u) / (u+v) );
}

vector<point> cutPolygon(point a, point b, const vector<point> &Q)
{
    vector<point> P;
    for(unsigned int i=0; i<Q.size(); i++)
    {
        double left1= cross(toVec(a,b), toVec(a,Q[i])), left2=0;
        if(i!= Q.size()-1) left2 =cross(toVec(a,b), toVec(a,Q[i+1]));
        if(left1 >-EPS ) P.push_back(Q[i]);
        if(left1 * left2 < -EPS) P.push_back(lineIntersectSeg(Q[i],Q[i+1],a,b));

    }

    if(!P.empty() && !(P.front()==P.back())) P.push_back(P.front());
    return P;
}
