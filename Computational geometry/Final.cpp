/*    tommy_trash    */
#include<bits/stdc++.h>

#define ll long long
#define PI acos(-1)
#define EPS 1e-9
#define GEOMETRY


using namespace std;

//////////////////////////          POINTS     /////////////////////////////////
struct point_i
{
    int x,y;
    point_i () { x=y=0; }
    point_i (int _x, int _y) : x(_x), y(_y) {}


};

struct point
{
    double x,y;
    point() { x=y=0.0; }
    point (double _x, double _y) : x(_x), y(_y) {}
    bool operator < (point other) const
    {   if (fabs(x - other.x) > EPS)
      return x < other.x;
      return y < other.y; }

    bool operator == (point other) const
    {   if (fabs(x - other.x) < EPS && fabs(y - other.y) < EPS)
                return true;
                return false;
    }


	point operator + (const point &p) const
	{
		return point(x + p.x, y + p.y);
	}
	point operator - (const point &p) const
	{
		return point(x - p.x, y - p.y);
	}

};

double dist(point p1, point p2) {return hypot(p1.x - p2.x, p1.y- p2.y);}

point rotate_theta(point p, double theta)
{
    double rad =(theta*PI)/180;
    return point (p.x * cos(rad) - p.y * sin(rad), p.x * sin(rad) + p.y * cos(rad));
}

/////////////////////////////            LINES         ///////////////////////////////

struct line { double a, b, c; };
struct line_m{ double m,c; };

void pointsToLine(point p1, point p2, line &l )
{
    if (fabs(p1.x -p2.x) < EPS) { l.a = 1.0; l.b = 0.0; l.c = -p1.x;}
    else
    {
        l.a = -(double)(p1.y - p2.y)/ (p1.x -p2.x);
        l.b= 1.0;
        l.c = -(double)(l.a *p1.x) - p1.y;
    }
}

bool areParallel(line l1, line l2) { return (fabs(l1.a-l2.a) < EPS) && (fabs(l1.b-l2.b) < EPS); }
bool areSame(line l1, line l2) { return areParallel( l1, l2) && (fabs(l1.c - l2.c) < EPS); }

bool areIntersect( line l1, line l2, point &p)
{
    if(areParallel(l1,l2)) return false;
    p.x= (l2.b * l1.c -l1.b * l2.c) / (l2.a * l1.b - l1.a * l2.b);
    if(fabs(l1.b) > EPS) p.y= -(l1.a * p.x + l1.c);
    else                     p.y= -(l2.a * p.x + l2.c);
    return true;
}

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

/////////////////////////////////          VECTOR -3D        /////////////////////////////////////////////
/////////////////////////////////          VECTOR -3D        /////////////////////////////////////////////

struct point3d
{
    double x,y,z;
    point3d () { x=y=z=0.0; }
    point3d (double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
};

struct vec3d
{
    double x,y,z;
    vec3d () { x=y=z=0.0; }
    vec3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
};

vec3d toVec3d(point3d a, point3d b) { return vec3d(b.x - a.x, b.y - a.y, b.z - a.z ); }


///////////////////////            CIRCLES             //////////////////////////////////////////
//integer version
struct circle{ point_i c; int r;};

int insideCircle ( point p, point c, double r )
{
    double dx=p.x - c.x, dy = p.y - c.y;
    double Euc = dx*dx + dy*dy, rSq = r*r;
    return Euc <= rSq ? 1 : 0;
}

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


point pivot(0,0);
bool angleCmp(point a, point b)
{
    if(collinear(pivot,a,b)) return dist(pivot,a) < dist(pivot,b);
    double d1x= a.x-pivot.x, d2x= b.x-pivot.x;
    double d1y= a.y-pivot.y, d2y= b.y-pivot.y;
    return (atan2(d1y, d1x) - atan2(d2y, d2x)) < 0;
}

vector<point> CH(vector<point> P)
{
    //pivot -figure out
    int i,j,n=(int)P.size();

    if(n<=3)
    {
        if(!(P[0]==P[n-1]))P.push_back(P[0]);
        return P;
    }

    int P0=0;
    for( i=1;i<n ;i++)
    {
        if(P[i].y < P[P0].y  || ( P[i].y==P[P0].y && P[i].x < P[P0].x ) ) P0=i;
    }

    point temp = P[0]; P[0]=P[P0]; P[P0]=temp;
    pivot = P[0];
    sort(++P.begin(), P.end(), angleCmp);

    vector<point> S;
    S.push_back(P[n-1]);S.push_back(P[0]);S.push_back(P[1]);
    i=2;
    while (i<n)
    {
        j=(int)S.size()-1;
        if(ccw(S[j-1], S[j], P[i]))
            S.push_back(P[i++]);
        else
        S.pop_back();
    }
    //if(!(S.front()==S.back())) S.push_back(S[0]);
    return S;
}

/*vector<point> andrew(vector<point> P)
{
    sort(P.begin(), P.end());
    int i,j,n=P.size();
    point pminmin=P[0], pmaxmax=P[n-1];

    point pminmax, pmaxmin;

    //calculating pminmax
    for(i=1;i<n;i++) if(P[i].x !=pminmin.x) break;
    pminmax= point(pminmin.x, P[i-1].y);

    //calculating pmaxmax
    for(i=n-2; i>=0; i--) if(P[i].x!= pmaxmin.x) break;
    pmaxmax= point(pmaxmin.x, P[i+1].y);

    vector<point> S;
    S.push_back(pminmin);

    //lower hull
    i=1;
    if(!(pminmin ==pminmax))
    {
        while(i<n)
        {
            if(ccw(pminmin, pminmax, P[i])) {i++; continue;}
            j=S.size();
            if(j<2) S.push_back(P[i++]);
            else if (ccw(S[j-2],S[j-1],P[i])) S.push_back(P[i++]);
            else S.pop_back();
        }
    }

    if(!(pmaxmin == pmaxmax)) S.pb()


}
*/

////////////////////////      CLOSEST PAIR  BY LINE SWEEP ALGORITHM        ////////////////////////////////
////////////////////////      CLOSEST PAIR  BY LINE SWEEP ALGORITHM        ////////////////////////////////

struct sortBYCoordinate {
    bool operator() (const point &a, const point &b) const
    {	

        if(fabs( a.y - b.y )>EPS )
        return a.y < b.y;
        return a.x < b.x;
    }
};

double CPlinesweep(vector<point> P)
{
    sort(P.begin(),P.end());

    int l=P.size();
    if(l<2) return DBL_MAX;

    set<point, sortBYCoordinate> box;
    double smldis=dist(P[0],P[1]);

    box.insert(P[0]);
    int left=0;

    for(int i=1;i<l;i++)
    {
      //Maintain only points to the left of the current point whose distance is less than bestDist
      while( left<i && P[i].x - P[left].x >smldis )
      {
          box.erase(P[left]);
          left++;
      }

      //Consider only points within bestDist of the current point
      for( typeof(box.begin()) it =box.lower_bound( point(P[i].x -smldis, P[i].y - smldis )); it!= box.end() && (P[i].y + smldis) >=it->y; it++)
      {
            smldis= min(smldis, dist(*it,P[i]));
      }
      box.insert(P[i]);
    }

    return smldis;
}


///////////////////////////       CLOSEST PAIR BY DIVIDE AND CONQUER            /////////////////////////////
///////////////////////////       CLOSEST PAIR BY DIVIDE AND CONQUER            /////////////////////////////

double smldis=DBL_MAX;

bool compareY(point &a, point &b)
{
    if( fabs(a.y - b.y)> EPS) return a.y < b.y;
    return a.x < b.x;
}

double CPDivandConq(int l,int r, vector<point> &P)
{
    int s=(r-l+1);
    if(l>r ) return DBL_MAX;
    if(s<2) return DBL_MAX;
    //if(s<=3)-----------------no need
    int m=(l+r)/2;

    smldis=min(smldis,min( CPDivandConq(l,m,P), CPDivandConq(m+1,r,P ) ) );

    vector<point> strip;
    strip.clear();
    for(int i=m-1;i>=l;i--)
    {
        if( (P[m].x - P[i].x) < smldis )
           { if( fabs(P[m].y - P[i].y) < smldis ) strip.push_back(point(P[i].x,P[i].y)); }
        else break;
    }

    for(int i=m+1;i<=r;i++)
    {
        if( (P[i].x - P[m].x) < smldis )
            { if( fabs(P[m].y - P[i].y) < smldis ) strip.push_back(point(P[i].x,P[i].y)); }
        else break;
    }

    for(vector<point>::iterator it=strip.begin(); it!=strip.end(); it++)
    {
            smldis=min(smldis, dist(P[m], *it ) );
    }

    return smldis;
}


///////////////////////////////////////////////////     GCDISTANCE      //////////////////////////////////
double gcDistance(double pLat, double pLong,
double qLat, double qLong, double radius) {
pLat *= PI / 180; pLong *= PI / 180;     // convert degree to radian
qLat *= PI / 180; qLong *= PI / 180;
return radius * acos(cos(pLat)*cos(pLong)*cos(qLat)*cos(qLong) +
cos(pLat)*sin(pLong)*cos(qLat)*sin(qLong) +
sin(pLat)*sin(qLat)); }



//haversine formula
//for short distances

#define R 6378.0
double spherical_distance(double lat1,double lon1,double lat2,double lon2) {

       lat1 *= PI / 180; lon1 *= PI / 180;     // convert degree to radian
        lat2*= PI / 180; lon2 *= PI / 180;
       double dlon = lon2 - lon1;
       double dlat = lat2 - lat1;
       double a = pow((sin(dlat/2)),2) + cos(lat1) * cos(lat2) * pow(sin(dlon/2), 2);
       double c = 2 * atan2(sqrt(a), sqrt(1-a));
       double d = R * c;
       return d;
}

///////////////////////////////////////   LINE SEGMENT INTERSECTIONS         //////////////////////////////////
///////////////////////////////////////   LINE SEGMENT INTERSECTIONS        //////////////////////////////////


//CASE OF COINCIDENT ENDPOINTS TO INTERSECT DIDN'T TAKEN INTO ACCOUNT- 6 DIFF CASES CAN BE MADE------
//IN THAT CASE POINTS SHOULD BE SORTED ON THE BASIS OF TYPE ---0-2-1
struct sortBYCoordinate {
    bool operator() (const point &a, const point &b) const
    {
        //if(fabs( a.y - b.y ) > EPS )
        return a.y < b.y;
        //return a.x < b.x;
    }
};

struct event
{
    point p1,p2;
    int type;
    event() {};
    event(point p1,point p2, int type) : p1(p1), p2(p2),type(type) {};
};

bool compare(event a,event b)
{
    if(a.p1.x==b.p1.x)
    {
        if(a.type<=b.type)
            {if(a.p1.y < b.p2.y) return 1;
            else return 0;}
        else
            return 0;
    }
    return a.p1.x < b.p1.x;
}
int e;
event events[1000005];
set<point,sortBYCoordinate> s;

void ls_intersections()
{
    for(int i=0;i<e;i++)
    {
        event c=events[i];

        if(c.type==0)
        {
            s.insert(c.p1);
        }
        else if(c.type==2)
        {
            s.erase(c.p1);
        }
        else
        {
            for( set<point>::iterator it= s.lower_bound(point(-1,c.p1.y)); it!=s.end() && it->y<=c.p2.y; it++ )

                printf("%.0f %.0f\n",c.p1.x, it->y);
        }
    }
}

void Line_Segments_Inter()
{
    int n;
    scanf("%d",&n);

    point p1,p2,p3;

    while(n--)
    {
            //so that sorting will be in y coordinate format
            scanf("%lf%lf%lf%lf",&p1.x,&p1.y,&p2.x,&p2.y);
            //check for order of points
           /* if(p1.x>p2.x)
            {
                p3=p2;
                p2=p1;
                p1=p3;
            }

            if(p1.y>p2.y)
            {
                p3=p2;
                p2=p1;
                p1=p3;
            }

            */
            if(p1.x==p2.x)
            {
                events[e++]=event(p1,p2,1);
            }
            else
            {
                events[e++]=event(p1,p2,0);
                events[e++]=event(p2,p1,2);
            }
    }

    sort(events,events + e,compare);
    //for(int i=0;i<e;i++)
      //  cout<<events[i].p1.x<<" "<<events[i].p1.y<<" "<<events[i].p2.x<<" "<<events[i].p2.y<<endl;
    ls_intersections();
}

////////////////////////////////////  UNION RECTANGLES   //////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////   UNION RECTANGLES   ////////////////////////////////////////////////////////////////////////////////

///// can e improved using bst for in_set instead of boolarray


#define MAX 1000
struct event 
{
    int ind; // Index of rectangle in rects
    bool type; // Type of event: 0 = Lower-left ; 1 = Upper-right
    event() {};
    event(int ind, int type) : ind(ind), type(type) {};
};
struct point 
{
    int x, y;
};
point rects [MAX][2]; // Each rectangle consists of 2 points: [0] = lower-left ; [1] = upper-right
bool compare_x(event a, event b) { return rects[a.ind][a.type].x<rects[b.ind][b.type].x; }
bool compare_y(event a, event b) { return rects[a.ind][a.type].y<rects[b.ind][b.type].y; }
int union_area(event events_v[],event events_h[],int n,int e)
{
       //n is the number of rectangles, e=2*n , e is the number of points (each rectangle has two points as described in declaration of rects)
        bool in_set[MAX]={0};int area=0;
        sort(events_v, events_v+e, compare_x);  //Pre-sort of vertical edges
        sort(events_h, events_h+e, compare_y); // Pre-sort set of horizontal edges
        in_set[events_v[0].ind] = 1;
        for (int i=1;i<e;++i) 
        { // Vertical sweep line
                event c = events_v[i];
                int cnt = 0; // Counter to indicate how many rectangles are currently overlapping
                // Delta_x: Distance between current sweep line and previous sweep line
                int delta_x = rects[c.ind][c.type].x - rects[events_v[i-1].ind][events_v[i-1].type].x;
                int begin_y;
                if (delta_x==0){
                        in_set[c.ind] = (c.type==0);
                        continue;
                }
                for (int j=0;j<e;++j)
                        if (in_set[events_h[j].ind]==1)                 //Horizontal sweep line for active rectangle
                        {
                                if (events_h[j].type==0)                //If it is a bottom edge of rectangle
                                {
                                        if (cnt==0) begin_y = rects[events_h[j].ind][0].y; // Block starts
                                        ++cnt;                          //incrementing number of overlapping rectangles
                                }
                                else                                    //If it is a top edge
                                {
                                        --cnt;                          //the rectangle is no more overlapping, so remove it
                                        if (cnt==0)                     //Block ends
                                        {
                                                int delta_y = (rects[events_h[j].ind][1].y-begin_y);//length of the vertical sweep line cut by rectangles
                                                area+=delta_x * delta_y;
                                        }
                                }
                        }
                in_set[c.ind] = (c.type==0);//If it is a left edge, the rectangle is in the active set else not
        }
    return area;
}



/////////////////////           SOLUTION STARTS HERE          ////////////////////////////////////////////
/////////////////////           SOLUTION STARTS HERE          ////////////////////////////////////////////

nt main()
{
    #ifdef GEOMETRY
		freopen("input.txt", "r", stdin);
		freopen("output.txt", "w", stdout);
	#endif

    int t=1,c;
    cin>>c;

    while(c--)
    {

    }
    return 0;
}
