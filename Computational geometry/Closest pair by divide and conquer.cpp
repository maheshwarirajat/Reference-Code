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
