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
