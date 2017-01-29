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
