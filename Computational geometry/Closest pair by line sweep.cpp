////////////////////////      CLOSEST PAIR  BY LINE SWEEP ALGORITHM        ////////////////////////////////
////////////////////////      CLOSEST PAIR  BY LINE SWEEP ALGORITHM        ////////////////////////////////


//vector<point> p contains all the points which are under consideration

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
