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

