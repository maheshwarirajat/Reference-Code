
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
