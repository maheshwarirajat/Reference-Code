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
