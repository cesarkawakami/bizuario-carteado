struct point
{
  huge x, y;
  point(huge x=0, huge y=0) : x(x), y(y) {}
  inline point operator+(const point &p) const {return point(x+p.x, y+p.y);}
  inline point operator-(const point &p) const {return point(x-p.x, y-p.y);}
  inline bool operator==(const point &p) const {return x==p.x&&y==p.y;}
};

struct line
{
  point p1, p2;
  line(point p1, point p2) : p1(p1), p2(p2) {}
};

inline huge dot(const point &a, const point &b) {return a.x*b.x + a.y*b.y;}
inline huge cross(const point &a, const point &b) {return a.x*b.y - a.y*b.x;}

huge ccw(const point &a, const point &b, const point &c)
{
  huge k;
  if( (k=cross(b-a,c-a)) > 0 ) return 1;
  if( k < 0 ) return -1;
  if( dot(c-a,b-a) < 0 ) return -1;
  if( dot(a-b,c-b) < 0 ) return 1;
  return 0;
}

// Falha se line a ou line b for pontual, ver se pode ocorrer!
bool intersect(const line &a, const line &b)
{
  return ( ccw(a.p1, a.p2, b.p1) * ccw(a.p1, a.p2, b.p2) <= 0 ) &&
    ( ccw(b.p1, b.p2, a.p1) * ccw(b.p1, b.p2, a.p2) <= 0 );
}

inline bool between(point a, point b, point p)
{
  return cross(b-a,p-a)==0 && dot(p-a,p-b) <= 0;
}
