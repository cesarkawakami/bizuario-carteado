//////////////////////////////////////////////////////////////////////////////
// GEOMETRIA 3D
//////////////////////////////////////////////////////////////////////////////

////
// PAPERWEIGHT: convex hull face, distances
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
using namespace std;

const long double eps = 1e-9;

struct point
{
  long double x, y, z;
};

point cross(point a, point b)
{
  return (point){
    a.y*b.z - a.z*b.y,
      a.z*b.x - a.x*b.z,
      a.x*b.y - a.y*b.x};
}

long double dot(point a, point b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

long double mixedproduct(point a, point b, point c)
{
  return dot(a, cross(b, c));
}

point operator-(const point &a, const point &b)
{
  return (point){a.x - b.x, a.y - b.y, a.z - b.z};
}

point operator+(const point &a, const point &b)
{
  return (point){a.x + b.x, a.y + b.y, a.z + b.z};
}

point operator/(const point &a, long double k)
{
  return (point){a.x / k, a.y / k, a.z / k};
}

point operator*(const point &a, long double k)
{
  return (point){a.x * k, a.y * k, a.z * k};
}

long double norm(point p) { return sqrt(dot(p, p)); }

bool ishullface(const vector<point> &poly, int a, int b, int c)
{
  long double mp[6];
  int pos=0,neg=0;
  for (int i=0; i<5; ++i)
    mp[i]=mixedproduct(poly[b] - poly[a], poly[c] - poly[a], poly[i] - poly[a]);
  for (int i=0;i<5;i++)
    {
      if (mp[i]<-eps)
	neg++;
      else if (mp[i]>eps)
	pos++;
    }
  if (neg*pos == 0)
    {
      //printf("ishullface(%d %d %d) = true\n", a, b, c);
      return true;
    }
  else
    {
      //printf("ishullface(%d %d %d) = false\n", a, b, c);
      return false;
    }
}

long double volume(point a, point b, point c, point d)
{
  return fabs(mixedproduct(b - a, c - a, d - a));
}

point cmass(point a, point b, point c, point d)
{
  return (a + b + c + d) / 4;
}

point totalcmass(const vector<point> &poly)
{
  point CM1 = cmass(poly[0], poly[1], poly[2], poly[3]);
  long double VM1 = volume(poly[0], poly[1], poly[2], poly[3]);
  point CM2 = cmass(poly[0], poly[1], poly[2], poly[4]);
  long double VM2 = volume(poly[0], poly[1], poly[2], poly[4]);
  return (CM1 * VM1 + CM2 * VM2) / (VM1 + VM2);
}

bool isinside(point p, point a, point b, point c)
{
  if (fabs(norm(cross(p - a, b - a)) + norm(cross(p - b, c - b)) 
	+ norm(cross(p - c, a - c)) - norm(cross(b - a, c - a))) < eps)
    return true;
  else
    return false;
}

point project(point p, point a, point b, point c)
{
  point v = cross(b - a, c - a);
  return p - v * dot(p-a, v) / dot(v, v);
}

long double pointlinedist(point p, point a, point b)
{
  return norm(cross(p - a, p - b)) / norm(b - a);
}

bool isreallyinside(point p, point a, point b, point c)
{
  if (pointlinedist(p, a, b) > .2 && pointlinedist(p, b, c) > .2 
      && pointlinedist(p, a, c) > .2)
    return true;
  else
    return false;
}

bool isreallyinside(point p, point a, point b, point c, point d)
{
  if (pointlinedist(p, a, c) > .2 && pointlinedist(p, b, c) > .2 
      && pointlinedist(p, a, d) > .2 && pointlinedist(p, b, d) > .2)
    return true;
  else
    return false;
}

bool isstable(const vector<point> &poly, int a, int b, int c)
{
  point tcmass = totalcmass(poly);
  point P = project(tcmass, poly[a], poly[b], poly[c]);
  if (isinside(P, poly[a], poly[b], poly[c]) 
      && isreallyinside(P, poly[a], poly[b], poly[c]))
    return true;
  else
    return false;
}

/// novo codigo
bool isstable(const vector<point> &poly, int a, int b, int c, int d)
{
  int other = 1 + 2 + 3 + 4 - a - b - c - d;
  point tcmass = totalcmass(poly);
  point P = project(tcmass, poly[a], poly[b], poly[c]);
  if ((isinside(P, poly[a], poly[c], poly[d]) 
       || isinside(P, poly[b], poly[c], poly[d]))
      && isreallyinside(P, poly[a], poly[b], poly[c], poly[d]))
    return true;
  else
    return false;
}

long double pointplanedist(point p, point a, point b, point c)
{
  return norm(p - project(p, a, b, c));
}

int main()
{
  vector<point> points;
  point F;
  
  for (int ncase = 1;;++ncase)
    {
      points.clear();
      for (int i=0; i<5; ++i)
	{
	  double x, y, z;
	  if (scanf(" %lf %lf %lf", &x, &y, &z) != 3)
	    return 0;
	  points.push_back((point){x, y, z});
	}

      scanf(" %Lf %Lf %Lf", &F.x, &F.y, &F.z);

      long double mindist = 1e100, maxdist = -1;

      for (int i=0; i<5; ++i)
	for (int j=0; j<i; ++j)		       
	  for (int k=0; k<j; ++k)
	    if (ishullface(points, i, j, k) && isstable(points, i, j, k))
	      {
		//printf("(%d %d %d)\n", i, j, k);
		mindist = min(mindist, 
		    pointplanedist(F, points[i], points[j], points[k]));
		maxdist = max(maxdist, 
		    pointplanedist(F, points[i], points[j], points[k]));
	      }

      for (int i=0; i<5; ++i)
	for (int j=0; j<i; ++j)
	  for (int k=0; k<j; ++k)
	    for (int l=0; l<k; ++l)
	      if (volume(points[i], points[j], points[k], points[l]) < eps 
	          && isstable(points, i, j, k, l))
	      {
		//printf("(%d %d %d)\n", i, j, k);
		mindist = min(mindist, 
		    pointplanedist(F, points[i], points[j], points[k]));
		maxdist = max(maxdist, 
		    pointplanedist(F, points[i], points[j], points[k]));
	      }

      printf("Case %d: %.5Lf %.5Lf\n", ncase, mindist, maxdist);
    }

  return 0;
}

////
// GUERRA
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <deque>
#include <queue>
#include <stack>
#include <functional>
#include <sstream>
#include <iostream>
#include <ctime>
#include <algorithm>
using namespace std;

#define all(v) (v).begin(), (v).end()
#define _foreach(it, b, e) for (__typeof__(b) it = (b); it != (e); ++it)
#define foreach(x...) _foreach(x)

typedef long long huge;

const int inf = 0x3f3f3f3f;
const huge hugeinf = 0x3f3f3f3f3f3f3f3fll;
const double eps = 1e-9;
struct point
{
  long double x, y, z;
  point(long double x=0, long double y=0, long double z=0) :x(x), y(y), z(z) {}
  inline point operator+(const point &p) const {return point(x+p.x, y+p.y, z+p.z);}
  inline point operator-(const point &p) const {return point(x-p.x, y-p.y, z-p.z);}
  inline point operator*(long double p) const {return point(x*p, y*p, z*p); }
  inline bool operator==(const point &p) const 
  {return fabsl(x - p.x) < eps && fabsl(y - p.y) < eps && fabsl(z - p.z) < eps; }
};
inline long double dot(const point &a, const point &b)
{
  return a.x*b.x+a.y*b.y+a.z*b.z;
}
inline point cross(const point &a, const point &b)
{
  return point(a.y*b.z-a.z*b.y,
               a.z*b.x-a.x*b.z,
               a.x*b.y-a.y*b.x);
} 

point proj(point R, point P) // projeta P sobre R
{
  return R * (dot(P, R) / dot(R, R));
}

long double abs(point P)
{
  return sqrt(dot(P, P));
}

bool between(point A, point B, point P) // inclusive, assume colinear
{
  return dot(P - A, P - B) < eps;
}

long double linedist(point A, point B, point P)
{
  point Pproj = proj(B - A, P - A) + A;
  if (between(A, B, Pproj))
    return abs(P - Pproj);
  else
    return min(abs(P - A), abs(P - B));
}

long double pointtri(point A, point B, point C, point P)
{
  point X = B - A;
  point Y = C - A;
  point Pori = P;
  P = P - A;
  point PP = P - proj(cross(X, Y), P);
  point PPP = PP + A;
  point R1 = cross(B - A, PPP - A);
  point R2 = cross(C - B, PPP - B);
  point R3 = cross(A - C, PPP - C);
  if (dot(R1, R2) > -eps && dot(R2, R3) > -eps && dot(R1, R3) > -eps)
    return abs(Pori - PPP);
  else
    return min(linedist(A, B, Pori), min(linedist(B, C, Pori), 
                                         linedist(C, A, Pori)));
}

int ccw(point A, point B, point C, point Ref)
{
  double k = dot(cross(B - A, C - A), Ref);
  if (k > eps) return 1;
  if (k < -eps) return -1;
  if (dot(C - A, B - A) < -eps) return -1;
  if (dot(A - B, C - B) < -eps) return 1;
  return 0;
}

bool intersect(point A, point B, point C, point D)
{
  point Ref;
  point X[4] = {A, B, C, D};
  for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
      for (int k=0; k<4; ++k)
        {
          point RC = cross(X[j] - X[i], X[k] - X[i]);
          if (abs(RC) > eps)
            Ref = RC;
        }
  assert(abs(Ref) > eps);
  return (ccw(A, B, C, Ref) * ccw(A, B, D, Ref) <= 0) 
         && (ccw(C, D, A, Ref) * ccw(C, D, B, Ref) <= 0);
}

long double segseg(point A, point B, point C, point D)
{
  point E = proj(cross(D - C, B - A), A - D);
  point Cl = C + E;
  point Dl = D + E;
  if (abs(cross(D - C, B - A)) > eps && intersect(A, B, Cl, Dl))
    return abs(E);
  else
    return min(min(abs(A - C), abs(A - D)), min(abs(B - C), abs(B - D)));
}

void read(point &P)
{
  long double x, y, z;
  scanf(" %Lf %Lf %Lf", &x, &y, &z);
  P = point(x, y, z);
}

int main()
{
  int ntests;
  scanf(" %d", &ntests);
  while (ntests--)
    {
      point T[2][4];
      for (int i=0; i<2; ++i)
        for (int j=0; j<4; ++j)
          read(T[i][j]);
      long double bestdist = 0x3f3f3f3f3f3f3f3fll;
      for (int i=0; i<2; ++i)
        {
          for (int j=0; j<4; ++j)
            {                                   
              point P = T[0][j];
              bestdist = min(bestdist,
                             pointtri(T[1][0], T[1][1], T[1][2], P));
              bestdist = min(bestdist,
                             pointtri(T[1][0], T[1][1], T[1][3], P));
              bestdist = min(bestdist,
                             pointtri(T[1][0], T[1][2], T[1][3], P));
              bestdist = min(bestdist,
                             pointtri(T[1][1], T[1][2], T[1][3], P));
            }
          for (int j=0; j<4; ++j)
            for (int l=j+1; l<4; ++l)
              for (int m=0; m<4; ++m)
                for (int n=m+1; n<4; ++n)
                  bestdist = min(bestdist, 
                                 segseg(T[0][j], T[0][l], T[1][m], T[1][n]));
          for (int j=0; j<4; ++j)
            swap(T[0][j], T[1][j]);
        }
      printf("%.2Lf\n", bestdist);
    }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// AHO-CORASICK
//////////////////////////////////////////////////////////////////////////////
#define maxc 1000000
#define maxn 10000
struct link {
  int v;
  link *n;
} out[maxn], *l;
struct node {
  node *s, *b, *f;
  link *o;
  char i;
} rt[maxc], *p, *q, *r, *s;
char buf[maxc+maxn], mtx[maxn], *str;
int word[maxn], nn, nm, n;
void push_str(int i)
{
  for(str=buf+word[i], p=rt; *str; ++str) {
    for(q=p->s; q; q=q->b)
      if (q->i==*str)
	break;
    if (q) p=q;
    else {
      q=p->s;
      p=(p->s=rt+nn++);
      p->i=*str; p->b=q; p->s=0; p->o=0; p->f=rt;
    }
  }
  (p->o=out+nm++)->v=i;
}
void fix_tree()
{
  queue<node*> td;
  for(r=rt->s; r; r=r->b)
    td.push(r);
  while(!td.empty()) {
    p=td.front();td.pop();
    for(r=p->s; r; r=r->b)
      for(q=p, td.push(r); q!=rt;)
	for(q=q->f, s=q->s; s; s=s->b)
	  if (s->i==r->i) { r->f=s; q=rt; break; }
    if (p->o) p->o->n=p->f->o;
    else p->o=p->f->o;
  }
}
void match_str(int i)
{
  for(str=buf+word[i], p=rt; *str; ++str) {
    for(q=p->s; p!=rt||q; q=q->b) {
      for(; !q; p=p->f, q=p->s);
      if (q->i==*str) {	p=q; break; }
    }
    for(l=p->o; l&&!mtx[l->v]; l=l->n)
      mtx[l->v]=1;
  }
}

/////////////////////////////////////////////////
// Colorfull Spanning Tree
/////////////////////////////////////////////////
const int maxn = 210;
list<int> mtx[maxn][maxn];
list<pair<int,int> > adj[maxn];
bool expande[maxn];
pair<int,int> aresta[maxn];
bool e_marcada[maxn];
bool v_marcado[maxn];
bool temp[maxn];
bool pode_adicionar;
int n, k;
void adiciona_aresta(int x)
{
  v_marcado[aresta[x].second]=true;
  adj[aresta[x].first].push_back(make_pair(aresta[x].second,x));
  adj[aresta[x].second].push_back(make_pair(aresta[x].first,x));
}
void remove_aresta(int a, int b)
{
  foreach(it, all(adj[a]))
    if (it->first==b)
      {
	adj[a].erase(it);
	break;
      }
  foreach(it, all(adj[b]))
    if (it->first==a)
      {
	adj[b].erase(it);
	break;
      }
}
bool dfs(int v, int root)
{
  temp[v]=true;
  if (pode_adicionar)
    foreach(it, all(mtx[v][root]))
      if (!e_marcada[*it])
	{
	  e_marcada[*it]=true;
	  adj[root].push_back(make_pair(v,*it));
	  adj[v].push_back(make_pair(root,*it));
	  temp[v]=false;
	  return true;
	}
  foreach (it, all(adj[v]))
    if (!temp[it->first])
      {
	if (pode_adicionar||!expande[it->second])
	  {
	    if (dfs(it->first,root))
	      {
		temp[v]=false;
		return true;
	      }
	  }
	else
	  {
	    pode_adicionar=true;
	    if (dfs(it->first,root))
	      {
		adiciona_aresta(it->second);
		remove_aresta(v,it->first);
		temp[v]=false;
		pode_adicionar=false;
		return true;
	      }
	    pode_adicionar=false;
	  }
      }
  temp[v]=false;
  return false;
}
bool augment()
{
  for(int i=0; i<k; ++i)
    expande[i]=false;
  for(int i=0; i<n; ++i)
    if (v_marcado[i])
      for(int j=0; j<n; ++j)
	if (!v_marcado[j])
	  foreach (it, all(mtx[i][j]))
	    {
	      if (!e_marcada[*it])
		{
		  e_marcada[*it]=true;
		  adj[i].push_back(make_pair(j,*it));
		  adj[j].push_back(make_pair(i,*it));
		  v_marcado[j]=true;
		  return true;
		}
	      if (!expande[*it])
		{
		  expande[*it]=true;
		  aresta[*it]=make_pair(i,j);
		}
	    }
  for(int i=0; i<n; ++i)
    if (v_marcado[i])
      if (dfs(i,i))
	return true;
  return false;
}
void print_selected()
{
  for(int i=0; i<n; ++i)
    {
      printf("%d:", i+1);
      foreach(it, all(adj[i]))
	printf(" (%d,%d)", it->first+1, it->second+1);
      printf("\n");
    }
  printf("-\n");
}
int main()
{
  int a, b, c, m;
  bool ok;
  srand(time(NULL));
  for(int i=0; i<maxn; ++i)
    temp[i]=false;
  pode_adicionar=false;
  for(int teste=1;;++teste)
    {
      if (scanf(" %d %d %d", &n, &m, &k)!=3)
	return 0;
      for(int i=0; i<n; ++i)
	{
	  for(int j=0; j<n; ++j)
	    mtx[i][j].clear();
	  v_marcado[i]=false;
	  adj[i].clear();
	}
      for(int i=0; i<k; ++i)
	e_marcada[i]=false;
      for(int i=0; i<m; ++i)
	{
	  scanf(" %d %d %d", &a, &b, &c);
	  --a;--b;--c;
	  if (a!=b)
	    {
	      mtx[a][b].push_back(c);
	      mtx[b][a].push_back(c);
	    }
	}
      printf("Instancia %d\n", teste);
      v_marcado[0]=true;
      ok=1;
      for(int i=1; ok&&i<n; ++i)
	{
	  if (!augment())
	    ok=0;
	  //print_selected();
	}
      if (ok)
	printf("sim\n\n");
      else
	printf("nao\n\n");
    }
}

//////////////////////////////////////////////
// Multiplicação de Nímeros O(n^2)
//////////////////////////////////////////////
const int maxn=200;
const int lim=16;
int tab[maxn][maxn];
bool used[maxn];
void multiplica(int a, int b)
{
  for(int i=0; i<maxn; ++i)
    used[i]=0;
  for(int i=0; i<a; ++i)
    for(int j=0; j<b; ++j)
      used[tab[i][b]^tab[a][j]^tab[i][j]]=true;
  for(int i=0; ;++i)
    if (!used[i]) {
      tab[a][b]=i;
      return;
    }
}

////////////////////////////////
/// Closest - Pair problem 2D
////////////////
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<cmath>
#include<vector>
#include<set>
#include<map>
#include<list>
#include<deque>
#include<queue>
#include<stack>
#include<functional>
#include<sstream>
#include<iostream>
#include<ctime>
#include<algorithm>
using namespace std;

#define DEBUG(x...) printf(x)
#define all(v) (v).begin(),(v).end()
#define rall(v) (v).begin(),(v).rend()
#define _foreach(it,b,e) for(__typeof__(b) it=(b); it!=(e);++it)
#define foreach(x...) _foreach(x)

typedef long long int huge;

const int inf = 0x3f3f3f3f;
const huge hugeinf = 0x3f3f3f3f3f3f3f3fll;//dois L's
const double eps = 1e-9;

int n;
pair<double,double> point[15000];
pair<double,double> tmp[15]; 
bool mcomp(const pair<double,double> &x,const pair<double,double> &y)
{
    return x.second<y.second;   
}
double go(const int a,const int b)
{
  if(b<a+3)
    {
      double ret = hugeinf;
      for(int i=a;i<b;i++)
	for(int j=i+1;j<b;j++)
	  ret=min(ret,hypot(point[i].first-point[j].first,point[i].second-point[j].second));
      return ret;
    }
  double da = go(a,(a+b)/2),db = go((a+b)/2,b);
  double ret = min(da,db);
      
    vector<pair<double,double> > left;
    for(int i=(a+b)/2;i<b;i++)
        if(point[i].first-point[(a+b)/2].first< ret)
            left.push_back(point[i]);
    sort(all(left),mcomp);
    for(int i=a;i<(a+b)/2;i++)
    {
        int x = lower_bound(all(left),make_pair(0,point[i].second-ret-0.01),mcomp)-left.begin();
        int y = upper_bound(all(left),make_pair(0,point[i].second+ret+0.01),mcomp)-left.begin();
        for(int j=x;j<y;j++)
        {
            ret=min(ret,hypot(point[i].first-left[j].first,point[i].second-left[j].second));
        }
    }
  return ret;
}
int main ()
{
  while(scanf("%d",&n)==1)
    {
      if(!n)break;
      for(int i=0;i<n;i++)
	scanf("%lf %lf",&point[i].first,&point[i].second);
      sort(point,point+n);
      double x = go(0,n);

      /* double y=hugeinf;
      for(int i=0;i<n;i++)
	for(int j=i+1;j<n;j++)
	  y=min(y,hypot(point[i].first-point[j].first,point[i].second-point[j].second));
	  printf("%lf ",y);*/
    
      if(x>=10000)
	printf("INFINITY\n");
      else
	printf("%.4lf\n",x);
    }
  return 0;
}
