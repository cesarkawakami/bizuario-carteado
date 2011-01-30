// testando com problema 10080

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

#define DEBUG(x...) //printf(x)
#define all(v) (v).begin(),(v).end()
#define rall(v) (v).rbegin(),(v).rend()
#define _foreach(it, b, e) for(typeof(b) it = (b); it != (e); ++it)
#define foreach(x...) _foreach(x)

typedef long long int huge;

const int inf = 0x3f3f3f3f;
const huge hugeinf = 0x3f3f3f3f3f3f3f3fll;
const double eps = 1e-9;

struct point
{
  double x, y;
  point(double x=0, double y=0) : x(x), y(y) {}
  inline point operator-(const point &p) {return point(x-p.x, y-p.y);}
  inline point operator+(const point &p) {return point(x+p.x, y+p.y);}
  inline point operator*(const double &c) {return point(x*c, y*c);}
  inline point operator/(const double &c) {return point(x/c, y/c);}
  
  inline bool operator==(const point &p) {return fabs(x-p.x)<eps && fabs(y-p.y)<eps;}
  inline bool operator!=(const point &p) {return !(*this==p);}
};
inline double dot(const point &a, const point &b) {return a.x*b.x+a.y*b.y;}
inline double norm(const point &p) {return sqrt(dot(p,p));}

////////////////////////////////////////////////////////////////////////////////
// UNWEIGHTED BIPARTITE MATCHING
////////////////////////////////////////////////////////////////////////////////
// Matching resultante: para matching[i]!=-1, (i, matching[i]).
// tab[i][j] -> de i para j. m linhas e n colunas.
const int nmax = 200;

struct UBMatching
{
  int m, n, size;
  int matching[nmax], seen[nmax];
  int (*tab)[nmax];
  bool augment(int v)
  {
    for(int i=0;i<m;i++)
      if(tab[i][v] && !seen[i])
	{
	  seen[i]=1;
	  if(matching[i]<0 || augment(matching[i]))
 	    {
	      matching[i]=v;
	      return 1;
	    }
	}
    return 0;
  }
  void init(int nm, int nn, int mtx[nmax][nmax])
  {
    m = nm; n = nn; size = 0; tab = mtx;
    memset(matching, -1, sizeof(matching));
  }
  void match()
  {
    for(int i=0;i<n;i++)
      {
	memset(seen, 0, sizeof(seen));
	size += augment(i);
      }
  }
} ubmatching;

int s, v;

int tab[nmax][nmax];

int main()
{
  int n, m;
  while( scanf(" %d %d %d %d", &n, &m, &s, &v)==4 )
    {
      double x, y;
      vector<point> gophers, holes;
      for(int i=0; i<n; ++i)
	{
	  scanf(" %lf %lf", &x, &y);
	  gophers.push_back(point(x, y));
	}
      for(int i=0; i<m; ++i)
	{
	  scanf(" %lf %lf", &x, &y);
	  holes.push_back(point(x,y));
	}
      memset(tab, 0, sizeof(tab));
      for(int i=0; i<n; ++i)
	for(int j=0; j<m; ++j)
	  if( norm(gophers[i] - holes[j]) < s*v + eps)
	    tab[i][j] = 1;
	  else
	    tab[i][j] = 0;
      ubmatching.init(n, m, tab);
      ubmatching.match();
      printf("%d\n", n - ubmatching.size);
    }
  return 0;
}

