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
#include <stack>
#include <functional>
#include <sstream>
#include <iostream>
#include <ctime>
#include <algorithm>
using namespace std;

#define DEBUG(x...) printf(x)
#define all(v) (v).begin(),(v).end()
#define rall(v) (v).rbegin(),(v).rend()
//#define _foreach(it, b, e) for(typeof(b) it = (b); it != (e); ++it)
//#define foreach(x...) _foreach(x)
#define foreach(it, b, e) for(typeof(b) it = (b); it != (e); ++it)
typedef long long int huge;

const int inf = 0x3f3f3f3f;
const huge hugeinf = 0x3f3f3f3f3f3f3f3fll;
const double eps = 1e-9;

////////////////////////////////////////////////////////////////////////////////
// UNWEIGHTED BIPARTITE MATCHING
////////////////////////////////////////////////////////////////////////////////
// Matching resultante: para matching[i]!=-1, (i, matching[i]).
struct UBMatching
{
  int m, n, size;
  vector<int> matching, seen;
  const vector<vector<int> > &tab;
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
  UBMatching(const vector<vector<int> > &table) // table[i][j] = i pra j
    : m(table.size()), n(table[0].size()), size(0),
      matching(m,-1), seen(m,0), tab(table)
  {
    for(int i=0;i<n;i++)
      {
	fill(all(seen),0);
	size += augment(i);
      }
  }
};

struct slide
{
  int xmin, xmax, ymin, ymax;
};

struct point
{
  int x, y;
};

int main()
{
  int n, test = 1;
  while( scanf(" %d", &n) == 1 && n!=0 )
    {
      printf("Heap %d\n", test);
      vector<slide> slides;
      vector<point> numbers;
      for(int i=0; i<n; ++i)
	{
	  int a, b, c, d;
	  scanf(" %d %d %d %d", &a, &b, &c, &d);
	  slides.push_back( (slide){a, b, c, d} );
	}
      for(int i=0; i<n; ++i)
	{
	  int x, y;
	  scanf(" %d %d", &x, &y);
	  numbers.push_back( (point){x,y} );
	}
      vector<vector<int> > tab(n, vector<int>(n, 0));
      for(int i=0; i<n; ++i)
	for(int j=0; j<n; ++j)
	  if( slides[i].xmin < numbers[j].x && numbers[j].x < slides[i].xmax 
	      && slides[i].ymin < numbers[j].y && numbers[j].y < slides[i].ymax )
	    tab[i][j] = 1;
      UBMatching m(tab);
      vector<pair<char, int> > uniquepairs;
      for(int i=0; i<n; ++i)
	{
	  // removing edge
	  tab[i][m.matching[i]] = 0;
	  UBMatching m2(tab);
	  if( m2.size < n )
	    uniquepairs.push_back(make_pair('A'+i,1+m.matching[i]));
	  tab[i][m.matching[i]] = 1;
	}
      if( uniquepairs.size() == 0 )
	printf("none");
      else
	{
	  char *spc = "";
	  foreach(it, all(uniquepairs))
	    printf("%s(%c,%d)", spc, it->first, it->second), spc=" ";
	}
      printf("\n\n");
      ++test;
    }

  return 0;
}
