#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
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

////////////////////////////////////////////////////////////////////////////////
// UNWEIGHTED BIPARTITE MATCHING
////////////////////////////////////////////////////////////////////////////////
// Matching resultante: para matching[i]!=-1, (i, matching[i]).
// tab[i][j] -> de i para j. m linhas e n colunas. (mtx no sentido matematico)
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

vector<string> hero, com;
int tab[nmax][nmax];
int tmp[nmax], tmp2[nmax];

int main()
{
  int n, e;
  cin >> n >> e;
  hero.resize(n);
  com.resize(n);
  for(int i=0; i<n; ++i)
    cin >> hero[i];
  for(int i=0; i<n; ++i)
    cin >> com[i];
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      tab[i][j] = 1;
  for(int k=0; k<e; ++k)
    {
      int m;
      cin >> m;
      vector<int> h, c;
      while(m--)
	{
	  string s;
	  cin >> s;
	  vector<string>::iterator it;
	  if((it = find(all(hero), s)) != hero.end())
	    h.push_back(it - hero.begin());
	  if((it = find(all(com), s)) != com.end())
	    c.push_back(it - com.begin());
	}
      foreach(it, all(h))
	foreach(jt, all(c))
	  tab[*it][*jt] = 0;
    }

//   for(int i=0; i<n; ++i)
//     {
//       for(int j=0; j<n; ++j)
// 	cout << tab[i][j] << " ";
//       cout << endl;
//     }
  
  ubmatching.init(n, n, tab);
  ubmatching.match();

  if(ubmatching.size < n)
    {
      cout << "IMPOSSIVEL" << endl;
      exit(0);
    }
  
  for(int i=0; i<n; ++i)
    {
      cout << hero[i] << ":";
      for(int j=0; j<n; ++j)
	{
	  if(ubmatching.matching[i] == j)
	    cout << " " << com[j];
	  else if(tab[i][j])
	    {
	      // backing up current state
	      memcpy(tmp, ubmatching.matching, sizeof(ubmatching.matching));
	      memcpy(tmp2, tab[i], sizeof(tab[i]));

	      // forcing i-j matching
	      int matchj = find(ubmatching.matching, ubmatching.matching+n, j) - ubmatching.matching;
	      int matchi = ubmatching.matching[i];
	      ubmatching.matching[i] = j;
	      ubmatching.matching[matchj] = -1;
	      memset(tab[i], 0, sizeof(tab[i]));
	      tab[i][j] = 1;

	      // trying to match matchi
	      memset(ubmatching.seen, 0, sizeof(ubmatching.seen));
	      if(ubmatching.augment(matchi))
		cout << " " << com[j];

	      // reversing changes
	      memcpy(ubmatching.matching, tmp, sizeof(tmp));
	      memcpy(tab[i], tmp2, sizeof(tmp2));
	    }
	}
      cout << endl;
    }
  return 0;
}
