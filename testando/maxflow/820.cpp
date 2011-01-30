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

#define DEBUG(x...) printf(x)
#define all(v) (v).begin(),(v).end()
#define rall(v) (v).rbegin(),(v).rend()
#define _foreach(it, b, e) for(__typeof__(b) it = (b); it != (e); ++it)
#define foreach(x...) _foreach(x)

typedef long long int huge;

const int inf = 0x3f3f3f3f;
const huge hugeinf = 0x3f3f3f3f3f3f3f3fll; // sao dois L's!!!
const double eps = 1e-9;

////////////////////////////////////////////////////////////////////////////////
// EDMONDS-KARP MAXFLOW - O(max(MaxFlow*V, VE^2))
// Muda o m para ser grafo residual.
////////////////////////////////////////////////////////////////////////////////
const int nmax = 200;

struct EKMaxFlow
{
  int (*m)[nmax];
  vector<int> *l;
  int marked[nmax], prev[nmax];
  int n, s, t;
  void init(int mtx[nmax][nmax], vector<int> *list, int nn, int ns, int nt)
  { m = mtx; l = list; n = nn; s = ns; t = nt; }
  int augment()
  {
    int ans = 0;
    while(true)
      {
	memset(marked, 0, sizeof(marked));
	memset(prev, 0xff, sizeof(prev));
	queue<int> q;
	q.push(s);
	marked[s] = true;
	while(!q.empty())
	  {
	    int v = q.front(); q.pop();
	    if(v == t) break;
	    foreach(it, all(l[v])) // se for mudar pra mtx, eh aqui
	      if(!marked[*it] && m[v][*it] > 0)
		{
		  marked[*it] = true;
		  q.push(*it);
		  prev[*it] = v;
		}
	  }
	if(prev[t] == -1) break;
	int cap = inf;
	for(int i=t; i!=s; i=prev[i])
	  cap = min(cap, m[prev[i]][i]);
	for(int i=t; i!=s; i=prev[i])
	  {
	    m[i][prev[i]] += cap;
	    m[prev[i]][i] -= cap;
	  }
	ans += cap;
      }
    return ans;
  }
} ekmf;

int m[nmax][nmax];
vector<int> l[nmax];

int main()
{
  int n, test = 1;
  while(scanf(" %d", &n)==1 && n!=0)
    {
      memset(m, 0, sizeof(m));
      for(int i=0; i<nmax; ++i)
	l[i].clear();
      int s, t, e;
      scanf(" %d %d %d", &s, &t, &e);
      --s; --t;
      for(int i=0; i<e; ++i)
	{
	  int a, b, c;
	  scanf(" %d %d %d", &a, &b, &c);
	  --a, --b;
	  if(m[a][b] == 0)
	    {
	      l[a].push_back(b);
	      l[b].push_back(a);
	    }
	  m[a][b] = (m[b][a] += c);
	}
      
      printf("Network %d\n", test++);
      ekmf.init(m, l, n, s, t);
      printf("The bandwidth is %d.\n\n", ekmf.augment());
    }
  return 0;
}

