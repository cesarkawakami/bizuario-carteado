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
// Na lista de adj, tem de colocar aresta nos dois sentidos, mesmo em directed
////////////////////////////////////////////////////////////////////////////////
const int nmax = 5050;

huge m[nmax][nmax];
vector<int> l[nmax];

struct EKMaxFlow
{
  int marked[nmax], prev[nmax];
  int n, s, t;
  void init(int nn, int ns, int nt)
  { n = nn; s = ns; t = nt; }
  huge augment()
  {
    huge ans = 0;
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
	huge cap = inf;
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

set<int> edges[nmax];

int main()
{
  int V, E;
  scanf(" %d %d", &V, &E);
  for(int i=0; i<E; ++i) { int a, b, c; scanf(" %d %d %d", &a, &b, &c); if(a!=b) { --a, --b; edges[a].insert(b); edges[b].insert(a); m[a][b] = m[b][a] += c; } }
  for(int i=0; i<V; ++i) l[i] = vector<int>(all(edges[i]));
  ekmf.init(V, 0, V-1);
  printf("%lld\n", ekmf.augment());
}
