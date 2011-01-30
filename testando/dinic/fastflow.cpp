#include <cstdio>
#include <vector>
#include <set>
#include <queue>
#include <cstring>
#include <map>
#include <cassert>
#include <ext/hash_set>
using namespace std;
using namespace __gnu_cxx;

#define all(v) (v).begin(), (v).end()
#define _foreach(it, b, e) for(__typeof__(b) it = (b); it != (e); ++it)
#define foreach(x...) _foreach(x)

typedef long long huge;

const int inf = 0x3f3f3f3f;
const huge hugeinf = 0x3f3f3f3f3f3f3f3fll;

////////////////////////////////////////////////////////////////////////////////
// Dinic's blocking flow algorithm (Ahuja & Orlin) - O(V^2 E)
////////////////////////////////////////////////////////////////////////////////

// graph structure
const int maxn = 5050;
vector<int> graph[maxn];
huge cap[maxn][maxn];
int V;

// algorithm data
int dist[maxn], father[maxn], ndst[maxn], curarc[maxn];
huge dinic(int s, int t)
{
  huge ans = 0;
  // reverse bfs
  fill(ndst, ndst+V+1, 0); ndst[V] = V;
  fill(dist, dist+V, V); --ndst[V]; dist[t] = 0; ++ndst[0];
  fill(curarc, curarc+V, 0);
  queue<int> q; q.push(t);
  while(!q.empty()) {
    int v = q.front(); q.pop();
    foreach(it, all(graph[v])) if(dist[*it] == V && cap[*it][v] > 0)
      { --ndst[V]; dist[*it] = dist[v]+1; ++ndst[dist[*it]]; q.push(*it); }
  }
  int i = s;
  father[i] = -1;
  while(dist[s] < V) {
    if(i == t) { // augment
      huge d = hugeinf;
      for(int j=t; j!=s; j=father[j]) d = min(d, cap[father[j]][j]);
      for(int j=t; j!=s; j=father[j]) cap[father[j]][j] -= d, cap[j][father[j]] += d;
      ans += d; i = s; continue;
    }
    bool found = false;
    foreach(it, curarc[i] + all(graph[i]))
      if(dist[i] == dist[*it]+1 && cap[i][*it] > 0)
	{ found = true; father[*it] = i; curarc[i] = it-graph[i].begin(); i = *it; break; }
    if(!found) { // retreat
      curarc[i] = 0;
      int tmp = dist[i];
      --ndst[dist[i]]; dist[i] = V;
      foreach(it, all(graph[i])) if(cap[i][*it] > 0) dist[i] = min(dist[i], dist[*it]+1);
      ++ndst[dist[i]];
      if(ndst[tmp] == 0) break;
      if(i != s) i = father[i];
    }
  }
  return ans;
}

hash_set<int> edges[maxn];

int main()
{
  int E;
  scanf(" %d %d", &V, &E);

  for(int i=0; i<E; ++i)
    {
      int a, b, c;
      scanf(" %d %d %d", &a, &b, &c);
      --a, --b;
      if(a != b)
	{
	  edges[a].insert(b);
	  edges[b].insert(a);
	  cap[a][b] = cap[b][a] += c;
	}
    }
  for(int i=0; i<V; ++i) graph[i] = vector<int>(all(edges[i]));
  printf("%lld\n", dinic(0, V-1));

  return 0;
}
