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

/*Fluxo de custo mínimo em O(n^2f) (successive shortest path, com potenciais*/
/*Se desejado o menor custo para um certo fluxo c, criar um novo nó com aresta
  para s com custo 0 e capacidade c*/
/*Observação: no final não há ciclo de custo negativo na matriz residual (condição de optimalidade)*/
/*vértices de 0 a n, 0 eh t, n eh s*/

#define MAXN 200
#define INF 10000000000000000LL

int n,m;
/*cap[0][1][0] = fluxo desejado, cst[0][1][0] = 0*/
long long cap[MAXN][MAXN][2];/*cap[i][j][0] = capacidades dadas, cap[i][j][1] = 0*/
long long cst[MAXN][MAXN][2];/*cst[i][j][0] = custos dados, cst[i][j][1] = -cst[j][i][0]*/
int par[MAXN],pr[MAXN];
long long p[MAXN];

void bellman(int s,int t){
    
  int sw;
  int i,j,k,h;
    
  for(i=0;i<=n;i++){
    p[i] = INF;
    pr[i] = -1;   
  }
  p[s] = 0;
  sw = 1;
  for(h=0;(h<=n)&&(sw);h++){
    sw = 0;
        
    for(i=0;i<=n;i++){
      for(j=0;j<=n;j++){
	for(k=0;k<=1;k++){
	  if(cap[i][j][k] <= 0)
	    continue;
                    
	  if(p[j] > p[i] + cst[i][j][k]){
	    p[j] = p[i] + cst[i][j][k];
	    sw++;
	  }
	}
      }
    }
        
  }
        
}


bool dijkstra(int s,int t){
    
  int v,i,k;
  bool intree[MAXN];
  long long dst[MAXN],men,c;
    
  for(i=0;i<=n;i++){
    dst[i] = INF;
    par[i] = -1;
    pr[i] = -1;
    intree[i] = false;
  }
  dst[s] = 0;
  v = s;
  while(!intree[v]){
    intree[v] = true;
        
    for(i=0;i<=n;i++){
      if(intree[i])
	continue;
            
      for(k=0;k<=1;k++){
	if(cap[v][i][k] <= 0)
	  continue;
                
	c = cst[v][i][k] + p[v] - p[i];
	if(dst[i] > dst[v] + c){
	  dst[i] = dst[v] + c;
	  par[i] = v;
	  pr[i] = k;
	}
      } 
    }
        
    men = INF;
    for(i=0;i<=n;i++){
      if(intree[i])
	continue;
            
      if(dst[i] < men){
	men = dst[i];
	v = i;
      }    
    }
  }
    
  for(i=0;i<=n;i++)
    p[i] += dst[i];   
    
  return intree[t];
}

long long mincost(int s,int t){
    
  int i,j,v;
  long long aug,ret;
    
  /*somente se houver aresta com custo negativo, senao comece com p[i] = 0*/
  bellman(s,t);
    
  while(1){

    if(!dijkstra(s,t))
      break;
        
    v = t;
    aug = INF;
    while(v != s){
      if(aug > cap[par[v]][v][pr[v]]){
	aug = cap[par[v]][v][pr[v]];
      }
      v = par[v];   
    }
        
    v = t;
    while(v != s){
      cap[par[v]][v][pr[v]] -= aug;
      cap[v][par[v]][1-pr[v]] += aug;
      v = par[v];
    }
  }
    
  ret = 0;
  for(i=0;i<=n;i++){
    for(j=0;j<=n;j++){
      ret += cap[j][i][1] * cst[i][j][0];
    }
  }
  return ret;     
}

int main()
{
  int teste = 1;
  while(scanf(" %d %d", &n, &m) == 2)
    {
      memset(cap, 0, sizeof(cap));
      memset(cst, 0, sizeof(cst));
      while(m--)
	{
	  int a, b, c;
	  scanf(" %d %d %d", &a, &b, &c);
	  cap[a][b][0] = cap[b][a][0] = 1;
	  cst[a][b][0] = cst[b][a][0] = c;
	  cst[a][b][1] = cst[b][a][1] = -cst[a][b][0];
	}
      int d, k;
      scanf(" %d %d", &d, &k);
      for(int i=0; i<=n; ++i)
	for(int j=0; j<=n; ++j)
	  if(cap[i][j][0] == 1)
	    cap[i][j][0] = cap[j][i][0] = k;
      cap[0][1][0] = d;

      printf("Instancia %d\n", teste++);
      long long mc = mincost(0, n);
      if(cap[0][1][0] != 0)
	printf("impossivel\n");
      else
	printf("%lld\n", mc);
      printf("\n");
    }
}
