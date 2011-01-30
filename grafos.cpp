////////////////////////////////////////////////////////////////////////////////
// DIJKSTRA COM SET E ADJACENCY LIST (EH FACIL MODIFICAR PRA ADJ MTX)
////////////////////////////////////////////////////////////////////////////////
vector<int> dist(tamanho do grafo);
void dijkstra(int ori)
{
  set<pair<int, int> > td;
  td.insert(make_pair(0, ori));
  dist[ori] = 0;
  while(!td.empty())
    {
      int v = td.begin()->second;
      td.erase(td.begin());
      foreach(it, all(graph[v])) // lista de adj
	if( dist[it->first] > dist[v] + it->second )
	  {
	    if( dist[it->first] != INF )
	      td.erase(td.find(make_pair(dist[it->first], it->first)));
	    dist[it->first] = dist[v] + it->second;
	    td.insert(make_pair(dist[it->first], it->first));
	  }
    }
}

////////////////////////////////////////////////////////////////////////////////
// UNWEIGHTED BIPARTITE MATCHING
////////////////////////////////////////////////////////////////////////////////
// Matching resultante: para matching[i]!=-1, (i, matching[i]).
// tab[i][j] -> de i para j. m linhas e n colunas. (mtx no sentido matematico)
////////////////////////////////////////////////////////////////////////////////
// Ideias / extensoes:
// - Para determinar se a aresta encontra-se em algum perfect matching, force
//   ela estar no matching e de um augment() no vertice que sobra (nao esquecer
//   de limpar seen)
////////////////////////////////////////////////////////////////////////////////
const int nmax = 200;

struct UBMatching
{
  int m, n, size, matching[nmax], seen[nmax];
  int (*tab)[nmax];
  bool augment(int v)
  {
    for(int i=0;i<m;i++)
      if(tab[i][v] && !seen[i])
	{
	  seen[i]=1;
	  if(matching[i]<0 || augment(matching[i]))
 	    { matching[i]=v; return 1; }
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
      { memset(seen, 0, sizeof(seen)); size += augment(i); }
  }
} ubm;

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

////////////////////////////////////////////////////////////////////////////////
// EDMONDS-KARP MAXFLOW - O(max(MaxFlow*V, VE^2))
// Muda o m para ser grafo residual.
// Na lista de adj, tem de colocar aresta nos dois sentidos, mesmo em directed
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

/////////////////////////////////////////////////////////////////////
//Testa se existe um perfect-matching num grafo generico.
/////////////////////////////////////////////////////////////////////
struct perfectmatching //depende de maxflow
{
  int floodfill(int v, int cor)
  {
    int r=1;
    comp[v]=cor;
    foreach(it, all(adj[v]))
      if (!comp[*it])
	r+=floodfill(*it,cor);
    return r;
  }
  bool perfectmatch()
  {
    int c=1, nk, fk;
    for(int i=1; i<=n; ++i)
      if (!comp[i])
	{
	  nk=floodfill(i,c);
	  if (nk&1)
	    return false;
	  for(int j=1; j<=n; ++j)
	    if (comp[j]==c)
	      {
		lista[0].push_back(j);
		lista[j].push_back(0);
		mtx[0][j]=1;
	      }
	  fk=ekmf.augment();
	  ++c;
	  if (fk!=nk)
	    return false;
	}
    return true;
  }
} pfmt;

//////////////////////////////////////////////////////////
// Stable Marriage Problem.
//////////////////////////////////////////////////////////
int girl_order[MAXN][MAXN]; //ordem de preferencia das garotas
int boy_rate[MAXN][MAXN];   //nota dada pelo garoto para a garota
int girl_pos[MAXN];         //posicao atual da garota - inicializar com zero.
int boy_pair[MAXN];         //par do garoto - inicializar com infinito.
int n;                      //o algoritmo é arroz!

bool smp()
{
  for(int i=0, k=0, j=0; k<n; i=++k)
    while(1)
      {
	if (girl_pos[i]<n)
	  j=girl_order[i][girl_pos[i]];
	else
	  return 0;
	if (boy_pair[j]>n)
	  {
	    boy_pair[j]=i;
	    break;
	  }
	else if (boy_rate[j][i]<boy_rate[j][boy_pair[j]])
	  {
	    girl_pos[boy_pair[j]]++;
	    swap(boy_pair[j], i);
	  }
	else
	  girl_pos[i]++;
      }
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
// HUNGARIAN ALGORITHM (MAXIMUM WEIGHTED BIPARTITE MATCHING)
////////////////////////////////////////////////////////////////////////////////
#define N 55             //max number of vertices in one part
#define INF 100000000    //just infinity

int cost[N][N];          //cost matrix
int n, max_match;        //n workers and n jobs
int lx[N], ly[N];        //labels of X and Y parts
int xy[N];               //xy[x] - vertex that is matched with x,
int yx[N];               //yx[y] - vertex that is matched with y
bool S[N], T[N];         //sets S and T in algorithm
int slack[N];            //as in the algorithm description
int slackx[N];           //slackx[y] such a vertex, that
                         // l(slackx[y]) + l(y) - w(slackx[y],y) = slack[y]
int prev[N];             //array for memorizing alternating paths

void init_labels()
{
  memset(lx, 0, sizeof(lx));
  memset(ly, 0, sizeof(ly));
  for (int x = 0; x < n; x++)
    for (int y = 0; y < n; y++)
      lx[x] = max(lx[x], cost[x][y]);
}

void update_labels()
{
  int x, y, delta = INF;             //init delta as infinity
  for (y = 0; y < n; y++)            //calculate delta using slack
    if (!T[y])
      delta = min(delta, slack[y]);
  for (x = 0; x < n; x++)            //update X labels
    if (S[x]) lx[x] -= delta;
  for (y = 0; y < n; y++)            //update Y labels
    if (T[y]) ly[y] += delta;
  for (y = 0; y < n; y++)            //update slack array
    if (!T[y])
      slack[y] -= delta;
}

void add_to_tree(int x, int prevx) 
{
  S[x] = true;                    //add x to S
  prev[x] = prevx;                //we need this when augmenting
  for (int y = 0; y < n; y++)    //update slacks, because we add new vertex to S
    if (lx[x] + ly[y] - cost[x][y] < slack[y])
      {
	slack[y] = lx[x] + ly[y] - cost[x][y];
	slackx[y] = x;
      }
}

void augment()                         //main function of the algorithm
{
  if (max_match == n) return;        //check wether matching is already perfect
  int x, y, root;                    //just counters and root vertex
  int q[N], wr = 0, rd = 0;          //q - queue for bfs, wr,rd - write and read
  memset(S, false, sizeof(S));       //init set S
  memset(T, false, sizeof(T));       //init set T
  memset(prev, -1, sizeof(prev));    //init set prev - for the alternating tree
  for (x = 0; x < n; x++)            //finding root of the tree
    if (xy[x] == -1)
      {
	q[wr++] = root = x;
	prev[x] = -2;
	S[x] = true;
	break;
      }

  for (y = 0; y < n; y++)            //initializing slack array
    {
      slack[y] = lx[root] + ly[y] - cost[root][y];
      slackx[y] = root;
    }
  while (true)
    {
      while (rd < wr)
        {
	  x = q[rd++];
	  for (y = 0; y < n; y++)
	    if (cost[x][y] == lx[x] + ly[y] &&  !T[y])
	      {
		if (yx[y] == -1) break;
		//augmenting path exists!
		T[y] = true;
		q[wr++] = yx[y];
		//with y, to the queue
		add_to_tree(yx[y], x);
	      }
	  if (y < n) break;
        }
      if (y < n) break;
      
      update_labels();
      wr = rd = 0;                
      for (y = 0; y < n; y++)        
	if (!T[y] &&  slack[y] == 0)
	  {
	    if (yx[y] == -1)
	      {
		x = slackx[y];
		break;
	      }
	    else
	      {
		T[y] = true;
		if (!S[yx[y]])    
		  {
		    q[wr++] = yx[y];
		    add_to_tree(yx[y], slackx[y]);
		  }
	      }
	  }
      if (y < n) break;
    }

  if (y < n)
    {
      max_match++;
      for (int cx = x, cy = y, ty; cx != -2; cx = prev[cx], cy = ty)
        {
	  ty = xy[cx];
	  yx[cy] = cx;
	  xy[cx] = cy;
        }
      augment();
    }
}

int hungarian()
{
  int ret = 0;                      //weight of the optimal matching
  max_match = 0;                    //number of vertices in current matching
  memset(xy, -1, sizeof(xy));    
  memset(yx, -1, sizeof(yx));
  init_labels();                    //step 0
  augment();                        //steps 1-3
  for (int x = 0; x < n; x++)       //forming answer there
    ret += cost[x][xy[x]];
  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// MIN COST MAX FLOW
////////////////////////////////////////////////////////////////////////////////
/*Fluxo de custo mínimo em O(n^2f) (successive shortest path, com potenciais*/
/*Se desejado o menor custo para um certo fluxo c, criar um novo nó com aresta
  para s com custo 0 e capacidade c*/
/*Observação: no final não há ciclo de custo negativo na matriz residual 
  (condição de optimalidade)*/
/*vértices de 0 a n, 0 eh t, n eh s*/

#define MAXN 200
#define INF 10000000000000000LL

int n,m;
/*cap[0][1][0] = fluxo desejado, cst[0][1][0] = 0*/
long long cap[MAXN][MAXN][2];/*cap[i][j][0] = capacidades dadas, cap[i][j][1] = 0*/
long long cst[MAXN][MAXN][2];/*cst[i][j][0] = custos dados, cst[i][j][1] = -cst[j][i][0]*/
int par[MAXN],pr[MAXN];
long long p[MAXN];

void bellman(int s,int t)
{
  int sw;
  int i,j,k,h;
    
  for(i=0;i<=n;i++) {
    p[i] = INF;
    pr[i] = -1;   
  }
  p[s] = 0;
  sw = 1;
  for(h=0;(h<=n)&&(sw);h++) {
    sw = 0;
        
    for(i=0;i<=n;i++) {
      for(j=0;j<=n;j++) {
	for(k=0;k<=1;k++) {
	  if(cap[i][j][k] <= 0)
	    continue;
                    
	  if(p[j] > p[i] + cst[i][j][k]) {
	    p[j] = p[i] + cst[i][j][k];
	    sw++;
	  }
	}
      }
    }
  }
}

bool dijkstra(int s,int t)
{
  int v,i,k;
  bool intree[MAXN];
  long long dst[MAXN],men,c;
    
  for(i=0;i<=n;i++)
    {
      dst[i] = INF;
      par[i] = -1;
      pr[i] = -1;
      intree[i] = false;
    }
  dst[s] = 0;
  v = s;
  while(!intree[v])
    {
      intree[v] = true;
        
      for(i=0;i<=n;i++)
	{
	  if(intree[i])
	    continue;
            
	  for(k=0;k<=1;k++)
	    {
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
      for(i=0;i<=n;i++)
	{
	  if(intree[i])
	    continue;
            
	  if(dst[i] < men)
	    {
	      men = dst[i];
	      v = i;
	    }    
	}
    }
    
  for(i=0;i<=n;i++)
    p[i] += dst[i];   
    
  return intree[t];
}

long long mincost(int s,int t)
{
  int i,j,k,v;
  long long aug,ret;
    
  /*somente se houver aresta com custo negativo, senao comece com p[i] = 0*/
  bellman(s,t);
    
  while(1)
    {
      if(!dijkstra(s,t))
	break;
        
      v = t;
      aug = INF;
      while(v != s)
	{
	  if(aug > cap[par[v]][v][pr[v]])
	    {
	      aug = cap[par[v]][v][pr[v]];
	    }
	  v = par[v];   
	}
        
      v = t;
      while(v != s)
	{
	  cap[par[v]][v][pr[v]] -= aug;
	  cap[v][par[v]][1-pr[v]] += aug;
	  v = par[v];
	}
    }
    
  ret = 0;
  for(i=0;i<=n;i++)
    {
      for(j=0;j<=n;j++)
	{
	  ret += cap[j][i][1] * cst[i][j][0];
	}
    }
  return ret;     
}

////////////////////////////////////////////////////////////////////////////////
// STOER-WAGNER GLOBAL MIN-CUT O(n^3)
////////////////////////////////////////////////////////////////////////////////
const int nmax = 200;

int graph[nmax][nmax];
bool valid[nmax];
bool marked[nmax];
int tightness[nmax];
int order[nmax];
int n, nvalid;

_inline(int minimumCutPhase)()
{
  memset(marked, 0, sizeof(marked));
  memset(tightness, 0, sizeof(tightness));
  int v;
  for(v=0; v<n; ++v)
    if(valid[v]) break;
  tightness[v] = 1;
  for(int tam=0; tam<nvalid; ++tam)
    {
      int c = max_element(tightness, tightness+n) - tightness;
      marked[c] = true;
      order[tam] = c;
      for(int i=0; i<n; ++i)
	if(!marked[i] && valid[i])
	  tightness[i] += graph[c][i];
      tightness[c] = 0;
    }
  int ans = 0;
  for(int i=0; i<n; ++i)
    if(i!=order[nvalid-1] && valid[i])
      ans += graph[i][order[nvalid-1]];
  int &a = order[nvalid-2], &b = order[nvalid-1];
  for(int i=0; i<n; ++i)
    graph[a][i] = graph[i][a] = graph[a][i] + graph[b][i];
  valid[b] = false;
  --nvalid;
  return ans;
}

_inline(int minimumCut)()
{
  int curmin = inf;
  fill(valid, valid+n, true);
  nvalid = n;
  while(nvalid>1)
    curmin = min(curmin, minimumCutPhase());
  return curmin;
}

//////////////////////////////////////////////////////////
// Eulerian Tour
//////////////////////////////////////////////////////////
list<int> euleriantour(int start) //para multigrafos não direcionados.
{
  int v;
  bool viz;
  list<int> res;
  stack<int> dfs;
  dfs.push(start);
  while(!dfs.empty())
    {
      v=dfs.top();
      viz=0;
      for(int i=0; !viz&&i<maxn; ++i)
	if (mtx[v][i])
	  {
	    dfs.push(i);
	    mtx[v][i]--;
	    mtx[i][v]--;
	    viz=1;
	  }
      if(!viz)
	{
	  dfs.pop();
	  res.push_front(v);
	}
    }
  return res;
}

