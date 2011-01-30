#define MAXN //N+1, pois os indices comecam em 1.
int bit1[MAXN]; //update em um elemento, query no intervalo.

void update_elem(int pos, int diff)
{
  if (pos)
    for(;pos<=N; pos+=pos& -pos)
      bit1[pos]+=diff;
}

int query_int(int pos) // [1,pos]
{
  int res=0;
  for(;pos>0; pos-=pos& -pos)
    res+=bit1[pos];
  return res;
}

int bit2[MAXN]; //update no intervalo, query no elemento

void update_int(int pos, int diff) // [1,pos]
{
  for(;pos>0;pos-=pos& -pos)
    bit2[pos]+=diff;
}

int query_pos(int pos)
{
  int res=0;
  if (pos)
    for(;pos<=N; pos+=pos& -pos)
      res+=bit2[pos];
  return res;
}

////////////////////////////////////////////////////////////////////////////////
// 2D
////////////////////////////////////////////////////////////////////////////////
#define MAXN 32
const int N = MAXN/2;
const int root = MAXN-1;
inline int pai(int v) {return N+(v>>1);}
int BIT[MAXN][MAXN]; //[y][x]
void update1(int p, int lx, int ux, int v)
{
  for(;lx<=ux;lx=pai(lx),ux=pai(ux))
    {
      if (lx&1)
	BIT[p][lx++]+=v;
      if (!(ux&1))
	BIT[p][ux--]+=v;
    }
}
void update2(int ly, int uy, int lx, int ux, int v)
{
  for(;ly<=uy;ly=pai(ly),uy=pai(uy))
    {
      if (ly&1)
	update1(ly++, lx, ux, v);
      if (!(uy&1))
	update1(uy--, lx, ux, v);
    }
}
int query(int y, int x)
{
  int r=0;
  for(;x<root;x=pai(x))
    for(int t=y; t<root; t=pai(t))
      r+=BIT[t][x];
  return r;
}
