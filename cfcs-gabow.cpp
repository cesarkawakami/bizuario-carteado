const int maxn = 161010;
vector<int> ladj[maxn];
int preorder[maxn],cfc[maxn];
int num_order,num_cfc,n,m;
stack<int> P,S;
void dfs_gabow(int v)
{
  preorder[v]=num_order++;
  P.push(v),S.push(v);
  foreach(it,all(ladj[v]))
    {
      if(preorder[*it]==-1)dfs_gabow(*it);
      else if(cfc[*it]==-1)while(preorder[P.top()]>preorder[*it]) P.pop();
    }
  if(P.top()==v)
    {
      while(S.top()!=v){cfc[S.top()]=num_cfc;S.pop();}
      S.pop(),P.pop();
      cfc[v]=num_cfc++;
    }
}
void gabow()
{
  num_cfc=num_order=0;
  for(int i=0;i<n;i++)
    cfc[i]=preorder[i]=-1;
  for(int i=0;i<n;i++)
    if(preorder[i]==-1)dfs_gabow(i);
}

