int addtree[4*MAXN+1]; //MAXN eh o tamanho maximo do vetor de interesse.
int N; // N eh 2*(a menor potencia de dois maior que o tamanho do vetor).
int valor(int p)
{//computa o valor do p-esimo elemento do vetor.
  int r=0, i=1, w=N/2;
  while(i<N&&w>0)
    {
      r+=addtree[i];
      i<<=1; w>>=1;
      if (p>=w) { p-=w; i|=1; }
    }
  return r;
}
void addint(int v, int lb, int ub, int p=1, int w=N/2)
{//addint(v,lb,ub) adiciona v no intervalo FECHADO [lb,ub].
  if (lb==0&&ub==(w-1))
    addtree[p]+=v;
  else if (lb<=ub)
    {
      w>>=1; p<<=1;
      addint(v, lb, ub<?(w-1), p, w);
      addint(v, (lb-w)>?0, ub-w, p|1, w);
    }
}
