////////////////////////////////////////////////////////////////////////////////
// KMP O(t+p)
////////////////////////////////////////////////////////////////////////////////
vector <const char *> v;
const char *kmp_search(const char *texto, const char *p)
{
  int T[2000];
  int i,j;
  const char  *result = NULL;
  if (p[0] == 0) return texto;
  T[0] = -1;
  for (i=0; p[i] != 0; i++)
    {
      T[i+1] = T[i] + 1;
      while (T[i+1] > 0 && p[i] != p[T[i+1]-1])
	T[i+1] = T[T[i+1]-1] + 1;
    }
  for (i=j=0; texto[i] != 0; )
    {
      if (j < 0 || texto[i] == p[j])
        {
	  ++i, ++j;
	  if (p[j] == 0)
            {
	      result=texto+i-j;
	      v.push_back(result);
            }
        }
      else j = T[j];
    }
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// Suffix Tree
////////////////////////////////////////////////////////////////////////////////
const int maxn = 100100; //10^5
char str[maxn];
struct node
{
  int st, ed;
  map<char,int> pt;
  int link;
};
node vet[2*maxn];
int total;
_inline(void canonize)(int &s, int &k, int p)
{
  if (k<=p)
    {
      int ss=vet[s].pt.find(str[k])->second,
	kk=vet[ss].st, pp=vet[ss].ed;
      while(pp-kk<=p-k)
	{
	  k+=pp-kk+1;
	  s=ss;
	  if (k<=p)
	    ss=vet[s].pt.find(str[k])->second,
	      kk=vet[ss].st, pp=vet[ss].ed;
	}
    }
}
_inline(void init)(node &x, int st, int ed) {x.st=st; x.ed=ed; x.pt.clear();}
void make_tree()
{
  int s=0, oldr, rr, r, ss, k=0, kk;
  map<char,int>::iterator mi;
  for(int i=0; (!i)||str[i-1]; ++i)
    {
      oldr=0;
      while(1)
	{
	  if (k<i)
	    {
	      ss=(mi=vet[s].pt.find(str[k]))->second;
	      kk=vet[ss].st;
	      if (str[i]==str[kk+i-k])
		break;
	      else
		{
		  init(vet[r=total++], kk, kk+i-k-1);
		  vet[r].pt.insert(make_pair(str[kk+i-k],ss));
		  mi->second=r;
		  vet[ss].st+=i-k;
		}
	    }
	  else if (vet[s].pt.find(str[i])==vet[s].pt.end())
	    r=s;
	  else
	    break;
	  init(vet[rr=total++], i,inf);
	  vet[r].pt.insert(make_pair(str[i],rr));
	  if (oldr)
	    vet[oldr].link=r;
	  oldr=r;
	  if (!s)
	    ++k;
	  s=vet[s].link;
	  canonize(s, k, i-1);
	}
      if (oldr)
	vet[oldr].link=s;
      canonize(s,k,i);
    }
}
