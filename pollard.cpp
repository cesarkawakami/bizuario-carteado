#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include <cctype>
#include <iostream>
#include <sstream>
#include <deque>
#include <stack>
#include <queue>
#include <list>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <utility>
#include <functional>
#include <algorithm>
using namespace std;
using namespace rel_ops;

#define _foreach(it, b, e) for(__typeof(b) it=(b); it!=(e); ++it)
#define foreach(x...) _foreach(x...)
#define all(v) (v).begin(), (v).end()
#define rall(v) (v).rbegin(), (v).rend()

typedef long long int huge;
const int inf = 0x3fff3fff;
const huge hugeinf = 0x3fff3fff3fff3fffLL;
const double eps = 1e-9;

typedef long long int huge;
typedef unsigned long long int u_huge;

int accur=9;
int magic_number32[]={2,7,61};
int magic_number64[]={2,3,5,7,11,13,17,23,29};
huge random_vet[]={3,49,1337,5,171,};

inline huge mult(u_huge a,u_huge b,u_huge n)
{
  u_huge y = (u_huge)((long double)a*(long double)b/n + (long double)1/2);//floor(a*b/m)
  y=y*n; //m*floor(a*b/m) %z
  u_huge r = a*b-y;
  if( (huge) r <0)//normalization needed?
    {
      r = r+n;
      // y = y-1; (a%b)/m quotient
    }
  return r;//(a*b)%m residue
}
huge exp(huge a,huge b,huge n)
{
  huge acc = 1;
  while(b)
    {
      if(b&1)
	acc=mult(acc,a,n);
      b>>=1;
      a=mult(a,a,n);
    }
  return acc;
}
bool test_composite(huge x,huge n,int s)
{
  if(x==1 || x==n-1)
    return false;
  for(int i=1;i<s;i++)
    {
      x=mult(x,x,n);
      if(x==1)return true;
      if(x==n-1)return false;
    }
  return true;
}
bool is_prime(huge n)
{
  if(n<0)
    n=-n;
  if(n==2) return true;
  if(n==0 || n==1 || !(n&1))return false;
  int s = __builtin_ctzll(n-1);
  for(int i=0;i<accur;i++)
    {
      if(n!=magic_number64[i] && test_composite(exp(magic_number64[i],(n-1)>>s,n),n,s))
	return false;
    }
  return true;
}
huge gcd(huge a, huge b)
{
  return a<b?gcd(b,a):b?gcd(b,a%b):a;
}
huge rseq(huge x, huge c, huge n)
{
  return (mult(x,x,n)+c)%n;
}
void pollard(huge n)
{
  if (n<0)
    n=-n;
  if (n<=1) { // 0 e 1 são casos especiais.
    printf (" %lld^1", n);
  }
  else if (is_prime(n)) { // sobrou só um único primo?
    printf(" %lld^1", n);
  }
  else { //finalmente o algoritmo de verdade!
    huge x=2, y, d;
    for(int c=0; c<5; ++c)
      {
	d=1;
	y=rseq(x,random_vet[c],n);
	for(int p=1, i=0; d==1; ++i)
	  {
	    d=gcd(llabs(x-y),n);
	    if (i&p)
	      {
		x=y;
		p<<=1;
	      }
	    y=rseq(y,random_vet[c],n);
	  }
	if (d!=n)
	  {
	    x=0;
	    while(n%d==0)
	      {
		n/=d;
		++x;
	      }
	    if (is_prime(d))
	      printf(" %lld^%lld", d, x);
	    else
	      {
		printf(" (");
		pollard(d);
		printf(" )^%lld", x);
	      }
	    if (n>1)
	      pollard(n);
	    return;
	  }
      }
    printf("FUUUUUUUU\n");
  }
}


int main()
{
  huge n;
  while(1)
    {
      scanf("%lld", &n);
      if (!n)
	break;
      printf("%lld:", n);
      pollard(n);
      printf("\n");
    }
  return 0;
}
