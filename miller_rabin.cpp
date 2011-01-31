#include<cstdio>
#include<ctime>
typedef long long int huge;
typedef unsigned long long int u_huge;


static const int accur=3;
const huge magic_number32[3]= {2,7,61};
const huge magic_number64[9]= {2,3,5,7,11,13,17,23,29};
/*inline huge mult(u_huge a,u_huge b,u_huge n)//para 64 bits
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
}*/
inline huge mult(huge a,huge b,huge n)
{
    return (a*b)%n;
}
inline huge exp(huge a,huge b,huge n)
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
    for(int i=1; i<s; i++)
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
    if(n==2 || n==7 || n==61) return true;
    if(n==0 || n==1)return false;
    for(int i=0; i<6; i++)//testando primos pequenos
        if(n==magic_number64[i])
            return true;
        else if(n%magic_number64[i]==0)
            return false;
    int s = __builtin_ctz(n-1);
    for(int i=0; i<accur; i++)
    {
        if(test_composite(exp(magic_number32[i],(n-1)>>s,n),n,s))
            return false;
    }
    return true;
}

int main ()
{
    printf("%d\n",is_prime(561));
    printf("%d\n",is_prime(1300031));
    return 0;
}
