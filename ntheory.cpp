int gcd(int x, int y) { return y?gcd(y,x%y):abs(x); }

int lcm(int x, int y)
{
  if(x&&y) return abs(x)/gcd(x,y)*abs(y);
  return abs(x|y);
}

bool is_prime(int n)
{
  if(n<0) return is_prime(-n);
  if(n<5||n%2==0||n%3==0) return (n==2||n==3);
  int maxn = sqrt(n)+2;
  for(int i=5; i<maxn; i+=6)
    if(n%i==0||n%(i+2)==0) return false;
  return true;
}

void squeeze(map<int,int> &f, int &n, int &p) {for(;n%p==0;n/=p) ++f[p];}
map<int,int> factor(int n)
{
  map<int,int> ans;
  if(n<0) return factor(-n);
  if(n<2) return ans;
  squeeze(ans,n,2); squeeze(ans,n,3);
  int maxn = sqrt(n)+2;
  for(int i=5; i<maxn; i+=6)
    squeeze(ans,n,i), squeeze(ans,n,i+2);
  if(n>1) ++ans[n];
  return ans;
}

typedef struct{int d, a, b;} bezout;

// retorna (d, a, b), onde d = ax + by.
bezout extgcd(int x, int y)
{
  if(y==0) return (bezout){x, 1, 0};
  bezout s = extgcd(y, x%y);
  return (bezout){s.d, s.b, s.a-x/y*s.b};
}
