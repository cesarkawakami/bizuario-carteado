// ou usar __gcd
int gcd(int a, int b) { return a?gcd(b,a%b):abs(a); }

int lcm(int x, int y)
{
  if(x&&y) return abs(x)/gcd(x,y)*abs(y);
  return abs(x|y);
}

// retorna (d, a, b), onde d = ax + by.
struct bezout {int d, a, b;};
bezout extgcd(int x, int y)
{
  if(y==0) return (bezout){x, 1, 0};
  bezout s = extgcd(y, x%y);
  return (bezout){s.d, s.b, s.a-x/y*s.b};
}
