////////////////////////////////////////////////////////////////////////////////
// BigNum Library
////////////////////////////////////////////////////////////////////////////////
// + e -: in-place; * e /: out-of-place
////////////////////////////////////////////////////////////////////////////////
const int EXP = 4;      // base = 10^EXP
const int BASE = 10000; // base. Cuidado por overflows em multiplicacao.
                        // 4 eh um bom valor.
const int TAM = 2048;   // Maximo numero representavel = BASE^TAM = 10^(EXP*TAM)

struct bignum
{
  int n;                        // numero de digitos
  int d[TAM];
  bignum(int x = 0) : n(1) { memset(d, 0, sizeof(d)); d[n++] = x; fix(); }
  bignum(const char *s) : n(1) { // De um trim em s antes!
    memset(d, 0, sizeof(d));
    char sign = 1;
    if(s[0] == '-') { sign = -1; ++s; }
    char *b = strdup(s), *e = b + strlen(b);
    while(e > b)
      { *e = 0; e = max(b, e-EXP); sscanf(e, "%d", d+n); d[n++] *= sign; }
    free(b); fix();
  }
  bignum &fix(int m = 0) {
    n = max(m, n);
    char sign = 0;
    for(int i=1, car=0; i <= n || car && (n = i); ++i)
      { d[i]+=car; car=d[i]/BASE; d[i]%=BASE; if(d[i]) sign=(d[i]>0)?1:-1; }
    for(int i=n-1; i>0; --i)
      if(d[i] * sign < 0) { d[i] += BASE * sign; d[i+1] -= sign; }
    while(n && !d[n]) --n;
    return *this;
  }
  char compare(const bignum &x = 0) const {
    for(int i=max(n, x.n); i>0; --i) {
      if(d[i] < x.d[i]) return -1;
      else if(d[i] > x.d[i]) return 1;
    }
    return 0;
  }
  bool operator<(const bignum &x) const { return compare(x) < 0; }
  bool operator==(const bignum &x) const { return compare(x) == 0; }
  bool operator!=(const bignum &x) const { return compare(x) != 0; }
  char *c_str() {               // dar free depois!
    char *s = (char*)malloc(EXP*n+10);
    for(int k=n-1, i=sprintf(s, "%d", d[n]); k>0; i+=sprintf(s+i, "%04d", abs(d[k--])));
    return s;
  }
  bignum &operator+=(const bignum &x)
  { for(int i=1; i<=x.n; ++i) d[i] += x.d[i]; return fix(x.n); }
  bignum &operator-=(const bignum &x)
  { for(int i=1; i<=x.n; ++i) d[i] -= x.d[i]; return fix(x.n); }
  bignum operator+(const bignum &x) { return bignum(*this) += x; }
  bignum operator-(const bignum &x) { return bignum(*this) -= x; }
  bignum operator-() { bignum b(0); return b -= *this; }
  void ams(const bignum &x, const int &m, const int &b) {
    for(int i=1, car=0; (i <= x.n || car) && (n = i+b); ++i)
      { d[i+b] += x.d[i]*m + car; car = d[i+b]/BASE; d[i+b] %= BASE; }
  }
  bignum operator*(const bignum &x) const {
    bignum r(0);
    for(int i=1; i<=n; ++i) r.ams(x, d[i], i-1);
    return r;
  }
  bignum &operator*=(const bignum &x) { return *this = *this * x; }
  bignum div(const bignum &x) {
    if(x == 0) throw "divisao por zero";
    bignum tr, q(0); q.n = max(n - x.n + 1, 0);
    for(int i=q.n; i>0; --i) {
      for(int l=0, u=BASE-1; u > l; ) {
	tr = 0;
	q.d[i] = (l+u+1)/2;
	tr.ams(x, q.d[i], i-1);
	if(*this < tr) u = q.d[i] - 1, q.d[i] = u;
	else l = q.d[i]; // tr <= *this
      }
      tr = 0; tr.ams(x, q.d[i], i-1); *this -= tr;
    }
    return q.fix();
  }
  bignum &operator/=(const bignum &x) { return *this = div(x); }
  bignum &operator%=(const bignum &x) { div(x); return *this; }
  bignum operator/(const bignum &x) { return bignum(*this).div(x); }
  bignum operator%(const bignum &x) { return bignum(*this) %= x; }
  bignum pow(int x) {        // Pode ser otimizado a tempo logaritmico,
                             // mas custaria mais espaco.
    if(x < 0) return (*this == 1 || *this == -1) ? pow(-x) : 0;
    bignum r = 1;
    for(int i=0; i<x; ++i) r *= *this;
    return r;
  }
  bignum root(int x) {
    if(compare() == 0 || compare() < 0 && x % 2 == 0) return 0;
    if(*this == 1 || x == 1) return *this;
    if(compare() < 0) return -(-*this).root(x);
    // melhorando o chute inicial
    bignum d(0);
    d.n = this->n/x + 2;
    for(int i=1; i<=d.n; ++i) d.d[i] = BASE-1;
    bignum a = 1; //, d = *this;
    while(d != 1) {
      bignum b = a + (d /= 2);
      if(compare(b.pow(x)) >= 0) { d += 1; a = b; }
    }
    return a;
  }
};
