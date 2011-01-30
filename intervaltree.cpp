template <int n>
class itree // (expt 2 14)16384 > 10^5
{
  int fim;
  int it[n<<1];
  int (*func)(int,int);
  int sentinel;
  inline int pai(int p) { return n|(p>>1); }
public:
  itree(int(*f)(int,int), int s)
  {
    func=f;
    sentinel=s;
    fim=(n<<1)-1;
  }
  void update(int pos, int val)
  {
    for(;pos<fim;pos=pai(pos))
      it[pos]=func(it[pos],val);
  }
  int query(int left, int right)
  {
    int res=sentinel;
    for(;left<=right;left=pai(left), right=pai(right))
      {
	if (left&1)
	  res=func(res, it[left++]);
	if (!(right&1))
	  res=func(res, it[right--]);
      }
    return res;
  }
  void clear()
  {
    for(int i=0; i<=fim; ++i)
      it[i]=sentinel;
  }
};
