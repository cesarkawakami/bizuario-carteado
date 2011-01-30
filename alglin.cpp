////////////////////////////////////////////////////////////////////////////////
// ELIMINACAO GAUSSIANA
////////////////////////////////////////////////////////////////////////////////
// obtem echelon form (singular) ou row-reduced echelon form (nonsingular)
// retorna true se nonsingular, false caso contrario
struct submult_op
{
  double p;
  submult_op(const double &a) : p(a) {}
  double operator()(const double &a, const double &b) {return a-p*b;}
};
bool GaussElim(vector<vector<double> > &A)
{
  int m=A.size(),n=A[0].size();
  int i,j;
  for(i=0,j=0; i<m&&j<n;)
    {
      int maxi = i;
      for(int k=i+1; k<m; ++k) if( fabs(A[k][j]) > fabs(A[maxi][j]) ) maxi = k;
      swap(A[i], A[maxi]);
      if( fabs(A[i][j]) > eps )
	{
	  transform(all(A[i]),A[i].begin(),bind2nd(divides<double>(),A[i][j]));
	  for(int k=i+1; k<m; ++k)
	    transform(all(A[k]),A[i].begin(),A[k].begin(),submult_op(A[k][j]));
	  ++i;
	}
      ++j;
    }
  if( i == m && j == m )
    {
      for(int i=0; i<m; ++i) for(int j=0; j<i; ++j)
	transform(all(A[j]),A[i].begin(),A[j].begin(),submult_op(A[j][i]));
      return true;
    }
  return false;
}
