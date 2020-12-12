#include <fstream>
#include <iostream>
#include <Dense>
#include <vector>
#include "Sparse"


using Eigen::MatrixXd;
using Eigen::VectorXd;



using namespace std;
using namespace Eigen;

int sgn(double x)
{
  if (x>0)
  return 1;
  else
  return -1;
}

/*vector<Eigen::Matrix<double, Dynamic, Dynamic>>  qr_decomposition(Eigen::MatrixXd A)
{
vector<Eigen::Matrix<double, Dynamic, Dynamic>> X,PK ;
Eigen::Matrix<double, Dynamic, Dynamic>  P,W,Q,R,S,I ;
Eigen::VectorXd r,z,p,e ;
int n = A.rows() ;
I = MatrixXd::Identity(n,n);
e.resize(n);
z.resize(n);
z.setZero() ;
W.setZero() ;
Q.setZero() ;
R.setZero() ;
P.setZero() ;
Q.resize(n,n) ;
R.resize(n,n) ;
W.resize(n,n) ;
P.resize(n,n) ;
int j(0) ;
double s(0.), beta(0.) ;
e.setZero() ;
e(0) = 1;
R.col(0).setZero() ;
R.col(0)(0) = -sgn(A(0,0))*sqrt(A.col(0).dot(A.col(0)));
for (int i = 0 ; i<n ; ++i)
{
if (i<0)
{
z(i)=0;
}
else if (i==0)
{
z(i)=-R.col(0)(0)+A(i,i);
}
else
{
z(i)=A(i,0);
}
}

W.col(0)=(1.0/(sqrt(z.dot(z))))*z;
P = I - 2*W.col(0)*(W.col(0).transpose());
Q.col(0) = P*e;
cout << Q.col(0) << endl;
PK.push_back(P);
for (int k = 1 ; k<n  ; ++k)
{
/*if (k>0)
{
S=PK[0];
for (int l = 1 ; l<k-1 ; ++l)
{
S=S*PK[l];
}
R.col(k)=S*A.col(k);
}








e.setZero() ;
e(k) = 1;
s=0;
for (int j=k ; j<n ; ++j )
{
s+=A(j,k)*A(j,k);
}
beta = sgn(A(k,k))*sqrt(s);
for (int i = 0 ; i<n ; ++i)
{
if (i<k)
{
z(i)=0;
}
else if (i==k)
{
z(i)=beta+A(i,i);
}
else
{
z(i)=A(i,k);
}



W.col(k)=(1.0/(sqrt(z.dot(z))))*z;
P = I - 2*W.col(k)*(W.col(k).transpose());

Q.col(k) = P*e;
cout << R.col(k) << endl ;

PK.push_back(P);
S=PK[0];
for (int l = 1 ; l<k-1 ; ++l)
{
S=S*PK[l];
}
R.col(k)=S*A.col(k);
Q.col(k)=P*S*A.col(k);



}
}


X.push_back(Q);
X.push_back(R);
return X;
}*/
vector<Eigen::Matrix<double, Dynamic, Dynamic>>  qr_decomposition(Eigen::MatrixXd A)
{
  int n = A.cols();
  double beta;
  Eigen::Matrix<double, Dynamic, Dynamic>  H,W,Q,R,S,B,D,I;
  VectorXd z,w,x,e;
  vector<Eigen::Matrix<double, Dynamic, Dynamic>> P,X ;
  R = A ;

  Q = MatrixXd::Identity(n,n);


  for (int k = 0 ; k<n ; ++k)
  {
    e.resize(n-k);
    e.setZero();
    e(0) = 1;
    H = MatrixXd::Identity(n,n);
    I = MatrixXd::Identity(n-k,n-k);

    z.resize(n-k) ;
    w.resize(n-k);
    z.setZero() ;


    x = R.col(k).segment(k,n-k);
    z = x + sqrt(x.dot(x))*e ;

    w = (1/(sqrt(z.dot(z))))*z ;

    D = I - 2*w*(w.transpose()) ;
    if (k==0)
    {
      H = D ;
    }
    H.block(k,k,n-k,n-k) = D;
    P.push_back(H);
    R = H*R;
  }
  for (int i = 0 ; i<n ; ++i)
  {
    Q=Q*P[i] ;
  }
  X.push_back(Q);
  X.push_back(R);
  return X;

}
int main()
{
  Matrix<double, 3, 3> A, B ;
  A << 2, -2, 18,                           // Initialize A. The elements can also be
  2, 1, 0,                                // matrices.
  1, 2, 0;
  vector<Eigen::Matrix<double, Dynamic, Dynamic>> X ;
  B = A ;
  X=qr_decomposition(A);
  cout << B << endl;
  cout << "***********************************************"<< endl;
  cout << X[1] << endl ;
  cout << "***********************************************"<< endl;
  cout << X[0] << endl;
  cout << "***********************************************"<< endl;




  return 0 ;
}
