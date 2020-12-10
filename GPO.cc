
#include <fstream>
#include <math.h>
#include <iostream>
#include <Dense>
#include "Sparse"

using Eigen::MatrixXd;
using Eigen::VectorXd;



using namespace std;
using namespace Eigen;
Eigen::VectorXd GPO(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double, Dynamic, Dynamic> A, double eps)
{
  int n = x0.size() ;
  A.resize(n,n) ;
  Eigen::VectorXd r, X ;
  r.resize(n);
  X.resize(n) ;
  X=x0 ;
  r = b - A*x0 ;
  int k = 0 ;
  while ((sqrt(r.dot(r)) > eps) && (k <= kmax))
  {
    Eigen::VectorXd z ;
    z.resize(n);
    z = A*r ;
    double alpha = (r.dot(r))/(r.dot(z));
    X += alpha*r ;
    r = r - alpha*z ;
    k += 1 ;
  }
  return X ;

  if (k > kmax)
  {
    cout << "Tolerance non atteinte : " << sqrt(r.dot(r)) << endl ;
  }
}
Eigen::VectorXd Residuminimum(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double, Dynamic, Dynamic> A, double eps)
{
  int n = x0.size() ;
  A.resize(n,n) ;
  Eigen::VectorXd r, X ;
  r.resize(n);
  X.resize(n) ;
  X=x0 ;
  r = b - A*x0 ;
  int k = 0 ;
  while ((sqrt(r.dot(r))>eps) && (k<=kmax))
  {
    Eigen::VectorXd z ;
    z.resize(n);
    z = A*r ;
    double alpha = (r.dot(z))/(z.dot(z)) ;
    X += alpha*r ;
    r = r - alpha*z ;
    k+=1;
  }
  return X ;
  if (k > kmax)
  {
    cout << "Tolerance non atteinte : " << sqrt(r.dot(r)) << endl ;
  }
}



Eigen::MatrixXd cholesky_decomposition(Eigen::MatrixXd A)
{


  int i,j,k,n;
  double s,p;
  n=A.cols();
  Eigen::SparseMatrix<double> Lm(n,n);
  Eigen::Matrix<double, Dynamic, Dynamic> L;
  L.resize(n,n);

  L = MatrixXd(Lm);
  for (i=0;i<n;i++)
  {
    s=0;
    for (k=0;k<i;k++)
    {
      s=s+pow(L(i,k),2);
    }
    p=A(i,i)-s;



    L(i,i)=sqrt(p);

    for(j=i+1;j<n;j++)
    {
      s=0;
      for (k=0;k<i;k++)
      {
        s=s+L(i,k)*L(j,k);
      }
      L(j,i)=(A(j,i)-s)/L(i,i);
    }
  }
  return L ;
}
Eigen::VectorXd cholesky_resolution(Eigen::MatrixXd A, Eigen::VectorXd b)
{
  Eigen::MatrixXd L;
  Eigen::VectorXd y, x;
  int n = A.cols();
  int i,j,k ;
  double s, p ;
  L = cholesky_decomposition(A) ;
  x.resize(n);
  y.resize(n);
  for(i=0;i<n;i++)
  {
    s=0;
    for(j=0;j<i;j++)
    {
      s=s+L(i,j)*y(j);
    }
    y(i)=(b(i)-s)/L(i,i);
  }

  for(i=n-1;i>=0;i--)
  {
    s=0;
    for(j=i+1;j<n;j++) s=s+L(j,i)*x(j);
    x(i)=(y(i)-s)/L(i,i);
  }
  return x;
}

vector<Eigen::Matrix<double, Dynamic, Dynamic>> Arnoldi(Eigen::VectorXd r, Eigen::Matrix<double, Dynamic, Dynamic> A)
{
  Eigen::Matrix<double, Dynamic, Dynamic> Vm, Hm ;
  Eigen::VectorXd v, s, z;
  vector<Eigen::Matrix<double, Dynamic, Dynamic>> X ;
  int m = r.size();
  Eigen::SparseMatrix<double> H(m,m) ;
  Vm.resize(m,m);
  v = (1./sqrt(r.dot(r)))*r;
  Vm.col(0) = v ;
  for (int j = 0 ; j < m ; ++j )
  {
    for (int i = 0 ; i < j ; ++i )
    {
      H.insert(i,j) = Vm.col(i).dot(A*Vm.col(j));
    }
    s.setZero();
    for (int p = 0 ; p < m ; ++p)
    {
      s+=H.coeffRef(p,j)*Vm.col(p);
    }
    z = A*Vm.col(j) - s ;
    H.coeffRef(j+1,j) = sqrt(z.dot(z));
    if (H.coeffRef(j+1,j) == 0)
    {
      break;
    }
    Vm.col(j+1) = (1/H.coeffRef(j+1,j))*z;
  }
  Hm = MatrixXd(H);
  X.push_back(Hm);
  X.push_back(Vm);

  return X;
}

Eigen::VectorXd GMRes(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double, Dynamic, Dynamic> A, double eps)
{
  int n = x0.size() ;
  A.resize(n,n) ;
  Eigen::VectorXd r, X, e1,y ;
  vector<Eigen::Matrix<double, Dynamic, Dynamic>> W;
  r.resize(n);
  X.resize(n) ;
  e1.resize(n) ;
  for (int i=1 ; i<n ; ++i)
  {
    e1(i)=0;
  }

  e1(0)=1;
  X=x0 ;
  r = b - A*x0 ;
  double Beta = sqrt(r.dot(r)) ;
  int k = 0;
  while ((Beta > eps) && (k<=kmax))
  {
    W=Arnoldi(r,A) ; //retourne un vecteur de deux elements Hm et Vm
    y = cholesky_resolution(W[0].transpose()*W[0],Beta*W[0].transpose()*e1); //Moindres carrés
    X = X + W[1]*y ;
    r = Beta*e1-W[0]*y;
    Beta = sqrt(r.dot(r)) ;
    k += 1 ;
  }
  return X;
  if (k > kmax)
  {
    cout << "Tolerance non atteinte : " << sqrt(r.dot(r)) << endl ;
  }




}









int main()
{
  Matrix<double, 4, 4> A, B ;
  A << 4, -6, 8, 2,                           // Initialize A. The elements can also be
  -6, 10, 15, -3,                                // matrices.
  8, -15, 26, -1,
  2, -3, -1, 62;
  B=cholesky_decomposition(A);
  cout << B << endl;
  cout << B*B.transpose() << endl;
  return 0 ;
}
