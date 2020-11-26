
#include <fstream>

#include <iostream>
#include <Dense>

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
      H.insert(i+1,j+1) = Vm.col(j).dot(A*Vm.col(i));
    }
    for (int p = 1 ; p < m+1 ; ++p)
    {
      s+=H.coeffRef(p,j+1)*Vm.col(p-1);
    }
    z = A*Vm.col(j) - s ;
    H.coeffRef(j+2,j+1) = sqrt(z.dot(z));
    if (H.coeffRef(j+2,j+1) == 0)
    {
      break;
    }
    Vm.col(j+1) = (1/H.coeffRef(j+2,j+1))*z;
  }
  Hm = MatrixXd(H);
  X.push_back(Hm);
  X.push_back(Vm);

  return X;
}








int main()
{
  return 0 ;
}
