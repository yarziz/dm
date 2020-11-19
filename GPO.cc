
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









int main()
{
  return 0 ;
}
