#include <fstream>
#include <iostream>
#include <Dense>
#include <vector>
#include "Sparse"
#include "methode.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;



using namespace std;
using namespace Eigen;



vector<Eigen::Matrix<double, Dynamic, Dynamic>>  qr_decomposition(Eigen::MatrixXd A)
{
  vector<Eigen::Matrix<double, Dynamic, Dynamic>> X ;
  Eigen::Matrix<double, Dynamic, Dynamic>  R ;
  int n = A.rows() ;
  int m = A.cols() ;
  double s(0.) ;
  for (int k = 0 ; k<n ; ++k)
  {
    for (int i = 0 ; i<k-1, ++i)
    {
      s=0 ;
      for (int j = 0 ; j<m, ++j)
      {
        s += A(j,i)*A(j,k);
      }
      R(i,k) = s;
      for (j = 0 ; j<m ; ++j)
      {
        A(j,k)=A(j,k) - A(j,i)*R(i,k) ;
      }
      s=0 ;
      for (j = 0 ; j<m ; ++j)
      {
        s += pow(A(j,k),2) ;
      }
      R(k,k)=sqrt(s) ;
      for (j = 0 ; j<m ; ++j)
      {
        A(j,k)=A(j,k)/R(k,k) ;
      }
    }
  }
  X.push_back(A);
  X.push_back(R);
  return X;
}
