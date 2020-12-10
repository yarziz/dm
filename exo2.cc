#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <Eigen>
#include <vector>
#include "methode.h"

using namespace std;
using namespace Eigen;




double somme(Matrix<double,Dynamic,Dynamic> a,int n){
  double ligne=0.;
  double max=0.;
  for(int i=0;i<n;i++){
    ligne=0.;
    for(int j=0;j<n;j++){
      ligne=ligne+abs(a(i,j));
    }
    if(ligne>max){
      max=ligne;
    }
  }
  return max;
}




Matrix<double,Dynamic,Dynamic> gene_matrice(int n){
  Matrix<double,Dynamic,Dynamic> b;
  Matrix<double,Dynamic,Dynamic> b_bt;
  Matrix<double,Dynamic,Dynamic> a;
  //double a;
  b.resize(n,n);
  b_bt.resize(n,n);
  a.resize(n,n);
  srand(time(0));
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      b(i,j)= ((double) rand())/(double) RAND_MAX;
    }
  }
  b_bt=b.transpose()*b;
  b_bt=somme(b_bt,n)*b_bt;
  a=b_bt;
  for(int i=0;i<n;i++){
    a(i,i)=a(i,i)+1;
  }
  return a;
}






int main(){
  int n;
  Eigen::Matrix<double, Dynamic, Dynamic>;
  Eigen::VectorXd b,X0;
  cout<<"donner n"<<endl;
  cin>>n;
  b.resize(n);
  X0.resize(n);
  X0(0)=1;
  
 
  GPO(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double, Dynamic, Dynamic> A, double eps)
  cout<<gene_matrice(n)<<endl;
  
 
  return 0;
}
