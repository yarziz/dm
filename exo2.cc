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
  int kmax=10000;
  double eps=0.00001;
  Eigen::Matrix<double, Dynamic, Dynamic> a;
  Eigen::VectorXd b,x0,x;
  cout<<"donner n"<<endl;
  cin>>n;
  b.resize(n);
  x0.resize(n);
  for(int i=0;i<n;i++){
    b(i)=1;
  }
  x0(0)=1;
  a=gene_matrice(n);
  cout<<a<<endl;
  
  x=GPO(x0,b,kmax,a,eps);
  cout<<"--------------------------"<<endl;
  cout<<x<<endl;
  cout<<"--------------------------"<<endl;
  cout<<a*x<<endl;
  cout<<"--------------------------"<<endl;
  cout<<(a*x-b).maxCoeff()<<endl;
  
 
  return 0;
}
