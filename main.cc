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



//somme est une fonction qui calcule alpha de la question 2
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

//------------------------------------------------------------------------------------//

//gene_matrice_tri est une fonction qui donne une matrice symétrique tridiagonal définie positive
Matrix<double,Dynamic,Dynamic> gene_matrice_tri(Matrix<double,Dynamic,Dynamic> a,int n){
  Matrix<double,Dynamic,Dynamic> tri; 
  a.resize(n,n);
  tri.resize(n,n);
  tri.setZero();
  for(int j=0;j<n;j++){
    tri(j,j)=a(j,j)+pow(n,5);
  }

  for(int i=0;i<n-1;i++){
    tri(i+1,i)=a(i+1,i);
    tri(i,i+1)=a(i,i+1);
  }
  return tri;
}

//-----------------------------------------------------------------------------------------//

//gene_matrice est la fonction qui qui donne la matrice aléatoire de taille n*n
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



//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


int main(){
  int n;
  int kmax=100000;
  double eps=0.00000001;
  Eigen::Matrix<double, Dynamic, Dynamic> a;
  Eigen::VectorXd b,x0,x;
  cout<<"donner n"<<endl;
  cin>>n;
  b.resize(n);
  x0.resize(n);
  //vecteur b de la question 2
  for(int i=0;i<n;i++){
    b(i)=1;
  }
  //x initiale
  x0(0)=1;
  
  a=gene_matrice(n);
  x= GMRes(x0,b,kmax,a,eps);
  //x=GPO(x0,b,kmax,a,eps);
  //x=Residuminimum(x0,b,kmax,a,eps);
  //x=GradienConjugue(x0,b,a,kmax,eps);
  //x=FOM(x0,b,a,kmax,eps);
  //cout<<"--------------------------"<<endl;
  //cout<<b<<endl;
  //cout<<"--------------------------"<<endl;
  //cout<<romove_vector(b)<<endl;
  cout<<"--------------------------"<<endl;
  cout<<(a*x-b).maxCoeff()<<endl;
  
  
 
 
  return 0;
}
