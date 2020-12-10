#include <fstream>
#include <math.h>
#include <iostream>
#include <Dense>
#include "Sparse"

using Eigen::MatrixXd;
using Eigen::VectorXd;



using namespace std;
using namespace Eigen;


Eigen::VectorXd GradienConjugue(Eigen::VectorXd x0, Eigen::VectorXd b, Eigen::Matrix<double, Dynamic, Dynamic> A, int kmax){
  int n=x0.size();
  int k=0;
  double epsilon;
  double alpha;
  Eigen::VectorXd r,rn,p,z,x;
  double beta=0,gamma=0;
  A.resize(n,n);
  r.resize(n);
  rn.resize(n);
  p.resize(n);
  z.resize(n);
  x.resize(n);
  x=x0;
  r=b-A*x0;
  p=r;
  beta=sqrt(r.dot(r));
  while((beta>epsilon)&&(k<kmax+1)){
    z=A*p;
    alpha=r.dot(r)/z.dot(p);
    x=x+alpha*p;
    rn=r-alpha*z;
    gamma=(rn).dot(rn)/r.dot(r);
    p=rn+gamma*p;
    beta=sqrt(r.dot(r));
    k=k+1;
    r=rn;
  }
  if(k>kmax){
    std::cout<<"Tolerance non atteinte: "<<beta<<std::endl; 
  }
  return x;
}

/*
  Eigen::VectorXd FOM(Eigen::VectorXd x0, Eigen::VectorXd b, Eigen::MatrixXd A, int kmax){
  int n=x0.size();
  int k=0;
  double epsilon;
  double alpha;
  Eigen::VectorXd r,p,z,x;
  double beta=0,gamma=0;
  A.resize(n,n);
  r.resize(n);
  p.resize(n);
  z.resize(n);
  x.resize(n);
  x=x0;
  r=b-A*x0;
  p=r;
  beta=sqrt(r.dot(r));
  while((beta>epsilon)&&(k<kmax+1)){
  z=A*p;
  alpha=r.dot(r)/z.dot(p);
  x=x-alpha*p;
  beta=sqrt(r.dot(r));
  gamma=(r-alpha*z).dot(r-alpha*z)/r.dot(r);
  r=r-alpha*z;
  p=r+gamma*p;
  k=k+1;
  }
  if(beta>epsilon){
  std::cout<<"Tolerance non atteinte: "<<beta<<std::endl; 
  }else{
  return x;
  }
  }
*/

vector<Eigen::Matrix<double, Dynamic, Dynamic>> QR(Eigen::Matrix<double, Dynamic, Dynamic> a){
  n=a.rows();
  m=a.cols();
  I=MatrixXd::Identity(n,n);
  Eigen::Matrix<double, Dynamic, Dynamic> p;
  Eigen::VectorXd w, s, z;
  double beta=0;
  double s=0;
  
  for(int k=0;k<m;k++){
    if(k>0){
      r.col(k)=
	}
    s=0;
    for(int i=k;i<n;i++){
      s=s+a(i,k)*a(i,k);
    }
    beta=sign(a(k,k))*sqrt(s);
    for(int j=0;j<n;j++){
      if(j<k){
	z(j)=0;
      }
      if(j==k){
	z(j)=beta+a(j,j);
      }
      if(j>k){
	z(j)=a(j,k);
      }
    }
    w=(1/sqrt(z.dot(z)))*z;
    p=I-2*w*transpose(w);
    r.col(k+1)=P*r.
    


  }





}



int main(){
  return 0; 
}
