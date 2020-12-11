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


Eigen::Matrix<double, Dynamic, Dynamic> new_matrix(Eigen::Matrix<double, Dynamic, Dynamic> Hm){
  Eigen::Matrix<double, Dynamic, Dynamic> Hm_;
  int n(0),m(0);
  n=Hm.rows();
  m=Hm.cols();
  
  if(n==m+1){
    cout<<"ok"<<endl;
  }
  Hm_.resize(m,m);
  for(int k=0;k<n-1;k++){
    Hm_.row(k)=Hm.row(k);
  }
 
  
  for(int i=0;i<m-2;i++){
    for(int j=2+i;j<m;j++){
      Hm_(i,j)=0;
    }
  }
  return Hm_;
}

Eigen::Matrix<double, Dynamic, Dynamic> new_matrixx(Eigen::Matrix<double, Dynamic, Dynamic> Hm){
  Eigen::Matrix<double, Dynamic, Dynamic> Hm_;
  int n(0),m(0);
  n=Hm.rows();
  m=Hm.cols();
  
  if(n==m-1){
    cout<<"ok"<<endl;
  }
  Hm_.resize(n,n);
  for(int k=0;k<m-1;k++){
    Hm_.col(k)=Hm.col(k);
  }
  return Hm_;
}





Eigen::VectorXd GPO(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double, Dynamic, Dynamic> A, double eps)
{
  Eigen::VectorXd z ;
  double alpha=0;
  int n = x0.size() ;
  Eigen::VectorXd r, X ;
  A.resize(n,n) ;
  z.resize(n);
  r.resize(n);
  X.resize(n) ;
  X=x0 ;
  r = b - A*x0 ;
  int k = 0 ;
  while ((sqrt(r.dot(r)) > eps) && (k <= kmax))
    {
      z = A*r ;
      alpha = (r.dot(r))/(r.dot(z));
      X += alpha*r ;
      r = r - alpha*z ;
      k += 1 ;
    }

  if (k > kmax)
    {
      cout << "Tolerance non atteinte : " << sqrt(r.dot(r)) << endl ;
    }
  
  return X ;

 
}
Eigen::VectorXd Residuminimum(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double, Dynamic, Dynamic> A, double eps)
{
  Eigen::VectorXd z ;
  double alpha=0;
  int n = x0.size() ;
  A.resize(n,n) ;
  z.resize(n);  
  Eigen::VectorXd r, X ;
  r.resize(n);
  X.resize(n) ;
  X=x0 ;
  r = b - A*x0 ;
  int k = 0 ;
  while ((sqrt(r.dot(r))>eps) && (k<=kmax))
    {
      z = A*r ;
      alpha = (r.dot(z))/(z.dot(z)) ;
      X += alpha*r ;
      r = r - alpha*z ;
      k+=1;
    }

  if (k > kmax)
    {
      cout << "Tolerance non atteinte : " << sqrt(r.dot(r)) << endl ;
    }
  return X ;
  
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
  //Eigen::SparseMatrix<double> H(m,m) ;
  Vm.resize(m,m+1);
  Hm.resize(m+1,m);
  Hm.setZero();
  Vm.setZero();
  v.resize(m);
  s.resize(m);
  z.resize(m);
  v = (1./sqrt(r.dot(r)))*r;
  Vm.col(0) = v ;
  
  for (int j = 0 ; j < m ; ++j )
    {
      for (int i = 0 ; i < j+1 ; ++i )
	{
	  Hm(i,j) = (Vm.col(i)).dot(A*(Vm.col(j)));
	}
      s.setZero();
      for (int p = 0 ; p < j+1 ; ++p)
	{
	  s+=Hm(p,j)*(Vm.col(p));
	}
      z = A*(Vm.col(j)) - s ;
      Hm(j+1,j) = sqrt(z.dot(z));
      if (Hm(j+1,j) == 0)
	{
	  break;
	}
      Vm.col(j+1) = (1/Hm(j+1,j))*z;
    }
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
      y = cholesky_resolution((W[0].transpose())*W[0],Beta*W[0].transpose()*e1); //Moindres carrés
      X = X + W[1]*y ;
      r = Beta*e1-W[0]*y;
      Beta = sqrt(r.dot(r)) ;
      k += 1 ;
    }
  if (k > kmax)
    {
      cout << "Tolerance non atteinte : " << sqrt(r.dot(r)) << endl ;
    }
  return X;
}

Eigen::VectorXd GradienConjugue(Eigen::VectorXd x0, Eigen::VectorXd b, Eigen::Matrix<double, Dynamic, Dynamic> A, int kmax, double epsilon){
  int n=x0.size();
  int k=0;
 
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


Eigen::VectorXd FOM(Eigen::VectorXd x0, Eigen::VectorXd b, Eigen::Matrix<double, Dynamic, Dynamic> A, int kmax,double epsilon){
  int n=x0.size();
  int k=0;
  Eigen::VectorXd r,y,e1,x;
  vector<Eigen::Matrix<double, Dynamic, Dynamic>> W;
  Eigen::Matrix<double, Dynamic, Dynamic> W0,W1;
  double beta=0,gamma=0;
  
  A.resize(n,n);
  W0.resize(n,n);
  W1.resize(n,n);
  r.resize(n);
  y.resize(n);
  //p.resize(n);
  //z.resize(n);
  x.resize(n);
  x=x0;
  r=b-A*x0;
  //p=r;
  beta=sqrt(r.dot(r));
  e1.resize(n) ;
  for (int i=0 ; i<n ; ++i)
    {
      e1(i)=0;
    }

  e1(0)=1;


  while((beta>epsilon)&&(k<kmax+1)){
    W=Arnoldi(r,A);
    W0=new_matrix(W[0]);
    W1=new_matrixx(W[1]);
    //cout<<W0<<endl;
    y=cholesky_resolution(W0,beta*e1);
    cout<<y<<endl;
    x=x+W1*y;
    r=-W[0](n,n-1)*y(n-1)*(W[1].col(n));
    beta=sqrt(r.dot(r));
    k=k+1;
  }
  
  return x;
}


/*
  vector<Eigen::Matrix<double, Dynamic, Dynamic>> QR(Eigen::Matrix<double, Dynamic, Dynamic> a){
  /*
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
  }*/


