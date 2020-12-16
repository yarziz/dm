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

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();

  if( rowToRemove < numRows )
  matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

  matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;

  if( colToRemove < numCols )
  matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

  matrix.conservativeResize(numRows,numCols);
}

Eigen::VectorXd romove_vector(Eigen::VectorXd b){
  int n=b.rows();
  VectorXd z;
  z=b;
  z.resize(n-1);
  return z;
}



Eigen::Matrix<double, Dynamic, Dynamic> new_matrix(Eigen::Matrix<double, Dynamic, Dynamic> Hm){
  Eigen::Matrix<double, Dynamic, Dynamic> Hm_;
  int n(0),m(0);
  n=Hm.rows();
  m=Hm.cols();
  
  if(n==m+1){
    // cout<<"ok"<<endl;
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
    // cout<<"ok"<<endl;
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
  double alpha=0, norm = 0;
  int n = x0.size() ;
  Eigen::VectorXd r, X ;
  A.resize(n,n) ;
  z.resize(n);
  r.resize(n);
  X.resize(n) ;
  X=x0 ;
  r = b - A*x0 ;
  norm = sqrt(r.dot(r));
  int k = 0 ;
  ofstream mon_flux;
  string name_file("gporesidu.txt");
  mon_flux.open(name_file, ios::out);
  if(mon_flux)
  {
    mon_flux << k+1 << " " << norm << endl;
    while ((norm > eps) && (k <= kmax))
    {
      z = A*r ;
      alpha = (r.dot(r))/(r.dot(z));
      X += alpha*r ;
      r = r - alpha*z ;
      norm = sqrt(r.dot(r));
      k += 1 ;
      mon_flux << k+1 << " " << norm << endl;
      
    }
  }
  else // Renvoie un message d’erreur si ce n’est pas le cas
  {
    cout << "ERREUR: Impossible d’ouvrir le fichier." << endl;
  }
  mon_flux.close();
  if (k > kmax)
  {
    cout << "Tolerance non atteinte : " << sqrt(r.dot(r)) << endl ;
  }
  return X ;
}
 

Eigen::VectorXd Residuminimum(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double, Dynamic, Dynamic> A, double eps)
{
  Eigen::VectorXd z ;
  double alpha=0, norm(0.);
  int n = x0.size() ;
  A.resize(n,n) ;
  z.resize(n);
  Eigen::VectorXd r, X ;
  r.resize(n);
  X.resize(n) ;
  X=x0 ;
  r = b - A*x0 ;
  norm=sqrt(r.dot(r));
  int k = 0 ;
  ofstream mon_flux;
  string name_file("Residuminimumresidu.txt");
  mon_flux.open(name_file, ios::out);
  if(mon_flux)
  {
    mon_flux << k+1 << " " << norm << endl;
    while ((norm>eps) && (k<=kmax))
    {
      z = A*r ;
      alpha = (r.dot(z))/(z.dot(z)) ;
      X += alpha*r ;
      r = r - alpha*z ;
      norm=sqrt(r.dot(r));
      k+=1;
      mon_flux << k+1 << " " << norm << endl;
     
    }
  }
  else // Renvoie un message d’erreur si ce n’est pas le cas
  {
    cout << "ERREUR: Impossible d’ouvrir le fichier." << endl;
  }
  mon_flux.close();

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
  Eigen::VectorXd r, X, e1,y,q,qq ;
  vector<Eigen::Matrix<double, Dynamic, Dynamic>> W,WW;
  Eigen::Matrix<double, Dynamic, Dynamic> Q,R,Vm,Rm;
  double g=0;
  r.resize(n);
  X.resize(n) ;
  e1.resize(n+1) ;
  Q.resize(n+1,n+1);
  R.resize(n+1,n);
  Rm.resize(n,n);
  Vm.resize(n,n);
  qq.resize(n);
  q.resize(n+1);
  for (int i=1 ; i<n+1 ; ++i)
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
      WW=qr_decomposition(W[0]);
      Q=WW[0];
      R=WW[1];
      q=Beta*(Q.transpose()*e1);
      g=q(n);
      Rm=new_matrix(R);
      qq=romove_vector(q);
      y=resolution_Gm(Rm,qq);
      //removeRow(W[0],W[0].rows()-1);
      //removeColumn(W[1],W[1].cols()-1);
      //y = cholesky_resolution((W[0].transpose())*W[0],Beta*W[0].transpose()*e1); //Moindres carrés
      Vm=new_matrixx(W[1]);
      X = X + Vm*y ;
      r = Beta*W[1]*((Q.transpose()).col(n));
      Beta = sqrt(r.dot(r));
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
  ofstream mon_flux;
  string name_file("GradienConjugeResidu.txt");
  mon_flux.open(name_file, ios::out);
  if(mon_flux)
  {
    mon_flux << k+1 << " " << sqrt((b-A*x).dot((b-A*x))) << endl;
    while((beta>epsilon)&&(k<kmax+1)){
      z=A*p;
      alpha=r.dot(r)/z.dot(p);
      x=x+alpha*p;
      rn=r-alpha*z;
      gamma=(rn).dot(rn)/r.dot(r);
      p=rn+gamma*p;
      beta=sqrt(r.dot(r));
      k=k+1;
      mon_flux << k+1 << " " << sqrt((b-A*x).dot((b-A*x))) << endl;
      r=rn;
    }
  }
  else // Renvoie un message d’erreur si ce n’est pas le cas
  {
    cout << "ERREUR: Impossible d’ouvrir le fichier." << endl;
  }
  mon_flux.close();
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
    //y=cholesky_resolution(W0,beta*e1);
    y=qr_resolution_fom(W0,beta*e1);
      //cout<<y<<endl;
    x=x+W1*y;
    r=-W[0](n,n-1)*y(n-1)*(W[1].col(n));
    beta=sqrt(r.dot(r));
    k=k+1;
  }
  
  return x;
}







vector<Eigen::Matrix<double, Dynamic, Dynamic>>  qr_decomposition(Eigen::MatrixXd A)
{
  int n = A.rows();
  int m = A.cols();
  int c(0);
  double beta;
  Eigen::Matrix<double, Dynamic, Dynamic>  H,W,Q,Rbar,R,S,B,Hbar,I;
  VectorXd z,w,x,e;
  vector<Eigen::Matrix<double, Dynamic, Dynamic>> P, X ;
  R=A;
  if (n<m)
  {
    c=n;
  }
  else
  {
    c=m;
  }
  Q.resize(n,n);
  for (int k = 0 ; k<c ; ++k)
  {
    H = MatrixXd::Identity(n,n);
    e.resize(n-k);
    e.setZero();
    e(0) = 1;
    I = MatrixXd::Identity(n-k,n-k);
    z.resize(n-k) ;
    w.resize(n-k);
    z.setZero() ;
    x = R.col(k).segment(k,n-k);
    z = x + sqrt(x.dot(x))*e ;
    w = (1/(sqrt(z.dot(z))))*z ;
    Hbar = I - 2*w*(w.transpose()) ;
    Rbar = Hbar*R.block(k,k,n-k,m-k);
    R.block(k,k,n-k,m-k) = Rbar;
    H.block(k,k,n-k,n-k) = Hbar;
    P.push_back(H);
  }
  Q=P[0];
  for (int i = 1 ; i<c ; ++i)
  {
    Q=Q*P[i] ;
  }
  X.push_back(Q);
  X.push_back(R);
  return X;
}













Eigen::VectorXd qr_resolution_fom(Eigen::MatrixXd A, Eigen::VectorXd b)
{
  vector<Eigen::Matrix<double, Dynamic, Dynamic>> QR;
  //Eigen::Matrix<double, Dynamic, Dynamic> Q;
  Eigen::Matrix<double, Dynamic, Dynamic> R;
  Eigen::VectorXd x,q;
  int n = A.cols();
  int i,j,k ;
  double s, p;
  b.resize(n);
  A.resize(n,n);
  q.resize(n);
  R.resize(n,n);
  QR=qr_decomposition(A) ;
  q=(QR[0].transpose())*b;
  R=QR[1];
  x.resize(n);
  x(n-1)=q(n-1)/R(n-1,n-1);
  for(int i=n-2;i>=0;i--){
    s=0;
    for(int j=i+1;j<n;j++){
      s=s+R(i,j)*x(j);
    }
    x(i)=(1/R(i,i))*(q(i)-s);
  }
  return x;

}




Eigen::VectorXd resolution_Gm(Eigen::MatrixXd R, Eigen::VectorXd q){
  Eigen::VectorXd x;
  int n = R.cols();
  double s;
  x.resize(n);
  x(n-1)=q(n-1)/R(n-1,n-1);
  for(int i=n-2;i>=0;i--){
    s=0;
    for(int j=i+1;j<n;j++){
      s=s+R(i,j)*x(j);
    }
    x(i)=(1/R(i,i))*(q(i)-s);
  }
  return x;
}
