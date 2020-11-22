#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <random>

using namespace std;

vector<vector<double>> transpose(vector<vector<double>> a,int n){
   vector<vector<double>> b(n,vector<double>(n,0));
   for(int i=0;i<n;i++){
     for(int j=0;j<n;j++){
       b[i][j]=a[j][i];
     }
   }
   return b;
}

double somme(vector<vector<double>> a,int n){
  double alpha=0;
   for(int i=0;i<n;i++){
     for(int j=0;j<n;j++){
       alpha=alpha+abs(a[i][j]);
     }
   }
   return alpha;
}




vector<vector<double>> gene_matrice(int n){
  vector<vector<double>> b(n,vector<double>(n,0));
  vector<vector<double>> b_bt(n,vector<double>(n,0));
  vector<vector<double>> a(n,vector<double>(n,0));
  //double a;
  srand(time(0));
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      b[i][j]= ((double) rand())/(double) RAND_MAX;
    }
  }
  b_bt=transpose(b,n)*b;
  b_bt=somme(b_bt,n)*b_bt;
  a=b_bt;
  for(int i=0;i<n;i++){
    a[i][i]=a[i][i]+1;
  }
}






int main(){
   int n;
   cout<<"donner n"<<endl;
   cin>>n;
   affichage_matrice(gene_matrice(n));
 
  return 0;
}
