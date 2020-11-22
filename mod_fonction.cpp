#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
//#include "LU.h"
#include "mod_fonction.h"

using namespace std;

void affichage_matrice(vector<vector<double>> kk,int n){
  string a="";
  for(int j=0;j<n+1;j++){
    a="";
    for(int i=0;i<n+1;i++){
      a=a+" "+to_string(kk[j][i]);
    }
    cout<<a<<endl;
  }
}


void affichage_vector(vector<double> x,int n){
  for(int j=0;j<n+1;j++){
       
    cout<<x[j]<<endl;
  }

}
