




#include "Sparse"
#include "iostream"
#include <vector>
#include "fstream"
#include "algorithm"

using namespace std;
using namespace Eigen;

int main() {
  ifstream fin("s3rmt3m3.mtx");
  int M, N, L;
  //Ignore headers and comments
  while (fin.peek() == '%')
  fin.ignore(2048, '\n');

  fin >> M >> N >> L;
  // Eigen::setNbThreads(8);

  MatrixXd Matrix(M, N);
  Matrix.setZero();


  for (int i = 0; i < M; ++i) {
    for (int j = 0 ; j<N ; ++j){


      int m, n;
      double data;
      fin >> m >> n >> data;
      Matrix(m-1,n-1) = data;
    }
  }
  cout << Matrix.rows() << endl;
  fin.close();





  return 0;


}
