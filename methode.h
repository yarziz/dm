#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <Dense>
#include "Sparse"

using Eigen::MatrixXd;
using Eigen::VectorXd;


Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> new_matrixx(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Hm);
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> new_matrix(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Hm);
Eigen::VectorXd GPO(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A, double eps);
std::vector<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> Arnoldi(Eigen::VectorXd r, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A);
Eigen::VectorXd cholesky_resolution(Eigen::MatrixXd A, Eigen::VectorXd b);
Eigen::VectorXd GMRes(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A, double eps);
Eigen::VectorXd GradienConjugue(Eigen::VectorXd x0, Eigen::VectorXd b, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A, int kmax,double epsilon);
Eigen::VectorXd Residuminimum(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A, double eps);
Eigen::VectorXd FOM(Eigen::VectorXd x0, Eigen::VectorXd b, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A, int kmax,double epsilon);
std::vector<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>>  qr_decomposition(Eigen::MatrixXd A);
Eigen::VectorXd qr_resolution_fom(Eigen::MatrixXd A, Eigen::VectorXd b);
Eigen::VectorXd romove_vector(Eigen::VectorXd b);
Eigen::VectorXd resolution_Gm(Eigen::MatrixXd R, Eigen::VectorXd q);
