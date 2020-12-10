#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <Dense>
#include "Sparse"

using Eigen::MatrixXd;
using Eigen::VectorXd;

Eigen::VectorXd GPO(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A, double eps);
std::vector<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> Arnoldi(Eigen::VectorXd r, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A);
Eigen::VectorXd cholesky_resolution(Eigen::MatrixXd A, Eigen::VectorXd b);
Eigen::VectorXd GMRes(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A, double eps);
Eigen::VectorXd GradienConjugue(Eigen::VectorXd x0, Eigen::VectorXd b, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A, int kmax,double epsilon);
Eigen::VectorXd Residuminimum(Eigen::VectorXd x0, Eigen::VectorXd b, int kmax, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A, double eps);
//std::vector<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> QR(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> a);
