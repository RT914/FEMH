#ifndef __NEWTONRAPHSONMETHOD_H__
#define __NEWTONRAPHSONMETHOD_H__

#include "Square.h"
#include <Eigen/Dense>

Eigen::VectorXd Newton(Square square);
Eigen::VectorXd Newton_one(Square square);
void Newton_copy(Square square);
int Newton_loop(Square square, int seed);
Eigen::MatrixXd calMatrixQ(Square square, Eigen::VectorXd v);
void calMatrixQPreparation(Square square);
Eigen::VectorXd calVectorb(Square square, Eigen::VectorXd v);
void calVectorbPreparation(Square square);
Eigen::MatrixXd calMatrixR(Eigen::MatrixXd M);
Eigen::VectorXd calVectorc(Eigen::VectorXd V, Eigen::MatrixXd M);
int GridToFlat(Eigen::Vector3i grid_index);
Eigen::Vector3i FlatToGrid(int flat_index);
double calNormalize(Eigen::VectorXd v);
void showDetail(double x, double y, double z, Square s, Eigen::VectorXd v);
double N(double x);
double N_dash(double x);
void calCheckMatrixPreparation(Square square);
Eigen::MatrixXd calCheckMatrix(Square square);

void calVectorb_part2_Preparation(Square square);
Eigen::VectorXd calVectorb_part2(Square square, Eigen::VectorXd v);

#endif
