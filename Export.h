#ifndef __EXPORT_H__
#define __EXPORT_H__

#include <Eigen/Dense>

void export_CSV(double MatrixQForX[3][3][3][3][3][3], double MatrixQForY[3][3][3][3][3][3], double MatrixQForZ[3][3][3][3][3][3]);
void export_CSV_X(double MatrixQForX[3][3][3][3][3][3]);
void export_CSV_Y(double MatrixQForY[3][3][3][3][3][3]);
void export_CSV_Z(double MatrixQForZ[3][3][3][3][3][3]);
void export2_CSV(double VectorbPreparation[3][3][3]);
void export_Q(Eigen::MatrixXd m, std::string str);
void export_loop_times(std::vector<int> loop_times, std::string str, int num);
#endif