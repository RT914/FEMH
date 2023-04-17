#include <Eigen/Dense>
#include <stdio.h>
#include <iostream>

#include "Matrix.h"

// 3�~3�̍s��
double det(Eigen::Matrix3d M) {

	// �]���q�W�J�𗘗p�����s�񎮌v�Z
	// 
	// (0, 0)�̏��s��
	Eigen::Matrix2d m1; 
	m1 << M(1, 1), M(1, 2), 
		M(2,1), M(2, 2);
	//std::cout << "m1 : " << m1 << "\n";

	// (0, 1)�̏��s��
	Eigen::Matrix2d m2;
	m2 << M(1, 0), M(1, 2), 
		M(2, 0), M(2, 2);
	//std::cout << "m1 : " << m2 << "\n";

	// (0, 2)�̏��s��
	Eigen::Matrix2d m3;
	m3 << M(1, 0), M(1, 1), 
		M(2, 0), M(2, 1);
	//std::cout << "m1 : " << m3 << "\n";

	// ���s��̍s��
	double det_m1 = (m1(0, 0) * m1(1, 1) - m1(0, 1) * m1(1, 0));
	//std::cout << "determinant_m1 : " << det_m1 << "\n";
	double det_m2 = (m2(0, 0) * m2(1, 1) - m2(0, 1) * m2(1, 0)) * (-1);
	//std::cout << "determinant_m2 : " << det(M) << "\n";
	double det_m3 = (m3(0, 0) * m3(1, 1) - m3(0, 1) * m3(1, 0));
	//std::cout << "determinant_m3 : " << det(M) << "\n";

	double det = M(0, 0) * det_m1 + M(0, 1) * det_m2 + M(0, 2) * det_m3;

	return det;
}