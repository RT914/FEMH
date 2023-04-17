#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "Calculation.h"

Kernel_Array createKernelArray() {
	Eigen::MatrixXd array_i_x = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_i_y = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_i_z = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_j_x = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_j_y = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_j_z = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_k_x = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_k_y = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_k_z = Eigen::MatrixXd::Zero(4, 2);

	

	// original point grid
	Eigen::MatrixXd array_o = Eigen::MatrixXd::Zero(4, 2);

	// 原点固定のkernel
	// [1+x][1-x]
	array_o(1, 0) = 1;
	array_o(1, 1) = 1;
	array_o(2, 0) = 1;
	array_o(2, 1) = -1;

	Kernel_Array ka = Kernel_Array(array_i_x, array_i_y, array_i_z,
		array_j_x, array_j_y, array_j_z,
		array_k_x, array_k_y, array_k_z, array_o);

	return ka;
};

Kernel_Array setKernelArrayForX(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_ix, int gn_jx, int gn_kx) {

	// 値の初期化
	ka.array_i_x = Eigen::MatrixXd::Zero(4, 2);
	ka.array_j_x = Eigen::MatrixXd::Zero(4, 2);
	ka.array_k_x = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_ix == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_i_x(2, 0) = 1;
			ka.array_i_x(3, 0) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_x(2, 0) = 1 - gn_ix;
			ka.array_i_x(2, 1) = 1;
			ka.array_i_x(3, 0) = 1 + gn_ix;
			ka.array_i_x(3, 1) = -1;

		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_x(2, 0) = 1 - gn_ix;
			ka.array_i_x(2, 1) = 1;
			ka.array_i_x(3, 0) = 1 + gn_ix;
			ka.array_i_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ix == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_i_x(1, 0) = 1;
			ka.array_i_x(2, 0) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_x(1, 0) = 1 - gn_ix;
			ka.array_i_x(1, 1) = 1;
			ka.array_i_x(2, 0) = 1 + gn_ix;
			ka.array_i_x(2, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_x(1, 0) = 1 - gn_ix;
			ka.array_i_x(1, 1) = 1;
			ka.array_i_x(2, 0) = 1 + gn_ix;
			ka.array_i_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ix == -1)
	{
		if (frag_i == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_i_x(0, 0) = 1;
			ka.array_i_x(1, 0) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_x(0, 0) = 1 - gn_ix;
			ka.array_i_x(0, 1) = 1;
			ka.array_i_x(1, 0) = 1 + gn_ix;
			ka.array_i_x(1, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_x(0, 0) = 1 - gn_ix;
			ka.array_i_x(0, 1) = 1;
			ka.array_i_x(1, 0) = 1 + gn_ix;
			ka.array_i_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node ix value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jx == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_j_x(2, 0) = 1;
			ka.array_j_x(3, 0) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_x(2, 0) = 1 - gn_jx;
			ka.array_j_x(2, 1) = 1;
			ka.array_j_x(3, 0) = 1 + gn_jx;
			ka.array_j_x(3, 1) = -1;

		}
		else if (frag_j == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_x(2, 0) = 1 - gn_jx;
			ka.array_j_x(2, 1) = 1;
			ka.array_j_x(3, 0) = 1 + gn_jx;
			ka.array_j_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jx == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_j_x(1, 0) = 1;
			ka.array_j_x(2, 0) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_x(1, 0) = 1 - gn_jx;
			ka.array_j_x(1, 1) = 1;
			ka.array_j_x(2, 0) = 1 + gn_jx;
			ka.array_j_x(2, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_x(1, 0) = 1 - gn_jx;
			ka.array_j_x(1, 1) = 1;
			ka.array_j_x(2, 0) = 1 + gn_jx;
			ka.array_j_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jx == -1)
	{
		if (frag_j == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_j_x(0, 0) = 1;
			ka.array_j_x(1, 0) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_x(0, 0) = 1 - gn_jx;
			ka.array_j_x(0, 1) = 1;
			ka.array_j_x(1, 0) = 1 + gn_jx;
			ka.array_j_x(1, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_x(0, 0) = 1 - gn_jx;
			ka.array_j_x(0, 1) = 1;
			ka.array_j_x(1, 0) = 1 + gn_jx;
			ka.array_j_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node jx value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_kx == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_k_x(2, 0) = 1;
			ka.array_k_x(3, 0) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_x(2, 0) = 1 - gn_kx;
			ka.array_k_x(2, 1) = 1;
			ka.array_k_x(3, 0) = 1 + gn_kx;
			ka.array_k_x(3, 1) = -1;

		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_x(2, 0) = 1 - gn_kx;
			ka.array_k_x(2, 1) = 1;
			ka.array_k_x(3, 0) = 1 + gn_kx;
			ka.array_k_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kx == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_k_x(1, 0) = 1;
			ka.array_k_x(2, 0) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_x(1, 0) = 1 - gn_kx;
			ka.array_k_x(1, 1) = 1;
			ka.array_k_x(2, 0) = 1 + gn_kx;
			ka.array_k_x(2, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_x(1, 0) = 1 - gn_kx;
			ka.array_k_x(1, 1) = 1;
			ka.array_k_x(2, 0) = 1 + gn_kx;
			ka.array_k_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kx == -1)
	{
		if (frag_k == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_k_x(0, 0) = 1;
			ka.array_k_x(1, 0) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_x(0, 0) = 1 - gn_kx;
			ka.array_k_x(0, 1) = 1;
			ka.array_k_x(1, 0) = 1 + gn_kx;
			ka.array_k_x(1, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_x(0, 0) = 1 - gn_kx;
			ka.array_k_x(0, 1) = 1;
			ka.array_k_x(1, 0) = 1 + gn_kx;
			ka.array_k_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node kx value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setKernelArrayForY(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iy, int gn_jy, int gn_ky) {

	// 値の初期化
	ka.array_i_y = Eigen::MatrixXd::Zero(4, 2);
	ka.array_j_y = Eigen::MatrixXd::Zero(4, 2);
	ka.array_k_y = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_iy == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_y(2, 0) = 1 - gn_iy;
			ka.array_i_y(2, 1) = 1;
			ka.array_i_y(3, 0) = 1 + gn_iy;
			ka.array_i_y(3, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_i_y(2, 0) = 1;
			ka.array_i_y(3, 0) = -1;

		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_y(2, 0) = 1 - gn_iy;
			ka.array_i_y(2, 1) = 1;
			ka.array_i_y(3, 0) = 1 + gn_iy;
			ka.array_i_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iy == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_y(1, 0) = 1 - gn_iy;
			ka.array_i_y(1, 1) = 1;
			ka.array_i_y(2, 0) = 1 + gn_iy;
			ka.array_i_y(2, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_i_y(1, 0) = 1;
			ka.array_i_y(2, 0) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_y(1, 0) = 1 - gn_iy;
			ka.array_i_y(1, 1) = 1;
			ka.array_i_y(2, 0) = 1 + gn_iy;
			ka.array_i_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iy == -1)
	{
		if (frag_i == 1)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_y(0, 0) = 1 - gn_iy;
			ka.array_i_y(0, 1) = 1;
			ka.array_i_y(1, 0) = 1 + gn_iy;
			ka.array_i_y(1, 1) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_i_y(0, 0) = 1;
			ka.array_i_y(1, 0) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_y(0, 0) = 1 - gn_iy;
			ka.array_i_y(0, 1) = 1;
			ka.array_i_y(1, 0) = 1 + gn_iy;
			ka.array_i_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node iy value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jy == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_y(2, 0) = 1 - gn_jy;
			ka.array_j_y(2, 1) = 1;
			ka.array_j_y(3, 0) = 1 + gn_jy;
			ka.array_j_y(3, 1) = -1;

		}
		else if (frag_j == 2)
		{
			
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_j_y(2, 0) = 1;
			ka.array_j_y(3, 0) = -1;

		}
		else if (frag_j == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_y(2, 0) = 1 - gn_jy;
			ka.array_j_y(2, 1) = 1;
			ka.array_j_y(3, 0) = 1 + gn_jy;
			ka.array_j_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jy == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_y(1, 0) = 1 - gn_jy;
			ka.array_j_y(1, 1) = 1;
			ka.array_j_y(2, 0) = 1 + gn_jy;
			ka.array_j_y(2, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_j_y(1, 0) = 1;
			ka.array_j_y(2, 0) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_y(1, 0) = 1 - gn_jy;
			ka.array_j_y(1, 1) = 1;
			ka.array_j_y(2, 0) = 1 + gn_jy;
			ka.array_j_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jy == -1)
	{
		if (frag_j == 1)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_y(0, 0) = 1 - gn_jy;
			ka.array_j_y(0, 1) = 1;
			ka.array_j_y(1, 0) = 1 + gn_jy;
			ka.array_j_y(1, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_j_y(0, 0) = 1;
			ka.array_j_y(1, 0) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_y(0, 0) = 1 - gn_jy;
			ka.array_j_y(0, 1) = 1;
			ka.array_j_y(1, 0) = 1 + gn_jy;
			ka.array_j_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node jy value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_ky == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_y(2, 0) = 1 - gn_ky;
			ka.array_k_y(2, 1) = 1;
			ka.array_k_y(3, 0) = 1 + gn_ky;
			ka.array_k_y(3, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_k_y(2, 0) = 1;
			ka.array_k_y(3, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_y(2, 0) = 1 - gn_ky;
			ka.array_k_y(2, 1) = 1;
			ka.array_k_y(3, 0) = 1 + gn_ky;
			ka.array_k_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ky == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_y(1, 0) = 1 - gn_ky;
			ka.array_k_y(1, 1) = 1;
			ka.array_k_y(2, 0) = 1 + gn_ky;
			ka.array_k_y(2, 1) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_k_y(1, 0) = 1;
			ka.array_k_y(2, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_y(1, 0) = 1 - gn_ky;
			ka.array_k_y(1, 1) = 1;
			ka.array_k_y(2, 0) = 1 + gn_ky;
			ka.array_k_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ky == -1)
	{
		if (frag_k == 1)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_y(0, 0) = 1 - gn_ky;
			ka.array_k_y(0, 1) = 1;
			ka.array_k_y(1, 0) = 1 + gn_ky;
			ka.array_k_y(1, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_k_y(0, 0) = 1;
			ka.array_k_y(1, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_y(0, 0) = 1 - gn_ky;
			ka.array_k_y(0, 1) = 1;
			ka.array_k_y(1, 0) = 1 + gn_ky;
			ka.array_k_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node ky value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setKernelArrayForZ(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iz, int gn_jz, int gn_kz) {

	// 値の初期化
	ka.array_i_z = Eigen::MatrixXd::Zero(4, 2);
	ka.array_j_z = Eigen::MatrixXd::Zero(4, 2);
	ka.array_k_z = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_iz == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_z(2, 0) = 1 - gn_iz;
			ka.array_i_z(2, 1) = 1;
			ka.array_i_z(3, 0) = 1 + gn_iz;
			ka.array_i_z(3, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_z(2, 0) = 1 - gn_iz;
			ka.array_i_z(2, 1) = 1;
			ka.array_i_z(3, 0) = 1 + gn_iz;
			ka.array_i_z(3, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_i_z(2, 0) = 1;
			ka.array_i_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iz == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_z(1, 0) = 1 - gn_iz;
			ka.array_i_z(1, 1) = 1;
			ka.array_i_z(2, 0) = 1 + gn_iz;
			ka.array_i_z(2, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_z(1, 0) = 1 - gn_iz;
			ka.array_i_z(1, 1) = 1;
			ka.array_i_z(2, 0) = 1 + gn_iz;
			ka.array_i_z(2, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_i_z(1, 0) = 1;
			ka.array_i_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iz == -1)
	{
		if (frag_i == 1)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_z(0, 0) = 1 - gn_iz;
			ka.array_i_z(0, 1) = 1;
			ka.array_i_z(1, 0) = 1 + gn_iz;
			ka.array_i_z(1, 1) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_z(0, 0) = 1 - gn_iz;
			ka.array_i_z(0, 1) = 1;
			ka.array_i_z(1, 0) = 1 + gn_iz;
			ka.array_i_z(1, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_i_z(0, 0) = 1;
			ka.array_i_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node i value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jz == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_z(2, 0) = 1 - gn_jz;
			ka.array_j_z(2, 1) = 1;
			ka.array_j_z(3, 0) = 1 + gn_jz;
			ka.array_j_z(3, 1) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_z(2, 0) = 1 - gn_jz;
			ka.array_j_z(2, 1) = 1;
			ka.array_j_z(3, 0) = 1 + gn_jz;
			ka.array_j_z(3, 1) = -1;
		}
		else if (frag_j == 3)
		{
			
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_j_z(2, 0) = 1;
			ka.array_j_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jz == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_z(1, 0) = 1 - gn_jz;
			ka.array_j_z(1, 1) = 1;
			ka.array_j_z(2, 0) = 1 + gn_jz;
			ka.array_j_z(2, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_z(1, 0) = 1 - gn_jz;
			ka.array_j_z(1, 1) = 1;
			ka.array_j_z(2, 0) = 1 + gn_jz;
			ka.array_j_z(2, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_j_z(1, 0) = 1;
			ka.array_j_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jz == -1)
	{
		if (frag_j == 1)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_z(0, 0) = 1 - gn_jz;
			ka.array_j_z(0, 1) = 1;
			ka.array_j_z(1, 0) = 1 + gn_jz;
			ka.array_j_z(1, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_z(0, 0) = 1 - gn_jz;
			ka.array_j_z(0, 1) = 1;
			ka.array_j_z(1, 0) = 1 + gn_jz;
			ka.array_j_z(1, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_j_z(0, 0) = 1;
			ka.array_j_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node j value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_kz == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_z(2, 0) = 1 - gn_kz;
			ka.array_k_z(2, 1) = 1;
			ka.array_k_z(3, 0) = 1 + gn_kz;
			ka.array_k_z(3, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_z(2, 0) = 1 - gn_kz;
			ka.array_k_z(2, 1) = 1;
			ka.array_k_z(3, 0) = 1 + gn_kz;
			ka.array_k_z(3, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_k_z(2, 0) = 1;
			ka.array_k_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kz == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_z(1, 0) = 1 - gn_kz;
			ka.array_k_z(1, 1) = 1;
			ka.array_k_z(2, 0) = 1 + gn_kz;
			ka.array_k_z(2, 1) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_z(1, 0) = 1 - gn_kz;
			ka.array_k_z(1, 1) = 1;
			ka.array_k_z(2, 0) = 1 + gn_kz;
			ka.array_k_z(2, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_k_z(1, 0) = 1;
			ka.array_k_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kz == -1)
	{
		if (frag_k == 1)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_z(0, 0) = 1 - gn_kz;
			ka.array_k_z(0, 1) = 1;
			ka.array_k_z(1, 0) = 1 + gn_kz;
			ka.array_k_z(1, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_z(0, 0) = 1 - gn_kz;
			ka.array_k_z(0, 1) = 1;
			ka.array_k_z(1, 0) = 1 + gn_kz;
			ka.array_k_z(1, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_k_z(0, 0) = 1;
			ka.array_k_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node k value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setKernelArray2ForX(Kernel_Array ka, int gn_ix) {

	// 値の初期化
	ka.array_i_x = Eigen::MatrixXd::Zero(4, 2);

	// grid node ix について
	if (gn_ix == 1)
	{
		ka.array_i_x(2, 0) = 1 - gn_ix;
		ka.array_i_x(2, 1) = 1;
		ka.array_i_x(3, 0) = 1 + gn_ix;
		ka.array_i_x(3, 1) = -1;
	}
	else if (gn_ix == 0)
	{
		ka.array_i_x(1, 0) = 1 - gn_ix;
		ka.array_i_x(1, 1) = 1;
		ka.array_i_x(2, 0) = 1 + gn_ix;
		ka.array_i_x(2, 1) = -1;
	}
	else if (gn_ix == -1)
	{
		ka.array_i_x(0, 0) = 1 - gn_ix;
		ka.array_i_x(0, 1) = 1;
		ka.array_i_x(1, 0) = 1 + gn_ix;
		ka.array_i_x(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node x value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setKernelArray2ForY(Kernel_Array ka, int gn_iy) {

	// 値の初期化
	ka.array_i_y = Eigen::MatrixXd::Zero(4, 2);

	// grid node ix について
	if (gn_iy == 1)
	{
		ka.array_i_y(2, 0) = 1 - gn_iy;
		ka.array_i_y(2, 1) = 1;
		ka.array_i_y(3, 0) = 1 + gn_iy;
		ka.array_i_y(3, 1) = -1;
	}
	else if (gn_iy == 0)
	{
		ka.array_i_y(1, 0) = 1 - gn_iy;
		ka.array_i_y(1, 1) = 1;
		ka.array_i_y(2, 0) = 1 + gn_iy;
		ka.array_i_y(2, 1) = -1;
	}
	else if (gn_iy == -1)
	{
		ka.array_i_y(0, 0) = 1 - gn_iy;
		ka.array_i_y(0, 1) = 1;
		ka.array_i_y(1, 0) = 1 + gn_iy;
		ka.array_i_y(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node y value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setKernelArray2ForZ(Kernel_Array ka, int gn_iz) {

	// 値の初期化
	ka.array_i_z = Eigen::MatrixXd::Zero(4, 2);

	// grid node ix について
	if (gn_iz == 1)
	{
		ka.array_i_z(2, 0) = 1 - gn_iz;
		ka.array_i_z(2, 1) = 1;
		ka.array_i_z(3, 0) = 1 + gn_iz;
		ka.array_i_z(3, 1) = -1;
	}
	else if (gn_iz == 0)
	{
		ka.array_i_z(1, 0) = 1 - gn_iz;
		ka.array_i_z(1, 1) = 1;
		ka.array_i_z(2, 0) = 1 + gn_iz;
		ka.array_i_z(2, 1) = -1;
	}
	else if (gn_iz == -1)
	{
		ka.array_i_z(0, 0) = 1 - gn_iz;
		ka.array_i_z(0, 1) = 1;
		ka.array_i_z(1, 0) = 1 + gn_iz;
		ka.array_i_z(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node z value reference.\n";
		std::exit(1);
	}

	return ka;
}

std::vector<double> product(std::vector<double> a1, std::vector<double> a2)
{
	int length1 = a1.size();
	int length2 = a2.size();

	int lem = length1 + length2 - 1;

	std::vector<double> ans_array;

	for (int i = 0; i < lem; i++ )
	{
		ans_array.emplace_back(0.0);
	}

	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {

			ans_array[i + j] += a1[i] * a2[j];
		}
	}

	return ans_array;
}

double integral(int start_p, int stop_p, std::vector<double> array) 
{
	double ans = 0.0;
	for (int i = 0; i < array.size(); i++) {
		ans += array[i] * pow(stop_p, i + 1) / (i + 1) - array[i] * pow(start_p, i + 1) / (i + 1);
	}

	return ans;
}

double product_x(Kernel_Array ka) {
	double ans_x = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_o_vector;

	// std::cout << "X ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.array_i_x(a, 0));
		a_i_vector.emplace_back(ka.array_i_x(a, 1));

		a_j_vector.emplace_back(ka.array_j_x(a, 0));
		a_j_vector.emplace_back(ka.array_j_x(a, 1));

		a_k_vector.emplace_back(ka.array_k_x(a, 0));
		a_k_vector.emplace_back(ka.array_k_x(a, 1));

		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(product(a_i_vector, a_j_vector), product(a_k_vector, a_o_vector));
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "x dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/* std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_x += ans;
			///std::abs(ans);
	}

	return ans_x;
}

double product_y(Kernel_Array ka) {
	double ans_y = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Y ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.array_i_y(a, 0));
		a_i_vector.emplace_back(ka.array_i_y(a, 1));

		a_j_vector.emplace_back(ka.array_j_y(a, 0));
		a_j_vector.emplace_back(ka.array_j_y(a, 1));

		a_k_vector.emplace_back(ka.array_k_y(a, 0));
		a_k_vector.emplace_back(ka.array_k_y(a, 1));

		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(product(a_i_vector, a_j_vector), product(a_k_vector, a_o_vector));
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "y dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_y += ans;
			//std::abs(ans);
	}

	return ans_y;
}

double product_z(Kernel_Array ka) {
	double ans_z = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Z ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.array_i_z(a, 0));
		a_i_vector.emplace_back(ka.array_i_z(a, 1));

		a_j_vector.emplace_back(ka.array_j_z(a, 0));
		a_j_vector.emplace_back(ka.array_j_z(a, 1));

		a_k_vector.emplace_back(ka.array_k_z(a, 0));
		a_k_vector.emplace_back(ka.array_k_z(a, 1));

		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(product(a_i_vector, a_j_vector), product(a_k_vector, a_o_vector));
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "z dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_z += ans;
			//std::abs(ans);
	}

	return ans_z;
}

double product2_x(Kernel_Array ka) {
	double ans_x = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_o_vector;

	// std::cout << "X ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.array_i_x(a, 0));
		a_i_vector.emplace_back(ka.array_i_x(a, 1));
		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(a_i_vector, a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "x dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/ 

		ans_x += ans;
			//std::abs(ans);
	}

	return ans_x;
}

double product2_y(Kernel_Array ka) {
	double ans_y = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Y ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.array_i_y(a, 0));
		a_i_vector.emplace_back(ka.array_i_y(a, 1));
		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(a_i_vector, a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "y dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_y += ans;
			//std::abs(ans);
	}

	return ans_y;
}

double product2_z(Kernel_Array ka) {
	double ans_z = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Z ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.array_i_z(a, 0));
		a_i_vector.emplace_back(ka.array_i_z(a, 1));
		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(a_i_vector, a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "z dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/* std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_z += ans;
			//std::abs(ans);
	}

	return ans_z;
}

Kernel_Array createKernelArray3() {
	
	//計算確認用
	Eigen::MatrixXd array_kw_x = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_kw_y = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_kw_z = Eigen::MatrixXd::Zero(4, 2);

	// original point grid
	Eigen::MatrixXd array_o = Eigen::MatrixXd::Zero(4, 2);

	// 原点固定のkernel
	// [1+x][1-x]
	array_o(1, 0) = 1;
	array_o(1, 1) = 1;
	array_o(2, 0) = 1;
	array_o(2, 1) = -1;

	Kernel_Array ka = Kernel_Array(array_kw_x, array_kw_y, array_kw_z, array_o);

	return ka;
};

Kernel_Array setKernelArray3ForX(Kernel_Array ka, int gn_kx) {

	// 値の初期化
	ka.array_kw_x = Eigen::MatrixXd::Zero(4, 2);

	// grid node kx について
	if (gn_kx == 1)
	{
		ka.array_kw_x(2, 0) = 1 - gn_kx;
		ka.array_kw_x(2, 1) = 1;
		ka.array_kw_x(3, 0) = 1 + gn_kx;
		ka.array_kw_x(3, 1) = -1;
	}
	else if (gn_kx == 0)
	{
		ka.array_kw_x(1, 0) = 1 - gn_kx;
		ka.array_kw_x(1, 1) = 1;
		ka.array_kw_x(2, 0) = 1 + gn_kx;
		ka.array_kw_x(2, 1) = -1;
	}
	else if (gn_kx == -1)
	{
		ka.array_kw_x(0, 0) = 1 - gn_kx;
		ka.array_kw_x(0, 1) = 1;
		ka.array_kw_x(1, 0) = 1 + gn_kx;
		ka.array_kw_x(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node x value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setKernelArray3ForY(Kernel_Array ka, int gn_ky) {

	// 値の初期化
	ka.array_kw_y = Eigen::MatrixXd::Zero(4, 2);

	// grid node ky について
	if (gn_ky == 1)
	{
		ka.array_kw_y(2, 0) = 1 - gn_ky;
		ka.array_kw_y(2, 1) = 1;
		ka.array_kw_y(3, 0) = 1 + gn_ky;
		ka.array_kw_y(3, 1) = -1;
	}
	else if (gn_ky == 0)
	{
		ka.array_kw_y(1, 0) = 1 - gn_ky;
		ka.array_kw_y(1, 1) = 1;
		ka.array_kw_y(2, 0) = 1 + gn_ky;
		ka.array_kw_y(2, 1) = -1;
	}
	else if (gn_ky == -1)
	{
		ka.array_kw_y(0, 0) = 1 - gn_ky;
		ka.array_kw_y(0, 1) = 1;
		ka.array_kw_y(1, 0) = 1 + gn_ky;
		ka.array_kw_y(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node y value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setKernelArray3ForZ(Kernel_Array ka, int gn_kz) {

	// 値の初期化
	ka.array_kw_z = Eigen::MatrixXd::Zero(4, 2);

	// grid node kx について
	if (gn_kz == 1)
	{
		ka.array_kw_z(2, 0) = 1;
		ka.array_kw_z(3, 0) = -1;
	}
	else if (gn_kz == 0)
	{
		ka.array_kw_z(1, 0) = 1;
		ka.array_kw_z(2, 0) = -1;
	}
	else if (gn_kz == -1)
	{
		ka.array_kw_z(0, 0) = 1;
		ka.array_kw_z(1, 0) = -1;
	}
	else
	{
		std::cout << "Impossible grid node z value reference.\n";
		std::exit(1);
	}

	return ka;
}

double product3_x(Kernel_Array ka) {
	double ans_x = 0.0;
	std::vector<double> formula;
	std::vector<double> a_kw_vector;
	std::vector<double> a_o_vector;

	// std::cout << "X ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_kw_vector.clear();
		a_o_vector.clear();

		a_kw_vector.emplace_back(ka.array_kw_x(a, 0));
		a_kw_vector.emplace_back(ka.array_kw_x(a, 1));
		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(a_kw_vector, a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "x dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_x += ans;
		//std::abs(ans);
	}

	return ans_x;
}

double product3_y(Kernel_Array ka) {
	double ans_y = 0.0;
	std::vector<double> formula;
	std::vector<double> a_kw_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Y ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_kw_vector.clear();
		a_o_vector.clear();

		a_kw_vector.emplace_back(ka.array_kw_y(a, 0));
		a_kw_vector.emplace_back(ka.array_kw_y(a, 1));
		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(a_kw_vector, a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "y dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_y += ans;
		//std::abs(ans);
	}

	return ans_y;
}

double product3_z(Kernel_Array ka) {
	double ans_z = 0.0;
	std::vector<double> formula;
	std::vector<double> a_kw_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Z ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_kw_vector.clear();
		a_o_vector.clear();

		a_kw_vector.emplace_back(ka.array_kw_z(a, 0));
		a_kw_vector.emplace_back(ka.array_kw_z(a, 1));
		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(a_kw_vector, a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "z dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/* std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_z += ans;
		//std::abs(ans);
	}

	return ans_z;
}

Kernel_Array createKernelArray4() {

	//計算確認用
	Eigen::MatrixXd array_i_x = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_i_y = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_i_z = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_j_x = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_j_y = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_j_z = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_k_x = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_k_y = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_k_z = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_l_x = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_l_y = Eigen::MatrixXd::Zero(4, 2);
	Eigen::MatrixXd array_l_z = Eigen::MatrixXd::Zero(4, 2);

	/// original point grid
	Eigen::MatrixXd array_o = Eigen::MatrixXd::Zero(4, 2);

	// 原点固定のkernel
	// [1+x][1-x]
	array_o(1, 0) = 1;
	array_o(1, 1) = 1;
	array_o(2, 0) = 1;
	array_o(2, 1) = -1;

	Kernel_Array ka = Kernel_Array(array_i_x, array_i_y, array_i_z,
		array_j_x, array_j_y, array_j_z,
		array_k_x, array_k_y, array_k_z, 
		array_l_x, array_l_y, array_l_z,
		array_o);

	return ka;
};

double product4_x(Kernel_Array ka) {
	double ans_x = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_l_vector;
	std::vector<double> a_o_vector;

	// std::cout << "X ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_l_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.array_i_x(a, 0));
		a_i_vector.emplace_back(ka.array_i_x(a, 1));

		a_j_vector.emplace_back(ka.array_j_x(a, 0));
		a_j_vector.emplace_back(ka.array_j_x(a, 1));

		a_k_vector.emplace_back(ka.array_k_x(a, 0));
		a_k_vector.emplace_back(ka.array_k_x(a, 1));

		a_l_vector.emplace_back(ka.array_l_x(a, 0));
		a_l_vector.emplace_back(ka.array_l_x(a, 1));

		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(product(product(a_i_vector, a_j_vector), product(a_k_vector, a_l_vector)), a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "x dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/* std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_x += ans;
		///std::abs(ans);
	}

	return ans_x;
}

double product4_y(Kernel_Array ka) {
	double ans_y = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_l_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Y ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_l_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.array_i_y(a, 0));
		a_i_vector.emplace_back(ka.array_i_y(a, 1));

		a_j_vector.emplace_back(ka.array_j_y(a, 0));
		a_j_vector.emplace_back(ka.array_j_y(a, 1));

		a_k_vector.emplace_back(ka.array_k_y(a, 0));
		a_k_vector.emplace_back(ka.array_k_y(a, 1));

		a_l_vector.emplace_back(ka.array_l_y(a, 0));
		a_l_vector.emplace_back(ka.array_l_y(a, 1));

		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(product(product(a_i_vector, a_j_vector), product(a_k_vector, a_l_vector)), a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "y dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_y += ans;
		//std::abs(ans);
	}

	return ans_y;
}

double product4_z(Kernel_Array ka) {
	double ans_z = 0.0;
	std::vector<double> formula;
	std::vector<double> a_i_vector;
	std::vector<double> a_j_vector;
	std::vector<double> a_k_vector;
	std::vector<double> a_l_vector;
	std::vector<double> a_o_vector;

	// std::cout << "Z ----------------------------------------------------\n";
	for (int a = 0; a < 4; a++) {
		formula.clear();
		a_i_vector.clear();
		a_j_vector.clear();
		a_k_vector.clear();
		a_l_vector.clear();
		a_o_vector.clear();

		a_i_vector.emplace_back(ka.array_i_z(a, 0));
		a_i_vector.emplace_back(ka.array_i_z(a, 1));

		a_j_vector.emplace_back(ka.array_j_z(a, 0));
		a_j_vector.emplace_back(ka.array_j_z(a, 1));

		a_k_vector.emplace_back(ka.array_k_z(a, 0));
		a_k_vector.emplace_back(ka.array_k_z(a, 1));

		a_l_vector.emplace_back(ka.array_l_z(a, 0));
		a_l_vector.emplace_back(ka.array_l_z(a, 1));

		a_o_vector.emplace_back(ka.array_o(a, 0));
		a_o_vector.emplace_back(ka.array_o(a, 1));

		formula = product(product(product(a_i_vector, a_j_vector), product(a_k_vector, a_l_vector)), a_o_vector);
		double ans = integral(-2 + a, -1 + a, formula);
		// std::cout << "z dimension : " << "range=" << a + 1 << " integral : " << ans << "\n";

		// output formula
		/*
		std::cout << "formula : [ ";
		for (int v = 0; v < formula.size(); v++)
		{
			std::cout << formula[v] << " ";
		}
		std::cout << "]\n";
		*/

		ans_z += ans;
		//std::abs(ans);
	}

	return ans_z;
}

Kernel_Array setKernelArray4ForX(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_ix, int gn_jx, int gn_kx, int gn_lx) {

	// 値の初期化
	ka.array_i_x = Eigen::MatrixXd::Zero(4, 2);
	ka.array_j_x = Eigen::MatrixXd::Zero(4, 2);
	ka.array_k_x = Eigen::MatrixXd::Zero(4, 2);
	ka.array_l_x = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_ix == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_i_x(2, 0) = 1;
			ka.array_i_x(3, 0) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_x(2, 0) = 1 - gn_ix;
			ka.array_i_x(2, 1) = 1;
			ka.array_i_x(3, 0) = 1 + gn_ix;
			ka.array_i_x(3, 1) = -1;

		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_x(2, 0) = 1 - gn_ix;
			ka.array_i_x(2, 1) = 1;
			ka.array_i_x(3, 0) = 1 + gn_ix;
			ka.array_i_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ix == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_i_x(1, 0) = 1;
			ka.array_i_x(2, 0) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_x(1, 0) = 1 - gn_ix;
			ka.array_i_x(1, 1) = 1;
			ka.array_i_x(2, 0) = 1 + gn_ix;
			ka.array_i_x(2, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_x(1, 0) = 1 - gn_ix;
			ka.array_i_x(1, 1) = 1;
			ka.array_i_x(2, 0) = 1 + gn_ix;
			ka.array_i_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ix == -1)
	{
		if (frag_i == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_i_x(0, 0) = 1;
			ka.array_i_x(1, 0) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_x(0, 0) = 1 - gn_ix;
			ka.array_i_x(0, 1) = 1;
			ka.array_i_x(1, 0) = 1 + gn_ix;
			ka.array_i_x(1, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_x(0, 0) = 1 - gn_ix;
			ka.array_i_x(0, 1) = 1;
			ka.array_i_x(1, 0) = 1 + gn_ix;
			ka.array_i_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node i value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jx == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_j_x(2, 0) = 1;
			ka.array_j_x(3, 0) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_x(2, 0) = 1 - gn_jx;
			ka.array_j_x(2, 1) = 1;
			ka.array_j_x(3, 0) = 1 + gn_jx;
			ka.array_j_x(3, 1) = -1;

		}
		else if (frag_j == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_x(2, 0) = 1 - gn_jx;
			ka.array_j_x(2, 1) = 1;
			ka.array_j_x(3, 0) = 1 + gn_jx;
			ka.array_j_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jx == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_j_x(1, 0) = 1;
			ka.array_j_x(2, 0) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_x(1, 0) = 1 - gn_jx;
			ka.array_j_x(1, 1) = 1;
			ka.array_j_x(2, 0) = 1 + gn_jx;
			ka.array_j_x(2, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_x(1, 0) = 1 - gn_jx;
			ka.array_j_x(1, 1) = 1;
			ka.array_j_x(2, 0) = 1 + gn_jx;
			ka.array_j_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jx == -1)
	{
		if (frag_j == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_j_x(0, 0) = 1;
			ka.array_j_x(1, 0) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_x(0, 0) = 1 - gn_jx;
			ka.array_j_x(0, 1) = 1;
			ka.array_j_x(1, 0) = 1 + gn_jx;
			ka.array_j_x(1, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_x(0, 0) = 1 - gn_jx;
			ka.array_j_x(0, 1) = 1;
			ka.array_j_x(1, 0) = 1 + gn_jx;
			ka.array_j_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node j value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_kx == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_k_x(2, 0) = 1;
			ka.array_k_x(3, 0) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_x(2, 0) = 1 - gn_kx;
			ka.array_k_x(2, 1) = 1;
			ka.array_k_x(3, 0) = 1 + gn_kx;
			ka.array_k_x(3, 1) = -1;

		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_x(2, 0) = 1 - gn_kx;
			ka.array_k_x(2, 1) = 1;
			ka.array_k_x(3, 0) = 1 + gn_kx;
			ka.array_k_x(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kx == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_k_x(1, 0) = 1;
			ka.array_k_x(2, 0) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_x(1, 0) = 1 - gn_kx;
			ka.array_k_x(1, 1) = 1;
			ka.array_k_x(2, 0) = 1 + gn_kx;
			ka.array_k_x(2, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_x(1, 0) = 1 - gn_kx;
			ka.array_k_x(1, 1) = 1;
			ka.array_k_x(2, 0) = 1 + gn_kx;
			ka.array_k_x(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kx == -1)
	{
		if (frag_k == 1)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_k_x(0, 0) = 1;
			ka.array_k_x(1, 0) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_x(0, 0) = 1 - gn_kx;
			ka.array_k_x(0, 1) = 1;
			ka.array_k_x(1, 0) = 1 + gn_kx;
			ka.array_k_x(1, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_x(0, 0) = 1 - gn_kx;
			ka.array_k_x(0, 1) = 1;
			ka.array_k_x(1, 0) = 1 + gn_kx;
			ka.array_k_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node k value reference.\n";
		std::exit(1);
	}

	// grid node lx について
		if (gn_lx == 1)
		{
			ka.array_l_x(2, 0) = 1 - gn_lx;
			ka.array_l_x(2, 1) = 1;
			ka.array_l_x(3, 0) = 1 + gn_lx;
			ka.array_l_x(3, 1) = -1;
		}
		else if (gn_lx == 0)
		{
			ka.array_l_x(1, 0) = 1 - gn_lx;
			ka.array_l_x(1, 1) = 1;
			ka.array_l_x(2, 0) = 1 + gn_lx;
			ka.array_l_x(2, 1) = -1;
		}
		else if (gn_lx == -1)
		{
			ka.array_l_x(0, 0) = 1 - gn_lx;
			ka.array_l_x(0, 1) = 1;
			ka.array_l_x(1, 0) = 1 + gn_lx;
			ka.array_l_x(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible grid node x value reference.\n";
			std::exit(1);
		}

	return ka;
}

Kernel_Array setKernelArray4ForY(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iy, int gn_jy, int gn_ky, int gn_ly) {

	// 値の初期化
	ka.array_i_y = Eigen::MatrixXd::Zero(4, 2);
	ka.array_j_y = Eigen::MatrixXd::Zero(4, 2);
	ka.array_k_y = Eigen::MatrixXd::Zero(4, 2);
	ka.array_l_y = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_iy == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_y(2, 0) = 1 - gn_iy;
			ka.array_i_y(2, 1) = 1;
			ka.array_i_y(3, 0) = 1 + gn_iy;
			ka.array_i_y(3, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_i_y(2, 0) = 1;
			ka.array_i_y(3, 0) = -1;

		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_y(2, 0) = 1 - gn_iy;
			ka.array_i_y(2, 1) = 1;
			ka.array_i_y(3, 0) = 1 + gn_iy;
			ka.array_i_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iy == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_y(1, 0) = 1 - gn_iy;
			ka.array_i_y(1, 1) = 1;
			ka.array_i_y(2, 0) = 1 + gn_iy;
			ka.array_i_y(2, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_i_y(1, 0) = 1;
			ka.array_i_y(2, 0) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_y(1, 0) = 1 - gn_iy;
			ka.array_i_y(1, 1) = 1;
			ka.array_i_y(2, 0) = 1 + gn_iy;
			ka.array_i_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iy == -1)
	{
		if (frag_i == 1)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_y(0, 0) = 1 - gn_iy;
			ka.array_i_y(0, 1) = 1;
			ka.array_i_y(1, 0) = 1 + gn_iy;
			ka.array_i_y(1, 1) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_i_y(0, 0) = 1;
			ka.array_i_y(1, 0) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_y(0, 0) = 1 - gn_iy;
			ka.array_i_y(0, 1) = 1;
			ka.array_i_y(1, 0) = 1 + gn_iy;
			ka.array_i_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node i value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jy == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_y(2, 0) = 1 - gn_jy;
			ka.array_j_y(2, 1) = 1;
			ka.array_j_y(3, 0) = 1 + gn_jy;
			ka.array_j_y(3, 1) = -1;

		}
		else if (frag_j == 2)
		{

			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_j_y(2, 0) = 1;
			ka.array_j_y(3, 0) = -1;

		}
		else if (frag_j == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_y(2, 0) = 1 - gn_jy;
			ka.array_j_y(2, 1) = 1;
			ka.array_j_y(3, 0) = 1 + gn_jy;
			ka.array_j_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jy == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_y(1, 0) = 1 - gn_jy;
			ka.array_j_y(1, 1) = 1;
			ka.array_j_y(2, 0) = 1 + gn_jy;
			ka.array_j_y(2, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_j_y(1, 0) = 1;
			ka.array_j_y(2, 0) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_y(1, 0) = 1 - gn_jy;
			ka.array_j_y(1, 1) = 1;
			ka.array_j_y(2, 0) = 1 + gn_jy;
			ka.array_j_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jy == -1)
	{
		if (frag_j == 1)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_y(0, 0) = 1 - gn_jy;
			ka.array_j_y(0, 1) = 1;
			ka.array_j_y(1, 0) = 1 + gn_jy;
			ka.array_j_y(1, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_j_y(0, 0) = 1;
			ka.array_j_y(1, 0) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_y(0, 0) = 1 - gn_jy;
			ka.array_j_y(0, 1) = 1;
			ka.array_j_y(1, 0) = 1 + gn_jy;
			ka.array_j_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node j value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_ky == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_y(2, 0) = 1 - gn_ky;
			ka.array_k_y(2, 1) = 1;
			ka.array_k_y(3, 0) = 1 + gn_ky;
			ka.array_k_y(3, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_k_y(2, 0) = 1;
			ka.array_k_y(3, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_y(2, 0) = 1 - gn_ky;
			ka.array_k_y(2, 1) = 1;
			ka.array_k_y(3, 0) = 1 + gn_ky;
			ka.array_k_y(3, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ky == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_y(1, 0) = 1 - gn_ky;
			ka.array_k_y(1, 1) = 1;
			ka.array_k_y(2, 0) = 1 + gn_ky;
			ka.array_k_y(2, 1) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_k_y(1, 0) = 1;
			ka.array_k_y(2, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_y(1, 0) = 1 - gn_ky;
			ka.array_k_y(1, 1) = 1;
			ka.array_k_y(2, 0) = 1 + gn_ky;
			ka.array_k_y(2, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_ky == -1)
	{
		if (frag_k == 1)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_y(0, 0) = 1 - gn_ky;
			ka.array_k_y(0, 1) = 1;
			ka.array_k_y(1, 0) = 1 + gn_ky;
			ka.array_k_y(1, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_k_y(0, 0) = 1;
			ka.array_k_y(1, 0) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_y(0, 0) = 1 - gn_ky;
			ka.array_k_y(0, 1) = 1;
			ka.array_k_y(1, 0) = 1 + gn_ky;
			ka.array_k_y(1, 1) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node k value reference.\n";
		std::exit(1);
	}

	// grid node ly について
	if (gn_ly == 1)
	{
		ka.array_l_y(2, 0) = 1 - gn_ly;
		ka.array_l_y(2, 1) = 1;
		ka.array_l_y(3, 0) = 1 + gn_ly;
		ka.array_l_y(3, 1) = -1;
	}
	else if (gn_ly == 0)
	{
		ka.array_l_y(1, 0) = 1 - gn_ly;
		ka.array_l_y(1, 1) = 1;
		ka.array_l_y(2, 0) = 1 + gn_ly;
		ka.array_l_y(2, 1) = -1;
	}
	else if (gn_ly == -1)
	{
		ka.array_l_y(0, 0) = 1 - gn_ly;
		ka.array_l_y(0, 1) = 1;
		ka.array_l_y(1, 0) = 1 + gn_ly;
		ka.array_l_y(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node y value reference.\n";
		std::exit(1);
	}

	return ka;
}

Kernel_Array setKernelArray4ForZ(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iz, int gn_jz, int gn_kz, int gn_lz) {

	// 値の初期化
	ka.array_i_z = Eigen::MatrixXd::Zero(4, 2);
	ka.array_j_z = Eigen::MatrixXd::Zero(4, 2);
	ka.array_k_z = Eigen::MatrixXd::Zero(4, 2);
	ka.array_l_z = Eigen::MatrixXd::Zero(4, 2);

	// grid node i について
	if (gn_iz == 1) {
		if (frag_i == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_z(2, 0) = 1 - gn_iz;
			ka.array_i_z(2, 1) = 1;
			ka.array_i_z(3, 0) = 1 + gn_iz;
			ka.array_i_z(3, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_i, 1], [1 + gn_i, -1]]
			ka.array_i_z(2, 0) = 1 - gn_iz;
			ka.array_i_z(2, 1) = 1;
			ka.array_i_z(3, 0) = 1 + gn_iz;
			ka.array_i_z(3, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_i_z(2, 0) = 1;
			ka.array_i_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iz == 0) {
		if (frag_i == 1)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_z(1, 0) = 1 - gn_iz;
			ka.array_i_z(1, 1) = 1;
			ka.array_i_z(2, 0) = 1 + gn_iz;
			ka.array_i_z(2, 1) = -1;

		}
		else if (frag_i == 2)
		{
			// [[0, 0], [1 - gn_i, 1], [1 + gn_i, -1], [0, 0]]
			ka.array_i_z(1, 0) = 1 - gn_iz;
			ka.array_i_z(1, 1) = 1;
			ka.array_i_z(2, 0) = 1 + gn_iz;
			ka.array_i_z(2, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_i_z(1, 0) = 1;
			ka.array_i_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_iz == -1)
	{
		if (frag_i == 1)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_z(0, 0) = 1 - gn_iz;
			ka.array_i_z(0, 1) = 1;
			ka.array_i_z(1, 0) = 1 + gn_iz;
			ka.array_i_z(1, 1) = -1;
		}
		else if (frag_i == 2)
		{
			// [[1 - gn_i, 1], [1 + gn_i, -1], [0, 0], [0, 0]]
			ka.array_i_z(0, 0) = 1 - gn_iz;
			ka.array_i_z(0, 1) = 1;
			ka.array_i_z(1, 0) = 1 + gn_iz;
			ka.array_i_z(1, 1) = -1;
		}
		else if (frag_i == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_i_z(0, 0) = 1;
			ka.array_i_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node i value reference.\n";
		std::exit(1);
	}


	// grid node j について
	if (gn_jz == 1) {
		if (frag_j == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_z(2, 0) = 1 - gn_jz;
			ka.array_j_z(2, 1) = 1;
			ka.array_j_z(3, 0) = 1 + gn_jz;
			ka.array_j_z(3, 1) = -1;

		}
		else if (frag_j == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_j, 1], [1 + gn_j, -1]]
			ka.array_j_z(2, 0) = 1 - gn_jz;
			ka.array_j_z(2, 1) = 1;
			ka.array_j_z(3, 0) = 1 + gn_jz;
			ka.array_j_z(3, 1) = -1;
		}
		else if (frag_j == 3)
		{

			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_j_z(2, 0) = 1;
			ka.array_j_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jz == 0) {
		if (frag_j == 1)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_z(1, 0) = 1 - gn_jz;
			ka.array_j_z(1, 1) = 1;
			ka.array_j_z(2, 0) = 1 + gn_jz;
			ka.array_j_z(2, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[0, 0], [1 - gn_j, 1], [1 + gn_j, -1], [0, 0]]
			ka.array_j_z(1, 0) = 1 - gn_jz;
			ka.array_j_z(1, 1) = 1;
			ka.array_j_z(2, 0) = 1 + gn_jz;
			ka.array_j_z(2, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_j_z(1, 0) = 1;
			ka.array_j_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_jz == -1)
	{
		if (frag_j == 1)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_z(0, 0) = 1 - gn_jz;
			ka.array_j_z(0, 1) = 1;
			ka.array_j_z(1, 0) = 1 + gn_jz;
			ka.array_j_z(1, 1) = -1;
		}
		else if (frag_j == 2)
		{
			// [[1 - gn_j, 1], [1 + gn_j, -1], [0, 0], [0, 0]]
			ka.array_j_z(0, 0) = 1 - gn_jz;
			ka.array_j_z(0, 1) = 1;
			ka.array_j_z(1, 0) = 1 + gn_jz;
			ka.array_j_z(1, 1) = -1;
		}
		else if (frag_j == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_j_z(0, 0) = 1;
			ka.array_j_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node j value reference.\n";
		std::exit(1);
	}

	// grid node k について
	if (gn_kz == 1) {
		if (frag_k == 1)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_z(2, 0) = 1 - gn_kz;
			ka.array_k_z(2, 1) = 1;
			ka.array_k_z(3, 0) = 1 + gn_kz;
			ka.array_k_z(3, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[0, 0], [0, 0], [1 - gn_k, 1], [1 + gn_k, -1]]
			ka.array_k_z(2, 0) = 1 - gn_kz;
			ka.array_k_z(2, 1) = 1;
			ka.array_k_z(3, 0) = 1 + gn_kz;
			ka.array_k_z(3, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [0, 0], [1, 0], [-1, 0]]
			ka.array_k_z(2, 0) = 1;
			ka.array_k_z(3, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kz == 0) {
		if (frag_k == 1)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_z(1, 0) = 1 - gn_kz;
			ka.array_k_z(1, 1) = 1;
			ka.array_k_z(2, 0) = 1 + gn_kz;
			ka.array_k_z(2, 1) = -1;

		}
		else if (frag_k == 2)
		{
			// [[0, 0], [1 - gn_k, 1], [1 + gn_k, -1], [0, 0]]
			ka.array_k_z(1, 0) = 1 - gn_kz;
			ka.array_k_z(1, 1) = 1;
			ka.array_k_z(2, 0) = 1 + gn_kz;
			ka.array_k_z(2, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[0, 0], [1, 0], [-1, 0], [0, 0]]
			ka.array_k_z(1, 0) = 1;
			ka.array_k_z(2, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}
	}
	else if (gn_kz == -1)
	{
		if (frag_k == 1)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_z(0, 0) = 1 - gn_kz;
			ka.array_k_z(0, 1) = 1;
			ka.array_k_z(1, 0) = 1 + gn_kz;
			ka.array_k_z(1, 1) = -1;
		}
		else if (frag_k == 2)
		{
			// [[1 - gn_k, 1], [1 + gn_k, -1], [0, 0], [0, 0]]
			ka.array_k_z(0, 0) = 1 - gn_kz;
			ka.array_k_z(0, 1) = 1;
			ka.array_k_z(1, 0) = 1 + gn_kz;
			ka.array_k_z(1, 1) = -1;
		}
		else if (frag_k == 3)
		{
			// [[1, 0], [-1, 0], [0, 0], [0, 0]]
			ka.array_k_z(0, 0) = 1;
			ka.array_k_z(1, 0) = -1;
		}
		else
		{
			std::cout << "Impossible flag value reference.\n";
			std::exit(1);
		}

	}
	else
	{
		std::cout << "Impossible grid node k value reference.\n";
		std::exit(1);
	}

	// grid node lz について
	if (gn_lz == 1)
	{
		ka.array_l_z(2, 0) = 1 - gn_lz;
		ka.array_l_z(2, 1) = 1;
		ka.array_l_z(3, 0) = 1 + gn_lz;
		ka.array_l_z(3, 1) = -1;
	}
	else if (gn_lz == 0)
	{
		ka.array_l_z(1, 0) = 1 - gn_lz;
		ka.array_l_z(1, 1) = 1;
		ka.array_l_z(2, 0) = 1 + gn_lz;
		ka.array_l_z(2, 1) = -1;
	}
	else if (gn_lz == -1)
	{
		ka.array_l_z(0, 0) = 1 - gn_lz;
		ka.array_l_z(0, 1) = 1;
		ka.array_l_z(1, 0) = 1 + gn_lz;
		ka.array_l_z(1, 1) = -1;
	}
	else
	{
		std::cout << "Impossible grid node z value reference.\n";
		std::exit(1);
	}

	return ka;
}

