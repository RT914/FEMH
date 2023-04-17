#ifndef __CALCULATION_H__
#define __CALCULATION_H__

#include <Eigen/Dense>
#include <vector>

struct Kernel_Array
{
	Kernel_Array(Eigen::MatrixXd& aix, Eigen::MatrixXd& aiy, Eigen::MatrixXd& aiz,
		Eigen::MatrixXd& ajx, Eigen::MatrixXd& ajy, Eigen::MatrixXd& ajz,
		Eigen::MatrixXd& akx, Eigen::MatrixXd& aky, Eigen::MatrixXd& akz, Eigen::MatrixXd& ao) :
		array_i_x(aix), array_i_y(aiy), array_i_z(aiz),
		array_j_x(ajx), array_j_y(ajy), array_j_z(ajz),
		array_k_x(akx), array_k_y(aky), array_k_z(akz), array_o(ao){};

	/*
	Kernel_Array(Eigen::MatrixXd& aix, Eigen::MatrixXd& aiy, Eigen::MatrixXd& aiz, Eigen::MatrixXd& ao) :
		array_i_x(aix), array_i_y(aiy), array_i_z(aiz), array_o(ao) {};
		*/

	Kernel_Array(Eigen::MatrixXd& aix, Eigen::MatrixXd& aiy, Eigen::MatrixXd& aiz, Eigen::MatrixXd& ao) :
		array_kw_x(aix), array_kw_y(aiy), array_kw_z(aiz), array_o(ao) {};

	Kernel_Array(Eigen::MatrixXd& aix, Eigen::MatrixXd& aiy, Eigen::MatrixXd& aiz,
		Eigen::MatrixXd& ajx, Eigen::MatrixXd& ajy, Eigen::MatrixXd& ajz,
		Eigen::MatrixXd& akx, Eigen::MatrixXd& aky, Eigen::MatrixXd& akz,
		Eigen::MatrixXd& alx, Eigen::MatrixXd& aly, Eigen::MatrixXd& alz, 
		Eigen::MatrixXd& ao) :
		array_i_x(aix), array_i_y(aiy), array_i_z(aiz),
		array_j_x(ajx), array_j_y(ajy), array_j_z(ajz),
		array_k_x(akx), array_k_y(aky), array_k_z(akz),
		array_l_x(akx), array_l_y(aky), array_l_z(akz), array_o(ao) {};

	Eigen::MatrixXd array_i_x;
	Eigen::MatrixXd array_i_y;
	Eigen::MatrixXd array_i_z;
	Eigen::MatrixXd array_j_x;
	Eigen::MatrixXd array_j_y;
	Eigen::MatrixXd array_j_z;
	Eigen::MatrixXd array_k_x;
	Eigen::MatrixXd array_k_y;
	Eigen::MatrixXd array_k_z;
	Eigen::MatrixXd array_l_x;
	Eigen::MatrixXd array_l_y;
	Eigen::MatrixXd array_l_z;
	Eigen::MatrixXd array_o;

	Eigen::MatrixXd array_kw_x;
	Eigen::MatrixXd array_kw_y;
	Eigen::MatrixXd array_kw_z;

	bool cal_flag = true;
};

Kernel_Array createKernelArray();
Kernel_Array createKernelArray3();
Kernel_Array createKernelArray4();

Kernel_Array setKernelArrayForX(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_ix, int gn_jx, int gn_kx);
Kernel_Array setKernelArrayForY(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iy, int gn_jy, int gn_ky);
Kernel_Array setKernelArrayForZ(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iz, int gn_jz, int gn_kz);

// setKernelArray2 シリーズ
Kernel_Array setKernelArray2ForX(Kernel_Array ka, int gn_ix);
Kernel_Array setKernelArray2ForY(Kernel_Array ka, int gn_iy);
Kernel_Array setKernelArray2ForZ(Kernel_Array ka, int gn_iz);

// setKernelArrayk シリーズ
Kernel_Array setKernelArray3ForX(Kernel_Array ka, int gn_kx);
Kernel_Array setKernelArray3ForY(Kernel_Array ka, int gn_ky);
Kernel_Array setKernelArray3ForZ(Kernel_Array ka, int gn_kz);

Kernel_Array setKernelArray4ForX(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_ix, int gn_jx, int gn_kx, int gn_lx);
Kernel_Array setKernelArray4ForY(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iy, int gn_jy, int gn_ky, int gn_ly);
Kernel_Array setKernelArray4ForZ(Kernel_Array ka, int frag_i, int frag_j, int frag_k, int gn_iz, int gn_jz, int gn_kz, int gn_lz);

std::vector<double> product(std::vector<double> a1, std::vector<double> a2);

double integral(int start_p, int stop_p, std::vector<double> array);

double product_x(Kernel_Array ka);
double product_y(Kernel_Array ka);
double product_z(Kernel_Array ka);

double product2_x(Kernel_Array ka);
double product2_y(Kernel_Array ka);
double product2_z(Kernel_Array ka);

double product3_x(Kernel_Array ka);
double product3_y(Kernel_Array ka);
double product3_z(Kernel_Array ka);

double product4_x(Kernel_Array ka);
double product4_y(Kernel_Array ka);
double product4_z(Kernel_Array ka);


#endif
