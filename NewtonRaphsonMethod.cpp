#include <Eigen/Dense>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>

#include "NewtonRaphsonMethod.h"
#include "Square.h"
#include "fem.h"
#include "Export.h"
#include "Calculation.h"

// [\alpha][\beta][\gamma][i-xi][j-xi][k-xi]
double MatrixQForX[3][3][3][3][3][3] = {};
double MatrixQForY[3][3][3][3][3][3] = {};
double MatrixQForZ[3][3][3][3][3][3] = {};

// [i-xi][j-xi][k-xi]
double VectorbPreparation[3][3][3] = {};

double Vectorb_part2_PreparationX[3][3][3][3][3][3][3] = {};
double Vectorb_part2_PreparationY[3][3][3][3][3][3][3] = {};
double Vectorb_part2_PreparationZ[3][3][3][3][3][3][3] = {};

//確認用
double CheckMatrixPreparation[3][3][3] = {};
Eigen::MatrixXd CheckMatrix(NumberOfParticles, NumberOfParticles);


// 行列の大きさは [3×NumberOfParticles][NumberOfParticles]　とする．
Eigen::MatrixXd MatrixQ(3 * NumberOfParticles, NumberOfParticles);
Eigen::MatrixXd MatrixQ1(NumberOfParticles, NumberOfParticles);

// ベクトルの大きさは [NumberOfParticles]　とする．
Eigen::VectorXd Vectorb1(NumberOfParticles);
Eigen::VectorXd Vectorb2(NumberOfParticles);
Eigen::VectorXd Vectorb(NumberOfParticles);

// 行列の大きさは [3×NumberOfParticles][3×NumberOfParticles]　とする．
Eigen::MatrixXd MatrixR;

// ベクトルの大きさは [3×NumberOfParticles]　とする．
Eigen::VectorXd Vectorc;

// ベクトルの大きさは[3×NumberOfParticles]　とする．参照用
Eigen::VectorXd re_barphi(3 * NumberOfParticles);

//回数
int det_cal_time = 0.0;

//ニュートン反復回数
int newton_num = 0;

Eigen::VectorXd Newton(Square square) {

	double NormVectorDeltaPhi = 1.0;
	int num = 0;
	int SquarePointsNumber = square.points.size();
	// ベクトルの大きさは [3×NumberOfParticles]　とする．
	Eigen::VectorXd barphi(3 * NumberOfParticles);

	// 座標の取得
	for (int i = 0; i < SquarePointsNumber; i++) {
		barphi(3*i) = square.points[i].position[0];
		barphi(3*i+1) = square.points[i].position[1];
		barphi(3*i+2) = square.points[i].position[2];
	}

	// std::cout << barphi << "\n";
	
	//while (NormVectorDeltaPhi > 1.0e-6 && num <= 10000) {
	while (NormVectorDeltaPhi > 1.0e-6) {
		MatrixQ = calMatrixQ(square, barphi);
		Vectorb = calVectorb(square, barphi);
		

		MatrixR = calMatrixR(MatrixQ);
		Vectorc = calVectorc(Vectorb, MatrixQ);

		Eigen::FullPivLU<Eigen::MatrixXd> LU(MatrixR);
		Eigen::VectorXd VectorDeltaPhi = LU.solve(Vectorc);

		NormVectorDeltaPhi = VectorDeltaPhi.norm();

		std::cout << "num = " << num << " : " << NormVectorDeltaPhi << "\n";

		// std::cout << "barphi size\n" << "行" << barphi.rows() << "\n" << "列" << barphi.cols() << "\n";

		// std::cout << num << "------------------------------------------\n";
		//std::cout << VectorDeltaPhi << "\n";

		
		// 座標の取得
		barphi += VectorDeltaPhi;
		
		num++;
		
		
	}

	// std::cout << barphi << "\n";
	return barphi;
};

int Newton_loop(Square square, int seed) {

	//std::string output_csv_file_name = "DeltaVect.csv";
	//std::ofstream data_file(output_csv_file_name);

	double NormVectorDeltaPhi = 1.0;
	int num = 0;
	int SquarePointsNumber = square.points.size();
	srand(seed);
	// ベクトルの大きさは [3×NumberOfParticles]　とする．
	Eigen::VectorXd barphi(3 * NumberOfParticles);

	// ベクトルの大きさは [3×NumberOfParticles]　とする．
	Eigen::VectorXd barphi_prime(3 * NumberOfParticles);

	Eigen::VectorXd barphi_dis(3 * NumberOfParticles);

	//double c = pow(10.0, (1 - seed));
	double c = pow(10.0, 1 - seed);
	srand(12);

	// 座標の取得
	for (int i = 0; i < SquarePointsNumber; i++) {
		barphi(3 * i) = square.points[i].position[0] + c * ((rand() % 2) - 1);
		barphi(3 * i + 1) = square.points[i].position[1] + c * ((rand() % 2) - 1);
		barphi(3 * i + 2) = square.points[i].position[2] + c * ((rand() % 2) - 1);

	}

	for (int i = 0; i < SquarePointsNumber; i++) {
		re_barphi(3 * i) = square.points[i].reference_position[0];
		re_barphi(3 * i + 1) = square.points[i].reference_position[1];
		re_barphi(3 * i + 2) = square.points[i].reference_position[2];
	}

	// Header
	//data_file << "Newton Times" << "," << "Value" << ",";
	//data_file << std::endl;

	while (NormVectorDeltaPhi > 1.0e-6) {
		MatrixQ = calMatrixQ(square, barphi);
		Vectorb = calVectorb_part2(square, barphi);

		MatrixR = calMatrixR(MatrixQ);
		Vectorc = calVectorc(Vectorb, MatrixQ);

		Eigen::FullPivLU<Eigen::MatrixXd> LU(MatrixR);
		Eigen::VectorXd VectorDeltaPhi = LU.solve(Vectorc);

		NormVectorDeltaPhi = VectorDeltaPhi.norm();

		std::cout << "num = " << num << " : " << NormVectorDeltaPhi << "\n";

		//data_file << num << "," << NormVectorDeltaPhi << ",";
		//data_file << std::endl;

		//Line Search with Wolfe
		double sigma = 0.5;
		double eps1 = 0.2;
		double eps2 = 0.5;
		double lambda = 1.0;

		barphi_prime = VectorDeltaPhi * lambda + barphi;

		/*
		//Armijo
		Eigen::VectorXd f_new = calVectorb_part2(square, barphi_prime);
		Eigen::VectorXd f = Vectorb;
		Eigen::VectorXd nabla_f = eps1 * lambda * MatrixQ.transpose() * VectorDeltaPhi;
		Eigen::VectorXd right_f = f + nabla_f;

		//curvature
		Eigen::VectorXd nabla_f_new = calMatrixQ(square, barphi_prime).transpose() * VectorDeltaPhi;
		Eigen::VectorXd nabla_f_prime = eps2 * calMatrixQ(square, barphi).transpose() * VectorDeltaPhi;


		while (f_new.norm() > right_f.norm() || abs(nabla_f_prime.norm()) < abs(nabla_f_prime.norm())) {
			lambda = sigma * lambda;
			barphi_prime = VectorDeltaPhi * lambda + barphi;
			//f_new = calVectorb(square, barphi_prime);
			f_new = calVectorb_part2(square, barphi_prime);
			nabla_f = MatrixQ.transpose() * eps1 * lambda * VectorDeltaPhi;
			right_f = f + nabla_f;
			nabla_f_new = calMatrixQ(square, barphi_prime).transpose() * VectorDeltaPhi;
			nabla_f_prime = eps2 * calMatrixQ(square, barphi).transpose() * VectorDeltaPhi;
		}

		//std::cout << lambda << "\n";
		*/

		// 座標の取得
		barphi += lambda * VectorDeltaPhi;
		//barphi += VectorDeltaPhi;
		//std::cout << VectorDeltaPhi << "\n";

		num++;

	}
	//std::cout << barphi << "\n";
	
	//data_file.close();
	
	return num;
};

Eigen::VectorXd Newton_one(Square square) {

	double NormVectorDeltaPhi = 1.0;
	
	int SquarePointsNumber = square.points.size();

	// ベクトルの大きさは [3×NumberOfParticles]　とする．
	Eigen::VectorXd barphi(3 * NumberOfParticles);

	// ベクトルの大きさは [3×NumberOfParticles]　とする．
	Eigen::VectorXd barphi_prime(3 * NumberOfParticles);

	// 座標の取得
	for (int i = 0; i < SquarePointsNumber; i++) {
		barphi(3 * i) = square.points[i].position[0];
		barphi(3 * i + 1) = square.points[i].position[1];
		barphi(3 * i + 2) = square.points[i].position[2];
	}

	for (int i = 0; i < SquarePointsNumber; i++) {
		re_barphi(3 * i) = square.points[i].reference_position[0];
		re_barphi(3 * i + 1) = square.points[i].reference_position[1];
		re_barphi(3 * i + 2) = square.points[i].reference_position[2];
	}

	/*
	std::cout << "reference" << "\n";
	std::cout << re_barphi << "\n";
	
	std::cout << re_barphi(39) << "\n";
	std::cout << re_barphi(40) << "\n";
	std::cout << re_barphi(41) << "\n";
	*/

	/*
	std::cout << "before" << "\n";
	std::cout << (re_barphi - barphi).norm() << "\n";
	*/

	/*
	std::cout << barphi(39) << "\n";
	std::cout << barphi(40) << "\n";
	std::cout << barphi(41) << "\n";
	*/


	if (NormVectorDeltaPhi > 1.0e-6) {
		MatrixQ = calMatrixQ(square, barphi);
		//Vectorb = calVectorb(square, barphi);
		Vectorb = calVectorb_part2(square, barphi);

		//std::cout << MatrixQ.rows() << ", " << MatrixQ.cols() << "\n";
		//std::cout << Vectorb.rows() << ", " << Vectorb.cols() << "\n";

		/*
		std::string filename = "matrixQ1.csv";
		for (int i = 0; i < MatrixQ1.rows(); i++) {
			for (int j = 0; j < MatrixQ1.cols(); j++) {
				MatrixQ1(i, j) = MatrixQ(3 * i, j);
			}
		}
		export_Q(MatrixQ1, filename);
		*/


		// std::cout << Vectorb << "\n";


		MatrixR = calMatrixR(MatrixQ);
		Vectorc = calVectorc(Vectorb, MatrixQ);
		//std::cout << MatrixR.rows() << ", " << MatrixR.cols() << "\n";
		//std::cout << Vectorc.rows() << ", " << Vectorc.cols() << "\n";
		// std::cout << Vectorc.rows() << Vectorc.cols() << "\n";

		Eigen::FullPivLU<Eigen::MatrixXd> LU(MatrixR);
		Eigen::VectorXd VectorDeltaPhi = LU.solve(Vectorc);
		//std::cout << VectorDeltaPhi.rows() << ", " << VectorDeltaPhi.cols() << "\n";

		NormVectorDeltaPhi = VectorDeltaPhi.norm();

		std::cout << "num = " << newton_num << " : " << NormVectorDeltaPhi << "\n";
		newton_num++;


		// std::cout << "barphi size\n" << "行" << barphi.rows() << "\n" << "列" << barphi.cols() << "\n";

		// std::cout << num << "------------------------------------------\n";
		//std::cout << VectorDeltaPhi << "\n";

		
		//Line Search with Wolfe
		double sigma = 0.5;
		double eps1 = 0.2;
		double eps2 = 0.5;
		double lambda = 1.0;
		
		barphi_prime = VectorDeltaPhi * lambda + barphi;

		//Armijo
		//Eigen::VectorXd f_new = calVectorb(square, barphi_prime);
		Eigen::VectorXd f_new = calVectorb_part2(square, barphi_prime);
		Eigen::VectorXd f = Vectorb;
		Eigen::VectorXd nabla_f = eps1 * lambda * MatrixQ.transpose() * VectorDeltaPhi;
		Eigen::VectorXd right_f = f + nabla_f;

		//curvature
		Eigen::VectorXd nabla_f_new = calMatrixQ(square, barphi_prime).transpose() * VectorDeltaPhi;
		Eigen::VectorXd nabla_f_prime = eps2 * calMatrixQ(square, barphi).transpose() * VectorDeltaPhi;

		
		while (f_new.norm() > right_f.norm() || nabla_f_prime.norm() < nabla_f_prime.norm()) {
			lambda = sigma * lambda;
			barphi_prime = VectorDeltaPhi * lambda + barphi;
			//f_new = calVectorb(square, barphi_prime);
			f_new = calVectorb_part2(square, barphi_prime);
			nabla_f = MatrixQ.transpose() * eps1 * lambda * VectorDeltaPhi;
			right_f = f + nabla_f;
			nabla_f_new = calMatrixQ(square, barphi_prime).transpose() * VectorDeltaPhi;
			nabla_f_prime = eps2 * calMatrixQ(square, barphi).transpose() * VectorDeltaPhi;
		}
		

		//std::cout << lambda << "\n";
		
		// 座標の取得
		barphi += lambda * VectorDeltaPhi;

	
	}
	
	/*
	std::cout << "result" << "\n";
	std::cout << (re_barphi - barphi).norm() << "\n";
	*/

	return barphi;
};

void Newton_copy(Square square) {

	// ベクトルの大きさは [3×NumberOfParticles]　とする．
	Eigen::VectorXd barphi(3 * NumberOfParticles);


	// 座標の取得
	for (int i = 0; i < square.points.size(); i++) {
		barphi(3 * i) = square.points[i].position[0];
		barphi(3 * i + 1) = square.points[i].position[1];
		barphi(3 * i + 2) = square.points[i].position[2];
	}

	//std::cout << "barphi : \n";
	//std::cout << barphi << "\n";

	CheckMatrix = calCheckMatrix(square);
	//std::cout << CheckMatrix.cols() << "\n";
	std::string output_csv_file_name = "result_w3w.csv";
	export_Q(CheckMatrix, output_csv_file_name);

	
	MatrixQ = calMatrixQ(square, barphi);
	Vectorb = calVectorb(square, barphi);
	//std::cout << MatrixQ << "\n";
	//std::cout << Vectorb << "\n";

	for (int i = 0; i < MatrixQ1.rows(); i++) {
		for (int j = 0; j < MatrixQ1.cols(); j++) {
			MatrixQ1(i, j) = MatrixQ(3 * i + 2, j);
		}
	}
	// std::cout << MatrixQ1 << "\n";

	output_csv_file_name = "result_matrixQ3.csv";
	export_Q(MatrixQ1, output_csv_file_name);

	output_csv_file_name = "result_matrixQ.csv";
	//export_Q(MatrixQ, output_csv_file_name);

	MatrixR = calMatrixR(MatrixQ);
	Vectorc = calVectorc(Vectorb, MatrixQ);
	//std::cout << "MatrixR : " << MatrixR << "\n";
	//std::cout << "Vectorc : " << Vectorc << "\n";

	// 出力結果の確認
	/*
	for (int i = 0; i < pow(NumberOfOneDemensionParticles, 3); i++) {
		std::cout << "bar_phi_i : " << square.points[i].position[0] << " "
			<< square.points[i].position[1] << " "
			<< square.points[i].position[2] << "\n";
	}
	*/

	/*
	std::cout << "Matrix Q---------------------------------" << "\n";
	for (int k = 0; k < NumberOfParticles * 3; k++) {
		for (int xi = 0; xi < NumberOfParticles; xi++) {
			std::cout << MatrixQ(k, xi) << "   ";
		}
		std::cout << "\n";
	}
	*/
	
	/*
	Eigen::MatrixXd MatrixQ_trans = MatrixQ.transpose();
	Eigen::MatrixXd Vectorb_trans = Vectorb.transpose();
	std::cout << "Matrix Q size\n" << "行" << MatrixQ.rows() << "\n" << "列" << MatrixQ.cols() << "\n";
	std::cout << "Matrix Q trans size\n" << "行" << MatrixQ_trans.rows() << "\n" << "列" << MatrixQ_trans.cols() << "\n";
	std::cout << "Vector b size\n" << "行" << Vectorb.rows() << "\n" << "列" << Vectorb.cols() << "\n";
	std::cout << "Vector b trans size\n" << "行" << Vectorb_trans.rows() << "\n" << "列" << Vectorb_trans.cols() << "\n";
	*/

	/*
	std::cout << "Matrix Q.transpose---------------------------------" << "\n";
	for (int k = 0; k < NumberOfParticles; k++) {
		for (int xi = 0; xi < NumberOfParticles * 3; xi++) {
			std::cout << MatrixQ_trans(k, xi) << "   ";
		}
		std::cout << "\n";
	}
	*/
	


	export_CSV(MatrixQForX, MatrixQForY, MatrixQForZ);
	export_CSV_X(MatrixQForX);
	export_CSV_Y(MatrixQForY);
	export_CSV_Z(MatrixQForZ);
	export2_CSV(VectorbPreparation);

	// export2_CSV(VectorbPreparation);

	/*
	std::cout << "\n" << "Vector b" << "\n";
	for (int xi = 0; xi < Vectorb.size(); xi++) {
		std::cout << xi << " : " << Vectorb(xi) << "\n";
	}
	*/
	

	/*
	std::cout << "Matrix R---------------------------------" << "\n";
	for (int k = 0; k < MatrixR.rows(); k++) {
		//std::cout << k <<" : ";
		for (int l = 0; l < MatrixR.cols(); l++) {
			std::cout << MatrixR(k, l) << "   ";
		}
		std::cout << "\n";
	}
	*/
	
	


	/*
	std::cout << "\n" << "Vector c" << "\n";
	for (int l = 0; l < Vectorc.size(); l++) {
		std::cout << l << " : " << Vectorc(l) << "\n";
	}
	*/
	
	//std::cout << "Matrix R size\n" << "行" << MatrixR.rows() << "\n" << "列" << MatrixR.cols() << "\n";
	//std::cout << "Vector c size\n" << "行" << Vectorc.rows() << "\n" << "列" << Vectorc.cols() << "\n";
	
	
	Eigen::FullPivLU<Eigen::MatrixXd> LU(MatrixR);
	Eigen::VectorXd VectorDeltaPhi = LU.solve(Vectorc);

	// std::cout << "\n" << "VectorDeltaPhi" << "\n";
	for (int l = 0; l < VectorDeltaPhi.size(); l++) {
		//std::cout << l << " : " << VectorDeltaPhi(l) << "\n";
	}
	

	// std::cout << "VectorPhi size\n" << "行" << VectorPhi.rows() << "\n" << "列" << VectorPhi.cols() << "\n";

	//showDetail( 0.1, 0.1, 0.1, square, barphi);
	//std::cout << "norm : " << VectorDeltaPhi.norm() << "\n";
};

void calMatrixQPreparation(Square square) {
	Kernel_Array ka = createKernelArray();

	// preparation
	// kernel index
	for (int l = 0; l < dimensions; l++) {
		for (int m = 0; m < dimensions; m++) {
			for (int n = 0; n < dimensions; n++) {

				int ki_i = l + 1;  // grid i において微分するkernel
				int	ki_j = m + 1;  // grid j において微分するkernel
				int	ki_k = n + 1;  // grid k において微分するkernel

				for (int ix = 0; ix < dimensions; ix++) {
					for (int jx = 0; jx < dimensions; jx++) {
						for (int kx = 0; kx < dimensions; kx++) {

							// case x
							// -1 ～ 1
							int gn_i_x = ix - 1;
							int gn_j_x = jx - 1;
							int gn_k_x = kx - 1;
						
							ka = setKernelArrayForX(ka, ki_i, ki_j, ki_k, gn_i_x, gn_j_x, gn_k_x);
							MatrixQForX[l][m][n][ix][jx][kx] = product_x(ka);

						}
					}
				}

				for (int iy = 0; iy < dimensions; iy++) {
					for (int jy = 0; jy < dimensions; jy++) {
						for (int ky = 0; ky < dimensions; ky++) {

							// case y
							// -1 ～ 1
							int gn_i_y = iy - 1;
							int gn_j_y = jy - 1;
							int gn_k_y = ky - 1;

							ka = setKernelArrayForY(ka, ki_i, ki_j, ki_k, gn_i_y, gn_j_y, gn_k_y);
							MatrixQForY[l][m][n][iy][jy][ky] = product_y(ka);

						}
					}
				}

				for (int iz = 0; iz < dimensions; iz++) {
					for (int jz = 0; jz < dimensions; jz++) {
						for (int kz = 0; kz < dimensions; kz++) {

							// case z
							// -1 ～ 1
							int gn_i_z = iz - 1;
							int gn_j_z = jz - 1;
							int gn_k_z = kz - 1;

							ka = setKernelArrayForZ(ka, ki_i, ki_j, ki_k, gn_i_z, gn_j_z, gn_k_z);
							MatrixQForZ[l][m][n][iz][jz][kz] = product_z(ka);
							// std::cout << MatrixQForZ[l][m][n][iz][jz][kz] << "\n";

						}
					}
				}



			}
		}
	}
}

Eigen::MatrixXd calMatrixQ(Square square, Eigen::VectorXd v) {

	int las = 0;

	bool FlagCalGrid;

	Eigen::MatrixXd OutputM(3*square.points.size(), square.points.size());

	//calculate Q matrix
	for (int xi = 0; xi < NumberOfParticles; xi++) {
		for (int k = 0; k < NumberOfParticles; k++) {

			Eigen::Vector3i grid_k = FlatToGrid(k);
			Eigen::Vector3i grid_xi = FlatToGrid(xi);
			Eigen::Vector3i k_minus_xi = grid_k - grid_xi;

			OutputM(3 * k, xi) = 0.0;
			OutputM(3 * k + 1, xi) = 0.0;
			OutputM(3 * k + 2, xi) = 0.0;

			double m_x[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			double m_y[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
			double m_z[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

			double Q1 = 0.0;
			double Q2 = 0.0;
			double Q3 = 0.0;

			int n = 0;

			for (int i = 0; i < NumberOfParticles; i++) {
				for (int j = 0; j < NumberOfParticles; j++) {

					FlagCalGrid = false;

					Eigen::Vector3i grid_i = FlatToGrid(i);
					Eigen::Vector3i grid_j = FlatToGrid(j);

					Eigen::Vector3i i_minus_xi = grid_i - grid_xi;
					Eigen::Vector3i j_minus_xi = grid_j - grid_xi;

					
					// \Phi = \bar{\varphi2}^i * \bar{\varphi3}^j - \bar{\varphi3}^i * \bar{\varphi2}^j
					//double Phi1 = (barphi(3 * i + 1) * barphi(3 * j + 2) - barphi(3 * i + 2) * barphi(3 * j + 1));
					double Phi1_1 = (v(3 * i + 1) * v(3 * j + 2));
					double Phi1_2 = (v(3 * i + 2) * v(3 * j + 1));
					

					// \Phi = \bar{\varphi3}^i * \bar{\varphi1}^j - \bar{\varphi1}^i * \bar{\varphi3}^j
					//double Phi2 = (barphi(3 * i + 2) * barphi(3 * j) - barphi(3 * i) * barphi(3 * j + 2));
					double Phi2_1 = (v(3 * i + 2) * v(3 * j));
					double Phi2_2 = (v(3 * i) * v(3 * j + 2));
					

					// \Phi = \bar{\varphi1}^i * \bar{\varphi2}^j - \bar{\varphi2}^i * \bar{\varphi1}^j
					//double Phi3 = (barphi(3 * i) * barphi(3 * j + 1) - barphi(3 * i + 1) * barphi(3 * j));
					double Phi3_1 = (v(3 * i) * v(3 * j + 1));
					double Phi3_2 = (v(3 * i + 1) * v(3 * j));

					/*
					if (k == 0 && xi == 0) 
					{
							std::cout << i << " , " << j << " Phi1: " << Phi1 << "\n";
							std::cout << i << " , " << j << " Phi2: " << Phi2 << "\n";
							std::cout << i << " , " << j << " Phi3: " << Phi3 << "\n";
					}
					*/
					

					// i - xi, j - xi, k - xi が -1～1の範囲であるかの判定
					//i - xi
					if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {

						//j - xi
						if ((abs(j_minus_xi[0]) <= 1) && (abs(j_minus_xi[1]) <= 1) && (abs(j_minus_xi[2]) <= 1)) {

							//k - xi
							if ((abs(k_minus_xi[0]) <= 1) && (abs(k_minus_xi[1]) <= 1) && (abs(k_minus_xi[2]) <= 1)) {

								FlagCalGrid = true;
								//std::cout << "true\n";

							}
						}
					}

					// 判定〇なら計算
					if (FlagCalGrid) {

						/*
						double w2w3w1 = MatrixQForX[1][2][0][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[1][2][0][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[1][2][0][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w2w3w2 = MatrixQForX[1][2][1][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[1][2][1][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[1][2][1][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w2w3w3 = MatrixQForX[1][2][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[1][2][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[1][2][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];


						double w1w3w1 = MatrixQForX[0][2][0][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][2][0][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][2][0][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w1w3w2 = MatrixQForX[0][2][1][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][2][1][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][2][1][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w1w3w3 = MatrixQForX[0][2][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][2][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][2][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];


						double w1w2w1 = MatrixQForX[0][1][0][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][1][0][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][1][0][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w1w2w2 = MatrixQForX[0][1][1][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][1][1][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][1][1][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w1w2w3 = MatrixQForX[0][1][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][1][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][1][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];
						

						double Phi1 = barphi(3 * i + 1) * barphi(3 * j + 2) - barphi(3 * i + 2) * barphi(3 * j + 1);
						double Phi2 = barphi(3 * i + 2) * barphi(3 * j) - barphi(3 * i) * barphi(3 * j + 2);
						double Phi3 = barphi(3 * i) * barphi(3 * j + 1) - barphi(3 * i + 1) * barphi(3 * j);

						Q1 += Phi1 * w2w3w1 + Phi2 * w2w3w2 + Phi3 * w2w3w3;
						Q2 += - Phi1 * w1w3w1 - Phi2 * w1w3w2 - Phi3 * w1w3w3;
						Q3 += Phi1 * w1w2w1 + Phi2 * w1w2w2 + Phi3 * w1w2w3;
						*/

						
						double w2w3w1 = MatrixQForX[1][2][0][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[1][2][0][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[1][2][0][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w1w3w2 = MatrixQForX[0][2][1][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][2][1][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][2][1][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w1w2w3 = MatrixQForX[0][1][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][1][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][1][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];
						
						//m_x += Phi1 * w2w3w1;
						//m_y += -Phi2 * w1w3w2;
						//m_z += Phi3 * w1w2w3;
						
						m_x[0] += Phi1_1 * w2w3w1;
						m_x[1] += Phi1_2 * w2w3w1;
						m_x[2] += Phi1_2 * w1w3w2;
						m_x[3] += Phi1_1 * w1w3w2;
						m_x[4] += Phi1_1 * w1w2w3;
						m_x[5] += Phi1_2 * w1w2w3;

						m_y[0] += Phi2_1 * w2w3w1;
						m_y[1] += Phi2_2 * w2w3w1;
						m_y[2] += Phi2_2 * w1w3w2;
						m_y[3] += Phi2_1 * w1w3w2;
						m_y[4] += Phi2_1 * w1w2w3;
						m_y[5] += Phi2_2 * w1w2w3;

						m_z[0] += Phi3_1 * w2w3w1;
						m_z[1] += Phi3_2 * w2w3w1;
						m_z[2] += Phi3_2 * w1w3w2;
						m_z[3] += Phi3_1 * w1w3w2;
						m_z[4] += Phi3_1 * w1w2w3;
						m_z[5] += Phi3_2 * w1w2w3;
						

					

						/*
						if ( Phi2 * w2w3w2 != 0) {
							std::cout << "Phi1 * w2w3w1 : " << Phi2 * w2w3w2 << "\n";
						}
						*/

						
						if (xi == 1 && k == 0) {
							n++;
							/*
							Eigen::Vector3i grid_i = FlatToGrid(i);
							Eigen::Vector3i grid_j = FlatToGrid(j);
							Eigen::Vector3i grid_xi = FlatToGrid(xi);
							std::cout << "i-xi : " << grid_i(0) - grid_xi(0) << grid_i(1) - grid_xi(1) << grid_i(2) - grid_xi(2) 
								<< ", j : " << grid_j(0) - grid_xi(0) << grid_j(1) - grid_xi(1) << grid_j(2) - grid_xi(2) 
								<< ", Phi1_2 * w1w3w2 : " << Phi1_2 << "\n";
								*/
							//std::cout << "times : " << n << "\n";
							//std::cout << ", grid i-xi : " << i_minus_xi(0) << i_minus_xi(1) << i_minus_xi(2) << "\n";
							//std::cout << ", grid j-xi : " << j_minus_xi(0) << j_minus_xi(1) << j_minus_xi(2) << "\n";
							//std::cout << ", grid k-xi : " << k_minus_xi(0) << k_minus_xi(1) << k_minus_xi(2) << "\n";
							//std::cout << ", Phi1 : " << Phi1_1  <<", w2w3w1 : " << w2w3w1 << "\n";
							//std::cout << ", Phi1 * w2w3w1: " << Phi1_1 * w1w3w2 << "\n";

						}
						
						
					}
				}
			}

			

			
			if (abs(m_x[1]) > 1.0e-6) {
				//std::cout << "xi : " << xi << ", k : " << k << ", m_x[1] : " << m_x[1] << "\n";
			}
			
			/*
			OutputM(3 * k, xi) = Q1;
			OutputM(3 * k + 1, xi) = Q2;
			OutputM(3 * k + 2, xi) = Q3;
			*/
			
			
			OutputM(3 * k, xi) = m_x[0] - m_x[1] + m_x[2] - m_x[3] + m_x[4] - m_x[5];
			OutputM(3 * k + 1, xi) = m_y[0] - m_y[1] + m_y[2] - m_y[3] + m_y[4] - m_y[5];
			OutputM(3 * k + 2, xi) = m_z[0] - m_z[1] + m_z[2] - m_z[3] + m_z[4] - m_z[5];
			
			

			// index比較用
			if (k == 1 && xi == 0) {
				//std::cout << "Phi1 : " << Phi1 << "\n";
				//std::cout << "w2w3w1 : " << w2w3w1 << "\n";
				//std::cout << "Phi1 * w2w3w1 : " << - Phi1 * w1w3w2 << "\n";
				//std::cout << "Phi2 : " << Phi2 << "\n";
				//std::cout << "w1w2w2 : " << w2w3w2 << "\n";
				///std::cout << "Q3 : " << Phi1 * w1w2w1 + Phi2 * w1w2w2 + Phi3 * w1w2w3 << "\n";
				//std::cout << "Q1 : " << OutputM(3 * k, xi) << "\n";
				//std::cout << "m_x[0] : " << m_x[5] << "\n";
			}
			
		}
	}

	

	return OutputM;
}

void calVectorbPreparation(Square square) {
	Kernel_Array ka = createKernelArray();

	// grid node
	for (int ix = 0; ix < dimensions; ix++) {
		for (int iy = 0; iy < dimensions; iy++) {
			for (int iz = 0; iz < dimensions; iz++) {

				// -1 ～ 1
				int gn_ix = ix - 1;
				int gn_iy = iy - 1;
				int gn_iz = iz - 1;

				// std::cout << "grid node index : ";
				// std::cout << "x " << gn_ix << ",y " << gn_iy << ",z " << gn_iz << "\n";

				ka = setKernelArray2ForX(ka, gn_ix);
				ka = setKernelArray2ForY(ka, gn_iy);
				ka = setKernelArray2ForZ(ka, gn_iz);

				double el = product2_x(ka) * product2_y(ka) * product2_z(ka);

				// std::cout << "result : " << el << "\n" << "\n";

				VectorbPreparation[ix][iy][iz] = { el };
			}
		}
	}
}

Eigen::VectorXd calVectorb(Square square, Eigen::VectorXd v) {
	bool FlagCalGrid;
	Eigen::VectorXd OutputV(square.points.size());
	//calculate b vector
	for (int xi = 0; xi < NumberOfParticles; xi++) {

		Eigen::Vector3i grid_xi = FlatToGrid(xi);

		Vectorb1(xi) = 0.0;
		Vectorb2(xi) = 0.0;

		for (int i = 0; i < NumberOfParticles; i++) {

			FlagCalGrid = false;
			Eigen::Vector3i grid_i = FlatToGrid(i);
			//double Theta_i = det(square.points[xi].deformation_gradient);
			double Theta_i = square.points[i].theta;
			// std::cout << Theta_i << "\n";

			Eigen::Vector3i i_minus_xi = grid_i - grid_xi;
			// std::cout << i_minus_xi << "\n";


			// i - xiが -1～1の範囲であるかの判定
			//i - xi
			if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {

				FlagCalGrid = true;
				//std::cout << "true\n";
				//std::cout << abs(i_minus_xi[0]) << "\n";
				//std::cout << abs(i_minus_xi[1]) << "\n";
				//std::cout << abs(i_minus_xi[2]) << "\n";

			}
			else {
				//std::cout << abs(i_minus_xi[0]) << "\n";
				//std::cout << abs(i_minus_xi[1]) << "\n";
				//std::cout << abs(i_minus_xi[2]) << "\n";
			}

			if (FlagCalGrid == true) {
				double ww = VectorbPreparation[i_minus_xi[0] + 1][i_minus_xi[1] + 1][i_minus_xi[2] + 1];
				// std::cout << ww << "\n";
				Vectorb1(xi) += Theta_i * ww;
				//Vectorb1(xi) += ww;
				// std::cout << "ww : " << ww << ", Theta_i : " << Theta_i << "\n";
				// std::cout << "xi : " << xi << ", V : " << Vectorb1(xi) << "\n";
				// index確認用
				/*
				if (xi == 0) {
					std::cout << "i : " << i << ",ww : " << ww << "\n";
				}
				*/

			}

		}

		for (int i = 0; i < NumberOfParticles; i++) {

			Eigen::Vector3i grid_i = FlatToGrid(i);
			Eigen::Vector3i i_minus_xi = grid_i - grid_xi;


			double bar_phi_i_x = v(3 * i);
			double bar_phi_i_y = v(3 * i + 1);
			double bar_phi_i_z = v(3 * i + 2);

			/*
			if (xi == 13) {
				//std::cout << "phi : " << Phi << "\n";
				std::cout << "bar_phi_i : " << i << " : " << square.points[i].position[0] << " "
					<< square.points[i].position[1] << " "
					<< square.points[i].position[2] << "\n";
			}
			*/


			for (int j = 0; j < NumberOfParticles; j++) {

				Eigen::Vector3i grid_j = FlatToGrid(j);
				Eigen::Vector3i j_minus_xi = grid_j - grid_xi;

				double bar_phi_j_x = v(3 * j);
				double bar_phi_j_y = v(3 * j + 1);
				double bar_phi_j_z = v(3 * j + 2);

				/*
				if (xi == 0 && i == 0) {

					std::cout << "bar_phi_j : " << square.points[j].position[0] << " "
						<< square.points[j].position[1] << " "
						<< square.points[j].position[2] << "\n";
				}
				*/


				for (int k = 0; k < NumberOfParticles; k++) {

					FlagCalGrid = false;

					Eigen::Vector3i grid_k = FlatToGrid(k);
					Eigen::Vector3i k_minus_xi = grid_k - grid_xi;

					double bar_phi_k_x = v(3 * k);
					double bar_phi_k_y = v(3 * k + 1);
					double bar_phi_k_z = v(3 * k + 2);

					/*
					if (xi == 0 && i == 0 && j == 0) {

						std::cout << "bar_phi_k : " << square.points[k].position[0] << " "
							<< square.points[k].position[1] << " "
							<< square.points[k].position[2] << "\n";
					}
					*/


					// i - xi, j - xi, k - xi が -1～1の範囲であるかの判定
					//i - xi
					if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {
						//j - xi
						if ((abs(j_minus_xi[0]) <= 1) && (abs(j_minus_xi[1]) <= 1) && (abs(j_minus_xi[2]) <= 1)) {
							//k - xi
							if ((abs(k_minus_xi[0]) <= 1) && (abs(k_minus_xi[1]) <= 1) && (abs(k_minus_xi[2]) <= 1)) {

								FlagCalGrid = true;
								//std::cout << "true\n";

							}
						}
					}

					if (FlagCalGrid == true) {
						
						double Phi2 = bar_phi_i_x *  bar_phi_j_y * bar_phi_k_z;
						
						
						double Phi1 = 
							bar_phi_i_x * (bar_phi_j_y * bar_phi_k_z - bar_phi_j_z * bar_phi_k_y) 
							+ bar_phi_i_y * (bar_phi_j_z * bar_phi_k_x - bar_phi_j_x * bar_phi_k_z)
							+ bar_phi_i_z * (bar_phi_j_x * bar_phi_k_y - bar_phi_j_y * bar_phi_k_x);
							

						double w1w2w3 = MatrixQForX[0][1][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][1][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][1][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w2w3w1 = MatrixQForX[1][2][0][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[1][2][0][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[1][2][0][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w3w1w2 = MatrixQForX[2][0][1][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[2][0][1][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[2][0][1][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w3w2w1 = MatrixQForX[2][1][0][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[2][1][0][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[2][1][0][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w2w1w3 = MatrixQForX[1][0][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[1][0][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[1][0][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w1w3w2 = MatrixQForX[0][2][1][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][2][1][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][2][1][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						// std::cout << Phi << "\n";

						Vectorb2(xi) += Phi1 * w1w2w3;
						//Vectorb2(xi) += Phi2 * (w1w2w3 + w2w3w1 + w3w1w2 - w3w2w1 - w2w1w3 - w1w3w2);
						
						
						if (xi == 0 && k == 0) {
							//std::cout << "i : " << i << ", j : " << j << ", k : " << k << "  ";
							//std::cout << "No.  " << xi << "  Phi1 * w1w2w3 = " << Phi1 * w1w2w3 << "";
							//std::cout << "V2 : " << Vectorb2(xi) << "\n";

							/*
							std::cout << "i : " << bar_phi_i_x << bar_phi_i_y << bar_phi_i_z 
								<< ", j : " << bar_phi_j_x << bar_phi_j_y << bar_phi_j_z
								<< ", k : " << bar_phi_k_x << bar_phi_k_y << bar_phi_k_z << "\n";
							*/
							//std::cout << bar_phi_i_x * (bar_phi_j_y * bar_phi_k_z - bar_phi_j_z * bar_phi_k_y) << "\n";
							//std::cout << bar_phi_i_y * (bar_phi_j_z * bar_phi_k_x - bar_phi_j_x * bar_phi_k_z) << "\n";
							//std::cout << bar_phi_i_z * (bar_phi_j_x * bar_phi_k_y - bar_phi_j_y * bar_phi_k_x) << "\n";
							
								

							//std::cout << "V2 : " << Phi1 * w1w2w3 << "\n";
							
							
							//std::cout << i << "grid_i : " << grid_i << "\n";
							//std::cout << j << "grid_j : " << grid_j << "\n";
							//std::cout << k << "grid_k : " << grid_k << "\n";
							
							//std::cout << GridToFlat(k_minus_xi) << "times" << "\n";
							//std::cout << "i_minus_xi : " << i_minus_xi << "\n";
							//std::cout << "j_minus_xi : " << j_minus_xi << "\n";
							//std::cout << "k_minus_xi : " << k_minus_xi << "\n";
							

							//std::cout << bar_phi_j_z << "\n";

							//std::cout << "No.  " << xi << "  w1w2w3 : " << w1w2w3 << "\n";
							///std::cout << "No.  " << xi << "  Phi1 : " << Phi1 << " w1w2w3 : " << w1w2w3 << "\n";
							//std::cout << "No.  " << xi << "  V2 : " << Phi1 * w1w2w3 << "\n";
							
						}
						
						//std::cout << "No.  " << xi << "  V2 : " << Phi << "\n";
						//std::cout << "No.  " << xi << "  V2 : " << w1w2w3 << "\n";
						//std::cout << "No.  " << xi << "  V2 : " << Phi1 * w1w2w3 << "\n";
					}

				}
			}
		}
		
		OutputV(xi) = ( pow(square.dx, 3) * Vectorb1(xi) ) - Vectorb2(xi);

		//std::cout << "No.  " << xi << "  V1 : " << pow(square.dx, 3) * Vectorb1(xi) << "\n";
		
		//std::cout << "No.  " << xi << "  V2 : " << Vectorb2(xi) << "\n";
		//std::cout << "No.  " << xi << "  V : " << OutputV(xi) << "\n";
		
		
		

	}

	// std::cout << Vectorb1 << "\n";

	return OutputV;
}

void calVectorb_part2_Preparation(Square square) {
	Kernel_Array ka = createKernelArray4();
	

	// preparation
	// kernel index
	for (int l = 0; l < dimensions; l++) {
		for (int m = 0; m < dimensions; m++) {
			for (int n = 0; n < dimensions; n++) {

				int ki_i = l + 1;  // grid i において微分するkernel
				int	ki_j = m + 1;  // grid j において微分するkernel
				int	ki_k = n + 1;  // grid k において微分するkernel


				for (int ix = 0; ix < dimensions; ix++) {
					for (int jx = 0; jx < dimensions; jx++) {
						for (int kx = 0; kx < dimensions; kx++) {
							for (int lx = 0; lx < dimensions; lx++) {
								// case x
								// -1 ～ 1
								int gn_i_x = ix - 1;
								int gn_j_x = jx - 1;
								int gn_k_x = kx - 1;
								int gn_l_x = lx - 1;

								ka = setKernelArray4ForX(ka, ki_i, ki_j, ki_k, gn_i_x, gn_j_x, gn_k_x, gn_l_x);
								Vectorb_part2_PreparationX[l][m][n][ix][jx][kx][lx] = product4_x(ka);
							}
						}
					}
				}

				for (int iy = 0; iy < dimensions; iy++) {
					for (int jy = 0; jy < dimensions; jy++) {
						for (int ky = 0; ky < dimensions; ky++) {
							for (int ly = 0; ly < dimensions; ly++) {
								// case x
								// -1 ～ 1
								int gn_i_y = iy - 1;
								int gn_j_y = jy - 1;
								int gn_k_y = ky - 1;
								int gn_l_y = ly - 1;

								ka = setKernelArray4ForY(ka, ki_i, ki_j, ki_k, gn_i_y, gn_j_y, gn_k_y, gn_l_y);
								Vectorb_part2_PreparationY[l][m][n][iy][jy][ky][ly] = product4_y(ka);
							}
						}
					}
				}

				for (int iz = 0; iz < dimensions; iz++) {
					for (int jz = 0; jz < dimensions; jz++) {
						for (int kz = 0; kz < dimensions; kz++) {
							for (int lz = 0; lz < dimensions; lz++) {
								// case x
								// -1 ～ 1
								int gn_i_z = iz - 1;
								int gn_j_z = jz - 1;
								int gn_k_z = kz - 1;
								int gn_l_z = lz - 1;

								ka = setKernelArray4ForZ(ka, ki_i, ki_j, ki_k, gn_i_z, gn_j_z, gn_k_z, gn_l_z);
								Vectorb_part2_PreparationZ[l][m][n][iz][jz][kz][lz] = product4_z(ka);
							}
						}
					}
				}



			}
		}
	}
}

Eigen::VectorXd calVectorb_part2(Square square, Eigen::VectorXd v) {
	
	bool FlagCalGrid;
	det_cal_time += 1;

	Eigen::VectorXd OutputV(square.points.size());


	//calculate b vector
	for (int xi = 0; xi < NumberOfParticles; xi++) {

		Eigen::Vector3i grid_xi = FlatToGrid(xi);

		Vectorb2(xi) = 0.0;

		if (det_cal_time == 1) {
			Vectorb1(xi) = 0.0;
			for (int i = 0; i < NumberOfParticles; i++) {

				Eigen::Vector3i grid_i = FlatToGrid(i);
				Eigen::Vector3i i_minus_xi = grid_i - grid_xi;


				double re_bar_phi_i_x = re_barphi(3 * i);
				double re_bar_phi_i_y = re_barphi(3 * i + 1);
				double re_bar_phi_i_z = re_barphi(3 * i + 2);

				/*
				if (xi == 13) {
					//std::cout << "phi : " << Phi << "\n";
					std::cout << "bar_phi_i : " << i << " : " << square.points[i].position[0] << " "
						<< square.points[i].position[1] << " "
						<< square.points[i].position[2] << "\n";
				}
				*/


				for (int j = 0; j < NumberOfParticles; j++) {

					Eigen::Vector3i grid_j = FlatToGrid(j);
					Eigen::Vector3i j_minus_xi = grid_j - grid_xi;

					double re_bar_phi_j_x = re_barphi(3 * j);
					double re_bar_phi_j_y = re_barphi(3 * j + 1);
					double re_bar_phi_j_z = re_barphi(3 * j + 2);



					for (int k = 0; k < NumberOfParticles; k++) {

						FlagCalGrid = false;

						Eigen::Vector3i grid_k = FlatToGrid(k);
						Eigen::Vector3i k_minus_xi = grid_k - grid_xi;

						double re_bar_phi_k_x = re_barphi(3 * k);
						double re_bar_phi_k_y = re_barphi(3 * k + 1);
						double re_bar_phi_k_z = re_barphi(3 * k + 2);


						// i - xi, j - xi, k - xi が -1～1の範囲であるかの判定
						//i - xi
						if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {
							//j - xi
							if ((abs(j_minus_xi[0]) <= 1) && (abs(j_minus_xi[1]) <= 1) && (abs(j_minus_xi[2]) <= 1)) {
								//k - xi
								if ((abs(k_minus_xi[0]) <= 1) && (abs(k_minus_xi[1]) <= 1) && (abs(k_minus_xi[2]) <= 1)) {

									FlagCalGrid = true;
									//std::cout << "true\n";

								}
							}
						}

						if (FlagCalGrid == true) {

							double Phi1_b1 =
								re_bar_phi_i_x * (re_bar_phi_j_y * re_bar_phi_k_z - re_bar_phi_j_z * re_bar_phi_k_y)
								+ re_bar_phi_i_y * (re_bar_phi_j_z * re_bar_phi_k_x - re_bar_phi_j_x * re_bar_phi_k_z)
								+ re_bar_phi_i_z * (re_bar_phi_j_x * re_bar_phi_k_y - re_bar_phi_j_y * re_bar_phi_k_x);


							double w1w2w3w = MatrixQForX[0][1][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
								* MatrixQForY[0][1][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
								* MatrixQForZ[0][1][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];


							Vectorb1(xi) += Phi1_b1 * w1w2w3w;
							
						}

					}
				}
			}
		}
		

		for (int i = 0; i < NumberOfParticles; i++) {

			Eigen::Vector3i grid_i = FlatToGrid(i);
			Eigen::Vector3i i_minus_xi = grid_i - grid_xi;


			double bar_phi_i_x = v(3 * i);
			double bar_phi_i_y = v(3 * i + 1);
			double bar_phi_i_z = v(3 * i + 2);

			/*
			if (xi == 13) {
				//std::cout << "phi : " << Phi << "\n";
				std::cout << "bar_phi_i : " << i << " : " << square.points[i].position[0] << " "
					<< square.points[i].position[1] << " "
					<< square.points[i].position[2] << "\n";
			}
			*/


			for (int j = 0; j < NumberOfParticles; j++) {

				Eigen::Vector3i grid_j = FlatToGrid(j);
				Eigen::Vector3i j_minus_xi = grid_j - grid_xi;

				double bar_phi_j_x = v(3 * j);
				double bar_phi_j_y = v(3 * j + 1);
				double bar_phi_j_z = v(3 * j + 2);

				/*
				if (xi == 0 && i == 0) {

					std::cout << "bar_phi_j : " << square.points[j].position[0] << " "
						<< square.points[j].position[1] << " "
						<< square.points[j].position[2] << "\n";
				}
				*/


				for (int k = 0; k < NumberOfParticles; k++) {

					FlagCalGrid = false;

					Eigen::Vector3i grid_k = FlatToGrid(k);
					Eigen::Vector3i k_minus_xi = grid_k - grid_xi;

					double bar_phi_k_x = v(3 * k);
					double bar_phi_k_y = v(3 * k + 1);
					double bar_phi_k_z = v(3 * k + 2);

					/*
					if (xi == 0 && i == 0 && j == 0) {

						std::cout << "bar_phi_k : " << square.points[k].position[0] << " "
							<< square.points[k].position[1] << " "
							<< square.points[k].position[2] << "\n";
					}
					*/


					// i - xi, j - xi, k - xi が -1～1の範囲であるかの判定
					//i - xi
					if ((abs(i_minus_xi[0]) <= 1) && (abs(i_minus_xi[1]) <= 1) && (abs(i_minus_xi[2]) <= 1)) {
						//j - xi
						if ((abs(j_minus_xi[0]) <= 1) && (abs(j_minus_xi[1]) <= 1) && (abs(j_minus_xi[2]) <= 1)) {
							//k - xi
							if ((abs(k_minus_xi[0]) <= 1) && (abs(k_minus_xi[1]) <= 1) && (abs(k_minus_xi[2]) <= 1)) {

								FlagCalGrid = true;
								//std::cout << "true\n";

							}
						}
					}

					if (FlagCalGrid == true) {

						double Phi2 = bar_phi_i_x * bar_phi_j_y * bar_phi_k_z;


						double Phi1 =
							bar_phi_i_x * (bar_phi_j_y * bar_phi_k_z - bar_phi_j_z * bar_phi_k_y)
							+ bar_phi_i_y * (bar_phi_j_z * bar_phi_k_x - bar_phi_j_x * bar_phi_k_z)
							+ bar_phi_i_z * (bar_phi_j_x * bar_phi_k_y - bar_phi_j_y * bar_phi_k_x);


						double w1w2w3 = MatrixQForX[0][1][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][1][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][1][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w2w3w1 = MatrixQForX[1][2][0][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[1][2][0][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[1][2][0][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w3w1w2 = MatrixQForX[2][0][1][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[2][0][1][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[2][0][1][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w3w2w1 = MatrixQForX[2][1][0][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[2][1][0][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[2][1][0][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w2w1w3 = MatrixQForX[1][0][2][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[1][0][2][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[1][0][2][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						double w1w3w2 = MatrixQForX[0][2][1][i_minus_xi[0] + 1][j_minus_xi[0] + 1][k_minus_xi[0] + 1]
							* MatrixQForY[0][2][1][i_minus_xi[1] + 1][j_minus_xi[1] + 1][k_minus_xi[1] + 1]
							* MatrixQForZ[0][2][1][i_minus_xi[2] + 1][j_minus_xi[2] + 1][k_minus_xi[2] + 1];

						// std::cout << Phi << "\n";

						Vectorb2(xi) += Phi1 * w1w2w3;
						//Vectorb2(xi) += Phi2 * (w1w2w3 + w2w3w1 + w3w1w2 - w3w2w1 - w2w1w3 - w1w3w2);


						if (xi == 0 && k == 0) {
							//std::cout << "i : " << i << ", j : " << j << ", k : " << k << "  ";
							//std::cout << "No.  " << xi << "  Phi1 * w1w2w3 = " << Phi1 * w1w2w3 << "";
							//std::cout << "V2 : " << Vectorb2(xi) << "\n";

							/*
							std::cout << "i : " << bar_phi_i_x << bar_phi_i_y << bar_phi_i_z
								<< ", j : " << bar_phi_j_x << bar_phi_j_y << bar_phi_j_z
								<< ", k : " << bar_phi_k_x << bar_phi_k_y << bar_phi_k_z << "\n";
							*/
							//std::cout << bar_phi_i_x * (bar_phi_j_y * bar_phi_k_z - bar_phi_j_z * bar_phi_k_y) << "\n";
							//std::cout << bar_phi_i_y * (bar_phi_j_z * bar_phi_k_x - bar_phi_j_x * bar_phi_k_z) << "\n";
							//std::cout << bar_phi_i_z * (bar_phi_j_x * bar_phi_k_y - bar_phi_j_y * bar_phi_k_x) << "\n";



							//std::cout << "V2 : " << Phi1 * w1w2w3 << "\n";


							//std::cout << i << "grid_i : " << grid_i << "\n";
							//std::cout << j << "grid_j : " << grid_j << "\n";
							//std::cout << k << "grid_k : " << grid_k << "\n";

							//std::cout << GridToFlat(k_minus_xi) << "times" << "\n";
							//std::cout << "i_minus_xi : " << i_minus_xi << "\n";
							//std::cout << "j_minus_xi : " << j_minus_xi << "\n";
							//std::cout << "k_minus_xi : " << k_minus_xi << "\n";


							//std::cout << bar_phi_j_z << "\n";

							//std::cout << "No.  " << xi << "  w1w2w3 : " << w1w2w3 << "\n";
							///std::cout << "No.  " << xi << "  Phi1 : " << Phi1 << " w1w2w3 : " << w1w2w3 << "\n";
							//std::cout << "No.  " << xi << "  V2 : " << Phi1 * w1w2w3 << "\n";

						}

						//std::cout << "No.  " << xi << "  V2 : " << Phi << "\n";
						//std::cout << "No.  " << xi << "  V2 : " << w1w2w3 << "\n";
						//std::cout << "No.  " << xi << "  V2 : " << Phi1 * w1w2w3 << "\n";
					}

				}
			}
		}

		OutputV(xi) = Vectorb1(xi) - Vectorb2(xi);




		//std::cout << "No.  " << xi << "  V1 : " << Vectorb1(xi) << "\n";
		//std::cout << "No.  " << xi << "  V2 : " << Vectorb2(xi) << "\n";
		//std::cout << "No.  " << xi << "  V : " << OutputV(xi) << "\n";




	}

	// std::cout << Vectorb1 << "\n";

	return OutputV;
}

Eigen::MatrixXd calMatrixR(Eigen::MatrixXd M) {
	return (M) * (M.transpose());
}

Eigen::VectorXd calVectorc(Eigen::VectorXd V, Eigen::MatrixXd M) {
	// 本来，vector b は転置しないが，列ベクトルであるため，行ベクトルにした．
	return V.transpose() * (M.transpose());
}

int GridToFlat(Eigen::Vector3i grid_index)
{
	int flat_index;
	int grid_x = grid_index[0];
	int grid_y = grid_index[1];
	int grid_z = grid_index[2];

	flat_index = grid_x * int(pow(NumberOfOneDemensionParticles, 2)) + grid_y * NumberOfOneDemensionParticles + grid_z;

	return flat_index;
};

Eigen::Vector3i FlatToGrid(int flat_index)
{
	Eigen::Vector3i grid_index = {};

	grid_index[0] = flat_index / int(pow(NumberOfOneDemensionParticles, 2));
	grid_index[1] = (flat_index % int(pow(NumberOfOneDemensionParticles, 2))) / NumberOfOneDemensionParticles;
	grid_index[2] = ((flat_index % int(pow(NumberOfOneDemensionParticles, 2))) % NumberOfOneDemensionParticles);

	return grid_index;
};

double calNormalize(Eigen::VectorXd v) {

	double norm = 0.0;

	for (int i = 0; i < v.size(); i++) {
		norm += pow(v(i), 2);
	}

	return pow(norm, 0.5);

}

void showDetail(double x, double y, double z, Square s, Eigen::VectorXd v) {
	
	Eigen::MatrixXd m_add = Eigen::MatrixXd::Zero(3, 3);
	Eigen::VectorXd weight = Eigen::MatrixXd::Zero(3, 3);
	
	for (int i = 0; i < NumberOfParticles; i++) {

		// std::cout << "phi x : " << barphi(3 * i) << ", phi y : " << barphi(3 * i + 1) << ", phi z : " << barphi(3 * i + 2) << "\n";
		
		double x_dist = (x - v(3 * i)) / s.dx;
		double y_dist = (y - v(3 * i + 1)) / s.dx;
		double z_dist = (z - v(3 * i + 2)) / s.dx;

		// std::cout << "x dist : " << x_dist << ", y dist : " << y_dist << ", z dist : " << z_dist << "\n";

		// std::cout << "x dist : " << x_dist << ", y dist : " << y_dist << ", z dist : " << z_dist << "\n";

		weight(0) = N_dash(x_dist) * N(y_dist) * N(z_dist) / s.dx;

		weight(1) = N(x_dist) * N_dash(y_dist) * N(z_dist) / s.dx;

		weight(2) = N(x_dist) * N(y_dist) * N_dash(z_dist) / s.dx;

		// std::cout << "w1 : " << weight(0) << ", w2 : " << weight(1) << ", w3 : " << weight(2) << "\n";

		// std::cout << "x dist : " << weight(0) << ", y dist : " << weight(1) << ", z dist : " << weight(2) << "\n";

		if (abs(x_dist) <= 1.0 && abs(y_dist) <= 1.0 && abs(z_dist) <= 1.0) {

			Eigen::MatrixXd m = Eigen::MatrixXd::Zero(3, 3);

			for (int a = 0; a < dimensions; a++) {
				for (int b = 0; b < dimensions; b++) {

					m(a, b) = v(3 * i + a) * weight(b);
					// std::cout << "a : " << a <<  ", b : " << b  << ", m : " << m(a, b) << "\n";

				}
			}

			m_add += m;
			// std::cout << "m : " << m << "\n";

		}
	}
	

	// std::cout << "m_add : " << m_add << "\n";
	//std::cout << "determinant : " << m_add.determinant() << "\n";

	
	double W = 0.0;

	for (int i = 0; i < NumberOfParticles; i++) {

		double bar_phi_i_x = v(3 * i);
		double bar_phi_i_y = v(3 * i + 1);
		double bar_phi_i_z = v(3 * i + 2);

		double x_i_dist = (x - v(3 * i)) / s.dx;
		double y_i_dist = (y - v(3 * i + 1)) / s.dx;
		double z_i_dist = (z - v(3 * i + 2)) / s.dx;

		double weight_i_1 = N_dash(x_i_dist) * N(y_i_dist) * N(z_i_dist) / s.dx;

		for (int j = 0; j < NumberOfParticles; j++) {

			double bar_phi_j_x = v(3 * j);
			double bar_phi_j_y = v(3 * j + 1);
			double bar_phi_j_z = v(3 * j + 2);

			double x_j_dist = (x - v(3 * j)) / s.dx;
			double y_j_dist = (y - v(3 * j + 1)) / s.dx;
			double z_j_dist = (z - v(3 * j + 2)) / s.dx;

			double weight_j_2 = N(x_j_dist) * N_dash(y_j_dist) * N(z_j_dist) / s.dx;

			for (int k = 0; k < NumberOfParticles; k++) {

				double bar_phi_k_x = v(3 * k);
				double bar_phi_k_y = v(3 * k + 1);
				double bar_phi_k_z = v(3 * k + 2);

				double x_k_dist = (x - v(3 * k)) / s.dx;
				double y_k_dist = (y - v(3 * k + 1)) / s.dx;
				double z_k_dist = (z - v(3 * k + 2)) / s.dx;

				double weight_k_3 = N(x_k_dist) * N(y_k_dist) * N_dash(z_k_dist) / s.dx;

				double Phi = bar_phi_i_x * (bar_phi_j_y * bar_phi_k_z - bar_phi_j_z * bar_phi_k_y)
					+ bar_phi_i_y * (bar_phi_j_z * bar_phi_k_x - bar_phi_j_x * bar_phi_k_z)
					+ bar_phi_i_z * (bar_phi_j_x * bar_phi_k_y - bar_phi_j_y * bar_phi_k_x);

				double det = Phi * weight_i_1 * weight_j_2 *  weight_k_3;

				W += det;

				if (det != 0) {
					//std::cout << det << "\n";
				}

				// std::cout << W << "\n";
			}
		}
	}
	
	//std::cout << "determinant : " << W << "\n";

	double TW = 0.0;

	for (int i = 0; i < NumberOfParticles; i++) {

		double x_dist = (x - v(3 * i)) / s.dx;
		double y_dist = (y - v(3 * i + 1)) / s.dx;
		double z_dist = (z - v(3 * i + 2)) / s.dx;

		if (abs(x_dist) <= 1.0 && abs(y_dist) <= 1.0 && abs(z_dist) <= 1.0) {
			double Theta = s.points[i].theta;
			double w = N(x_dist) * N(y_dist) * N(z_dist);
			TW += Theta * w;
		}
	}

	//std::cout << "Theta * w : " << TW << "\n";
	

}

double N(double x) {
	double dist = abs(x);
	double kernel_value = 0.0;

	if (dist >= 1.0) {
		kernel_value = 0.0;
	}
	else if (dist >= 0.0 && dist < 1.0) {
		kernel_value = 1.0 - dist;
	}

	return kernel_value;
}

double N_dash(double x) {
	double sgn;
	if (x >= 0.0) {
		sgn = 1.0;
	}else {
		sgn = -1.0;
	}

	double dist = abs(x);
	double kernel_value = 0.0;
	if (dist >= 1.0) {
		kernel_value = 0.0;
	}else if(dist >= 0.0 && dist < 1.0) {
		kernel_value = - 1.0;
	}

	return kernel_value * sgn;
}


void calCheckMatrixPreparation(Square square) {
	Kernel_Array ka = createKernelArray3();

	// grid node
	for (int kx = 0; kx < dimensions; kx++) {
		for (int ky = 0; ky < dimensions; ky++) {
			for (int kz = 0; kz < dimensions; kz++) {

				// -1 ～ 1
				int gn_kx = kx - 1;
				int gn_ky = ky - 1;
				int gn_kz = kz - 1;

				// std::cout << "grid node index : ";
				// std::cout << "x " << gn_ix << ",y " << gn_iy << ",z " << gn_iz << "\n";

				ka = setKernelArray3ForX(ka, gn_kx);
				ka = setKernelArray3ForY(ka, gn_ky);
				ka = setKernelArray3ForZ(ka, gn_kz);

				double el = product3_x(ka) * product3_y(ka) * product3_z(ka);

				// std::cout << "result : " << el << "\n" << "\n";

				CheckMatrixPreparation[kx][ky][kz] = { el };
				if (kx == 2 && ky == 2 && kz == 1) {
					//std::cout << el << "\n";
				}
			}
		}
	}
}

Eigen::MatrixXd calCheckMatrix(Square square) {

	bool FlagCalGrid;

	Eigen::MatrixXd OutputM(square.points.size(), square.points.size());

	//calculate Q matrix
	for (int xi = 0; xi < NumberOfParticles; xi++) {
		Eigen::Vector3i grid_xi = FlatToGrid(xi);

		for (int k = 0; k < NumberOfParticles; k++) {
			Eigen::Vector3i grid_k = FlatToGrid(k);
			FlagCalGrid = false;
			OutputM(k, xi) = 0.0;
			
			Eigen::Vector3i k_minus_xi = grid_k - grid_xi;

			if ((abs(k_minus_xi[0]) <= 1) && (abs(k_minus_xi[1]) <= 1) && (abs(k_minus_xi[2]) <= 1)) {

				double w1w = CheckMatrixPreparation[k_minus_xi[0] + 1][k_minus_xi[1] + 1][k_minus_xi[2] + 1];
				OutputM(k, xi) = w1w;
				if (xi == 13 && k == 0) {
					//std::cout << OutputM(k, xi) << "\n";
				}

			}

		}
	}



	return OutputM;
}