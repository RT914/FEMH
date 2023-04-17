#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#include "fem.h"

void export_CSV(double MatrixQForX[3][3][3][3][3][3], double MatrixQForY[3][3][3][3][3][3], double MatrixQForZ[3][3][3][3][3][3])
{
	std::string output_csv_file_name = "result.csv";
	std::ofstream data_file(output_csv_file_name);

	// Header
	data_file << "kx-ƒÌx" << "," << "ky-ƒÌy" << "," << "kz-ƒÌz" << ",";
	data_file << "ix-ƒÌx" << "," << "iy-ƒÌy" << "," << "iz-ƒÌz" << ",";
	data_file << "jx-ƒÌx" << "," << "jy-ƒÌy" << "," << "jz-ƒÌz" << ",";

	// kernel index
	for (int wi = 0; wi < dimensions; wi++) {
		for (int wj = 0; wj < dimensions; wj++) {
			for (int wk = 0; wk < dimensions; wk++) {

				data_file << "w" << wi + 1 << "w" << wj + 1 << "w" << wk + 1 << ",";

			}
		}
	}

	data_file << std::endl;


	// Grid index & Value
	// grid node					
	for (int kx = 0; kx < dimensions; kx++) {
		for (int ky = 0; ky < dimensions; ky++) {
			for (int kz = 0; kz < dimensions; kz++) {

				// case i
				// -1 ` 1
				int gn_k_x = kx - 1;
				int gn_k_y = ky - 1;
				int gn_k_z = kz - 1;

				for (int ix = 0; ix < dimensions; ix++) {
					for (int iy = 0; iy < dimensions; iy++) {
						for (int iz = 0; iz < dimensions; iz++) {

							// case j
							// -1 ` 1
							int gn_i_x = ix - 1;
							int gn_i_y = iy - 1;
							int gn_i_z = iz - 1;

							for (int jx = 0; jx < dimensions; jx++) {
								for (int jy = 0; jy < dimensions; jy++) {
									for (int jz = 0; jz < dimensions; jz++) {

										// case k
										// -1 ` 1
										
										int gn_j_x = jx - 1;
										int gn_j_y = jy - 1;
										int gn_j_z = jz - 1;
					
										
										data_file << gn_k_x << "," << gn_k_y << "," << gn_k_z << ","
											<< gn_i_x << "," << gn_i_y << "," << gn_i_z << ","
											<< gn_j_x << "," << gn_j_y << "," << gn_j_z << ",";


										for (int wi = 0; wi < dimensions; wi++) {
											for (int wj = 0; wj < dimensions; wj++) {
												for (int wk = 0; wk < dimensions; wk++) {

													data_file << MatrixQForX[wi][wj][wk][ix][jx][kx] * MatrixQForY[wi][wj][wk][iy][jy][ky] * MatrixQForZ[wi][wj][wk][iz][jz][kz] << ",";

												}
											}
										}
										


										data_file << std::endl;

										
									}
								}
							}


						}
					}
				}



			}
		}
	}

	data_file.close();

}

/*
void export_CSV(double MatrixQForX[3][3][3][3][3][3], double MatrixQForY[3][3][3][3][3][3], double MatrixQForZ[3][3][3][3][3][3])
{
	std::string output_csv_file_name = "result_a.csv";

	std::ofstream data_file(output_csv_file_name);

	// Header

	data_file << "ix-ƒÌx" << "," << "iy-ƒÌy" << "," << "iz-ƒÌz" << ",";
	data_file << "jx-ƒÌx" << "," << "jy-ƒÌy" << "," << "jz-ƒÌz" << ",";
	data_file << "kx-ƒÌx" << "," << "ky-ƒÌy" << "," << "kz-ƒÌz" << ",";

	// kernel index
	for (int wi = 0; wi < dimensions; wi++) {
		for (int wj = 0; wj < dimensions; wj++) {
			for (int wk = 0; wk < dimensions; wk++) {

				data_file << "w" << wi + 1 << "w" << wj + 1 << "w" << wk + 1 << ",";

			}
		}
	}

	data_file << std::endl;


	// Grid index & Value
	// grid node					
	for (int ix = 0; ix < dimensions; ix++) {
		for (int iy = 0; iy < dimensions; iy++) {
			for (int iz = 0; iz < dimensions; iz++) {

				// case i
				// -1 ` 1
				int gn_i_x = ix - 1;
				int gn_i_y = iy - 1;
				int gn_i_z = iz - 1;

				for (int jx = 0; jx < dimensions; jx++) {
					for (int jy = 0; jy < dimensions; jy++) {
						for (int jz = 0; jz < dimensions; jz++) {

							// case j
							// -1 ` 1
							int gn_j_x = jx - 1;
							int gn_j_y = jy - 1;
							int gn_j_z = jz - 1;

							for (int kx = 0; kx < dimensions; kx++) {
								for (int ky = 0; ky < dimensions; ky++) {
									for (int kz = 0; kz < dimensions; kz++) {

										// case k
										// -1 ` 1

										int gn_k_x = kx - 1;
										int gn_k_y = ky - 1;
										int gn_k_z = kz - 1;



										
										if (gn_i_x != -1 && gn_i_y != -1 && gn_i_z != -1 && gn_j_x != -1 && gn_j_y != -1 && gn_j_z != -1 && gn_k_x != -1 && gn_k_y != -1 && gn_k_z != -1) {
											data_file << gn_i_x << "," << gn_i_y << "," << gn_i_z << ","
												<< gn_j_x << "," << gn_j_y << "," << gn_j_z << ","
												<< gn_k_x << "," << gn_k_y << "," << gn_k_z << ",";


											for (int wi = 0; wi < dimensions; wi++) {
												for (int wj = 0; wj < dimensions; wj++) {
													for (int wk = 0; wk < dimensions; wk++) {

														data_file << MatrixQForX[wi][wj][wk][ix][jx][kx] * MatrixQForY[wi][wj][wk][iy][jy][ky] * MatrixQForZ[wi][wj][wk][iz][jz][kz] << ",";

													}
												}
											}

											data_file << std::endl;
										}
										


										data_file << gn_i_x << "," << gn_i_y << "," << gn_i_z << ","
											<< gn_j_x << "," << gn_j_y << "," << gn_j_z << ","
											<< gn_k_x << "," << gn_k_y << "," << gn_k_z << ",";


										for (int wi = 0; wi < dimensions; wi++) {
											for (int wj = 0; wj < dimensions; wj++) {
												for (int wk = 0; wk < dimensions; wk++) {

													data_file << MatrixQForX[wi][wj][wk][ix][jx][kx] * MatrixQForY[wi][wj][wk][iy][jy][ky] * MatrixQForZ[wi][wj][wk][iz][jz][kz] << ",";

												}
											}
										}



										data_file << std::endl;


									}
								}
							}


						}
					}
				}



			}
		}
	}

	data_file.close();

}
*/

void export_CSV_X(double MatrixQForX[3][3][3][3][3][3])
{
	std::string output_csv_file_name = "result_x.csv";

	std::ofstream data_file(output_csv_file_name);

	// Header

	data_file << "ix-ƒÌx" << ",";
	data_file << "jx-ƒÌx" << ",";
	data_file << "kx-ƒÌx" << ",";

	// kernel index
	for (int wi = 0; wi < dimensions; wi++) {
		for (int wj = 0; wj < dimensions; wj++) {
			for (int wk = 0; wk < dimensions; wk++) {

				data_file << "w" << wi + 1 << "w" << wj + 1 << "w" << wk + 1 << ",";

			}
		}
	}

	data_file << std::endl;


	// Grid index & Value
	// grid node					
	for (int ix = 0; ix < dimensions; ix++) {
		// case i
		// -1 ` 1
		int gn_i_x = ix - 1;

		for (int jx = 0; jx < dimensions; jx++) {
			// case j
			// -1 ` 1
			int gn_j_x = jx - 1;

			for (int kx = 0; kx < dimensions; kx++) {
				// case k
				// -1 ` 1
				int gn_k_x = kx - 1;

				data_file << gn_i_x << "," << gn_j_x << "," << gn_k_x << ",";


				for (int wi = 0; wi < dimensions; wi++) {
					for (int wj = 0; wj < dimensions; wj++) {
						for (int wk = 0; wk < dimensions; wk++) {

							data_file << MatrixQForX[wi][wj][wk][ix][jx][kx] << ",";

						}
					}
				}

				data_file << std::endl;
			}
		}

	}

	data_file.close();

}

void export_CSV_Y(double MatrixQForY[3][3][3][3][3][3])
{
	std::string output_csv_file_name = "result_y.csv";

	std::ofstream data_file(output_csv_file_name);

	// Header

	data_file << "iy-ƒÌy" << ",";
	data_file << "jy-ƒÌy" << ",";
	data_file << "ky-ƒÌy" << ",";

	// kernel index
	for (int wi = 0; wi < dimensions; wi++) {
		for (int wj = 0; wj < dimensions; wj++) {
			for (int wk = 0; wk < dimensions; wk++) {

				data_file << "w" << wi + 1 << "w" << wj + 1 << "w" << wk + 1 << ",";

			}
		}
	}

	data_file << std::endl;


	// Grid index & Value
	// grid node					
	for (int iy = 0; iy < dimensions; iy++) {
		// case i
		// -1 ` 1
		int gn_i_y = iy - 1;

		for (int jy = 0; jy < dimensions; jy++) {
			// case j
			// -1 ` 1
			int gn_j_y = jy - 1;

			for (int ky = 0; ky < dimensions; ky++) {
				// case k
				// -1 ` 1
				int gn_k_y = ky - 1;

				data_file << gn_i_y << "," << gn_j_y << "," << gn_k_y << ",";


				for (int wi = 0; wi < dimensions; wi++) {
					for (int wj = 0; wj < dimensions; wj++) {
						for (int wk = 0; wk < dimensions; wk++) {

							data_file << MatrixQForY[wi][wj][wk][iy][jy][ky] << ",";

						}
					}
				}

				data_file << std::endl;
			}
		}

	}

	data_file.close();

}

void export_CSV_Z(double MatrixQForZ[3][3][3][3][3][3])
{
	std::string output_csv_file_name = "result_z.csv";

	std::ofstream data_file(output_csv_file_name);

	// Header

	data_file << "iz-ƒÌz" << ",";
	data_file << "jz-ƒÌz" << ",";
	data_file << "kz-ƒÌz" << ",";

	// kernel index
	for (int wi = 0; wi < dimensions; wi++) {
		for (int wj = 0; wj < dimensions; wj++) {
			for (int wk = 0; wk < dimensions; wk++) {

				data_file << "w" << wi + 1 << "w" << wj + 1 << "w" << wk + 1 << ",";

			}
		}
	}

	data_file << std::endl;


	// Grid index & Value
	// grid node					
	for (int iz = 0; iz < dimensions; iz++) {
		// case i
		// -1 ` 1
		int gn_i_z = iz - 1;

		for (int jz = 0; jz < dimensions; jz++) {
			// case j
			// -1 ` 1
			int gn_j_z = jz - 1;

			for (int kz = 0; kz < dimensions; kz++) {
				// case k
				// -1 ` 1
				int gn_k_z = kz - 1;

				data_file << gn_i_z << "," << gn_j_z << "," << gn_k_z << ",";


				for (int wi = 0; wi < dimensions; wi++) {
					for (int wj = 0; wj < dimensions; wj++) {
						for (int wk = 0; wk < dimensions; wk++) {

							data_file << MatrixQForZ[wi][wj][wk][iz][jz][kz] << ",";

						}
					}
				}

				data_file << std::endl;
			}
		}

	}

	data_file.close();

}

void export2_CSV(double VectorbPreparation[3][3][3])
{
	std::string output_csv_file_name = "result2.csv";

	std::ofstream data_file(output_csv_file_name);

	// Header

	data_file << "ix-ƒÌx" << "," << "iy-ƒÌy" << "," << "iz-ƒÌz" << ",";

	data_file << std::endl;


	// Grid index & Value
	// grid node
	for (int ix = 0; ix < dimensions; ix++) {
		for (int iy = 0; iy < dimensions; iy++) {
			for (int iz = 0; iz < dimensions; iz++) {

				// -1 ` 1
				int gn_ix = ix - 1;
				int gn_iy = iy - 1;
				int gn_iz = iz - 1;


				data_file << gn_ix << "," << gn_iy << "," << gn_iz << ",";
				data_file << VectorbPreparation[ix][iy][iz] << ",";

				data_file << std::endl;
			}
		}
	}

	data_file.close();

}

void export_Q(Eigen::MatrixXd m, std::string str)
{

	std::ofstream data_file(str);

	// Header
	data_file << "matrix" << ",";
	for (int i = 0; i < m.cols(); i++) {

		data_file << i << ",";

	}
	data_file << std::endl;

	for (int i = 0; i < m.rows(); i++) {

		data_file << i << ",";

		for (int j = 0; j < m.cols(); j++) {

			// Header
			data_file << m(i, j) << ",";

		}
		data_file << std::endl;
	}

	data_file.close();

}

void export_loop_times(std::vector<int> loop_times, std::string str, int num) {

	std::ofstream data_file(str);
	// Header
	data_file << "Loop Times" << ",";
	for (int i = 0; i < num; i++) {

		data_file << i << ",";

	}
	data_file << std::endl;

	// Header
	data_file << "Value" << ",";
	for (int i = 0; i < num; i++) {

		data_file << loop_times[i] << ",";
		
	}
	data_file << std::endl;

	data_file.close();
}