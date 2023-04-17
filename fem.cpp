#include <GL/freeglut.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>

#include "draw.h"
#include "square.h"
#include "fem.h"
#include "NewtonRaphsonMethod.h"
#include "export.h"

Square square = createSquare(NumberOfOneDemensionParticles);

// fem_for_key用
int calculation_times = 0;


double dt = 1.0e-3;
Eigen::Vector3d gravity{ 0.0, 0.0, -9.81 };



void calVelocity()
{
	for (int i = 0; i < pow(NumberOfOneDemensionParticles, 3); i++) {
		square.points[i].velocity = square.points[i].velocity + gravity * dt;

	}
}


Eigen::Vector3d calConflict(Eigen::Vector3d vel, Eigen::Vector3d pos)
{
	if (pos.z() <= 0.001) {
		return { vel.x(), vel.y(), 0.0 };
	}

	return vel;
};


void calPosition()
{
	for (int i = 0; i < pow(NumberOfOneDemensionParticles, 3); i++) {
		square.points[i].velocity = calConflict(square.points[i].velocity, square.points[i].position); //衝突判定・計算
		square.points[i].position = square.points[i].position + square.points[i].velocity * dt;
	}
};


void fem(int SimulationTime)
{

	if (SimulationTime == 1)
	{
		
		calMatrixQPreparation(square);
		calVectorbPreparation(square);
		//calCheckMatrixPreparation(square);

		
		int loop_n = 6;
		std::vector<int> loop_times;
		
		for (int i = 0; i < loop_n; i++) {
			std::cout << "Loop_times : " << i << "\n";
			int t = Newton_loop(square, i);
			loop_times.emplace_back(t);
		}
		std::cout << "All Clear!!" << "\n";
	
		std::string output_csv_file_name = "result_loop_times.csv";
		export_loop_times(loop_times, output_csv_file_name, loop_n);
		

		/*
		Eigen::VectorXd new_phi = Newton(square);

		
		for(int i = 0;i < NumberOfParticles;i++){
			square.points[i].position[0] = new_phi(3 * i);
			square.points[i].position[1] = new_phi(3 * i + 1);
			square.points[i].position[2] = new_phi(3 * i + 2);
		}
		*/
		
		
		
		// Newton_copy(square);
		
		
		/*
		Eigen::MatrixXd m(2, 2);
		m(0, 0) = 3;
		m(0, 1) = 0;
		m(1, 0) = 0;
		m(1, 1) = 4;
		Eigen::Vector2d v;
		v << 3, 8;
		Eigen::FullPivLU<Eigen::MatrixXd> LU(m);
		Eigen::VectorXd hi = LU.solve(v);
		// Eigen::VectorXd hi = m.inverse() * v;
		std::cout << hi << "\n";
		*/
		
		/*
		Eigen::Vector3d DeltaX;
		DeltaX << 0.0, 0.0, 0.0;
		Eigen::Vector3d b;
		b << 0.1, 0.0, 0.1;
		Eigen::Matrix3d Q;
		Q << 1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0;
		

		Eigen::FullPivLU<Eigen::MatrixXd> LU(Q);
		Eigen::VectorXd answer = LU.solve(b);
		std::cout << answer << "\n";
		*/

		/*
		//計算確認用
		Eigen::VectorXd bp(375);
		Eigen::VectorXd b(125);
		Eigen::VectorXd delta(375);
		for (int i = 0; i < 125; i++) {
			Eigen::Vector3i grid = FlatToGrid(i);
			// x
			bp(3 * i) = grid(0) * 2 - 2;
			bp(3 * i + 1) = grid(1) - 2;
			bp(3 * i + 2) = grid(2) - 2;
			//std::cout << i << "\n";
			//std::cout << Delta(3 * i) << "\n";
			//std::cout << Delta(3 * i + 1) << "\n";
			//std::cout << Delta(3 * i + 2) << "\n";
			b(i) = 1; //Theta - 1
			delta(3 * i) = double(int(i / 25));
			delta(3 * i + 1) = 0;
			delta(3 * i + 2) = 0;
			//std::cout << i << "\n";
			//std::cout << delta(3 * i) << "\n";
			//std::cout << delta(3 * i + 1) << "\n";
			//std::cout << delta(3 * i + 2) << "\n";
			if (i == 31) {
				std::cout << delta(3 * i) << "\n";
				std::cout << b(i) << "\n";
			}
			
		}
		*/
		
	}

	

	//calVelocity();
	//calPosition();
	
	glColor3f(0.5, 0.0, 0.0);
	drawSquare(square);
	Ground();
	
};

void fem_for_key(int SimulationTime, bool key)
{

	if (key == true)
	{
		
		if (calculation_times == 0) {
			calMatrixQPreparation(square);
			calVectorbPreparation(square);
			//calVectorb_part2_Preparation(square);
			//calCheckMatrixPreparation(square);
			calculation_times++;
		}

		Eigen::VectorXd new_phi = Newton_one(square);

		for (int i = 0; i < NumberOfParticles; i++) {
			square.points[i].position[0] = new_phi(3 * i);
			square.points[i].position[1] = new_phi(3 * i + 1);
			square.points[i].position[2] = new_phi(3 * i + 2);
		}
	}

	glColor3f(0.5, 0.0, 0.0);
	drawSquare(square);
	Ground();
};
