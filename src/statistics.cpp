#pragma once
#include <vector>
#include <Eigen/Dense>

struct SimResult {
	std::vector<double> time;
	std::vector<double> dt_history;
	std::vector<int> newton_iters;
	std::vector<int> rollbacks;
	std::vector<double> front_sharpness;
	double breakthrough_time = -1.0;
	int nx = 0, ny = 0;
	Eigen::VectorXd Sw_final;
	Eigen::VectorXd p_final;
};