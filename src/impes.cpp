#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Sparse>
#include "model.cpp"
#include <iostream>

inline Eigen::VectorXd impesPressure(const Model& m, const std::vector<double>& Sw) {
	const Grid& g = m.grid;
	const Fluid& fl = m.fluid;
	int N = g.numCells();
	std::vector<Eigen::Triplet<double>> trip;
	trip.reserve(N * 5);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(N);

	

	for (int c = 0;c < N;++c) {
		const Well* w = m.wellAt(c);
		if (w) {
			trip.emplace_back(c, c, 1.0);
			b(c) = w->p_bhp;
			continue;
		}
		double diag = 0.0;
		for (const Neighbor& f : g.neighbors(c)) {
			int j = f.id;
			double T = transmissibility(f, m.k(c), m.k(j));
			double lam = 0.5 * (fl.lambda_t(Sw[c]) + fl.lambda_t(Sw[j]));
			double Tl = T * lam;
			diag += Tl;
			trip.emplace_back(c, j, -Tl);
		}
		trip.emplace_back(c, c, diag);
	}
	
	Eigen::SparseMatrix<double> A(N, N);
	
	A.setFromTriplets(trip.begin(), trip.end());

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);

	return solver.solve(b);
}

inline double impesStep(const Model& m, std::vector<double>& Sw, double dt_try) {
	
	const Grid& g = m.grid;
	
	const Fluid& fl = m.fluid;
	
	int N = g.numCells();
	
	double V = g.cellVolume();

	Eigen::VectorXd p = impesPressure(m, Sw);

	double maxSpeed = 0.0;
	for (int c = 0; c < N; ++c) {
		for (const Neighbor& f : g.neighbors(c)) {
			int j = f.id;
			if (j < c) continue;
			double T = transmissibility(f, m.k(c), m.k(j));
			double dp = p(c) - p(j);
			int up = (dp >= 0.0) ? c : j;
			double F = T * fl.lambda_t(Sw[up]) * dp;
			double ds = Sw[c] - Sw[j];
			double dfw = (std::abs(ds) > 1e-12) ? (fl.fw(Sw[c]) - fl.fw(Sw[j])) / ds : fl.dfw_dSw(Sw[up]);
			double spd = std::abs(F) * std::abs(dfw);
			if (spd > maxSpeed) maxSpeed = spd;
		}
	}
	double dt = dt_try;
	if (maxSpeed > 0.0) {
		dt = std::min(dt_try, m.cfl_safety * m.phi * V / maxSpeed);
	}

	std::vector<double> dSw(N, 0.0), dTot(N, 0.0);
	for(int c=0; c<N; ++c){
		for (const Neighbor& f: g.neighbors(c)) {
			int j = f.id;
			if (j < c) continue;
			double T = transmissibility(f, m.k(c), m.k(j));
			double dp = p(c) - p(j);
			int up = (dp >= 0.0) ? c : j;
			double F = T * fl.lambda_t(Sw[up]) * dp;
			double Fw = fl.fw(Sw[up]) * F;
			dSw[c] -= Fw;
			dSw[j] += Fw;
			dTot[c] -= F;
			dTot[j] += F;
		}
	}
	for (int c = 0; c < N; ++c) {
		const Well* w = m.wellAt(c);
		if (w && w->type == Well::Injector) {
			Sw[c] = 1.0 - fl.kr.Sor;
			continue;
		}
		
		double src = (w && w->type == Well::Producer) ? fl.fw(Sw[c]) * dTot[c] : 0.0;
		Sw[c] += dt / (m.phi * V) * (dSw[c] - src);
		double lo = fl.kr.Swc;
		double hi = 1.0 - fl.kr.Sor;
		Sw[c] = std::min(hi, std::max(lo, Sw[c]));
	}
	return dt;
}