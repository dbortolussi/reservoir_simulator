#pragma once
#include <string>
#include <format>
#include <cmath>
#include <stdexcept>
#include <Eigen/Dense>
#include "model.cpp"
#include "statistics.cpp"
#include "impes.cpp"
#include "newton.cpp"
#include <iostream>
#include <optional>

class Simulator {
public :
	Model model;
	double dt_init;
	double dt_min;

	Simulator(int nx, int ny, double Lx, double Ly, double dz = 1.0):
		model(Grid(nx, ny, Lx, Ly, dz)), dt_init(8.64e4), dt_min(1.0){
	}


	void set_permeability(const Eigen::MatrixXd& K) {
		const Grid& g = model.grid;
		if (K.rows() != g.nx || K.cols() != g.ny) {
			throw std::runtime_error(std::format("permeability array must be {} by {}", g.nx, g.ny));
		}
		for (int i = 0; i < g.nx; ++i) {
			for (int j = 0; j < g.ny; ++j) {
				model.perm[g.idx(i, j)] = K(i, j);
			}
		}
	}

	void set_fluid(double mu_w, double mu_o) {
		model.fluid.mu_w = mu_w; model.fluid.mu_o = mu_o;
	}

	void set_relperm(double Swc, double Sor, double n, double krw0 = 1.0, double kro0 = 1.0) {
		model.fluid.kr = BrooksCorey{ krw0, kro0, Swc, Sor, n };
	}

	void set_porosity(double phi) { model.phi = phi; }

	void set_cfl(double c) { model.cfl_safety = c; }

	void add_well(int i, int j, const std::string& type, double p_bhp) {
		Well w;
		w.cell = model.grid.idx(i, j);
		w.p_bhp = p_bhp;
		w.type = (type == "injector") ? Well::Injector : Well::Producer;
		model.wells.push_back(w);
	}

	SimResult run(const std::string& solver, double t_end, std::optional<std::vector<double>> Sw_init = std::nullopt) {
		const Grid& g = model.grid;
		int N = g.numCells();
		SimResult res; res.nx = g.nx; res.ny = g.ny;
		double p_inj = 0, p_prd = 0;
		std::vector<double> Sw = Sw_init ? std::move(*Sw_init) : std::vector<double>(N, model.fluid.kr.Swc);
		Eigen::VectorXd p = Eigen::VectorXd::Zero(N);
		int injector_cell = -1, producer_cell = -1;
		for (const Well& w : model.wells) {
			if (w.type == Well::Injector) injector_cell = w.cell;
			else producer_cell = w.cell;
		}
		if (injector_cell >= 0 && producer_cell >= 0) {
			p_inj = model.wellAt(injector_cell)->p_bhp;
			p_prd = model.wellAt(producer_cell)->p_bhp;
			
			int i_inj, j_inj, i_prod, j_prod;
			g.ij(injector_cell, i_inj, j_inj);
			g.ij(producer_cell, i_prod, j_prod);
			double S_high = 1.0 - model.fluid.kr.Sor;
			double S_low = model.fluid.kr.Swc;
			for (int c = 0; c < N; ++c) {
				int i, j; g.ij(c, i, j);
				std::fill(Sw.begin(), Sw.end(), model.fluid.kr.Swc);

				double frac = double(i + j) / double((model.grid.nx - 1) + (model.grid.ny - 1));
				p(c) = p_inj + frac * (p_prd - p_inj);
			}
			for (const Well& w : model.wells)
				if (w.type == Well::Injector)
					Sw[w.cell] = 1.0 - model.fluid.kr.Sor;
		}
		else {
			std::fill(Sw.begin(), Sw.end(), model.fluid.kr.Swc);
		}
		for (const Well& w : model.wells) {
			if (w.type == Well::Injector) {
				Sw[w.cell] = 1.0 - model.fluid.kr.Sor;
			}
			else {
				producer_cell = w.cell;
			}
		}
		double Sbt = model.fluid.kr.Swc + 0.1 * (1.0 - model.fluid.kr.Sor - model.fluid.kr.Swc);
		double t = 0.0;
		double dt = dt_init;
		

		while (t < t_end) {
			int nit = 0; 
			int rolls = 0;
			double dt_used;

			if (solver == "impes") {
				dt_used = impesStep(model, Sw, std::min(dt, t_end - t));
			}
			else {
				std::vector<double> Sw_old = Sw;
				Eigen::VectorXd p_old = p;
				double dt_try = std::min(dt, t_end - t);
				Eigen::VectorXd x(2 * N);
				for (;;) {
					Sw = Sw_old;
					p = p_old;
					for (int c = 0;c < N;++c) {
						x(2 * c) = p(c);
						x(2 * c + 1) = Sw[c];
					}
					NewtonReport rep = newtonSolve(model, x, Sw_old, dt_try);
					
					nit = rep.iters;
					if (rep.converged) {
						for (int c = 0; c < N; ++c) {
							p(c) = x(2 * c);
							Sw[c] = x(2 * c + 1);
						}
						dt_used = dt_try;
						dt = std::min(dt_init, dt_try * 1.5);
						break;
					}
					dt_try *= 0.5;
					rolls++;
					if (dt_try < dt_min) {
						throw std::runtime_error(std::format("Timestep fell below minimum {}", dt_min));
					}
				}
			}
			t += dt_used;
			res.time.push_back(t);
			res.dt_history.push_back(dt_used);
			res.newton_iters.push_back(nit);
			res.rollbacks.push_back(rolls);
			res.front_sharpness.push_back(frontSharpness(Sw));
			if (res.breakthrough_time < 0.0 && producer_cell >= 0.0 && Sw[producer_cell] > Sbt) {
				res.breakthrough_time = t;
			}
		}
		if (solver == "impes") {
			p = impesPressure(model, Sw);
		}
		res.Sw_final = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Sw.data(), Sw.size());
		res.p_final = p;
		return res;
		
	}

private:
	double frontSharpness(const std::vector<double>& Sw) const {
		const Grid& g = model.grid;
		double mx = 0.0;
		for (int i = 1; i < g.nx - 1; ++i) {
			for (int j = 1; j < g.ny - 1;++j) {
				double gx = 0.5 * (Sw[g.idx(i + 1, j)] - Sw[g.idx(i - 1, j)]) / g.dx;
				double gy = 0.5 * (Sw[g.idx(i, j + 1)] - Sw[g.idx(i, j - 1)]) / g.dy;
				double mag = std::sqrt(gx * gx + gy * gy);
				mx = std::max(mx, mag);
			}
		}
		return mx;
	}

};
