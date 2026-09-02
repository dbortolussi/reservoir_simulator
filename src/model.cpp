#pragma once
#include <vector>
#include "grid.cpp"
#include "fluid.cpp"
#include "wells.cpp"

struct Model {
	Grid				grid;
	Fluid				fluid;
	std::vector<double>	perm;
	std::vector<Well>	wells;
	double				phi = 0.2;
	double				cfl_safety = 0.2;

	explicit Model(const Grid& g) :
		grid(g), perm(g.numCells(), 1.0e-3) {
	}

	double k(int cell) const { return perm[cell]; }

	const Well* wellAt(int cell) const {
		for (const auto& w : wells) {
			if (w.cell == cell) return &w;
		}
		return nullptr;
	}
};