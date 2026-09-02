#pragma once
#include <vector>
#include <cmath>
#include <iostream>

enum class Face { XMIN, XMAX, YMIN, YMAX, INTERIOR };

struct Neighbor {
	int id;
	double area;
	double dist;
};

class Grid {
public:
	const int nx, ny;
	const double dx, dy, dz;
	Grid(int nx_, int ny_, double Lx, double Ly, double dz_ = 1.0) :
		nx(nx_), ny(ny_), dx(Lx / nx_), dy(Ly / ny_), dz(dz_) {
	}

	int idx(int i, int j) const { return i * ny + j; }
	int numCells() const { return nx * ny; }
	void ij(int c, int& i, int& j) const { i = c / ny; j = c % ny; }

	std::pair<double, double> center(int id) const {
		int i, j; ij(id, i, j);
		return { (i + 0.5) * dx, (j + 0.5) * dy };
	}

	double cellVolume() const { return dx * dy * dz; }

	std::vector<Neighbor> neighbors(int id) const {
		int i, j; ij(id, i, j);
		std::vector<Neighbor> nb;
		if (i > 0)      nb.push_back({ idx(i - 1, j), dy * dz, dx });
		if (i < nx - 1) nb.push_back({ idx(i + 1, j), dy * dz, dx });
		if (j > 0)      nb.push_back({ idx(i, j - 1), dx * dz, dy });
		if (j < ny - 1) nb.push_back({ idx(i, j + 1), dx * dz, dy });
		return nb;
	}

	bool isInterior(const Neighbor& n) const { return n.id >= 0; }
};

inline double transmissibility(const Neighbor& f, double k_a, double k_b) {
	double kh = 2.0 * k_a * k_b / (k_a + k_b);
	return f.area * kh / f.dist;
}



