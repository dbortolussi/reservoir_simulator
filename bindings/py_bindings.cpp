#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "../src/simulator.cpp"

using namespace pybind11::literals;
namespace py = pybind11;

PYBIND11_MODULE(pyressim, m) {
	m.doc() = "Two-phase reservoir simulator";

	py::class_<SimResult>(m, "SimResult")
		.def_readonly("time", &SimResult::time)
		.def_readonly("dt_history", &SimResult::dt_history)
		.def_readonly("newton_iters", &SimResult::newton_iters)
		.def_readonly("rollbacks", &SimResult::rollbacks)
		.def_readonly("front_sharpness", &SimResult::front_sharpness)
		.def_readonly("breakthrough_time", &SimResult::breakthrough_time)
		.def_readonly("nx", &SimResult::nx)
		.def_readonly("ny", &SimResult::ny)
		.def_property_readonly("Sw_final", [](const SimResult& r) {return r.Sw_final;})
		.def_property_readonly("p_final", [](const SimResult& r) {return r.p_final;});

	py::class_<Simulator>(m, "Simulator")
		.def(py::init<int, int, double, double, double>(),
			"nx"_a, "ny"_a, "Lx"_a, "Ly"_a, "dz"_a = 1.0)
		//.def("set_permeability", &Simulator::set_permeability, "K"_a)
		.def("set_fluid", &Simulator::set_fluid, "mu_w"_a, "mu_o"_a)
		.def("set_relperm", &Simulator::set_relperm, "Swc"_a, "Sor"_a, "n"_a, "krw0"_a = 1.0, "kro0"_a = 1.0)
		.def("set_porosity", &Simulator::set_porosity, "phi"_a)
		.def("set_cfl", &Simulator::set_cfl, "c"_a)
		.def("add_well", &Simulator::add_well, "i"_a, "j"_a, "type"_a, "p_bhp"_a)
		.def("run", &Simulator::run, "solver"_a, "t_end"_a, py::arg("vec") = py::none())
		.def_readwrite("dt_init", &Simulator::dt_init)
		.def_readwrite("dt_min", &Simulator::dt_min)
		.def("set_permeability", [](Simulator& self, py::array_t<double, py::array::c_style | py::array::forcecast> K) {
		auto buf = K.request();
		if (buf.ndim != 2) {
			throw std::runtime_error("set_permeability expects a 2D array");
		}
		Eigen::MatrixXd M(buf.shape[0], buf.shape[1]);
		std::memcpy(M.data(), buf.ptr, sizeof(double) * buf.shape[0] * buf.shape[1]);
		self.set_permeability(M);
			}, "K"_a);
}