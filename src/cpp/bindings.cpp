#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "Vector2D.hpp"
#include "Particle.hpp"
#include "MDSimulation.hpp"
#include "PotentialsAndForces.hpp"

namespace py = pybind11;

PYBIND11_MODULE(simple_md, m) {
    m.doc() = "2D Molecular Dynamics simulation with modified Lennard-Jones potential";

    py::class_<Vector2D>(m, "Vector2D")
        .def(py::init<double, double>(), py::arg("x") = 0.0, py::arg("y") = 0.0)
        .def_readwrite("x", &Vector2D::x)
        .def_readwrite("y", &Vector2D::y)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / double())
        .def(-py::self)
        .def(py::self += py::self)
        .def(py::self -= py::self)
        .def(py::self *= double())
        .def(py::self /= double())
        .def("dot", &Vector2D::dot)
        .def("norm", &Vector2D::norm)
        .def("norm_sq", &Vector2D::norm_sq)
        .def("__repr__", [](const Vector2D& v) {
            return "Vector2D(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ")";
        });

    py::class_<Particle>(m, "Particle")
        .def_readwrite("position", &Particle::position)
        .def_readwrite("prev_position", &Particle::prev_position)
        .def_readwrite("velocity", &Particle::velocity)
        .def_readwrite("force", &Particle::force);

    py::class_<MDSimulation>(m, "MDSimulation")
        .def(py::init<const std::vector<Vector2D>&,
                       const std::vector<Vector2D>&,
                       double, double>(),
             py::arg("positions"), py::arg("velocities"),
             py::arg("L"), py::arg("dt"))
        .def_readwrite("L", &MDSimulation::L)
        .def_readwrite("dt", &MDSimulation::dt)
        .def_readonly("num_particles", &MDSimulation::num_particles)
        .def_readonly("step_count", &MDSimulation::step_count)
        .def("compute_forces", &MDSimulation::compute_forces)
        .def("step", &MDSimulation::step)
        .def("apply_periodic_bc", &MDSimulation::apply_periodic_bc)
        .def("shortest_distance", &MDSimulation::shortest_distance)
        .def("zero_velocities", &MDSimulation::zero_velocities)
        .def("get_mean_velocity", &MDSimulation::get_mean_velocity)
        .def("remove_mean_velocity", &MDSimulation::remove_mean_velocity)
        .def("get_positions", &MDSimulation::get_positions)
        .def("get_velocities", &MDSimulation::get_velocities)
        .def("get_kinetic_energy", &MDSimulation::get_kinetic_energy)
        .def("get_temperature", &MDSimulation::get_temperature)
        .def("get_potential_energy", &MDSimulation::get_potential_energy)
        .def("get_total_energy", &MDSimulation::get_total_energy)
        .def("get_pressure", &MDSimulation::get_pressure);

    m.def("lj_potential", &potentials::lj_potential, py::arg("x"));
    m.def("lj_force_magnitude_signed", &potentials::lj_force_magnitude_signed, py::arg("x"));
    m.def("pairwise_force", &potentials::pairwise_force, py::arg("r_ij"));
}
