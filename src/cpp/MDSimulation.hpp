#ifndef MD_SIMULATION_HPP
#define MD_SIMULATION_HPP

#include <vector>
#include "Vector2D.hpp"
#include "Particle.hpp"

class MDSimulation {
public:
    std::vector<Particle> particles;
    double L;
    double dt;
    std::size_t num_particles;
    std::size_t step_count;

    MDSimulation(const std::vector<Vector2D>& positions,
                 const std::vector<Vector2D>& velocities,
                 double L, double dt);

    void compute_forces();
    void step();
    void apply_periodic_bc();
    Vector2D shortest_distance(const Vector2D& r_ij) const;
    void zero_velocities();
    Vector2D get_mean_velocity() const;
    void remove_mean_velocity();

    std::vector<Vector2D> get_positions() const;
    std::vector<Vector2D> get_velocities() const;
    double get_kinetic_energy() const;
    double get_temperature() const;
    double get_potential_energy() const;
    double get_total_energy() const;
    double get_pressure() const;
};

#endif // MD_SIMULATION_HPP
