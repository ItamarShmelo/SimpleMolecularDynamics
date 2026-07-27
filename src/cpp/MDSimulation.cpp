#include "MDSimulation.hpp"
#include "Particle.hpp"
#include "PotentialsAndForces.hpp"
#include "Vector2D.hpp"
#include <cmath>

MDSimulation::MDSimulation(const std::vector<Vector2D>& positions,
                           const std::vector<Vector2D>& velocities,
                           double _L, double _dt)
                           :
                           particles(positions.size(), Particle{}),
                           L{_L}, dt{_dt},
                           num_particles{particles.size()},
                           step_count{0UL}
{
    for (std::size_t i=0; i < num_particles; ++i){
        particles[i] = Particle(positions[i], velocities[i], dt);
    }
}

void MDSimulation::compute_forces()
{
    for(auto& p : particles){
        p.force = {0.0, 0.0};
    }

    for (std::size_t i = 0; i < num_particles; ++i){
        auto& p_i = particles[i];
        for(std::size_t j = i+1; j < num_particles; ++j){
            auto& p_j = particles[j];
            auto const r_ij = shortest_distance(p_i.position - p_j.position);

            auto const force_ij = potentials::pairwise_force(r_ij);
            p_i.force += force_ij;
            p_j.force -= force_ij;
        }
    }
}

void MDSimulation::step()
{
    compute_forces();

    double dt2 = dt * dt;
    for (auto& p : particles) {
        Vector2D disp = shortest_distance(p.position - p.prev_position);
        Vector2D new_pos = p.position + disp + p.force * dt2;
        p.velocity = (disp * 2.0 + p.force * dt2) / (2.0 * dt);
        p.prev_position = p.position;
        p.position = new_pos;
    }

    apply_periodic_bc();
    ++step_count;
}

void MDSimulation::apply_periodic_bc()
{
    auto fmod_positive = [&](double x) {
        x = std::fmod(x, L);
        if (x < 0.0) x += L;
        return x;
    };
    for(auto& p : particles){
        p.position = {fmod_positive(p.position.x), fmod_positive(p.position.y)};
    }
}

Vector2D MDSimulation::shortest_distance(const Vector2D& r_ij) const
{
    auto shortest_ij{r_ij};
    
    if(r_ij.x > 0.5*L) shortest_ij.x -= L;
    else if(r_ij.x < -0.5*L) shortest_ij.x += L;
    
    if(r_ij.y > 0.5*L) shortest_ij.y -= L;
    else if(r_ij.y < -0.5*L) shortest_ij.y += L;

    return shortest_ij;
}

void MDSimulation::zero_velocities()
{
    for (auto& p : particles){
        p.prev_position = p.position;
        p.velocity = {0.0, 0.0};
    }
}

Vector2D MDSimulation::get_mean_velocity() const
{
    Vector2D sum(0.0, 0.0);
    for (const auto& p : particles) {
        sum += p.velocity;
    }
    return sum / static_cast<double>(num_particles);
}

void MDSimulation::remove_mean_velocity()
{
    Vector2D mean = get_mean_velocity();
    for (auto& p : particles) {
        p.velocity -= mean;
    }
}

std::vector<Vector2D> MDSimulation::get_positions() const
{
    std::vector<Vector2D> pos;
    pos.reserve(num_particles);
    for (const auto& p : particles) {
        pos.push_back(p.position);
    }
    return pos;
}

std::vector<Vector2D> MDSimulation::get_velocities() const
{
    std::vector<Vector2D> vel;
    vel.reserve(num_particles);
    for (const auto& p : particles) {
        vel.push_back(p.velocity);
    }
    return vel;
}

double MDSimulation::get_kinetic_energy() const
{
    double K = 0.0;
    for (const auto& p : particles) {
        K += 0.5 * p.velocity.norm_sq();
    }
    return K;
}

double MDSimulation::get_temperature() const
{
    Vector2D v_mean = get_mean_velocity();
    double K_thermal = 0.0;
    for (const auto& p : particles) {
        Vector2D v_rel = p.velocity - v_mean;
        K_thermal += 0.5 * v_rel.norm_sq();
    }
    return K_thermal / num_particles;
}

double MDSimulation::get_potential_energy() const
{
    double U = 0.0;
    for (std::size_t i = 0; i < num_particles; ++i) {
        for (std::size_t j = i + 1; j < num_particles; ++j) {
            Vector2D r_ij = shortest_distance(
                particles[i].position - particles[j].position);
            U += potentials::lj_potential(r_ij.norm());
        }
    }
    return U;
}

double MDSimulation::get_total_energy() const
{
    return get_kinetic_energy() + get_potential_energy();
}

double MDSimulation::get_pressure() const
{
    double forces_term = 0.0;
    for (std::size_t i = 0; i < num_particles; ++i) {
        for (std::size_t j = i + 1; j < num_particles; ++j) {
            Vector2D r_ij = shortest_distance(
                particles[i].position - particles[j].position);
            double x = r_ij.norm();
            forces_term += potentials::lj_force_magnitude_signed(x) * x;
        }
    }

    Vector2D v_mean = get_mean_velocity();
    double kinetic_sum = 0.0;
    for (const auto& p : particles) {
        Vector2D v_rel = p.velocity - v_mean;
        kinetic_sum += v_rel.norm_sq();
    }

    double V = L * L;
    return (forces_term + kinetic_sum) / (2.0 * V);
}
