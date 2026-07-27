#include "Particle.hpp"
#include "Vector2D.hpp"

Particle::Particle()
    :
    position{0.0, 0.0},
    prev_position(0.0, 0.0),
    velocity(0.0, 0.0),
    force(0.0, 0.0)
{}

Particle::Particle(const Vector2D& _position, const Vector2D& _velocity, double dt)
    : 
    position{_position}, 
    prev_position{_position - _velocity * dt},
    velocity{_velocity},
    force{0.0, 0.0}
{}
