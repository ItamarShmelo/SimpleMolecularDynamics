#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "Vector2D.hpp"

class Particle {
public:
    Vector2D position;
    Vector2D prev_position;
    Vector2D velocity;
    Vector2D force;

    Particle();
    Particle(const Vector2D& position, const Vector2D& velocity, double dt);
};

#endif // PARTICLE_HPP
