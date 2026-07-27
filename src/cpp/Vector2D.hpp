#ifndef VECTOR2D_HPP
#define VECTOR2D_HPP

#include <cmath>
struct Vector2D {
    double x;
    double y;

    Vector2D(double _x, double _y) : x{_x}, y{_y} {}

    Vector2D operator+(const Vector2D& other) const{
        return Vector2D{x + other.x, y + other.y};
    }
    
    Vector2D operator-(const Vector2D& other) const{
        return Vector2D{x - other.x, y - other.y};
    }
    
    Vector2D operator*(double scalar) const {
        return Vector2D{scalar*x, scalar*y};
    }

    Vector2D operator/(double scalar) const {
        return Vector2D{x / scalar, y / scalar};
    }

    Vector2D operator-() const {
        return Vector2D{-x, -y};
    }

    Vector2D& operator+=(const Vector2D& other){
        x += other.x;
        y += other.y;
        return *this;
    }

    Vector2D& operator-=(const Vector2D& other){
        x -= other.x;
        y -= other.y;
        return *this;
    }

    Vector2D& operator*=(double scalar){
        x *= scalar;
        y *= scalar;
        return *this;
    }

    Vector2D& operator/=(double scalar) {
        x /= scalar;
        y /= scalar;
        return *this;
    }

    double dot(const Vector2D& other) const {
        return x*other.x + y*other.y;
    }

    
    double norm_sq() const {
        return dot(*this);
    }

    double norm() const{
        return std::sqrt(norm_sq());
    }
};

inline Vector2D operator*(double scalar, const Vector2D& v){
    return v*scalar;
}

#endif // VECTOR2D_HPP
