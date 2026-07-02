#ifndef VECTOR_H
#define VECTOR_H

struct Vector2D {
	double x;
	double y;
};

Vector2D operator+(const Vector2D& lhs, const Vector2D& rhs) {
	Vector2D out;

	out.x = lhs.x + rhs.x;
	out.y = lhs.y + rhs.y;

	return out;
}

#endif