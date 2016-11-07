#ifndef _MATH_TYPES_H_
#define _MATH_TYPES_H_

#include <math.h>

template<class T>
class Vector3
{
public:
	union {
		struct { T x, y, z; };
		T data[3];
	};

	explicit Vector3(T init = T(0))
		{ x = y = z = init; }

	Vector3(T x, T y, T z)
		: x(x), y(y), z(z) { }

	template<class T2>
	Vector3(const Vector3<T2> &v)
	{
		x = (T)v.x;
		y = (T)v.y;
		z = (T)v.z;
	}

	T Length() const;
	T Square() const;
	void Normalize();

	T &operator[](int i) { return data[i]; }
	T operator[](int i) const { return data[i]; }

	bool operator==(const Vector3 &v) const;
	bool operator!=(const Vector3 &v) const;

	Vector3 operator+(const Vector3 &v) const;
	Vector3 operator-(const Vector3 &v) const;
	Vector3 operator*(const Vector3 &v) const;
	Vector3 operator/(const Vector3 &v) const;
	Vector3 operator-() const;
	Vector3 operator*(T scale) const;
	Vector3 operator/(T scale) const;

	Vector3 &operator+=(const Vector3 &v);
	Vector3 &operator-=(const Vector3 &v);
	Vector3 &operator*=(T scale);
	Vector3 &operator/=(T scale);
};

typedef Vector3<double> Vector3d;
typedef Vector3<float> Vector3f;
typedef Vector3<int> Vector3i;

#define OVERLOAD_V1(func) \
template<class T> \
Vector3<T> func(const Vector3<T> &v) { \
	return Vector3<T>(func(v.x), func(v.y), func(v.z)); \
}

#define OVERLOAD_V2(func) \
template<class T> \
Vector3<T> func(const Vector3<T> &v1, const Vector3<T> &v2) { \
	return Vector3<T>(func(v1.x, v2.x), func(v1.y, v2.y), func(v1.z, v2.z)); \
}

OVERLOAD_V1(abs)
OVERLOAD_V1(floor)
OVERLOAD_V1(ceil)
OVERLOAD_V2(min)
OVERLOAD_V2(max)

#undef OVERLOAD_V1
#undef OVERLOAD_V2

template<class T>
T Dot(const Vector3<T> &v1, const Vector3<T> &v2) {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

template<class T>
Vector3<T> Cross(const Vector3<T> &v1, const Vector3<T> &v2) {
	return Vector3<T>(
		v1.y*v2.z - v1.z*v2.y,
		v1.z*v2.x - v1.x*v2.z,
		v1.x*v2.y - v1.y*v2.x);
}

template<class T>
T Vector3<T>::Length() const {
	return sqrt(x*x + y*y + z*z);
}

template<class T>
T Vector3<T>::Square() const {
	return x*x + y*y + z*z;
}

template<class T>
void Vector3<T>::Normalize() {
	T f = T(1) / Length();
	x *= f; y *= f; z *= f;
}

template<class T>
bool Vector3<T>::operator==(const Vector3<T> &v) const {
	return x == v.x && y == v.y && z == v.z;
}

template<class T>
bool Vector3<T>::operator!=(const Vector3<T> &v) const {
	return !operator==(v);
}

template<class T>
Vector3<T> Vector3<T>::operator+(const Vector3<T> &v) const {
	return Vector3<T>(x + v.x, y + v.y, z + v.z);
}

template<class T>
Vector3<T> Vector3<T>::operator-(const Vector3<T> &v) const {
	return Vector3<T>(x - v.x, y - v.y, z - v.z);
}

template<class T>
Vector3<T> Vector3<T>::operator*(const Vector3<T> &v) const {
	return Vector3<T>(x * v.x, y * v.y, z * v.z);
}

template<class T>
Vector3<T> Vector3<T>::operator/(const Vector3<T> &v) const {
	return Vector3<T>(x / v.x, y / v.y, z / v.z);
}

template<class T>
Vector3<T> Vector3<T>::operator-() const {
	return Vector3<T>(-x, -y, -z);
}

template<class T>
Vector3<T> Vector3<T>::operator*(T scale) const {
	return Vector3<T>(x * scale, y * scale, z * scale);
}

template<class T>
Vector3<T> operator*(T scale, const Vector3<T> &v) {
	return v*scale;
}

template<class T>
Vector3<T> Vector3<T>::operator/(T scale) const {
	T f = T(1) / scale;
	return Vector3<T>(x * f, y * f, z * f);
}

template<class T>
Vector3<T> &Vector3<T>::operator+=(const Vector3<T> &v) {
	x += v.x; y += v.y; z += v.z;
	return *this;
}

template<class T>
Vector3<T> &Vector3<T>::operator-=(const Vector3<T> &v) {
	x -= v.x; y -= v.y; z -= v.z;
	return *this;
}

template<class T>
Vector3<T> &Vector3<T>::operator*=(T scale) {
	x *= scale; y *= scale; z *= scale;
	return *this;
}

template<class T>
Vector3<T> &Vector3<T>::operator/=(T scale) {
	x /= scale; y /= scale; z /= scale;
	return *this;
}

#endif // _MATH_TYPES_H_