#ifndef _ENSEMBLE_H_
#define _ENSEMBLE_H_

#include <vector>
#include "mathtypes.h"
#include "particle.h"

template<class T>
class Iterator
{
public:
	Iterator(int reserveSize = 0) : cur(0)
	{
		elems.reserve(reserveSize);
	}

	T &Current() { return elems[cur]; }
	const T &Current() const { return elems[cur]; }

	int GetSize() const { return elems.size(); }

	void Add(const T &elem) { elems.push_back(elem); }
	void Delete() {
		elems[cur--] = elems.back();
		elems.pop_back();
	}
	void Clear() {
		cur = 0;
		elems.clear();
	}

	void Begin() const { cur = 0; }
	void Next() const { cur++; }
	bool End() const { return cur >= elems.size(); };
private:
	mutable unsigned cur;
	std::vector<T> elems;
};

typedef Iterator<Particle> ParticleSystem;

class Ensemble
{
public:
	Ensemble() { }
	Ensemble(const Vector3d &vmin, const Vector3d &cellSize, const Vector3i &size) :
		vmin(vmin),
		vmax(vmin + cellSize * size),
		size(size),
		systems(size.x * size.y * size.z)
	{ }

	Vector3i GetSize() const { return size; }
	Vector3d GetMin() const { return vmin; }
	Vector3d GetMax() const { return vmax; }

	bool IsPointInside(const Vector3d &p)
	{
		return
			p.x >= vmin.x && p.x <= vmax.x &&
			p.y >= vmin.y && p.y <= vmax.y &&
			p.z >= vmin.z && p.z <= vmax.z;
	}

	ParticleSystem &operator()(int i, int j, int k) {
		return systems[(k * size.y + j) * size.x + i];
	}
	const ParticleSystem &operator()(int i, int j, int k) const {
		return systems[(k * size.y + j) * size.x + i];
	}
private:
	Vector3d vmin, vmax;
	Vector3i size;
	std::vector<ParticleSystem> systems;
};

#endif // _ENSEMBLE_H_