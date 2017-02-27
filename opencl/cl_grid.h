#ifndef _CL_GRID_H_
#define _CL_GRID_H_

#include <CL\cl.hpp>
#include "mathtypes.h"
#include "real_t.h"

class cl_Lattice
{
public:
	cl::Buffer buffer;

	cl_Lattice(const Vector3i &size, cl::CommandQueue &queue) :
		data(size.x * size.y * size.z),
		size(size)
	{
		cl::Context ctx = queue.getInfo<CL_QUEUE_CONTEXT>();
		std::size_t memSize = data.size() * sizeof(real_t);
		buffer = cl::Buffer(ctx, CL_MEM_READ_WRITE, memSize);
		queue.enqueueWriteBuffer(buffer, true, 0, memSize, data.data());
	}

	Vector3i GetSize() const { return size; }
	std::vector<real_t> &GetData() { return data; }
	const std::vector<real_t> &GetData() const { return data; }

	const real_t &operator()(int i, int j, int k) const {
		return data[(k * size.y + j) * size.x + i];
	}
	real_t &operator()(int i, int j, int k) {
		return data[(k * size.y + j) * size.x + i];
	}

	void ReadBuffer(cl::CommandQueue &queue)
	{
		queue.enqueueReadBuffer(buffer, true, 0, data.size() * sizeof(real_t), data.data());
	}
private:
	std::vector<real_t> data;
	Vector3i size;
};

class cl_Grid
{
private:
	Vector3i numCells;
public:
	cl_Grid(const Vector3f &vmin, const Vector3f &vmax, const Vector3i &numInnerCells, cl::CommandQueue &q) :
		vmin(vmin),
		vmax(vmax),
		numCells(numInnerCells + Vector3i(2)),
		cellSize((vmax - vmin) / numInnerCells),

		Ex(numCells + Vector3i(0, 1, 1), q),
		Ey(numCells + Vector3i(1, 0, 1), q),
		Ez(numCells + Vector3i(1, 1, 0), q),

		Bx(numCells + Vector3i(1, 0, 0), q),
		By(numCells + Vector3i(0, 1, 0), q),
		Bz(numCells + Vector3i(0, 0, 1), q),

		Jx(numCells + Vector3i(0, 1, 1), q),
		Jy(numCells + Vector3i(1, 0, 1), q),
		Jz(numCells + Vector3i(1, 1, 0), q)
	{ }

	Vector3f GetMin() const { return vmin; }
	Vector3f GetMax() const { return vmax; }
	Vector3i GetNumCells() const { return numCells; }
	Vector3i GetNumInnerCells() const { return numCells - Vector3i(2); }
	Vector3f GetCellSize() const { return cellSize; }

	cl_Lattice Ex, Ey, Ez;
	cl_Lattice Bx, By, Bz;
	cl_Lattice Jx, Jy, Jz;

private:
	Vector3f vmin, vmax;
	Vector3f cellSize;
};

#endif // _CL_GRID_H_