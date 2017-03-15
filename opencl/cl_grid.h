#ifndef _CL_GRID_H_
#define _CL_GRID_H_

#include "common.h"
#include "real_t.h"
#include "mathtypes.h"
#include "cl_descriptor.h"

class cl_Lattice
{
public:
	cl::Buffer buffer;

	cl_Lattice(cl_Descriptor &cld, const Vector3i &size) :
		q(cld.queue),
		data(size.x * size.y * size.z),
		size(size)
	{
		cl::Context ctx = q.getInfo<CL_QUEUE_CONTEXT>();
		std::size_t memSize = data.size() * sizeof(real_t);
		buffer = cl::Buffer(ctx, CL_MEM_READ_WRITE, memSize);
		q.enqueueWriteBuffer(buffer, true, 0, memSize, data.data());
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

	void ReadBuffer()
	{
		q.enqueueReadBuffer(buffer, true, 0, data.size() * sizeof(real_t), data.data());
	}
private:
	std::vector<real_t> data;
	Vector3i size;
	cl::CommandQueue q;
};

class cl_Grid
{
private:
	Vector3i numCells;
public:
	cl_Grid(cl_Descriptor &cld, const Vector3f &vmin, const Vector3f &vmax,
		    const Vector3i &numInnerCells, const Vector3i &groupNum) :
		vmin(vmin),
		vmax(vmax),
		numCells(numInnerCells + Vector3i(2)),
		cellSize((vmax - vmin) / numInnerCells),
		groupNum(groupNum),
		groupSize(numInnerCells / groupNum),

		Ex(cld, numCells + Vector3i(0, 1, 1)),
		Ey(cld, numCells + Vector3i(1, 0, 1)),
		Ez(cld, numCells + Vector3i(1, 1, 0)),

		Bx(cld, numCells + Vector3i(1, 0, 0)),
		By(cld, numCells + Vector3i(0, 1, 0)),
		Bz(cld, numCells + Vector3i(0, 0, 1)),

		Jx(cld, numCells + Vector3i(0, 1, 1)),
		Jy(cld, numCells + Vector3i(1, 0, 1)),
		Jz(cld, numCells + Vector3i(1, 1, 0))
	{ }

	Vector3f GetMin() const { return vmin; }
	Vector3f GetMax() const { return vmax; }
	Vector3i GetNumCells() const { return numCells; }
	Vector3i GetNumInnerCells() const { return numCells - Vector3i(2); }
	Vector3f GetCellSize() const { return cellSize; }
	Vector3i GetGroupNum() const { return groupNum; }
	Vector3i GetGroupSize() const { return groupSize; }

	cl_Lattice Ex, Ey, Ez;
	cl_Lattice Bx, By, Bz;
	cl_Lattice Jx, Jy, Jz;

private:
	Vector3f vmin, vmax;
	Vector3f cellSize;
	Vector3i groupNum;
	Vector3i groupSize;
};

#endif // _CL_GRID_H_