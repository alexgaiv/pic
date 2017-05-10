#ifndef _CL_GRID_H_
#define _CL_GRID_H_

#include "common.h"
#include "real_t.h"
#include "mathtypes.h"
#include "cl_descriptor.h"
#include "cl_buffer.h"

class cl_Lattice
{
public:
    cl_Buffer<real_t> cl_buffer;

    cl_Lattice(cl_Descriptor &cld, const Vector3i &size) :
        cl_buffer(cl_Buffer<real_t>(cld, CL_MEM_READ_WRITE, size.x * size.y * size.z)),
        size(size),
        data(cl_buffer.data)
    { }

    Vector3i GetSize() const { return size; }
    void ReadBuffer() { cl_buffer.Read(); }

    const real_t &operator()(int i, int j, int k) const {
        return data[(k * size.y + j) * size.x + i];
    }
    real_t &operator()(int i, int j, int k) {
        return data[(k * size.y + j) * size.x + i];
    }
private:
    Vector3i size;
    std::vector<real_t> &data;
};

class cl_Grid
{
private:
    Vector3i numCells;
public:
    cl_Grid(cl_Descriptor &cld, const Vector3f &vmin, const Vector3f &vmax,
        const Vector3i &numInnerCells, const Vector3i &numGroups) :
        vmin(vmin),
        vmax(vmax),
        numCells(numInnerCells + Vector3i(2)),
        cellSize((vmax - vmin) / numInnerCells),
        numGroups(numGroups),
        groupSize(numInnerCells / numGroups),

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
    Vector3i GetNumGroups() const { return numGroups; }
    Vector3i GetGroupSize() const { return groupSize; }
    Vector3f GetCellSize() const { return cellSize; }

    cl_Lattice Ex, Ey, Ez;
    cl_Lattice Bx, By, Bz;
    cl_Lattice Jx, Jy, Jz;

private:
    Vector3f vmin, vmax;
    Vector3f cellSize;
    Vector3i numGroups;
    Vector3i groupSize;
};

#endif // _CL_GRID_H_