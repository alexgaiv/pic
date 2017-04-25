#ifndef _PIC_KERNEL_H_
#define _PIC_KERNEL_H_

#include "common.h"
#include <iostream>
#include "real_t.h"
#include "utils.h"
#include "cl_buffer.h"
#include "cl_grid.h"

class PicKernel
{
public:
	PicKernel() : lastArgIdx(0) { }

	PicKernel(cl_Descriptor &cld, const char *filename, cl_Grid &grid, bool verbose = false)
		: cld(cld), lastArgIdx(0)
	{
		Init(cld, filename, grid, verbose);
	}

	cl::Kernel GetKernel() const { return kernel; }

	void Init(cl_Descriptor &cld, const char *filename, cl_Grid &grid, bool verbose);

	template<class T>
	void AddArg(T value) {
		kernel.setArg(++lastArgIdx, value);
	}
	
	void AddArg(::size_t size, void *argPtr) {
		kernel.setArg(++lastArgIdx, size, argPtr);
	}

	template<class DataType>
	void AddArg(const cl_Buffer<DataType> &cl_buffer) {
		kernel.setArg(++lastArgIdx, cl_buffer.buffer);
	}

	template<class T>
	void AddArgLocal(::size_t size) {
		kernel.setArg(++lastArgIdx, size * sizeof(T), NULL);
	}

	void Run() {
		cld.queue.enqueueNDRangeKernel(kernel, cl::NullRange, globalRange, localRange);
		cld.queue.flush();
		cld.queue.finish();
	}
private:
	cl_Descriptor cld;
	cl_uint lastArgIdx;
	cl::Kernel kernel;
	cl::NDRange globalRange, localRange;
};

#endif // _PIC_KERNEL_H_