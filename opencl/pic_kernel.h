#ifndef _PIC_KERNEL_H_
#define _PIC_KERNEL_H_

#include <iostream>
#include <vector>

#include "common.h"
#include "real_t.h"
#include "utils.h"
#include "cl_buffer.h"
#include "cl_grid.h"

class PicKernel
{
public:
    PicKernel() : lastArgIdx(0) { }
    PicKernel(cl_Descriptor &cld, const char *filename,
        const char *include_dir, cl_Grid &grid, bool verbose = false);

    cl::Kernel GetKernel() const { return kernel; }

    template<class T>
    void AddArg(T value) {
        kernel.setArg(lastArgIdx++, value);
    }
    
    void AddArg(::size_t size, void *argPtr) {
        kernel.setArg(lastArgIdx++, size, argPtr);
    }

    template<class DataType>
    void AddArg(const cl_Buffer<DataType> &cl_buffer) {
        kernel.setArg(lastArgIdx++, cl_buffer.buffer);
    }

    template<class T>
    void AddLocalBufferArg(::size_t size) {
        kernel.setArg(lastArgIdx++, size * sizeof(T), NULL);
    }

    template<class T>
    void AddGlobalBufferArg(::size_t size, cl_mem_flags flags) {
        cl::Buffer buffer(cld.ctx, flags, size * sizeof(T));
        globalBuffers.push_back(buffer);
        AddArg(buffer);
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
    std::vector<cl::Buffer> globalBuffers;
};

#endif // _PIC_KERNEL_H_