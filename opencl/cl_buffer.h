#ifndef _CL_BUFFER_H_
#define _CL_BUFFER_H_

#include "common.h"
#include "cl_descriptor.h"

template<class DataType>
class cl_Buffer
{
public:
    std::vector<DataType> data;
    cl::Buffer buffer;

    cl_Buffer() { }

    cl_Buffer(cl_Descriptor &cld, cl_mem_flags flags, ::size_t size) :
        queue(cld.queue),
        buffer(cld.ctx, flags, size * sizeof(DataType)),
        data(size)
    { }

    void Read() {
        queue.enqueueReadBuffer(buffer, true, 0, data.size() * sizeof(DataType), data.data());
    }

    void Write() {
        queue.enqueueWriteBuffer(buffer, true, 0, data.size() * sizeof(DataType), data.data());
    }
private:
    cl::CommandQueue queue;
};

#endif // _CL_BUFFER_H_