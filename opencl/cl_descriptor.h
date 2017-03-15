#ifndef _CL_DESCRIPTOR_H_
#define _CL_DESCRIPTOR_H_

#include "common.h"

class cl_Descriptor
{
public:
	cl::Context ctx;
	cl::Device device;
	cl::CommandQueue queue;

	cl_Descriptor() { }

	cl_Descriptor(const cl::Context &ctx) : ctx(ctx)
	{
		device = ctx.getInfo<CL_CONTEXT_DEVICES>()[0];
		queue = cl::CommandQueue(ctx, device);
	}
};

#endif // _CL_DESCRIPTOR_H_