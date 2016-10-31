#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#define __CL_ENABLE_EXCEPTIONS

#include <stdlib.h>
#include <iostream>
#include <cl\cl.hpp>
#include <vector>
#include <fstream>
#include "utils.h"

using namespace cl;
using namespace std;

const std::size_t vec_size = 1024;

int main()
{
	double A[vec_size];
	double B[vec_size];
	double C[vec_size];

	for (int i = 0; i < vec_size; i++)
	{
		A[i] = rand() % 100;
		B[i] = rand() % 100;
		C[i] = 0;
	}

	try
	{
		std::size_t mem_size = vec_size * sizeof(double);

		cout << "init context...";
		Context ctx(CL_DEVICE_TYPE_GPU);
		cout << "done\n";

		Device device = ctx.getInfo<CL_CONTEXT_DEVICES>()[0];
		CommandQueue commandQueue(ctx, device);
		Program program = CreateProgramFromFile("kernel.cl", ctx);

		try
		{
			cout << "building a program...";
			program.build();
			cout << "done\n";
		}
		catch (Error e)
		{
			string log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
			OutputDebugStringA(log.c_str());
			throw e;
		}

		Buffer A_mem(ctx, CL_MEM_READ_ONLY, mem_size);
		Buffer B_mem(ctx, CL_MEM_READ_ONLY, mem_size);
		Buffer C_mem(ctx, CL_MEM_WRITE_ONLY, mem_size);
		Kernel kernel(program, "main");

		kernel.setArg(0, 0.5);
		kernel.setArg(1, A_mem);
		kernel.setArg(2, B_mem);
		kernel.setArg(3, B_mem);

		commandQueue.enqueueWriteBuffer(A_mem, true, 0, mem_size, A);
		commandQueue.enqueueWriteBuffer(B_mem, true, 0, mem_size, B);
		//commandQueue.enqueueNDRangeKernel(kernel, )

		
	}
	catch (Error e)
	{
		cout << endl << e.what() << " failed";
	}
	

	getchar();
	return 0;
}