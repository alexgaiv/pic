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
	float *A = new float[vec_size];
	float *B = new float[vec_size];
	float *C = new float[vec_size];
	float alpha = 0.5f;

	for (int i = 0; i < vec_size; i++)
	{
		A[i] = float(rand() % 100);
		B[i] = float(rand() % 100);
		C[i] = -1.0f;
	}

	try
	{
		std::size_t mem_size = vec_size * sizeof(float);

		cout << "init context...";
		Context ctx(CL_DEVICE_TYPE_GPU);
		cout << "done\n";

		Device device = ctx.getInfo<CL_CONTEXT_DEVICES>()[0];
		CommandQueue queue(ctx, device);
		Program program = CreateProgramFromFile("kernel.cl", ctx);

		auto s = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();

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

		kernel.setArg(0, alpha);
		kernel.setArg(1, A_mem);
		kernel.setArg(2, B_mem);
		kernel.setArg(3, C_mem);

		queue.enqueueWriteBuffer(A_mem, true, 0, mem_size, A);
		queue.enqueueWriteBuffer(B_mem, true, 0, mem_size, B);
		queue.enqueueNDRangeKernel(kernel, NullRange, NDRange(vec_size), NDRange(64));
		queue.enqueueReadBuffer(C_mem, true, 0, mem_size, C);
		
		queue.flush();
		queue.finish();

		bool f = true;
		for (int i = 0; i < vec_size; i++)
		{
			if (C[i] != alpha * A[i] + B[i]) {
				f = false;
				break;
			}
		}
		cout << (f ? "all good" : "all bad");
	}
	catch (const Error &e)
	{
		cout << endl << e.what() << " failed";
	}
	
	delete[] A;
	delete[] B;
	delete[] C;

	getchar();
	return 0;
}