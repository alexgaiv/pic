#include "pic_kernel.h"

void PicKernel::Init(cl_Descriptor &cld, const char *filename, cl_Grid &grid, bool verbose)
{
	cl::Program program = CreateProgramFromFile(filename, cld.ctx);

	try
	{
		if (verbose) std::cout << "building a program...";
		program.build("-I cl/include");
		if (verbose) std::cout << "done\n";
		std::string log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(cld.device);
		if (log != "") OutputDebugStringA(log.c_str());
	}
	catch (const cl::Error &e)
	{
		std::string log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(cld.device);
		OutputDebugStringA(log.c_str());
		throw e;
	}

	Vector3i numCells = grid.GetNumInnerCells();
	Vector3i groupSize = grid.GetGroupSize();

	this->globalRange = cl::NDRange(numCells.x, numCells.y, numCells.z);
	this->localRange = cl::NDRange(groupSize.x, groupSize.y, groupSize.z);
	this->kernel = cl::Kernel(program, "main");

	AddArg(v2v(grid.GetMin()));
	AddArg(v2v(grid.GetMax()));
          
	AddArg(grid.Ex.cl_buffer);
	AddArg(grid.Ey.cl_buffer);
	AddArg(grid.Ez.cl_buffer);
	AddArg(grid.Bx.cl_buffer);
	AddArg(grid.By.cl_buffer);
	AddArg(grid.Bz.cl_buffer);
	AddArg(grid.Jx.cl_buffer);
	AddArg(grid.Jy.cl_buffer);
	AddArg(grid.Jz.cl_buffer);

	Vector3i s = groupSize + Vector3i(2);
	AddLocalBufferArg<real_t>(s.x * (s.y + 1) * (s.z + 1)); // Ex local
	AddLocalBufferArg<real_t>((s.x + 1) * s.y * (s.z + 1)); // Ey local
	AddLocalBufferArg<real_t>((s.x + 1) * (s.y + 1) * s.z); // Ez local

	AddLocalBufferArg<real_t>((s.x + 1) * s.y * s.z);       // Bx local
	AddLocalBufferArg<real_t>(s.x * (s.y + 1) * s.z);       // By local
	AddLocalBufferArg<real_t>(s.x * s.y * (s.z + 1));       // Bz local

	AddLocalBufferArg<real_t>(s.x * (s.y + 1) * (s.z + 1)); // Jx local
	AddLocalBufferArg<real_t>((s.x + 1) * s.y * (s.z + 1)); // Jy local
	AddLocalBufferArg<real_t>((s.x + 1) * (s.y + 1) * s.z); // Jz local

	lastArgIdx = 20;
}