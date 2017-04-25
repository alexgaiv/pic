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

	kernel.setArg(0, v2v(grid.GetMin()));
	kernel.setArg(1, v2v(grid.GetMax()));
	kernel.setArg(2, v2v(numCells));

	kernel.setArg(3, grid.Ex.cl_buffer.buffer);
	kernel.setArg(4, grid.Ey.cl_buffer.buffer);
	kernel.setArg(5, grid.Ez.cl_buffer.buffer);
	kernel.setArg(6, grid.Bx.cl_buffer.buffer);
	kernel.setArg(7, grid.By.cl_buffer.buffer);
	kernel.setArg(8, grid.Bz.cl_buffer.buffer);
	kernel.setArg(9, grid.Jx.cl_buffer.buffer);
	kernel.setArg(10, grid.Jy.cl_buffer.buffer);
	kernel.setArg(11, grid.Jz.cl_buffer.buffer);

	Vector3i s = groupSize + Vector3i(2);
	kernel.setArg(12, s.x * (s.y + 1) * (s.z + 1) * sizeof(real_t), NULL);
	kernel.setArg(13, (s.x + 1) * s.y * (s.z + 1) * sizeof(real_t), NULL);
	kernel.setArg(14, (s.x + 1) * (s.y + 1) * s.z * sizeof(real_t), NULL);

	kernel.setArg(15, (s.x + 1) * s.y * s.z * sizeof(real_t), NULL);
	kernel.setArg(16, s.x * (s.y + 1) * s.z * sizeof(real_t), NULL);
	kernel.setArg(17, s.x * s.y * (s.z + 1) * sizeof(real_t), NULL);

	kernel.setArg(18, s.x * (s.y + 1) * (s.z + 1) * sizeof(real_t), NULL);
	kernel.setArg(19, (s.x + 1) * s.y * (s.z + 1) * sizeof(real_t), NULL);
	kernel.setArg(20, (s.x + 1) * (s.y + 1) * s.z * sizeof(real_t), NULL);

	lastArgIdx = 20;
}