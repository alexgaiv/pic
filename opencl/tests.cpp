#include "tests.h"
#include "pic_kernel.h"
#include <iostream>

using namespace std;
using namespace cl;

void TestGrid(cl_Descriptor &cld)
{
	const Vector3f vmin(0);
	const Vector3f vmax(1, 0.125, 0.125);
	const Vector3i numCells(64, 8, 8);
	const Vector3i groupNum(4, 4, 4);

	cl_Grid grid(cld, vmin, vmax, numCells, groupNum);
	const Vector3i groupSize = grid.GetGroupSize();

	int testSize = numCells.x * numCells.y * numCells.z;
	int *testData = new int[testSize];
	cl::Buffer testMem(cld.ctx, CL_MEM_WRITE_ONLY, testSize * sizeof(int));

	PicKernel kernel(cld, "cl/test_grid.cl", grid, true);
	kernel.AddArg(testMem);

	int dt = timeGetTime();
	kernel.Run();

	grid.Ex.ReadBuffer();
	grid.Ey.ReadBuffer();
	grid.Ez.ReadBuffer();
	grid.Bx.ReadBuffer();
	grid.By.ReadBuffer();
	grid.Bz.ReadBuffer();
	cld.queue.enqueueReadBuffer(testMem, true, 0, testSize * sizeof(int), testData);

	dt = timeGetTime() - dt;
	cout << "time: " << dt << "ms" << endl;

	bool passed = true;
	FOR3(i, j, k, numCells)
	{
		int idx = (k * numCells.y + j) * numCells.x + i;
		if (testData[idx] != 1) {
			passed = false;
			break;
		}
	}
	delete[] testData;
	cout << "grid test: " << (passed ? "passed" : "failed");
}