#include <stdlib.h>
#include <cl\cl.h>

int main()
{
	cl_platform_id platformId;
	cl_device_id deviceId;

	clGetPlatformIDs(1, &platformId, NULL);
	clGetDeviceIDs(platformId, CL_DEVICE_TYPE_GPU, 1, &deviceId, NULL);

	return 0;
}