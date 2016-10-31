kernel void main(double alpha, global float *A, global float *B, global float *C)
{
	int i = get_global_id(0);
	C[i] = alpha * A[i] + B[i];
}