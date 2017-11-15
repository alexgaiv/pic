#include "tests.h"
#include "mgl_suppress_warnings.h"
#include <mgl2\mgl.h>
#include <iostream>

#include "pic_kernel.h"
#include "cl_buffer.h"

using std::cout;
using std::endl;

static void BuildErrSubPlot(mglGraph *gr, const char *title, const char *legend, const char *color,
    float *values, int num_values, float x_min, float x_max, float k_min, float k_max)
{
    mglData y;

    gr->Title(title);
    gr->ClearLegend();
    gr->AddLegend(legend, color);
    gr->SetRanges(x_min, x_max, values[num_values - 1] * k_min, values[0] * k_max);
    gr->Label('x', "steps");
    gr->Label('y', "error");
    gr->Axis();
    gr->Legend();
    y.Set(values, num_values);
    gr->Plot(y, color);
}

static void BuildErrPlot(const char *filename, const char *title, float *rel_values, float *abs_values,
    int num_values, float x_min, float x_max, float k_min, float k_max)
{
    mglGraph gr;

    gr.SetSize(1500, 600);
    gr.SubPlot(2, 1, 0);
    BuildErrSubPlot(&gr, title, "relative error", "b", rel_values, num_values, x_min, x_max, k_min, k_max);
    gr.SubPlot(2, 1, 1);
    BuildErrSubPlot(&gr, title, "absolute error", "g", abs_values, num_values, x_min, x_max, k_min, k_max);
    gr.WriteBMP(filename);
    gr.ClearFrame();
}

bool TestGrid(cl_Descriptor &cld)
{
    cl_Grid grid(cld, Vector3f(0), Vector3f(1, 0.125, 0.125), Vector3i(64, 8, 8), Vector3i(4, 4, 4));
    const Vector3i numCells = grid.GetNumInnerCells();
    const Vector3i groupSize = grid.GetGroupSize();

    int result_buffer_size = numCells.x * numCells.y * numCells.z;
    cl_Buffer<int> result_buffer(cld,  CL_MEM_WRITE_ONLY, result_buffer_size);

    PicKernel kernel(cld, "cl/test_grid.cl", "cl/include", grid, true);
    kernel.AddArg(result_buffer);

    int dt = timeGetTime();
    kernel.Run();

    grid.Ex.ReadBuffer();
    grid.Ey.ReadBuffer();
    grid.Ez.ReadBuffer();
    grid.Bx.ReadBuffer();
    grid.By.ReadBuffer();
    grid.Bz.ReadBuffer();
    result_buffer.Read();

    dt = timeGetTime() - dt;
    cout << "time: " << dt << "ms" << endl;

    bool passed = true;
    FOR3(i, j, k, numCells)
    {
        int idx = (k * numCells.y + j) * numCells.x + i;
        if (result_buffer.data[idx] != 1) {
            passed = false;
            break;
        }
    }
    
    return passed;
}

void TestBoris(cl_Descriptor &cld)
{
    cl_Grid grid(cld, Vector3f(0), Vector3f(1), Vector3i(16, 8, 8), Vector3i(4, 4, 4));
    PicKernel kernel(cld, "cl/test_boris.cl", "cl/include", grid, true);

    const int startN = 100;
    const int dN = 50;
    const int numIters = 50;

    cl_Buffer<float> buffers[8];

    kernel.AddArg(startN);
    kernel.AddArg(dN);
    kernel.AddArg(numIters);

    for (int i = 0; i < 8; i++) {
        buffers[i] = cl_Buffer<float>(cld, CL_MEM_WRITE_ONLY, numIters);
        kernel.AddArg(buffers[i]);
    }

    kernel.Run();

    for (int i = 0; i < 8; i++)
        buffers[i].Read();

    int maxN = dN * numIters;
    BuildErrPlot("plot/boris/test1-1.bmp", "Test 1 (Coords)",
        buffers[0].RawData(), buffers[1].RawData(), numIters, (float)startN, (float)maxN, 0.1f, 1.0f);
    BuildErrPlot("plot/boris/test1-2.bmp", "Test 1 (Momentum)",
        buffers[2].RawData(), buffers[3].RawData(), numIters, (float)startN, (float)maxN, 0.99f, 1.0f);

    BuildErrPlot("plot/boris/test2-1.bmp", "Test 2 (Coords)",
        buffers[4].RawData(), buffers[5].RawData(), numIters, (float)startN, (float)maxN, 1.0f, 1.0f);
    BuildErrPlot("plot/boris/test2-2.bmp", "Test 2 (Momentum)",
        buffers[6].RawData(), buffers[7].RawData(), numIters, (float)startN, (float)maxN, 1.0f, 1.0f);
}