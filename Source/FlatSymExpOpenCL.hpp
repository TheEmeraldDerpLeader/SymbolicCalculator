#pragma once

#include "FlatSymExp.hpp"
#include "Helpers.hpp"

#include <iostream>
#include <ostream>
#include <CL/opencl.hpp>

#define float float


//error codes: https://gist.github.com/bmount/4a7144ce801e5569a0b6

//OpenCL docs
//OpenCL C language specification https://registry.khronos.org/OpenCL/specs/3.0-unified/html/OpenCL_C.html
//Unified model specification https://registry.khronos.org/OpenCL/specs/3.0-unified/html/OpenCL_API.html
//SDK simple guide https://www.khronos.org/files/opencl30-reference-guide.pdf
//SDK specification https://registry.khronos.org/OpenCL/sdk/3.0/docs/man/html/
//C++ Bindings https://github.khronos.org/OpenCL-CLHPP/index.html

void QueryDevices(std::ostream& out = std::cout);

void QueryDevices(std::vector<cl::Platform>& platformsOut, std::vector<std::vector<cl::Device>>& devicesOut, std::ostream& out = std::cout);

class FlatKernelExecutor
{
public:
	cl::Device device;
	cl::Context context;
	cl::CommandQueue comQue;
	cl::Program prog;
	cl::Kernel kern;

	cl::Buffer flatSymData;
	cl::Buffer gradientData;
	cl::Buffer testCoordsBuf; //coords to perform nm on
	cl::Buffer rootsBuf; //roots to check final coords on
	cl::Buffer finalCoordsBuf; //the final coords of each point after applying nm
	cl::Buffer itersBuf; //the number of iterations that nm took before satisfying convergence condition (-1 if failed)
	cl::Buffer idsBuf; //the root index a final coord is close to (-1 if not in root set, -2 if failed to converge)
	cl::Buffer indexesBuf; //index of finalCoords which have new roots (id -1)
	cl::Buffer initialBuf; //initial coord for NMnTomRect
	cl::Buffer base1Buf; //basis vector for x coordinate change for NMnTomRect
	cl::Buffer base2Buf; //basis vector for y coordinate change for NMnTomRect
	cl::Buffer texBuf; //RGBA texture

	FlatKernelExecutor(cl::Device& deviceRef);

	void InitializeBuffers(int flatSize, int gradSize, int width, int height, int idC, int rootC, int checkCount);
	void InitializeRootBuffer(int idC, int rootC);

	void RunNMnTom(FlatSymExp& flat, FlatSymExp& grad, int coordsCount, float* testCoords, float* finalCoords, int* iters, float threshold = 0.00001f, int maxIter = 80);
	void RunNMnTomRect(FlatSymExp& flat, FlatSymExp& grad, int coordsX, int cooordsY, float* initial, float* base1, float* base2, float* finalCoords, int* iters, float threshold = 0.00001f, int maxIter = 80);
	void RunAssociateCoords(int coordsCount, int idC, int rootC, float* finalCoords, float* roots, int* iters, int* ids);
	void RunFindNewRoots(int checkCount, int* ids, int* indexes, int checkSpan, int checkOffset);
	void RunColorTex(int coordsCount, int* iters, unsigned char* tex);
	void RunColorTex(int coordsCount, int* iters, unsigned char* tex, int* ids);

};

void FindAndCollectRootsCPU(int idC, int coordsCount, std::vector<float>& roots, float* outCoords, int* iters, int* ids);
void CollectNewRoots(int idC, int checkCount, std::vector<float>& roots, float* coordsOut, int* indexes);

#undef float