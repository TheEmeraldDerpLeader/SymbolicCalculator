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
	cl::Buffer inCoords;
	cl::Buffer outCoords;

	size_t workGroupSize = 0;

	FlatKernelExecutor(cl::Device& deviceRef);

	void Execute(FlatSymExp& flat, FlatSymExp& grad, int coordsCount, float* in, float* out, float threshold = 0.00001f, int maxIter = 80);
};


#undef float