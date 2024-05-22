#include "FlatSymExpOpenCL.hpp"

#define float float

static struct ErrorHandle
{
public:
	cl_int errorVal = 0;
	ErrorHandle() {}

	ErrorHandle& operator=(cl_int error) 
	{
		errorVal = error;
		Eval();

		return *this;
	}
	void Eval()
	{
		if (errorVal != 0)
			std::cout << "OpenCL Error Code: " << errorVal << '\n';
	}
	cl_int* errorValP() { return &errorVal; }
};

static ErrorHandle error;

void QueryDevices(std::ostream& out)
{
	cl::Device device;
	bool gotDevice = false;

	ErrorHandle error;
	std::vector<cl::Platform> platforms;
	error = cl::Platform::get(&platforms);

	std::string info;
	cl_uint infoI;
	cl_ulong infoL;
	for (int i = 0; i < platforms.size(); i++)
	{
		cl::Platform& platform = platforms[i];

		out << "#" << i << ":\n";
		error = platform.getInfo(CL_PLATFORM_NAME, &info);
		out << "Platform name: " << info << '\n';
		error = platform.getInfo(CL_PLATFORM_VERSION, &info);
		out << "Version: " << info << '\n';

		std::vector<cl::Device> devices;
		error = platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
		for (int j = 0; j < devices.size(); j++)
		{
			out << "##" << j << ":\n";
			if (i == 0 && j == 0)
			{
				device = devices[0];
				gotDevice = true;
			}
			error = devices[j].getInfo(CL_DEVICE_VENDOR, &info);
			out << "Vendor: " << info << '\n';
			error = devices[j].getInfo(CL_DEVICE_NAME, &info);
			out << "Device name: " << info << '\n';
			error = devices[j].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &infoL);
			out << "GLobal memory size: " << infoL/(1024*1024) << " MB" << '\n';
			error = devices[j].getInfo(CL_DEVICE_MAX_CLOCK_FREQUENCY, &infoI);
			out << "Max clock speed: " << infoI << " MHz" << '\n';
		}
		out << '\n';
	}
}

void QueryDevices(std::vector<cl::Platform>& platformsOut, std::vector<std::vector<cl::Device>>& devicesOut, std::ostream& out)
{
	cl::Device device;
	bool gotDevice = false;

	ErrorHandle error;
	error = cl::Platform::get(&platformsOut);
	devicesOut.resize(platformsOut.size());

	std::string info;
	cl_uint infoI;
	cl_ulong infoL;
	for (int i = 0; i < platformsOut.size(); i++)
	{
		cl::Platform& platform = platformsOut[i];

		out << "#" << i << ":\n";
		error = platform.getInfo(CL_PLATFORM_NAME, &info);
		out << "Platform name: " << info << '\n';
		error = platform.getInfo(CL_PLATFORM_VERSION, &info);
		out << "Version: " << info << '\n';

		std::vector<cl::Device>& devices = devicesOut[i];
		error = platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
		for (int j = 0; j < devices.size(); j++)
		{
			out << "##" << j << ":\n";
			if (i == 0 && j == 0)
			{
				device = devices[0];
				gotDevice = true;
			}
			error = devices[j].getInfo(CL_DEVICE_VENDOR, &info);
			out << "Vendor: " << info << '\n';
			error = devices[j].getInfo(CL_DEVICE_NAME, &info);
			out << "Device name: " << info << '\n';
			error = devices[j].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &infoL);
			out << "Global memory size: " << infoL/(1024*1024) << " MB" << '\n';
			//error = devices[j].getInfo(CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, &infoL);
			//out << "Max constant buffer size: " << infoL/(1024) << " KB" << '\n';
			error = devices[j].getInfo(CL_DEVICE_MAX_CLOCK_FREQUENCY, &infoI);
			out << "Max clock speed: " << infoI << " MHz" << '\n';
		}
		out << '\n';
	}
}

FlatKernelExecutor::FlatKernelExecutor(cl::Device& deviceRef)
{
	device = deviceRef;
	context = cl::Context(device);

	comQue = cl::CommandQueue(context, device, NULL, error.errorValP());
	error.Eval();

	std::vector<std::string> sources;
	sources.push_back(std::string());


	std::string filePath = "Assets/Kernels/NMnTom.cl";

	StreamBuffer sb(filePath);
	while (sb.CanPop())
		sources[0].push_back(sb.Pop());

	prog = cl::Program(context, sources, error.errorValP());
	error.Eval();
	error = prog.build("");

	if (error.errorVal != 0)
	{
		std::string buildLog;
		error = prog.getBuildInfo(device, CL_PROGRAM_BUILD_LOG, &buildLog);
		std::cout << buildLog << '\n';
	}

	kern = cl::Kernel(prog, "NMnTom", error.errorValP());
	error.Eval();

}

void FlatKernelExecutor::Execute(FlatSymExp& flat, FlatSymExp& grad, int coordsCount, float* in, float* out, float threshold, int maxIter)
{
	//flatSymData = cl::Buffer(context, CL_MEM_READ_ONLY, , nullptr, error.errorValP()); error.Eval();
	//gradientData = cl::Buffer(context, CL_MEM_READ_ONLY, , nullptr, error.errorValP()); error.Eval();
	//inCoords = cl::Buffer(context, CL_MEM_READ_WRITE, , nullptr, error.errorValP()); error.Eval();
	//outCoords = cl::Buffer(context, CL_MEM_READ_WRITE, , nullptr, error.errorValP()); error.Eval();

	int expC = flat.expCount();
	int idC = flat.idCount();

	error = kern.getWorkGroupInfo(device, CL_KERNEL_WORK_GROUP_SIZE, &workGroupSize);

	if (idC == 0 || expC == 0 || workGroupSize == 0 || coordsCount == 0)
		return;

	flatSymData = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*flat.size, nullptr, error.errorValP()); error.Eval();
	gradientData = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*grad.size, nullptr, error.errorValP()); error.Eval();
	inCoords = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float)*coordsCount*idC, nullptr, error.errorValP()); error.Eval();
	outCoords = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float)*coordsCount*idC, nullptr, error.errorValP()); error.Eval();

	error = kern.setArg(0,flatSymData);
	error = kern.setArg(1,gradientData);
	error = kern.setArg(2,inCoords);
	error = kern.setArg(3,outCoords);
	error = kern.setArg(4,idC*sizeof(float)*workGroupSize,nullptr); //testValBuf
	error = kern.setArg(5,idC*sizeof(float)*workGroupSize,nullptr); //lastValBuf
	error = kern.setArg(6,expC*sizeof(float)*workGroupSize,nullptr); //evalBuf
	error = kern.setArg(7,expC*idC*sizeof(float)*workGroupSize,nullptr); //gradEvalBuf
	error = kern.setArg(8,idC*idC*sizeof(float)*workGroupSize,nullptr); //vBasisBuf
	error = kern.setArg(9,expC*idC*sizeof(float)*workGroupSize,nullptr); //wBasisBuf
	error = kern.setArg(10,idC*sizeof(float)*workGroupSize,nullptr); //wDivForVBuf
	error = kern.setArg(11,threshold);
	error = kern.setArg(12,maxIter);
	error = kern.setArg(13,expC);
	error = kern.setArg(14,idC);

	error = comQue.enqueueWriteBuffer(flatSymData, false, 0, sizeof(unsigned char)*flat.size, flat.data);
	error = comQue.enqueueWriteBuffer(gradientData, false, 0, sizeof(unsigned char)*grad.size, grad.data);
	error = comQue.enqueueWriteBuffer(inCoords, false, 0, sizeof(float)*coordsCount*idC, in);

	error = comQue.enqueueNDRangeKernel(kern,0,coordsCount,workGroupSize);

	error = comQue.enqueueReadBuffer(outCoords, true, 0, coordsCount*idC, out);
	//error = comQue.enqueueNDRangeKernel(kern,0,64);	

	//error = comQue.enqueueReadBuffer(testOut, true, 0, 64*sizeof(float), fOut.data()); //W

	//User provides FlatSym, Gradient, inCoords, outCoords, threshold, and maxCount
}



#undef float