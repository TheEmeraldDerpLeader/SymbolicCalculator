#include "FlatSymExpOpenCL.hpp"

#define float float

struct ErrorHandle
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

	comQue = cl::CommandQueue(context, device, 0, error.errorValP());
	error.Eval();

	std::vector<std::string> sources;
	sources.push_back(std::string());
	sources.push_back(std::string());


	std::string filePath = "Assets/Kernels/NMnTom.cl";
	StreamBuffer sb(filePath);
	while (sb.CanPop())
		sources[0].push_back(sb.Pop());

	filePath = "Assets/Kernels/RootOrganizing.cl";
	sb = StreamBuffer(filePath);
	while (sb.CanPop())
		sources[1].push_back(sb.Pop());

	prog = cl::Program(context, sources, error.errorValP());
	error.Eval();
	error = prog.build("");

	if (error.errorVal != 0)
	{
		std::string buildLog;
		error = prog.getBuildInfo(device, CL_PROGRAM_BUILD_LOG, &buildLog);
		std::cout << buildLog << '\n';
	}
}

void FlatKernelExecutor::InitializeBuffers(int flatSize, int gradSize, int width, int height, int idC, int rootC, int checkCount)
{
	int coordsCount = width*height;

	if (coordsCount == 0 || flatSize == 0 || gradSize == 0 || idC == 0)
	{
		std::cout << "Invalid inputs for initial buffer sizes (essential buffers)\n";
		abort();
	}

	flatSymData = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*flatSize, nullptr, error.errorValP()); error.Eval();
	gradientData = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(unsigned char)*gradSize, nullptr, error.errorValP()); error.Eval();
	testCoordsBuf = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float)*coordsCount*idC, nullptr, error.errorValP()); error.Eval();
	finalCoordsBuf = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float)*coordsCount*idC, nullptr, error.errorValP()); error.Eval();
	itersBuf = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(int)*coordsCount, nullptr, error.errorValP()); error.Eval();
	initialBuf = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float)*idC, nullptr, error.errorValP()); error.Eval();
	base1Buf = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float)*idC, nullptr, error.errorValP()); error.Eval();
	base2Buf = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float)*idC, nullptr, error.errorValP()); error.Eval();
	idsBuf = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(int)*coordsCount, nullptr, error.errorValP()); error.Eval();
	texBuf = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(unsigned char)*4*coordsCount, nullptr, error.errorValP()); error.Eval();
	
	InitializeRootBuffer(idC, rootC);

	if (checkCount != 0)
		indexesBuf = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(int)*checkCount, nullptr, error.errorValP()); error.Eval(); //checkCount of 0 causes RunFindNewRoots to skip 
}

void FlatKernelExecutor::InitializeRootBuffer(int idC, int rootC)
{
	if (rootC == 0)
	{
		rootsBuf = cl::Buffer(context, CL_MEM_READ_ONLY, 1, nullptr, error.errorValP()); error.Eval(); //having no roots is a valid input
	}
	else
	{
		rootsBuf = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float)*rootC*idC, nullptr, error.errorValP()); error.Eval();
	}
}

void FlatKernelExecutor::RunNMnTom(FlatSymExp& flat, FlatSymExp& grad, int coordsCount, float* testCoords, float* finalCoords, int* iters, float threshold, int maxIter)
{
	kern = cl::Kernel(prog, "NMnTom", error.errorValP());
	error.Eval();

	int expC = flat.expCount();
	int idC = flat.idCount();

	size_t workGroupSize;
	error = kern.getWorkGroupInfo(device, CL_KERNEL_WORK_GROUP_SIZE, &workGroupSize);

	if (idC == 0 || expC == 0 || workGroupSize == 0 || coordsCount == 0)
		return;

	error = kern.setArg(0,flatSymData);
	error = kern.setArg(1,gradientData);
	error = kern.setArg(2,testCoordsBuf);
	error = kern.setArg(3,finalCoordsBuf);
	error = kern.setArg(4,itersBuf);
	error = kern.setArg(5,idC*sizeof(float)*workGroupSize,nullptr); //testValBuf
	error = kern.setArg(6,idC*sizeof(float)*workGroupSize,nullptr); //lastValBuf
	error = kern.setArg(7,expC*sizeof(float)*workGroupSize,nullptr); //evalBuf
	error = kern.setArg(8,expC*idC*sizeof(float)*workGroupSize,nullptr); //gradEvalBuf
	error = kern.setArg(9,idC*idC*sizeof(float)*workGroupSize,nullptr); //vBasisBuf
	error = kern.setArg(10,expC*expC*sizeof(float)*workGroupSize,nullptr); //wBasisBuf
	error = kern.setArg(11,idC*sizeof(float)*workGroupSize,nullptr); //wDivForVBuf
	error = kern.setArg(12,threshold);
	error = kern.setArg(13,maxIter);
	error = kern.setArg(14,expC);
	error = kern.setArg(15,idC);

	error = comQue.enqueueWriteBuffer(flatSymData, false, 0, sizeof(unsigned char)*flat.size, flat.data);
	error = comQue.enqueueWriteBuffer(gradientData, false, 0, sizeof(unsigned char)*grad.size, grad.data);
	error = comQue.enqueueWriteBuffer(testCoordsBuf, false, 0, sizeof(float)*coordsCount*idC, testCoords);

	int remLocal = coordsCount%workGroupSize;
	error = comQue.enqueueNDRangeKernel(kern,0,coordsCount-remLocal, workGroupSize);
	if (remLocal != 0)
		error = comQue.enqueueNDRangeKernel(kern, coordsCount-remLocal, remLocal, remLocal);


	error = comQue.enqueueReadBuffer(finalCoordsBuf, false, 0, sizeof(float)*coordsCount*idC, finalCoords);
	error = comQue.enqueueReadBuffer(itersBuf, true, 0, sizeof(int)*coordsCount, iters);

	//error = comQue.enqueueNDRangeKernel(kern,0,64);	

	//error = comQue.enqueueReadBuffer(testOut, true, 0, 64*sizeof(float), fOut.data()); //W

	//User provides FlatSym, Gradient, inCoords, outCoords, threshold, and maxCount
}

void FlatKernelExecutor::RunNMnTomRect(FlatSymExp& flat, FlatSymExp& grad, int coordsX, int coordsY, float* initial, float* base1, float* base2, float* finalCoords, int* iters, float threshold, int maxIter)
{
	kern = cl::Kernel(prog, "NMnTomRect", error.errorValP());
	error.Eval();

	int expC = flat.expCount();
	int idC = flat.idCount();

	size_t workGroupSize;
	error = kern.getWorkGroupInfo(device, CL_KERNEL_WORK_GROUP_SIZE, &workGroupSize);

	if (idC == 0 || expC == 0 || workGroupSize == 0 || coordsX*coordsY == 0)
		return;

	error = kern.setArg(0,flatSymData);
	error = kern.setArg(1,gradientData);
	error = kern.setArg(2,initialBuf);
	error = kern.setArg(3,base1Buf);
	error = kern.setArg(4,base2Buf);
	error = kern.setArg(5,finalCoordsBuf);
	error = kern.setArg(6,itersBuf);
	error = kern.setArg(7,idC*sizeof(float)*workGroupSize,nullptr); //testValBuf
	error = kern.setArg(8,idC*sizeof(float)*workGroupSize,nullptr); //lastValBuf
	error = kern.setArg(9,expC*sizeof(float)*workGroupSize,nullptr); //evalBuf
	error = kern.setArg(10,expC*idC*sizeof(float)*workGroupSize,nullptr); //gradEvalBuf
	error = kern.setArg(11,idC*idC*sizeof(float)*workGroupSize,nullptr); //vBasisBuf
	error = kern.setArg(12,expC*expC*sizeof(float)*workGroupSize,nullptr); //wBasisBuf
	error = kern.setArg(13,idC*sizeof(float)*workGroupSize,nullptr); //wDivForVBuf
	error = kern.setArg(14,threshold);
	error = kern.setArg(15,maxIter);
	error = kern.setArg(16,expC);
	error = kern.setArg(17,idC);
	error = kern.setArg(18,coordsX);
	error = kern.setArg(19,coordsY);

	error = comQue.enqueueWriteBuffer(flatSymData, false, 0, sizeof(unsigned char)*flat.size, flat.data);
	error = comQue.enqueueWriteBuffer(gradientData, false, 0, sizeof(unsigned char)*grad.size, grad.data);
	error = comQue.enqueueWriteBuffer(initialBuf, false, 0, sizeof(float)*idC, initial);
	error = comQue.enqueueWriteBuffer(base1Buf, false, 0, sizeof(float)*idC, base1);
	error = comQue.enqueueWriteBuffer(base2Buf, false, 0, sizeof(float)*idC, base2);



	int remLocal = (coordsX*coordsY)%workGroupSize;
	error = comQue.enqueueNDRangeKernel(kern,0,(coordsX*coordsY)-remLocal, workGroupSize);
	if (remLocal != 0)
		error = comQue.enqueueNDRangeKernel(kern, (coordsX*coordsY)-remLocal, remLocal, remLocal);

	error = comQue.enqueueReadBuffer(finalCoordsBuf, false, 0, sizeof(float)*coordsX*coordsY*idC, finalCoords);
	error = comQue.enqueueReadBuffer(itersBuf, true, 0, sizeof(int)*coordsX*coordsY, iters);
	//error = comQue.enqueueNDRangeKernel(kern,0,64);	

	//error = comQue.enqueueReadBuffer(testOut, true, 0, 64*sizeof(float), fOut.data()); //W

	//User provides FlatSym, Gradient, inCoords, outCoords, threshold, and maxCount
}

void FlatKernelExecutor::RunAssociateCoords(int coordsCount, int idC, int rootC, float* finalCoords, float* roots, int* iters, int* ids)
{
	kern = cl::Kernel(prog, "AssociateCoords", error.errorValP());
	error.Eval();

	size_t workGroupSize;
	error = kern.getWorkGroupInfo(device, CL_KERNEL_WORK_GROUP_SIZE, &workGroupSize);

	if (workGroupSize == 0 || coordsCount == 0)
		return;

	if (rootC == 0 && false)
	{
		for (int i = 0; i < coordsCount; i++)
		{
			if (iters[i] == -1)
				ids[i] = -2;
			else
				ids[i] = -1;
		}
		idsBuf = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(int)*coordsCount, nullptr, error.errorValP()); error.Eval();
		comQue.enqueueWriteBuffer(idsBuf, false, 0, sizeof(int)*coordsCount, ids);
		return;
	}

	error = kern.setArg(0,finalCoordsBuf);
	error = kern.setArg(1,rootsBuf);
	error = kern.setArg(2,itersBuf);
	error = kern.setArg(3,idsBuf);
	error = kern.setArg(4,sizeof(float)*idC*workGroupSize, nullptr); //coordBuf
	error = kern.setArg(5,idC);
	error = kern.setArg(6,rootC);

	//error = comQue.enqueueWriteBuffer(outCoordsBuf, false, 0, sizeof(float)*idC*coordsCount, outCoords);
	if (rootC != 0)
		error = comQue.enqueueWriteBuffer(rootsBuf, false, 0, sizeof(float)*idC*rootC, roots);
	//error = comQue.enqueueWriteBuffer(itersBuf, false, 0, sizeof(int)*coordsCount, iters);

	int remLocal = coordsCount%workGroupSize;
	error = comQue.enqueueNDRangeKernel(kern,0,coordsCount-remLocal, workGroupSize);
	if (remLocal != 0)
		error = comQue.enqueueNDRangeKernel(kern, coordsCount-remLocal, remLocal, remLocal);

	error = comQue.enqueueReadBuffer(idsBuf, true, 0, sizeof(int)*coordsCount, ids);
}

void FlatKernelExecutor::RunFindNewRoots(int checkCount, int* ids, int* indexes, int checkSpan, int checkOffset)
{
	kern = cl::Kernel(prog, "FindNewRoots", error.errorValP());
	error.Eval();

	size_t workGroupSize;
	error = kern.getWorkGroupInfo(device, CL_KERNEL_WORK_GROUP_SIZE, &workGroupSize);

	if (workGroupSize == 0 || checkCount == 0)
		return;

	error = kern.setArg(0,idsBuf);
	error = kern.setArg(1,indexesBuf);
	error = kern.setArg(2,checkSpan);
	error = kern.setArg(3,checkOffset);

	int remLocal = checkCount%workGroupSize;
	error = comQue.enqueueNDRangeKernel(kern,0,checkCount-remLocal, workGroupSize);
	if (remLocal != 0)
		error = comQue.enqueueNDRangeKernel(kern, checkCount-remLocal, remLocal, remLocal);

	error = comQue.enqueueReadBuffer(indexesBuf, true, 0, sizeof(int)*checkCount, indexes);
}

void FlatKernelExecutor::RunColorTex(int coordsCount, int* iters, unsigned char* tex)
{
	kern = cl::Kernel(prog, "ColorTex", error.errorValP());
	error.Eval();

	size_t workGroupSize;
	error = kern.getWorkGroupInfo(device, CL_KERNEL_WORK_GROUP_SIZE, &workGroupSize);

	if (workGroupSize == 0 || coordsCount == 0)
		return;

	error = kern.setArg(0,idsBuf);
	error = kern.setArg(1,itersBuf);
	error = kern.setArg(2,texBuf);

	int remLocal = coordsCount%workGroupSize;
	error = comQue.enqueueNDRangeKernel(kern,0,coordsCount-remLocal, workGroupSize);
	if (remLocal != 0)
		error = comQue.enqueueNDRangeKernel(kern, coordsCount-remLocal, remLocal, remLocal);

	error = comQue.enqueueReadBuffer(texBuf, true, 0, sizeof(unsigned char)*4*coordsCount, tex);
}

void FlatKernelExecutor::RunColorTex(int coordsCount, int* iters, unsigned char* tex, int* ids)
{
	kern = cl::Kernel(prog, "ColorTex", error.errorValP());
	error.Eval();

	size_t workGroupSize;
	error = kern.getWorkGroupInfo(device, CL_KERNEL_WORK_GROUP_SIZE, &workGroupSize);

	if (workGroupSize == 0 || coordsCount == 0)
		return;

	error = kern.setArg(0,idsBuf);
	error = kern.setArg(1,itersBuf);
	error = kern.setArg(2,texBuf);

	error = comQue.enqueueWriteBuffer(idsBuf, false, 0, sizeof(int)*coordsCount, ids);

	int remLocal = coordsCount%workGroupSize;
	error = comQue.enqueueNDRangeKernel(kern,0,coordsCount-remLocal, workGroupSize);
	if (remLocal != 0)
		error = comQue.enqueueNDRangeKernel(kern, coordsCount-remLocal, remLocal, remLocal);

	error = comQue.enqueueReadBuffer(texBuf, true, 0, sizeof(unsigned char)*4*coordsCount, tex);
}

void FindAndCollectRootsCPU(int idC, int coordsCount, std::vector<float>& roots, float* outCoords, int* iters, int* ids)
{
	roots.clear();
	std::vector<int> rootIndexes;
	for (int i = 0; i < coordsCount; i++)
	{
		if (iters[i] == -1)
		{
			ids[i] = -2;
			continue;
		}
		float* coord = outCoords+(idC*i);

		//check if equal to roots
		int id = -1;
		for (int i = 0; i < roots.size()/idC; i++)
		{
			float dif = 0;
			int rootIndex = rootIndexes[i];
			for (int j = 0; j < idC; j++)
				dif += std::abs(coord[j] - roots[(rootIndex*idC)+j]);
			if (dif <= 0.1f)
			{
				id = i;
				break;
			}
		}

		if (id != -1)
		{
			ids[i] = id;
			continue;
		}
		ids[i] = rootIndexes.size();

		//add to roots and sort
		rootIndexes.push_back(roots.size()/idC);
		for (int i = 0; i < idC; i++)
			roots.push_back(coord[i]);
		for (int i = (roots.size()/idC)-1; i > 0; i--)
		{
			bool isLess = false;
			for (int j = idC-1; j >= 0; j--)
			{
				// dif =      current       -    previous    
				float dif = roots[(i*idC)+j]-roots[((i-1)*idC)+j];
				if (j == idC-1 || dif >= 0.001f || dif <= -0.001f)
				{
					isLess = (dif < 0.0f);
					break;
				}
			}
			if (isLess == true)
			{
				for (int j = 0; j < idC; j++)
				{
					float hold = roots[(i*idC)+j];
					roots[(i*idC)+j] = roots[((i-1)*idC)+j];
					roots[((i-1)*idC)+j] = hold;
				}
				for (int j = 0; j < rootIndexes.size()-1; j++)
					if (rootIndexes[j] == i-1)
					{
						rootIndexes[j]++;
						break;
					}
				rootIndexes[rootIndexes.size()-1]--;
			}
			else
				break;
		}
	}

	//assigning to id
	for (int i = 0; i < coordsCount; i++)
	{
		if (ids[i] == -2)
			continue;
		ids[i] = rootIndexes[ids[i]];
	}
}

void CollectNewRoots(int idC, int checkCount, std::vector<float>& roots, float* coordsOut, int* indexes)
{
	std::vector<float> newRoots; newRoots.reserve(64);
	//sort through new roots to get unique values
	for (int i = 0; i < checkCount; i++)
	{
		int ind = indexes[i];
		if (ind == -1)
			continue;

		float* coord = coordsOut+(ind*idC);
		bool isNew = true;
		for (int i = 0; i < newRoots.size()/idC; i++)
		{
			float sum = 0;
			for (int j = 0; j < idC; j++)
				sum += std::abs(coord[j]-newRoots[(i*idC)+j]);
			if (sum <= 0.1f)
			{
				isNew = false;
				break;
			}
		}
		if (isNew)
			for (int i = 0; i < idC; i++)
				newRoots.push_back(coord[i]);
	}
	//sort new roots
	for (int i = 0; i < newRoots.size()/idC; i++)
	{
		for (int j = i; j > 0; j--)
		{
			//if this condition is true, current root is "less" (in this case, less y, or less x if ys are close)
			bool isLess = false;
			for (int k = idC-1; k >= 0; k--)
			{
				// dif =      current       -    previous    
				float dif = newRoots[(j*idC)+k]-newRoots[((j-1)*idC)+k];
				if (k == idC-1 || dif >= 0.001f || dif <= -0.001f)
				{
					isLess = (dif < 0.0f);
					break;
				}
			}
			if (isLess)
			{
				for (int k = 0; k < idC; k++)
				{
					float hold = newRoots[(j*idC)+k];
					newRoots[(j*idC)+k] = newRoots[((j-1)*idC)+k];
					newRoots[((j-1)*idC)+k] = hold;
				}
			}
			else
				break;
		}
	}
	{
		int rootIndex = 0;
		int i = 0;
		while (rootIndex < roots.size()/idC && i < newRoots.size()/idC)
		{
			bool isLess = false;
			for (int j = idC-1; j >= 0; j--)
			{
				// dif =      current      -     root    
				float dif = newRoots[(i*idC)+j]-roots[(rootIndex*idC)+j];
				if (j == idC-1 || dif >= 0.001f || dif <= -0.001f)
				{
					isLess = (dif < 0.0f);
					break;
				}
			}
			if (isLess)
			{
				for (int j = idC-1; j >= 0; j--)
					roots.insert(roots.begin()+rootIndex, newRoots[(i*idC)+j]);
				i++;
			}
			rootIndex++;
		}
		for (; i < newRoots.size()/idC; i++)
		{
			for (int j = 0; j < idC; j++)
				roots.push_back(newRoots[(i*idC)+j]);
		}
	}
	if (roots.size() > 2048*2)
	{
		std::cout << "Abnormally high number of roots, either root set has redundancy or root set is changing\n";
		roots.clear();
		//abort();
	}
}



#undef float