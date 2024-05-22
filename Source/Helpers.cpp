#include "Helpers.hpp"

#include <random>


StreamBuffer::StreamBuffer()
{
	for (int i = 0; i < 512; i++)
		buf[i] = 0;
}
StreamBuffer::StreamBuffer(std::string filename)
{
	stream = std::fstream(filename, std::ios::binary | std::ios::in | std::ios::out);
	Seek(0);
}

unsigned char StreamBuffer::Pop()
{
	if (bufIndex>=validBuf)
	{
		Seek(fileIndex+bufIndex);
		bufIndex = 1;
		return buf[0];
	}
	else
	{
		bufIndex++;
		return buf[bufIndex-1];
	}
}

unsigned char StreamBuffer::Read(int offset)
{
	if (offset >= 512)
	{
		std::cout << "Cannot read past end of buffer\n";
		abort();
	}
	if (bufIndex+offset >= validBuf)
		Seek(fileIndex+bufIndex);
	if (bufIndex+offset >= validBuf)
		return 'X';
	else
		return buf[bufIndex+offset];
}

void StreamBuffer::Seek(int location)
{
	fileIndex = location;
	bufIndex = 0;

	stream.clear();
	stream.seekg(location);
	stream.read(reinterpret_cast<char*>(buf),512);
	validBuf = stream.gcount();
}

unsigned char StreamBuffer::operator[](int offset)
{
	return Read(offset);
}

StreamBuffer& StreamBuffer::operator+=(int offset)
{
	bufIndex += offset;
	return *this;
}

std::mt19937 rng;

struct __SetupRand_ { public: __SetupRand_() { randSeed(time(NULL)); } }; __SetupRand_ __sr_;

void randSeed(unsigned int seed)
{
	rng.seed(seed);
}
int randI(int firstInt, int secondInt) ///inclusive
{
	if (firstInt > secondInt)
		std::swap(firstInt, secondInt);
	std::uniform_int_distribution<> dist(firstInt, secondInt);
	return dist(rng);
}
float randF(double firstFloat, double secondFloat)
{
	if (firstFloat > secondFloat)
		std::swap(firstFloat, secondFloat);
	std::uniform_real_distribution<> dist(firstFloat, secondFloat);
	return dist(rng);
}