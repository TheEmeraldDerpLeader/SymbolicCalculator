#pragma once

#include <fstream>
#include <iostream>

class StreamBuffer
{
public:
	std::fstream stream;
	unsigned char buf[512];
	int fileIndex = 0;
	int bufIndex = 0;
	int validBuf = 0;

	StreamBuffer();
	StreamBuffer(std::string filename);

	unsigned char Pop();
	unsigned char Read(int offset);

	void Seek(int location);
	inline bool CanPop() { return !stream.eof() || bufIndex < validBuf; }

	unsigned char operator[](int offset);
	StreamBuffer& operator+=(int offset);
};

void randSeed(unsigned int seed);
int randI(int firstInt, int secondInt);
float randF(double firstFloat, double secondFloat);