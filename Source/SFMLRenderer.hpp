#pragma once

class SFMLRenderer;

class SFMLRendererAPI
{
private:
	SFMLRenderer* sfr; //unoptimized PIMPL, shouldn't be creating these often anyway 
public:
	SFMLRendererAPI();
	int RunLoop(); //0 if good, 1 if closed
	void SetTex(unsigned char* tex, int width, int height);
	~SFMLRendererAPI();
};