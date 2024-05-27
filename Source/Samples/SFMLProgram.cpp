#include "Samples.hpp"

#include <SymExp.hpp>
#include <FlatSymExp.hpp>
#include <Helpers.hpp>

#include <FlatSymExpOpenCL.hpp>

#include <SFML/Graphics.hpp>

#include <iostream>
#include <vector>

//To Do:
//Fixup var names for this file and FlatSymExpOpenCL
//Investigate static on redundant example
//Test cross sections of three var functions
//Makefile, building, and api


static void FractalGenImage(FlatSymExp& flat, FlatSymExp& grad, unsigned char* rootTex, FlatKernelExecutor& fke, int width, int height, float centerX, float centerY, float scale, std::vector<float>& roots);

static void SetTex(unsigned char* tex, int width, int height);

static sf::Texture* sfTexPtr;
static sf::Sprite* drawRectanglePtr;

int SFMLMovement()
{
	int width = 800; int height = 800;

	sf::Texture sfTex; sfTexPtr = &sfTex;
	sf::Sprite drawRectangle; drawRectanglePtr = &drawRectangle;


	sf::RenderWindow window;
	window.create(sf::VideoMode(width, height), "Display");
	drawRectangle.setTexture(sfTex, true);
	window.setFramerateLimit(60);


	randSeed(0);
	std::vector<SymExp> rendExp; 



	SymExp hold;

	//SimpleQuad
	//rendExp.push_back(SymExp(1));
	//hold = SymExp(-1); hold.terms.push_back(Product(1, 0, 1)); hold.terms.push_back(Product(1, 1, 1)); rendExp[0] *= hold;
	//hold = SymExp(-1); hold.terms.push_back(Product(-1, 0, 1)); hold.terms.push_back(Product(1, 1, 1)); rendExp[0] *= hold;
	//rendExp.push_back(SymExp(1));
	//hold = SymExp(1); hold.terms.push_back(Product(1, 1, 1)); rendExp[1] *= hold;
	//hold = SymExp(0); hold.terms.push_back(Product(1, 0, 1)); rendExp[1] *= hold;
	//rendExp[0].scalar -= -0.25f;
	//rendExp[1].scalar -= -0.25f;

	/*float m = glm::sqrt(0.5f);
	float n = glm::sqrt(3.0f)/2.0f;
	rendExp.push_back(SymExp(1)); rendExp.push_back(SymExp(1));
	hold = SymExp(0.5); hold.terms.push_back(Product(1, 1, 1)); rendExp[0] *= hold;
	hold = SymExp(-m); hold.terms.push_back(Product(1.5f*m/n, 0, 1)); hold.terms.push_back(Product(m, 1, 1)); rendExp[0] *= hold;
	hold = SymExp(-m); hold.terms.push_back(Product(-1.5f*m/n, 0, 1)); hold.terms.push_back(Product(m, 1, 1)); rendExp[0] *= hold;
	hold = SymExp(0); hold.terms.push_back(Product(1, 0, 1)); rendExp[1] *= hold;
	hold = SymExp(0); hold.terms.push_back(Product(1.0f/n, 0, 1)); hold.terms.push_back(Product(-2, 1, 1)); rendExp[1] *= hold;
	hold = SymExp(0); hold.terms.push_back(Product(-1.0f/n, 0, 1)); hold.terms.push_back(Product(-2, 1, 1)); rendExp[1] *= hold;
	rendExp[0].scalar += 0.05;
	rendExp[1].scalar += 0.05;*/

	//simple redundant, needs custom id associating
	rendExp.push_back(SymExp(1)); rendExp.push_back(SymExp(0)); 
	hold = SymExp(-1); hold.terms.push_back(Product(1, 0, 2)); hold.terms.push_back(Product(1, 1, 2)); rendExp[0] *= hold;
	hold = SymExp(-4); hold.terms.push_back(Product(1, 0, 2)); hold.terms.push_back(Product(1.5f, 1, 2)); rendExp[0] *= hold;
	hold = SymExp(-9); hold.terms.push_back(Product(1.5f, 0, 2)); hold.terms.push_back(Product(1, 1, 2)); rendExp[0] *= hold;
	hold = SymExp(-32); hold.terms.push_back(Product(1, 0, 3)); hold.terms.push_back(Product(2.5f, 1, 2)); rendExp[0] *= hold;


	//rendExp.push_back(SymExp(-1)); rendExp[0].terms.push_back(Product(1,0)); rendExp[0].terms.push_back(Product(1,1));
	//rendExp.push_back(SymExp(0));

	//Random
	//rendExp.push_back(SymExp(GenRandomPoly(2, 4, 4))); rendExp.push_back(GenRandomPoly(2, 4, 4));
	//rendExp[0].scalar = 0; rendExp[1].scalar = 0;

	//NM
	//rendExp.push_back(SymExp(-1)); rendExp.push_back(SymExp(0));
	//rendExp[0].terms.push_back(Product(1, 0, 3));	rendExp[0].terms.push_back(Product(-3, 0, 1)); rendExp[0].terms[1].MultId(1,2);
	//rendExp[1].terms.push_back(Product(-1, 1, 3)); rendExp[1].terms.push_back(Product(3, 0, 2)); rendExp[1].terms[1].MultId(1,1);

	unsigned char* rootTex = new unsigned char[width*height*4];

	std::cout << rendExp[0].ToString() << '\n';
	std::cout << rendExp[1].ToString() << '\n';

	//Start using ExecuteRect
	//Finish RootOrganizing.cl

	FlatSymExp rendFlat(rendExp);
	FlatSymExp rendGradFlat(rendFlat.Gradient());

	std::vector<cl::Platform> platforms;
	std::vector<std::vector<cl::Device>> devices;
	QueryDevices(platforms,devices);
	FlatKernelExecutor fke(devices[0][0]);

	std::vector<float> roots;

	//Kernel for not using coords input
	//Kernel for associating coords out to given roots

	//Kernel for getting roots from coords in? (probs for real time, tho would not sure how to parallelize)
	//Could have work item size == global item size, out for each work item which is the root found, each work item goes through coords until it finds one that isn't already a root, if it doesn't find any new\
	coords it returns nan or sets a bool out buffer (do every few frame, probs not computationally expensive, compound results each frame) probably works fine

	float x = 0;
	float y = 0;
	float scale = 4.0f;


    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

		FractalGenImage(rendFlat, rendGradFlat, rootTex, fke, width, height, x, y, scale, roots);
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::W))
			y += 0.05f*scale;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
			y -= 0.05f*scale;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::A))
			x -= 0.05f*scale;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::D))
			x += 0.05f*scale;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Q))
			scale *= 1.1f;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::E))
			scale /= 1.1f;
		
		SetTex(rootTex, width, height);

        window.clear();
        window.draw(drawRectangle);
        window.display();
    }

	return 0;
}

static void FractalGenImage(FlatSymExp& flat, FlatSymExp& grad, unsigned char* rootTex, FlatKernelExecutor& fke, int width, int height, float centerX, float centerY, float scale, std::vector<float>& roots)
{
	thread_local std::vector<float> coordsOut;
	thread_local std::vector<int> iters;
	thread_local std::vector<int> ids;
	thread_local std::vector<int> indexes;
	coordsOut.resize(width*height*2);
	iters.resize(width*height);
	ids.resize(width*height);
	indexes.resize(width);

	std::vector<float> initial; initial.push_back(centerX-(scale*((float)width/(float)height))); initial.push_back(centerY-scale);
	std::vector<float> base1; base1.push_back(2*scale*(1.0f/(float)height)); base1.push_back(0);
	std::vector<float> base2; base2.push_back(0); base2.push_back(2*scale*(1.0f/(float)height));

	fke.InitializeBuffers(flat.size, grad.size, width, height, 2, roots.size()/2, width);

	fke.RunNMnTomRect(flat, grad, width, height, initial.data(), base1.data(), base2.data(), coordsOut.data(), iters.data(), 0.00001, 4);
	
	
	/*
	//getting ids using convergence points
	fke.RunAssociateCoords(width*height, flat.idCount(), roots.size()/2, coordsOut.data(), roots.data(), iters.data(), ids.data());
	fke.RunFindNewRoots(width, ids.data(), indexes.data(), height, height);
	std::vector<float> newRoots; newRoots.reserve(64);
	//sort through new roots to get unique values
	for (int i = 0; i < width; i++)
	{
		int ind = indexes[i];
		if (ind == -1)
			continue;
		float f1 = coordsOut[(ind*2)+0];
		float f2 = coordsOut[(ind*2)+1];
		bool isNew = true;
		for (int j = 0; j < newRoots.size()/2; j++)
		{
			if (std::abs(f1-newRoots[(j*2)+0]) + std::abs(f2-newRoots[(j*2)+1]) <= 0.1f)
			{
				isNew = false;
				break;
			}
		}
		if (isNew)
		{
			newRoots.push_back(f1);
			newRoots.push_back(f2);
		}
	}
	//sort new roots
	for (int i = 0; i < newRoots.size()/2; i++)
	{
		for (int j = i-1; j >= 0; j--)
		{
			int jm = j+1;
			//if this condition is true, current root is "less" (in this case, less y, or less x if ys are same)
			if (newRoots[(jm*2)+1] < newRoots[(j*2)+1] || (std::abs(newRoots[(jm*2)+1] - newRoots[(j*2)+1]) < 0.01f && newRoots[(jm*2)+0] < newRoots[(j*2)+0]))
			{
				float hold = newRoots[(jm*2)+0];
				newRoots[(jm*2)+0] = newRoots[(j*2)+0];
				newRoots[(j*2)+0] = hold;
				hold = newRoots[(jm*2)+1];
				newRoots[(jm*2)+1] = newRoots[(j*2)+1];
				newRoots[(j*2)+1] = hold;
			}
		}
	}
	{
		int rootIndex = 0;
		int i = 0;
		while (rootIndex < roots.size()/2 && i < newRoots.size()/2)
		{
			if (newRoots[(i*2)+1] < roots[(rootIndex*2)+1] || (std::abs(newRoots[(i*2)+1] - roots[(rootIndex*2)+1]) < 0.01f && newRoots[(i*2)+0] < roots[(rootIndex*2)+0]))
			{
				roots.insert(roots.begin()+rootIndex, newRoots[(i*2)+1]);
				roots.insert(roots.begin()+rootIndex, newRoots[(i*2)+0]);
				i++;
			}
			rootIndex++;
		}
		for (; i < newRoots.size()/2; i++)
		{
			roots.push_back(newRoots[(i*2)+0]);
			roots.push_back(newRoots[(i*2)+1]);
		}
	}
	if (roots.size() > 1024*2)
	{
		std::cout << "Abnormally high number of roots, likely fork bomb\n";
		roots.clear();
		abort();
	}
	*/
	
	//Id association for redundant example
	for (int i = 0; i < width*height; i++)
	{
		float rendVal[2]; rendVal[0] = coordsOut[(i*2)+0]; rendVal[1] = coordsOut[(i*2)+1];
		if (rendVal[0]*rendVal[0] + rendVal[1]*rendVal[1] < 1.5f)
			ids[i] = 0;
		else if (rendVal[0]*rendVal[0] + rendVal[1]*rendVal[1] < 4.5f)
			ids[i] = 1;
		else if (rendVal[0]*rendVal[0] + rendVal[1]*rendVal[1] < 9.5f)
			ids[i] = 2;
		else
			ids[i] = 3;
	}
	fke.comQue.enqueueWriteBuffer(fke.idsBuf, false, 0, sizeof(int)*width*height, ids.data());
	

	//coloring
	fke.RunColorTex(width*height, ids.data(), iters.data(), rootTex);
	/*for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			switch (ids[(height*x)+(height-1-y)])
			{
			default: rootTex[(width*4*y)+(4*x)+0] = 128; rootTex[(width*4*y)+(4*x)+1] = 128; rootTex[(width*4*y)+(4*x)+2] = 128; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 0: rootTex[(width*4*y)+(4*x)+0] = 255; rootTex[(width*4*y)+(4*x)+1] = 0; rootTex[(width*4*y)+(4*x)+2] = 0; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 1: rootTex[(width*4*y)+(4*x)+0] = 0; rootTex[(width*4*y)+(4*x)+1] = 255; rootTex[(width*4*y)+(4*x)+2] = 0; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 2: rootTex[(width*4*y)+(4*x)+0] = 0; rootTex[(width*4*y)+(4*x)+1] = 0; rootTex[(width*4*y)+(4*x)+2] = 255; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 3: rootTex[(width*4*y)+(4*x)+0] = 0; rootTex[(width*4*y)+(4*x)+1] = 255; rootTex[(width*4*y)+(4*x)+2] = 255; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 4: rootTex[(width*4*y)+(4*x)+0] = 255; rootTex[(width*4*y)+(4*x)+1] = 0; rootTex[(width*4*y)+(4*x)+2] = 255; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 5: rootTex[(width*4*y)+(4*x)+0] = 255; rootTex[(width*4*y)+(4*x)+1] = 255; rootTex[(width*4*y)+(4*x)+2] = 0; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 6: rootTex[(width*4*y)+(4*x)+0] = 128; rootTex[(width*4*y)+(4*x)+1] = 0; rootTex[(width*4*y)+(4*x)+2] = 0; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 7: rootTex[(width*4*y)+(4*x)+0] = 0; rootTex[(width*4*y)+(4*x)+1] = 128; rootTex[(width*4*y)+(4*x)+2] = 0; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 8: rootTex[(width*4*y)+(4*x)+0] = 0; rootTex[(width*4*y)+(4*x)+1] = 0; rootTex[(width*4*y)+(4*x)+2] = 128; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 9: rootTex[(width*4*y)+(4*x)+0] = 0; rootTex[(width*4*y)+(4*x)+1] = 128; rootTex[(width*4*y)+(4*x)+2] = 128; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 10: rootTex[(width*4*y)+(4*x)+0] = 128; rootTex[(width*4*y)+(4*x)+1] = 0; rootTex[(width*4*y)+(4*x)+2] = 128; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 11: rootTex[(width*4*y)+(4*x)+0] = 128; rootTex[(width*4*y)+(4*x)+1] = 128; rootTex[(width*4*y)+(4*x)+2] = 0; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 12: rootTex[(width*4*y)+(4*x)+0] = 180; rootTex[(width*4*y)+(4*x)+1] = 60; rootTex[(width*4*y)+(4*x)+2] = 60; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 13: rootTex[(width*4*y)+(4*x)+0] = 60; rootTex[(width*4*y)+(4*x)+1] = 180; rootTex[(width*4*y)+(4*x)+2] = 60; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 14: rootTex[(width*4*y)+(4*x)+0] = 60; rootTex[(width*4*y)+(4*x)+1] = 60; rootTex[(width*4*y)+(4*x)+2] = 180; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 15: rootTex[(width*4*y)+(4*x)+0] = 60; rootTex[(width*4*y)+(4*x)+1] = 180; rootTex[(width*4*y)+(4*x)+2] = 180; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 16: rootTex[(width*4*y)+(4*x)+0] = 180; rootTex[(width*4*y)+(4*x)+1] = 60; rootTex[(width*4*y)+(4*x)+2] = 180; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			case 17: rootTex[(width*4*y)+(4*x)+0] = 180; rootTex[(width*4*y)+(4*x)+1] = 180; rootTex[(width*4*y)+(4*x)+2] = 60; rootTex[(width*4*y)+(4*x)+3] = 255; break;
			}
		}
	}*/
}

static void SetTex(unsigned char* tex, int width, int height)
{
    sf::Image sfImg; sfImg.create(width, height, tex);
    sfTexPtr->loadFromImage(sfImg);
    drawRectanglePtr->setTexture(*sfTexPtr, true);
    drawRectanglePtr->setScale(width/drawRectanglePtr->getTextureRect().width, height/drawRectanglePtr->getTextureRect().height);
}