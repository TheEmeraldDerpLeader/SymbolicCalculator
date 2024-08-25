#include "Samples.hpp"

#include <SymExp.hpp>
#include <FlatSymExp.hpp>
#include <Helpers.hpp>

#include <FlatSymExpOpenCL.hpp>

#include <SFML/Graphics.hpp>

#include <iostream>
#include <vector>
#include <chrono>

//To Do:
//Test cross sections of three var functions (done, kinda meh)
//makefile, building


static void FractalGenImage2D(FlatSymExp& flat, FlatSymExp& grad, unsigned char* tex, FlatKernelExecutor& fke, int width, int height, float centerX, float centerY, float scale, std::vector<float>& roots);
static void FractalGenImage(FlatSymExp& flat, FlatSymExp& grad, unsigned char* tex, FlatKernelExecutor& fke, int width, int height, std::vector<float>& center, std::vector<float>& base1, std::vector<float>& base2, float scale, std::vector<float>& roots);
static void FractalGenImageChangingRoots(FlatSymExp& flat, FlatSymExp& grad, unsigned char* tex, FlatKernelExecutor& fke, int width, int height, std::vector<float>& center, std::vector<float>& base1, std::vector<float>& base2, float scale, std::vector<float>& roots);

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
	//rendExp.push_back(SymExp(1)); rendExp.push_back(SymExp(0)); 
	//hold = SymExp(-1); hold.terms.push_back(Product(1, 0, 2)); hold.terms.push_back(Product(1, 1, 2)); rendExp[0] *= hold;
	//hold = SymExp(-4); hold.terms.push_back(Product(1, 0, 2)); hold.terms.push_back(Product(1.5f, 1, 2)); rendExp[0] *= hold;
	//hold = SymExp(-9); hold.terms.push_back(Product(1.5f, 0, 2)); hold.terms.push_back(Product(1, 1, 2)); rendExp[0] *= hold;
	//hold = SymExp(-32); hold.terms.push_back(Product(1, 0, 3)); hold.terms.push_back(Product(2.5f, 1, 2)); rendExp[0] *= hold;

	//rendExp.push_back(SymExp(-1)); rendExp[0].terms.push_back(Product(1,0)); rendExp[0].terms.push_back(Product(1,1));
	//rendExp.push_back(SymExp(0));

	//Random
	//rendExp.push_back(SymExp(GenRandomPoly(2, 4, 4))); rendExp.push_back(GenRandomPoly(2, 4, 4));
	//rendExp[0].scalar = 0; rendExp[1].scalar = 0;
	//Random Triple
	//rendExp.push_back(SymExp(GenRandomPoly(3, 2, 4))); rendExp.push_back(GenRandomPoly(2, 2, 4)); rendExp.push_back(GenRandomPoly(2, 2, 4));
	//rendExp[0].scalar = 0; rendExp[1].scalar = 0; rendExp[2].scalar = 0;

	//NM
	//rendExp.push_back(SymExp(-1)); rendExp.push_back(SymExp(0));
	//rendExp[0].terms.push_back(Product(1, 0, 3));	rendExp[0].terms.push_back(Product(-3, 0, 1)); rendExp[0].terms[1].MultId(1,2);
	//rendExp[1].terms.push_back(Product(-1, 1, 3)); rendExp[1].terms.push_back(Product(3, 0, 2)); rendExp[1].terms[1].MultId(1,1);

	//Alt NM
	SymExp testExpNM(-1);
	testExpNM.terms.push_back(Product(1, 0, 3)); testExpNM.terms.push_back(Product(0, 0, 2)); testExpNM.terms.push_back(Product(0, 0, 1)); //testExpNM.terms[2].MultId(2, 1);
	rendExp = ComplexPoly(testExpNM,2);

	unsigned char* rootTex = new unsigned char[width*height*4];

	std::cout << rendExp[0].ToString() << '\n';
	std::cout << rendExp[1].ToString() << '\n';

	FlatSymExp rendFlat(rendExp);
	FlatSymExp rendGradFlat(rendFlat.Gradient());

	std::vector<cl::Platform> platforms;
	std::vector<std::vector<cl::Device>> devices;
	QueryDevices(platforms,devices);
	FlatKernelExecutor fke(devices[0][0]);

	std::vector<float> roots;

	std::vector<float> center; center.push_back(0); center.push_back(0); center.push_back(0);
	std::vector<float> base1; base1.push_back(1); base1.push_back(0); base1.push_back(0);
	std::vector<float> base2; base2.push_back(0); base2.push_back(1); base2.push_back(0);
	float scale = 4.0f;

	sf::Text fpsText;
	fpsText.setPosition(20, 750);
	fpsText.setCharacterSize(30);
	fpsText.setString("");

	sf::Font font;
	font.loadFromFile("Assets\\Fonts\\LeagueMono-2.300\\static\\TTF\\LeagueMono-Regular.ttf");
	fpsText.setFont(font);

	sf::Clock fpsClock;
	fpsClock.restart();

	float animTime = 0;

    while (window.isOpen())
    {
		sf::Time dif = fpsClock.getElapsedTime();
		if (dif != sf::Time::Zero)
			fpsText.setString(std::to_string(int(1.0f/dif.asSeconds())));
		fpsClock.restart();
		animTime += dif.asSeconds();
		if (animTime > 10)
			animTime -= 10;

		testExpNM.terms[2].coeff = 0.01f*(9.57f-30.0f);
		rendExp = ComplexPoly(testExpNM);
		rendExp[1].terms[1].coeff = -3.0f+(animTime-5.0f)*-0.18f;
		rendFlat = FlatSymExp(rendExp);
		rendGradFlat = rendFlat.Gradient();


        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

		//FractalGenImage2D(rendFlat, rendGradFlat, rootTex, fke, width, height, center[0], center[1], scale, roots);
		//FractalGenImage(rendFlat, rendGradFlat, rootTex, fke, width, height, center, base1, base2, scale, roots);
		FractalGenImageChangingRoots(rendFlat, rendGradFlat, rootTex, fke, width, height, center, base1, base2, scale, roots);


		if (sf::Keyboard::isKeyPressed(sf::Keyboard::W))
			center[1] += 0.05f*scale;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
			center[1] -= 0.05f*scale;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::A))
			center[0] -= 0.05f*scale;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::D))
			center[0] += 0.05f*scale;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Q))
			scale *= 1.1f;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::E))
			scale /= 1.1f;
		if (center.size() > 2)
		{
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Num1))
				center[2] -= 0.05f*scale;
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Num3))
				center[2] += 0.05f*scale;
		}
		
		SetTex(rootTex, width, height);

        window.clear();
        window.draw(drawRectangle);
		window.draw(fpsText);
        window.display();
    }

	return 0;
}

thread_local std::vector<float> coordsOut;
thread_local std::vector<int> iters;
thread_local std::vector<int> ids;
thread_local std::vector<int> indexes;

static void FractalGenImage2D(FlatSymExp& flat, FlatSymExp& grad, unsigned char* tex, FlatKernelExecutor& fke, int width, int height, float centerX, float centerY, float scale, std::vector<float>& roots)
{
	coordsOut.resize(width*height*2);
	iters.resize(width*height);
	ids.resize(width*height);
	indexes.resize(width);

	std::vector<float> initial; initial.push_back(centerX-(scale*((float)width/(float)height))); initial.push_back(centerY-scale);
	std::vector<float> base1; base1.push_back(2*scale*(1.0f/(float)height)); base1.push_back(0);
	std::vector<float> base2; base2.push_back(0); base2.push_back(2*scale*(1.0f/(float)height));

	//getting ids using convergence points
	fke.InitializeBuffers(flat.size, grad.size, width, height, 2, roots.size()/2, width);
	fke.RunNMnTomRect(flat, grad, width, height, initial.data(), base1.data(), base2.data(), coordsOut.data(), iters.data(), 0.00001, 60);	
	fke.RunAssociateCoords(width*height, flat.idCount(), roots.size()/2, coordsOut.data(), roots.data(), iters.data(), ids.data());
	fke.RunFindNewRoots(width, ids.data(), indexes.data(), height, height);
	CollectNewRoots(2, indexes.size(), roots, coordsOut.data(), indexes.data());
	fke.RunColorTex(width*height, iters.data(), tex);
	
	//getting ids from potentially changing roots (use for animations) uncomment second RunAssociateCoords and initBuffers
	//fke.InitializeBuffers(flat.size, grad.size, width, height, 2, roots.size()/2, width*height);
	//fke.RunNMnTomRect(flat, grad, width, height, initial.data(), base1.data(), base2.data(), coordsOut.data(), iters.data(), 0.00001, 60);	
	//indexes.resize(width*height);
	//FindAndCollectRootsCPU(2, coordsOut.size()/2, roots, coordsOut.data(), iters.data(), ids.data());
	//fke.RunColorTex(width*height, iters.data(), tex, ids.data());
	
	/*
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
	fke.RunColorTex(width*height, iters.data(), tex);
	*/
}

static void FractalGenImage(FlatSymExp& flat, FlatSymExp& grad, unsigned char* tex, FlatKernelExecutor& fke, int width, int height, std::vector<float>& center, std::vector<float>& base1, std::vector<float>& base2, float scale, std::vector<float>& roots)
{
	int idC = flat.idCount();

	coordsOut.resize(width*height*idC);
	iters.resize(width*height);
	ids.resize(width*height);
	indexes.resize(width);

	std::vector<float> initial; initial.resize(idC);
	std::vector<float> mBase1; mBase1.resize(idC);
	std::vector<float> mBase2; mBase2.resize(idC);
	for (int i = 0; i < idC; i++)
	{
		initial[i] = center[i] - scale*((base1[i]*float(width)/float(height)) + base2[i]);
		mBase1[i] = 2.0f*scale*base1[i]/float(height);
		mBase2[i] = 2.0f*scale*base2[i]/float(width);
	}

	fke.InitializeBuffers(flat.size, grad.size, width, height, idC, roots.size()/idC, width);
	fke.RunNMnTomRect(flat, grad, width, height, initial.data(), mBase1.data(), mBase2.data(), coordsOut.data(), iters.data(), 0.00001, 60);	
	fke.RunAssociateCoords(width*height, flat.idCount(), roots.size()/idC, coordsOut.data(), roots.data(), iters.data(), ids.data());
	fke.RunFindNewRoots(width, ids.data(), indexes.data(), height, height);
	CollectNewRoots(idC, indexes.size(), roots, coordsOut.data(), indexes.data());
	fke.RunColorTex(width*height, iters.data(), tex);
}
static void FractalGenImageChangingRoots(FlatSymExp& flat, FlatSymExp& grad, unsigned char* tex, FlatKernelExecutor& fke, int width, int height, std::vector<float>& center, std::vector<float>& base1, std::vector<float>& base2, float scale, std::vector<float>& roots)
{
	int idC = flat.idCount();

	coordsOut.resize(width*height*idC);
	iters.resize(width*height);
	ids.resize(width*height);
	indexes.resize(width*height);
	roots.clear();

	std::vector<float> initial; initial.resize(idC);
	std::vector<float> mBase1; mBase1.resize(idC);
	std::vector<float> mBase2; mBase2.resize(idC);
	for (int i = 0; i < idC; i++)
	{
		initial[i] = center[i] - scale*((base1[i]*float(width)/float(height)) + base2[i]);
		mBase1[i] = 2.0f*scale*base1[i]/float(height);
		mBase2[i] = 2.0f*scale*base2[i]/float(width);
	}

	fke.InitializeBuffers(flat.size, grad.size, width, height, idC, roots.size()/idC, width*height);
	fke.RunNMnTomRect(flat, grad, width, height, initial.data(), mBase1.data(), mBase2.data(), coordsOut.data(), iters.data(), 0.00001, 60);	
	FindAndCollectRootsCPU(2, width*height, roots, coordsOut.data(), iters.data(), ids.data());
	fke.RunColorTex(width*height, iters.data(), tex, ids.data());

	roots.clear();
}

static void SetTex(unsigned char* tex, int width, int height)
{
    sf::Image sfImg; sfImg.create(width, height, tex);
    sfTexPtr->loadFromImage(sfImg);
    drawRectanglePtr->setTexture(*sfTexPtr, true);
    drawRectanglePtr->setScale(width/drawRectanglePtr->getTextureRect().width, height/drawRectanglePtr->getTextureRect().height);
}