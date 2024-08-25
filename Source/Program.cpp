
#include "SymExp.hpp"
#include "FlatSymExp.hpp"
#include "Helpers.hpp"
#include "SFMLRenderer.hpp"

#include "FlatSymExpOpenCL.hpp"

#include "Samples/Samples.hpp"

#include <iostream>
#include <vector>

std::vector<float> NewtonFractalTest(std::vector<float> initial, float threshold = 0.00001);

void FractalGenImage(FlatSymExp& flat, FlatSymExp& grad, int* rootAssoc, unsigned char* rootTex, FlatKernelExecutor& fke, int width, int height, float centerX, float centerY, float scale, std::vector<float>* roots);

int main()
{
	return SFMLMovement();
	//return CPPBindingTest();
	/*SymParser sp; sp.Add(0,"a1"); sp.Add(1,"a2"); sp.Add(2,"a3"); sp.Add(3,"a4"); sp.Add(4,"a5"); sp.Add(5,"b1"); sp.Add(6,"b2"); sp.Add(7,"b3"); sp.Add(8,"b4"); sp.Add(9,"b5");
	sp.Add(10,"c1"); sp.Add(11,"c2"); sp.Add(12,"c3"); sp.Add(13,"c4"); sp.Add(14,"c5"); sp.Add(15,"d1"); sp.Add(16,"d2"); sp.Add(17,"d3"); sp.Add(18,"d4"); sp.Add(19,"d5");
	sp.Add(20,"f1"); sp.Add(21,"f2"); sp.Add(22,"f3"); sp.Add(23,"f4"); sp.Add(24,"f5");
	
	SymExpTable tab;
	for (int i = 0; i < 25; i++)
		tab.Add(i, SymExp(0));

	Product ps[25];
	for (int i = 0; i < 25; i++)
	{
		ps[i].coeff = 1;
		ps[i].MultId(i);
	}
	SymExp ses[5];
	for (int j = 0; j < 5; j++)
	{
		for (int i = 0; i < 5; i++)
			ses[j].terms.push_back(ps[(j*5)+i]);
	}
	SymExp sum = ses[0]+ses[1]+ses[2]+ses[3]+ses[4];
	std::cout << (sum*sum).ToString(sp) << '\n';*/
	
	/*SymParser sp; sp.Add(0, "x"); sp.Add(1, "y"); sp.Add(2, "z");
	SymExp se; se.terms.push_back(Product(1, 0)); se.terms.push_back(Product(1, 1)); se.terms.push_back(Product(1, 2));
	std::cout << se.ToString(sp) << '\n';
	se = se*se;
	std::cout << se.ToString(sp) << '\n';
	SymExpTable set; set.Add(0, SymExp(2)); set.Add(1, SymExp(3));
	std::cout << se.Eval(set).ToString(sp) << '\n';
	set.Add(2, -1);
	std::cout << se.Eval(set).ToString(sp) << ' ' << se.SclEval(set) << '\n';*/
	
	SymExp test = GenRandomPoly(4, 6, 4);
	SymExpTable evalt; evalt.Add(0, 0); evalt.Add(1, 0); evalt.Add(2, 0); evalt.Add(3, 0);
	std::cout << test.ToString() << '\n';
	std::cout << test.SclEval(evalt) << ' ' << test.SclEval(SymExpTable()) << '\n';
	std::cout << "\n\n\n";

	SymParser sp;
	std::vector<SymExp> prodExps;
	int p = 2;
	int ind = 0;
	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j < p; j++)
		{
			prodExps.push_back(SymExp());
			for (int k = 0; k < p+1; k++)
			{
				sp.Add(ind, LettersFromNum((i*p)+j)+std::to_string(k+1));
				prodExps.back().terms.push_back(Product(1,ind));
				ind++;
			}
		}
	}
	for (int i = 0; i < p*p; i++)
		std::cout << prodExps[i].ToString(sp) << '\n';

	SymExp poly = SymExp(2); poly.terms.push_back(Product(-2,0,1)); poly.terms.push_back(Product(0.25,0,3));
	std::vector<int> ids; ids.push_back(0); ids.push_back(1);
	std::cout << "P(x):\n" << poly.ToString() << '\n';
	std::cout << "delP(x):\n" << poly.Gradient(ids)[0].ToString() << '\n';
	std::cout << "del2P(x):\n" << poly.Gradient(ids)[0].Gradient(ids)[0].ToString() << '\n';
	std::cout << "del3P(x):\n" << poly.Gradient(ids)[0].Gradient(ids)[0].Gradient(ids)[0].ToString() << '\n';

	poly = GenRandomPoly(2,2,3);
	std::cout << "P(x):\n" << poly.ToString() << '\n';
	std::cout << "delP(x):\n< " << poly.Gradient(ids)[0].ToString() << ", " << poly.Gradient(ids)[1].ToString() << ">\n";

	poly = SymExp(2); poly.terms.push_back(Product(-2,0,1)); poly.terms.push_back(Product(0.25,0,3));
	std::vector<float> initTest; initTest.push_back(2.5);
	ids = std::vector<int>(); ids.push_back(0);
	initTest = NMnTo1(poly,ids,initTest);
	for (int i = 0; i < initTest.size(); i++)
		std::cout << initTest[i] << ' ';
	std::cout << '\n' << poly.SclEval(ids,initTest) << '\n';

	poly = GenRandomPoly(3, 2, 4);
	initTest = std::vector<float>(); initTest.push_back(0); initTest.push_back(0); initTest.push_back(0);
	ids = std::vector<int>(); ids.push_back(0); ids.push_back(1); ids.push_back(2);
	initTest = NMnTo1(poly,ids,initTest);
	std::cout << "P(x):\n" << poly.ToString() << '\n';
	for (int i = 0; i < initTest.size(); i++)
		std::cout << initTest[i] << ' ';
	std::cout << '\n' << poly.SclEval(ids,initTest) << '\n';

	std::vector<SymExp> testPolys; testPolys.push_back(SymExp(6)); testPolys.push_back(SymExp(18)); testPolys.push_back(SymExp(-6));
	testPolys[0].terms.push_back(Product(6, 0)); testPolys[1].terms.push_back(Product(6, 0));
	testPolys[1].terms.push_back(Product(6, 1)); testPolys[2].terms.push_back(Product(6, 1));
	testPolys[2].terms.push_back(Product(6, 2)); testPolys[0].terms.push_back(Product(6, 2));
	initTest = std::vector<float>(); initTest.push_back(1); initTest.push_back(2); initTest.push_back(3);
	ids = std::vector<int>(); ids.push_back(0); ids.push_back(1); ids.push_back(2);
	initTest = NMnTom(testPolys,ids,initTest);
	std::cout << "P(x):\n" << testPolys[0].ToString() << '\n' << testPolys[1].ToString() << '\n' << testPolys[2].ToString() << "\n\n";
	for (int i = 0; i < initTest.size(); i++)
		std::cout << initTest[i] << ' ';
	std::vector<float> calcT = SclEvalVec(testPolys, ids, initTest);
	std::cout << '\n' << calcT[0] << ' ' << calcT[1] << ' ' << calcT[2] << '\n';

	std::vector<SymExp> polys; polys.push_back(GenRandomPoly(3,3,4)); polys.push_back(GenRandomPoly(3,3,4)); polys.push_back(GenRandomPoly(3,3,4));
	initTest = std::vector<float>(); initTest.push_back(0); initTest.push_back(0); initTest.push_back(0);
	ids = std::vector<int>(); ids.push_back(0); ids.push_back(1); ids.push_back(2);
	initTest = NMnTom(polys,ids,initTest);
	std::cout << "P(x):\n" << polys[0].ToString() << '\n' << polys[1].ToString() << '\n' << polys[2].ToString() << "\n\n";
	for (int i = 0; i < initTest.size(); i++)
		std::cout << initTest[i] << ' ';
	std::vector<float> calc = SclEvalVec(polys, ids, initTest);
	std::cout << '\n' << calc[0] << ' ' << calc[1] << ' ' << calc[2] << '\n';

	polys = std::vector<SymExp>(); polys.push_back(GenRandomPoly(1,3,4)); polys.push_back(GenRandomPoly(1,3,4)); polys.push_back(GenRandomPoly(1,3,4));
	SymExpTable coordChange; coordChange.Add(0, SymExp(1)); coordChange.exps[0].terms.push_back(Product(1, 0));
	polys[0].scalar = 0; polys[1].scalar = 0; polys[2].scalar = 0;
	polys[0] = polys[0].Eval(coordChange); polys[1] = polys[1].Eval(coordChange); polys[2] = polys[2].Eval(coordChange);
	initTest = std::vector<float>(); initTest.push_back(4);
	ids = std::vector<int>(); ids.push_back(0);
	initTest = NMnTom(polys,ids,initTest);
	std::cout << "P(x):\n" << polys[0].ToString() << '\n' << polys[1].ToString() << '\n' << polys[2].ToString() << "\n\n";
	for (int i = 0; i < initTest.size(); i++)
		std::cout << initTest[i] << ' ';
	calc = SclEvalVec(polys, ids, initTest);
	std::cout << '\n' << calc[0] << ' ' << calc[1] << ' ' << calc[2] << '\n';

	unsigned char tex[16*16*4];
	for (int x = 0; x < 16; x++)
	{
		for (int y = 0; y < 16; y++)
		{
			tex[(16*4*x)+(4*y)+0] = x*16;
			tex[(16*4*x)+(4*y)+1] = y*16;
			tex[(16*4*x)+(4*y)+2] = 0;
			tex[(16*4*x)+(4*y)+3] = 255;
		}
	}

	std::vector<SymExp> testM;

	testM.push_back(SymExp()); testM[0].scalar = 1;
	testM[0].terms.push_back(Product(2, 3, 4));
	testM.push_back(GenRandomPoly(3,3,3));
	/*testM.push_back(SymExp());
	testM[1].terms.push_back(Product(2,0,1)); testM[1].terms[0].MultId(1,1);
	testM[1].terms.push_back(Product(2,0,1)); testM[1].terms[1].MultId(3,2);
	testM[1].terms.push_back(Product(2,0,1)); testM[1].terms[2].MultId(1,1); testM[1].terms[2].MultId(3,1);*/
	FlatSymExp flatTest = FlatSymExp(testM);

	std::vector<int> flatIds; flatIds.push_back(0); flatIds.push_back(1); flatIds.push_back(2); flatIds.push_back(3);
	std::vector<float> flatCoords; flatCoords.push_back(1); flatCoords.push_back(2); flatCoords.push_back(3); flatCoords.push_back(4);

	std::vector<float> f1 = SclEvalVec(testM, flatIds, flatCoords);
	std::vector<float> f2 = flatTest.SclEval(flatIds, flatCoords);

	FlatSymExp flatTestGrad = flatTest.Gradient();

	Vector2D<SymExp> testMGrad = Gradient(testM);

	f1 = SclEvalVec(testMGrad, flatIds, flatCoords);
	f2 = flatTestGrad.SclEval(flatIds, flatCoords);

	SFMLRendererAPI sfmlapi;
	sfmlapi.SetTex(tex, 16, 16);

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
	//rendExp[0].scalar -= 0.25f;
	//rendExp[1].scalar -= 0.25f;

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

	//simple redundant
	//rendExp.push_back(SymExp(1)); rendExp.push_back(SymExp(0)); 
	//hold = SymExp(-1); hold.terms.push_back(Product(1, 0, 2)); hold.terms.push_back(Product(1, 1, 2)); rendExp[0] *= hold;
	//hold = SymExp(-4); hold.terms.push_back(Product(1, 0, 2)); hold.terms.push_back(Product(1.5f, 1, 2)); rendExp[0] *= hold;
	//hold = SymExp(-9); hold.terms.push_back(Product(1.5f, 0, 2)); hold.terms.push_back(Product(1, 1, 2)); rendExp[0] *= hold;
	//hold = SymExp(-32); hold.terms.push_back(Product(1, 0, 3)); hold.terms.push_back(Product(2.5f, 1, 2)); rendExp[0] *= hold;


	//rendExp.push_back(SymExp(-1)); rendExp[0].terms.push_back(Product(1,0)); rendExp[0].terms.push_back(Product(1,1));
	//rendExp.push_back(SymExp(0));

	//Random
	rendExp.push_back(SymExp(GenRandomPoly(2, 7, 4))); rendExp.push_back(GenRandomPoly(2, 7, 4));
	rendExp[0].scalar = 0; rendExp[1].scalar = 0;
	
	//NM
	//rendExp.push_back(SymExp(-1)); rendExp.push_back(SymExp(0));
	//rendExp[0].terms.push_back(Product(1, 0, 3));	rendExp[0].terms.push_back(Product(-3, 0, 1)); rendExp[0].terms[1].MultId(1,2);
	//rendExp[1].terms.push_back(Product(-1, 1, 3)); rendExp[1].terms.push_back(Product(3, 0, 2)); rendExp[1].terms[1].MultId(1,1);
	int* rootAssoc = new int[800*800];
	unsigned char* rootTex = new unsigned char[800*800*4];

	std::cout << rendExp[0].ToString() << '\n';
	std::cout << rendExp[1].ToString() << '\n';

	FlatSymExp rendFlat(rendExp);
	FlatSymExp rendGradFlat(rendFlat.Gradient());

	std::vector<cl::Platform> platforms;
	std::vector<std::vector<cl::Device>> devices;
	QueryDevices(platforms,devices);
	FlatKernelExecutor fke(devices[0][0]);

	std::vector<float> roots[2]; roots[0].clear(); roots[1].clear();

	float x = 0;

	while (true)
	{
		FractalGenImage(rendFlat, rendGradFlat, rootAssoc, rootTex, fke, 800, 800, x, 0, 4.0f, roots);
		x += 0.25f;
		sfmlapi.SetTex(rootTex, 800, 800);
		if (sfmlapi.RunLoop() == 1)
			break;

	}
	
	delete[] rootAssoc;
	delete[] rootTex;

	return 0;
}

void FractalGenImage(FlatSymExp& flat, FlatSymExp& grad, int* rootAssoc, unsigned char* rootTex, FlatKernelExecutor& fke, int width, int height, float centerX, float centerY, float scale, std::vector<float>* roots)
{
	thread_local std::vector<int> rendIds; rendIds.push_back(0); rendIds.push_back(1);
	thread_local std::vector<float> rendVal; rendVal.resize(2);

	bool resetRoots = (roots[0].size() == 0);

	thread_local std::vector<float> coordsIn;
	thread_local std::vector<float> coordsOut;
	thread_local std::vector<int> itersOut;
	coordsIn.resize(width*height*2);
	coordsOut.resize(width*height*2);
	itersOut.resize(width*height);

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			coordsIn[(x*height*2)+(y*2)+0] = (2.0f*((x/(float)width)-0.5f)*scale*((float)width/(float)height))+centerX;
			coordsIn[(x*height*2)+(y*2)+1] = (2.0f*((y/(float)height)-0.5f)*scale)+centerY;
		}
	}

	fke.InitializeBuffers(flat.size, grad.size, width, height, 2, 0, 0);

	fke.RunNMnTom(flat, grad, width*height, coordsIn.data(), coordsOut.data(), itersOut.data());

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			rendVal[0] = (2.0f*((x/(float)width)-0.5f)*scale*((float)width/(float)height))+centerX;
			rendVal[1] = (2.0f*((y/(float)height)-0.5f)*scale)+centerY;

			//rendVal = NMnTom(rendExp, rendIds, rendVal);
			///rendVal = NewtonFractalTest(rendVal);
			//Vector<float> eval = SclEvalVec(rendExp, rendIds, rendVal);

			//rendVal = flat.NewtonsMethodSolve(grad, rendVal);
			rendVal[0] = coordsOut[(x*800*2)+(y*2)+0];
			rendVal[1] = coordsOut[(x*800*2)+(y*2)+1];

			Vector<float> eval = flat.SclEval(rendVal);

			//color by final position
			//if (std::abs(eval[0]) + std::abs(eval[1]) > 0.1 || (std::isnan(eval[0]) || std::isnan(eval[1]) ))
			//	rootAssoc[(800*x)+y] = -1;
			/*
			if (rendVal[0]*rendVal[0] + rendVal[1]*rendVal[1] < 1.5f)
			rootAssoc[(800*x)+y] = 0;
			else if (rendVal[0]*rendVal[0] + rendVal[1]*rendVal[1] < 4.5f)
			rootAssoc[(800*x)+y] = 1;
			else if (rendVal[0]*rendVal[0] + rendVal[1]*rendVal[1] < 9.5f)
			rootAssoc[(800*x)+y] = 2;
			else
			rootAssoc[(800*x)+y] = 3;
			*/

			//color by root

			int id = roots[0].size();
			for (int i = 0; i < roots[0].size(); i++)
			{
				if (std::abs(rendVal[0]-roots[0][i]) + std::abs(rendVal[1] - roots[1][i]) < 0.1)
				{
					id = i;
					break;
				}
			}
			if (std::abs(eval[0]) + std::abs(eval[1]) > 0.1 || (std::isnan(eval[0]) || std::isnan(eval[1]) ))
				id = -1;
			if (id == roots[0].size())
			{
				if (resetRoots)
				{
					roots[0].push_back(rendVal[0]);
					roots[1].push_back(rendVal[1]);
				}
				else
					id == -1;
			}
			rootAssoc[(height*x)+y] = id;


			//color by eval (comment NM)
			/*float evalD = Dot(eval, eval); if (evalD < 0) evalD *= -1;
			rootAssoc[(800*x)+y] = 0;
			while (evalD > 0.001f)
			{
			rootAssoc[(800*x)+y]++;
			evalD /= 1.5f;
			}*/


		}
	}

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			switch (rootAssoc[(height*x)+(height-1-y)])
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
	}
}

std::vector<float> NewtonFractalTest(std::vector<float> initial, float threshold)
{
	std::vector<int> valIds; valIds.push_back(0); valIds.push_back(1);\

	std::vector<SymExp> polys;
	std::vector<SymExp> derPolys;
	
	polys.push_back(SymExp(-1)); polys.push_back(SymExp(0));
	polys[0].terms.push_back(Product(1, 0, 3));	polys[0].terms.push_back(Product(-3, 0, 1)); polys[0].terms[1].MultId(1,2);
	polys[1].terms.push_back(Product(-1, 1, 3)); polys[1].terms.push_back(Product(3, 0, 2)); polys[1].terms[1].MultId(1,1);
	
	derPolys.push_back(SymExp(0)); derPolys.push_back(SymExp(0));
	derPolys[0].terms.push_back(Product(3, 0, 2));	derPolys[0].terms.push_back(Product(-3, 1, 2));
	derPolys[1].terms.push_back(Product(6, 0, 1)); derPolys[1].terms[0].MultId(1,1);

	Vector<float> testValue = initial;
	Vector<float> lastValue = initial;
	Vector<float> difVec = initial;
	float dif = threshold+1;

	Vector2D<SymExp> grad = Gradient(polys,valIds);

	int count = 0;
	Vector2D<float> gradEval(grad.width, grad.height);
	Vector2D<float> orthoVBasis(initial.size(), initial.size());
	Vector2D<float> gradEvalOrthoBasis(grad.width, grad.height);
	Vector<float> eval; eval.resize(polys.size());
	Vector<float> derEval; derEval.resize(polys.size());
	Vector<float> change = initial; //should staticify these variables and add descriptions
	Vector<float> sum; sum.resize(polys.size());
	while (dif > threshold && dif < 250 && count < 80)
	{
		eval = SclEvalVec(polys, valIds, testValue);
		derEval = SclEvalVec(derPolys, valIds, testValue);

		float denom = (derEval[0]*derEval[0])+(derEval[1]*derEval[1]);
		if (denom == 0)
			break;
		derEval[0] /= denom;
		derEval[1] /= -denom;
		
		change[0] = (eval[0]*derEval[0])-(eval[1]*derEval[1]);
		change[1] = (eval[1]*derEval[0])+(eval[0]*derEval[1]);

		lastValue = testValue;

		SubEq(testValue, change);
		//testValue -= change;

		Sub(testValue, lastValue, difVec);
		//difVec = testValue-lastValue;
		dif = SumAbs(difVec);

		count++;
	}

	return reinterpret_cast<std::vector<float>&&>(testValue);
}