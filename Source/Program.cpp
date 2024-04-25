
#include "SymExp.hpp"
#include "Helpers.hpp"
#include "SFMLRenderer.hpp"

#include <iostream>
#include <vector>

//ToDo:
//string to SymExp + SymParser
//Find roots for associativity equations for of RAlg3
//Visualize root fractal (done)
//Function to convert real valued expression to complex (meh)
//Implement traditional method for finding complex roots of polynomial (for og Newtons fractal, which can be compared with full NM) (done, images are nearly identical,\
	gradient makes an orthogonal basis at the roots, so maybe that's why they are the same)
//Flatten SymExp: ProductRef, make SymExp some (one?) vectors for sum of products
//Port evaluation, gradient, and NM to OpenCL
//Figure out makefile (maybe use non-VS sfml libs?)
std::vector<SymExp> GenRule(std::vector<SymExp>& prodExps, int p, int i1, int i2, int i3)
{
	std::vector<SymExp> resL; for (int i = 0; i < p+1; i++) resL.push_back(SymExp());
	for (int i = 0; i < p+1; i++) //terms in ii
	{
		if (i == 0)
			resL[i3+1] += SymExp(prodExps[(i1*p)+i2].terms[i]);
		else
			for (int j = 0; j < p+1; j++) //components in result
				//           coeff of ii_i                  symbol of ii_i * i                          e.g. i = 2 gives j component of ii * i, then add each compoent to res
				resL[j] += SymExp(prodExps[(i1*p)+i2].terms[i])*prodExps[((i-1)*p)+i3].terms[j];

	}

	std::vector<SymExp> resR; for (int i = 0; i < p+1; i++) resR.push_back(SymExp());
	for (int i = 0; i < p+1; i++) //terms in ii
	{
		if (i == 0)
			resR[i1+1] += SymExp(prodExps[(i2*p)+i3].terms[i]);
		else
			for (int j = 0; j < p+1; j++) //components in result
				//           coeff of ii_i                  symbol of ii_i * i                          e.g. i = 2 gives j component of ii * i, then add each compoent to res
				resR[j] += SymExp(prodExps[(i2*p)+i3].terms[i])*prodExps[(i1*p)+(i-1)].terms[j];

	}

	std::vector<SymExp> rule;
	for (int i = 0; i < p+1; i++)
	{
		rule.push_back(resL[i]-resR[i]);
		rule[i].Simplify();
	}
	return rule;
}

SymExp GenRandomPoly(int varCount, int maxP, float range)
{
	SymExp hold = SymExp(randF(-range, range));
	if (varCount <= 0)
		return hold;
	std::vector<int> pows; for (int i = 0; i < varCount; i++) pows.push_back(0);
	for (int p = 1; p <= maxP; p++)
	{
		pows[0] = p;
		while (true)
		{
			Product prod(randF(-range, range));
			for (int i = 0; i < varCount; i++)
			{
				if (pows[i] == 0)
					continue;
				prod.MultId(i, pows[i]);
			}
			hold.terms.push_back(prod);

			int grabma = pows[varCount-1]+1;
			pows[varCount-1] = 0;
			for (int i = varCount-2; i >= 0; i--)
			{
				if (pows[i] != 0)
				{
					pows[i]--;
					pows[i+1] += grabma;
					grabma = -1;
					break;
				}
			}
			if (grabma != -1)
				break;
		}
	}
	return hold;
	
}

std::vector<float> NewtonFractalTest(const std::vector<float> initial, const float threshold = 0.00001);

int main()
{
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
	
	//find a way to represent iii
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

	/*int i1 = 0;
	int i2 = 0;
	int i3 = 0;

	std::cout << "\n(i" << i1 << "*i" << i2 << ")*i" << i3 << ":\n";
	std::vector<SymExp> resL; for (int i = 0; i < p+1; i++) resL.push_back(SymExp());
	for (int i = 0; i < p+1; i++) //terms in ii
	{
		if (i == 0)
			resL[i3+1] += SymExp(prodExps[(i1*p)+i2].terms[i]);
		else
			for (int j = 0; j < p+1; j++) //components in result
				//           coeff of ii_i                  symbol of ii_i * i                          e.g. i = 2 gives j component of ii * i, then add each compoent to res
				resL[j] += SymExp(prodExps[(i1*p)+i2].terms[i])*prodExps[((i-1)*p)+i3].terms[j];
		
	}
	std::cout << "s: " << resL[0].ToString(sp) << '\n';
	for (int i = 1; i < p+1; i++)
		std::cout << "i" << i << ": " << resL[i].ToString(sp) << '\n';
	
	std::cout << '\n';
	std::cout << "\ni" << i1 << "*(i" << i2 << "*i" << i3 << "):\n";
	std::vector<SymExp> resR; for (int i = 0; i < p+1; i++) resR.push_back(SymExp());
	for (int i = 0; i < p+1; i++) //terms in ii
	{
		if (i == 0)
			resR[i1+1] += SymExp(prodExps[(i2*p)+i3].terms[i]);
		else
			for (int j = 0; j < p+1; j++) //components in result
				//           coeff of ii_i                  symbol of ii_i * i                          e.g. i = 2 gives j component of ii * i, then add each compoent to res
				resR[j] += SymExp(prodExps[(i2*p)+i3].terms[i])*prodExps[(i1*p)+(i-1)].terms[j];

	}
	std::cout << "s: " << resR[0].ToString(sp) << '\n';
	for (int i = 1; i < p+1; i++)
		std::cout << "i" << i << ": " << resR[i].ToString(sp) << '\n';


	std::cout << '\n';
	std::vector<SymExp> rule;
	for (int i = 0; i < p+1; i++)
	{
		rule.push_back(resL[i]-resR[i]);
		rule[i].Simplify();
	}
	std::cout << "Associativity rule for i" << i1 << "*i" << i2 << "*i" << i3 << '\n';
	std::cout << "s: " << rule[0].ToString(sp) << " = 0\n";
	for (int i = 1; i < p+1; i++)
		std::cout << "i" << i << ": " << rule[i].ToString(sp) << " = 0\n";*/

	std::vector<SymExp> rules;
	for (int i1 = 0; i1 < p; i1++)
	{
		for (int i2 = 0; i2 < p; i2++)
		{
			for (int i3 = 0; i3 < p; i3++)
			{
				std::vector<SymExp> rule = GenRule(prodExps, p, i1, i2, i3);
				rules.insert(rules.end(),rule.begin(),rule.end());
			}
		}
	}
	for (int i = 0; i < p*p*p*(p+1); i++)
	{
		std::cout << rules[i].ToString(sp) << " = 0\n";
		if (i%(p+1) == p)
			std::cout << '\n';
	}

	SymExpTable eval;
	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j < p; j++)
		{
			for (int k = 0; k < p+1; k++)
			{
				eval.Add((i*p*(p+1))+(j*(p+1))+k, SymExp(0));
			}
		}
	}

	eval.Add(2,SymExp(1));
	eval.Add(3,SymExp(1));
	eval.Add(6,SymExp(1));
	eval.Add(10,SymExp(1));

	for (int i = 0; i < p*p*p*(p+1); i++)
	{
		std::cout << rules[i].SclEval(eval) << '\n';
		if (i%(p+1) == p)
			std::cout << '\n';
	}

	std::cout << '\n';
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

	SFMLRendererAPI sfmlapi;
	sfmlapi.SetTex(tex, 16, 16);

	randSeed(0);
	//err //pass reference for Gradient, reset vectors, change vector assignment to resize. 1:15 time, bubbles going left, lots of cyan bottom right, boot shape at center (like louisiana)
	std::vector<SymExp> rendExp; rendExp.push_back(SymExp(GenRandomPoly(2, 4, 4))); rendExp.push_back(GenRandomPoly(2, 4, 4));
	rendExp[0].scalar = 0; rendExp[1].scalar = 0;
	//rendExp.push_back(SymExp(-1)); rendExp.push_back(SymExp(0));
	//rendExp[0].terms.push_back(Product(1, 0, 3));	rendExp[0].terms.push_back(Product(-3, 0, 1)); rendExp[0].terms[1].MultId(1,2);
	//rendExp[1].terms.push_back(Product(-1, 1, 3)); rendExp[1].terms.push_back(Product(3, 0, 2)); rendExp[1].terms[1].MultId(1,1);
	std::vector<int> rendIds; rendIds.push_back(0); rendIds.push_back(1);
	std::vector<float> rendVal; rendVal.resize(2);
	std::vector<float> roots[2];
	int* rootAssoc = new int[800*800];
	unsigned char* rootTex = new unsigned char[800*800*4];

	for (int x = 0; x < 800; x++)
	{
		for (int y = 0; y < 800; y++)
		{
			rendVal[0] = ((x/800.0f)-0.5f)*8.0f;
			rendVal[1] = ((y/800.0f)-0.5f)*8.0f;
			rendVal = NMnTom(rendExp, rendIds, rendVal);
			//rendVal = NewtonFractalTest(rendVal);
			std::vector<float> eval = SclEvalVec(rendExp, rendIds, rendVal);
			int id = roots[0].size();
			for (int i = 0; i < roots[0].size(); i++)
			{
				if (std::abs(rendVal[0]-roots[0][i]) + std::abs(rendVal[1] - roots[1][i]) < 0.1)
				{
					id = i;
					break;
				}
			}
			if (std::abs(eval[0]) + std::abs(eval[1]) > 0.1)
				id = -1;
			if (id == roots[0].size())
			{
				roots[0].push_back(rendVal[0]);
				roots[1].push_back(rendVal[1]);
			}
			rootAssoc[(800*x)+y] = id;
		}
	}
	for (int x = 0; x < 800; x++)
	{
		for (int y = 0; y < 800; y++)
		{
			switch (rootAssoc[(800*x)+y])
			{
			default: rootTex[(800*4*x)+(y*4)+0] = 128; rootTex[(800*4*x)+(y*4)+1] = 128; rootTex[(800*4*x)+(y*4)+2] = 128; rootTex[(800*4*x)+(y*4)+3] = 255; break;
			case 0: rootTex[(800*4*x)+(y*4)+0] = 255; rootTex[(800*4*x)+(y*4)+1] = 0; rootTex[(800*4*x)+(y*4)+2] = 0; rootTex[(800*4*x)+(y*4)+3] = 255; break;
			case 1: rootTex[(800*4*x)+(y*4)+0] = 0; rootTex[(800*4*x)+(y*4)+1] = 255; rootTex[(800*4*x)+(y*4)+2] = 0; rootTex[(800*4*x)+(y*4)+3] = 255; break;
			case 2: rootTex[(800*4*x)+(y*4)+0] = 0; rootTex[(800*4*x)+(y*4)+1] = 0; rootTex[(800*4*x)+(y*4)+2] = 255; rootTex[(800*4*x)+(y*4)+3] = 255; break;
			case 3: rootTex[(800*4*x)+(y*4)+0] = 0; rootTex[(800*4*x)+(y*4)+1] = 255; rootTex[(800*4*x)+(y*4)+2] = 255; rootTex[(800*4*x)+(y*4)+3] = 255; break;
			case 4: rootTex[(800*4*x)+(y*4)+0] = 255; rootTex[(800*4*x)+(y*4)+1] = 0; rootTex[(800*4*x)+(y*4)+2] = 255; rootTex[(800*4*x)+(y*4)+3] = 255; break;
			case 5: rootTex[(800*4*x)+(y*4)+0] = 255; rootTex[(800*4*x)+(y*4)+1] = 255; rootTex[(800*4*x)+(y*4)+2] = 0; rootTex[(800*4*x)+(y*4)+3] = 255; break;
			}
		}
	}
	
	sfmlapi.SetTex(rootTex, 800, 800);

	while (true)
	{
		if (sfmlapi.RunLoop() == 1)
			break;

	}
	
	delete[] rootAssoc;
	delete[] rootTex;

	return 0;
}

std::vector<float> NewtonFractalTest(const std::vector<float> initial, const float threshold)
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