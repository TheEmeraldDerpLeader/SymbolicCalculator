#include "SymExp.hpp"

#include <sstream>
#include <iomanip>

#define float float

template<typename T>
inline T&& mov(T& t) { return reinterpret_cast<T&&>(t); }

thread_local Vector<float> _testValue;
thread_local Vector<float> _lastValue;
thread_local Vector<float> _difVec;
thread_local Vector<SymExp> _grad1D;
thread_local Vector2D<SymExp> _grad2D;
thread_local Vector<float> _gradEval1D;
thread_local Vector2D<float> _gradEval2D;
thread_local Vector2D<float> _orthoVBasis;
thread_local Vector2D<float> _gradEvalOrthoBasis;
thread_local Vector<float> _eval;
thread_local Vector<float> _change;
thread_local Vector2D<float> _hold;
thread_local Vector2D<float> _newMap;
thread_local std::vector<Product> _prodGrad;
thread_local std::vector<SymExp> _expGrad;


int BinSearch(const std::vector<int>& vec, int val) //value at index is greater than or equal to val, can return index of vec.size()
{
    int upper = vec.size();
    int lower = 0;
    int ind = upper/2;
    while (upper != lower)
    {
        int v = vec[ind];
        if (val > v)
        {
            lower = ind+1;
            ind = (upper+lower)/2;
        }
        else if (val < v)
        {
            upper = ind;
            ind = (upper+lower)/2;
        }
        else
        {
            return ind;
        }
    }
    return ind;
}
int BinSearch(const std::vector<int>& vec, int lower, int upper, int val)
{
    int ind = upper/2;
    while (upper != lower)
    {
        int v = vec[ind];
        if (val > v)
        {
            lower = ind+1;
            ind = (upper+lower)/2;
        }
        else if (val < v)
        {
            upper = ind;
            ind = (upper+lower)/2;
        }
        else
        {
            return ind;
        }
    }
    return ind;
}

//Product ==============================================================
void Product::ClearEmpty()
{
    int dif = 0;
    for (int i = 0; i < ids.size()-dif; i++)
    {
        if (pows[i] == 0)
            dif++;
        else
        {
            ids[i-dif] = ids[i];
            pows[i-dif] = pows[i];
        }
    }
    ids.resize(ids.size()-dif);
}

Product& Product::MultId(int id, int pow)
{
    int ind = BinSearch(ids, id);
    if (ind == ids.size())
    {
        ids.push_back(id);
        pows.push_back(pow);
        return *this;
    }
    if (ids[ind] == id)
        pows[ind] += pow;
    else
    {
        ids.insert(ids.begin()+ind, id);
        pows.insert(pows.begin()+ind, pow);
    }
    
    return *this;
}

SymExp Product::operator+(const Product& prod)
{
    SymExp hold;
    hold.terms.push_back(*this);
    hold.terms.push_back(prod);
    return hold;
}

SymExp Product::operator-(const Product& prod)
{
    SymExp hold;
    hold.terms.push_back(*this);
    hold.terms.push_back(prod);
    hold.terms.back().coeff *= -1;
    return hold;
}

Product Product::operator*(const Product& prod)
{
    int i1 = 0;
    int i2 = 0;
    std::vector<int> idsT; std::vector<int> powsT;
    while (i1 < ids.size() && i2 < prod.ids.size())
    {
        int v1 = ids[i1]; int v2 = prod.ids[i2];
        if (v1 < v2)
        {
            idsT.push_back(v1);
            powsT.push_back(pows[i1]);
            i1++;
        }
        else if (v2 < v1)
        {
            idsT.push_back(v2);
            powsT.push_back(prod.pows[i2]);
            i2++;
        }
        else
        {
            int newPow = pows[i1]+prod.pows[i2];
            if (newPow < 0)
                continue;
            idsT.push_back(v1);
            powsT.push_back(newPow);
            i1++; i2++;
        }
    }
    while (i1 < ids.size())
    {
        idsT.push_back(ids[i1]);
        powsT.push_back(pows[i1]);
        i1++;
    }
    while (i2 < prod.ids.size())
    {
        idsT.push_back(prod.ids[i2]);
        powsT.push_back(prod.pows[i2]);
        i2++;
    }

    return Product(coeff*prod.coeff,mov(idsT),mov(powsT));
}

Product Product::operator/(const Product& prod)
{
    int i1 = 0;
    int i2 = 0;
    std::vector<int> idsT; std::vector<int> powsT;
    while (i1 < ids.size() && i2 < prod.ids.size())
    {
        int v1 = ids[i1]; int v2 = prod.ids[i2];
        if (v1 < v2)
        {
            idsT.push_back(v1);
            powsT.push_back(pows[i1]);
            i1++;
        }
        else if (v2 < v1)
        {
            idsT.push_back(v2);
            powsT.push_back(-prod.pows[i2]);
            i2++;
        }
        else
        {
            int newPow = pows[i1]-prod.pows[i2];
            if (newPow < 0)
                continue;
            idsT.push_back(v1);
            powsT.push_back(newPow);
            i1++; i2++;
        }
    }
    while (i1 < ids.size())
    {
        idsT.push_back(ids[i1]);
        powsT.push_back(pows[i1]);
        i1++;
    }
    while (i2 < prod.ids.size())
    {
        idsT.push_back(prod.ids[i2]);
        powsT.push_back(-prod.pows[i2]);
        i2++;
    }

    return Product(coeff*prod.coeff,mov(idsT),mov(powsT));
}

bool Product::CheckType(Product& prod) const
{
    if (ids.size() != prod.ids.size())
        return false;
    for (int i = 0; i < ids.size(); i++)
        if (ids[i] != prod.ids[i] || pows[i] != prod.pows[i])
            return false;
    return true;
}

SymExp Product::Eval(const SymExpTable& table) const
{
    SymExp exp(coeff);
    int j = 0;
    for (int i = 0; i < ids.size(); i++)
    {
        bool subbed = false;
        while (j < table.lookup.size())
        {
            if (table.lookup[j] == ids[i])
            {
                for (int k = 0; k < pows[i]; k++)
                    exp *= table.exps[j];
                subbed = true;
            }
            if (table.lookup[j] > ids[i])
                break;
            j++;
        }
        if (subbed == false)
        {
            SymExp t; t.terms.push_back(Product(1, ids[i], pows[i]));
            exp *= t;
        }
    }
    return exp;
}

float Product::SclEval(const SymExpTable& table) const
{
    float hold = coeff;
    int j = 0;
    for (int i = 0; i < ids.size(); i++)
    {
        while (j < table.lookup.size())
        {
            if (table.lookup[j] == ids[i])
            {
                for (int k = 0; k < pows[i]; k++)
                    hold *= table.exps[j].scalar;
            }
            if (table.lookup[j] > ids[i])
                break;
            j++;
        }
    }
    return hold;
}

SymExp Product::Eval(const std::vector<int> valIds, const std::vector<float> values) const
{
    SymExp exp; exp.terms.push_back(Product(coeff));
    Product& prod = exp.terms[0];
    int j = 0;
    for (int i = 0; i < ids.size(); i++)
    {
        bool subbed = false;
        while (j < valIds.size())
        {
            if (valIds[j] == ids[i])
            {
                for (int k = 0; k < pows[i]; k++)
                    prod.coeff *= values[j];
                subbed = true;
            }
            if (valIds[j] > ids[i])
                break;
            j++;
        }
        if (subbed == false)
        {
            prod.ids.push_back(ids[i]);
            prod.pows.push_back(pows[i]);
        }
    }
    return exp;
}

float Product::SclEval(const std::vector<int> valIds, const std::vector<float> values) const
{
    float hold = coeff;
    int j = 0;
    for (int i = 0; i < ids.size(); i++)
    {
        while (j < valIds.size())
        {
            if (valIds[j] == ids[i])
            {
                for (int k = 0; k < pows[i]; k++)
                    hold *= values[j];
            }
            if (valIds[j] > ids[i])
                break;
            j++;
        }
    }
    return hold;
}


std::string Product::ToString() const { SymParser p; return ToString(p, true); }

std::string Product::ToString(SymParser& parser, bool genDef) const
{
    if (genDef)
        parser.GenerateDefaults(*this);
    std::string hold;
    if (coeff != 1)
        hold += FloatToString(coeff) + "*";
    for (int i = 0; i < ids.size(); i++)
    {
        if (pows[i] == 1)
            hold += parser.SearchId(ids[i]);
        else if (pows[i] != 0)
            hold += parser.SearchId(ids[i]) + "^" + std::to_string(pows[i]);
        if (i < ids.size()-1)
            hold += "*";
    }
    return hold;
}

std::vector<int> Product::GetIds() const
{
    return std::vector<int>(ids);
}

std::vector<Product> Product::Gradient() const{ return Gradient(GetIds()); }
void Product::Gradient(std::vector<Product>& grad) const{ Gradient(GetIds(), grad); }

std::vector<Product> Product::Gradient(const std::vector<int>& gradIds) const
{
    std::vector<Product> grad; grad.resize(gradIds.size());
    for (int i = 0; i < gradIds.size(); i++)
    {
        int id = gradIds[i];
        int index = BinSearch(ids, id);
        if (index == ids.size() || ids[index] != id)
        {
            grad[i].coeff = 0;
            continue;
        }
        grad[i] = *this;
        grad[i].coeff *= grad[i].pows[index];
        grad[i].pows[index]--;
        if (grad[i].pows[index] == 0)
        {
            grad[i].ids.erase(grad[i].ids.begin()+index);
            grad[i].pows.erase(grad[i].pows.begin()+index);
        }
    }

    return mov(grad);
}
void Product::Gradient(const std::vector<int>& gradIds, std::vector<Product>& grad) const
{
    grad.resize(0);
    grad.resize(gradIds.size());
    for (int i = 0; i < gradIds.size(); i++)
    {
        int id = gradIds[i];
        int index = BinSearch(ids, id);
        if (index == ids.size() || ids[index] != id)
        {
            grad[i].coeff = 0;
            continue;
        }
        grad[i] = *this;
        grad[i].coeff *= grad[i].pows[index];
        grad[i].pows[index]--;
        if (grad[i].pows[index] == 0)
        {
            grad[i].ids.erase(grad[i].ids.begin()+index);
            grad[i].pows.erase(grad[i].pows.begin()+index);
        }
    }
}


//SymExp =====================================================
SymExp SymExp::operator+(const SymExp& symExp)
{
    std::vector<Product> termsP;
    termsP.insert(termsP.end(),terms.begin(),terms.end());
    termsP.insert(termsP.end(),symExp.terms.begin(),symExp.terms.end());
    return SymExp(scalar+symExp.scalar, termsP);
}

SymExp SymExp::operator-(const SymExp& symExp)
{
    std::vector<Product> termsP;
    termsP.insert(termsP.end(),terms.begin(),terms.end());
    int ind = termsP.size();
    termsP.insert(termsP.end(),symExp.terms.begin(),symExp.terms.end());
    for (; ind < termsP.size(); ind++)
        termsP[ind].coeff *= -1;
    return SymExp(scalar+symExp.scalar, termsP);
}

SymExp SymExp::operator*(const SymExp& symExp)
{
    SymExp hold;
    std::vector<Product>& termsP = hold.terms;
    for (int i = 0; i < terms.size(); i++)
        for (int j = 0; j < symExp.terms.size(); j++)
            termsP.push_back(terms[i]*symExp.terms[j]);
    if (symExp.scalar != 0)
        for (int i = 0; i < terms.size(); i++)
            termsP.push_back(terms[i]*symExp.scalar);
    if (scalar != 0)
        for (int j = 0; j < symExp.terms.size(); j++)
            termsP.push_back(Product(scalar)*symExp.terms[j]);
    
    hold.scalar = scalar*symExp.scalar;
    hold.Simplify();
    return mov(hold);
}

SymExp SymExp::operator/(const Product& prod)
{
    SymExp hold;
    Product p;
    p.coeff = scalar;
    p /= prod; 
    if (p.ids.size() == 0)
        hold.scalar += p.coeff;
    else
        hold.terms.push_back(p);

    for (int i = 0; i < terms.size(); i++)
    {
        p = terms[i]/prod;
        if (p.ids.size() == 0)
            hold.scalar += p.coeff;
        else
            hold.terms.push_back(p);
    }

    hold.Simplify();
    return mov(hold);
}

void SymExp::Simplify()
{
    for (int i = 0; i < terms.size(); i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (terms[i].CheckType(terms[j]))
            {
                terms[j].coeff += terms[i].coeff;
                terms[i].coeff = 0;
                if (terms[j].coeff == 0)
                {
                    terms.erase(terms.begin()+j);
                    i--;
                }
            }
        }
        if (terms[i].ids.size() == 0)
        {
            scalar += terms[i].coeff;
            terms[i].coeff = 0;
        }
        if (terms[i].coeff == 0)
        {
            terms.erase(terms.begin()+i);
            i--;
        }
    }
}

std::vector<int> SymExp::GetIds() const
{
    std::vector<int> ids;
    for (int i = 0; i < terms.size(); i++)
    {
        const std::vector<int>& termIds = terms[i].ids;
        for (int j = 0; j < termIds.size(); j++)
        {
            int ind = BinSearch(ids, termIds[j]);
            if (ids.size() == ind || ids[ind] != ids[j])
                ids.insert(ids.begin()+ind, termIds[j]);
        }
    }
    return mov(ids);
}

std::string SymExp::ToString() const { SymParser p; return ToString(p, true); }

std::string SymExp::ToString(SymParser& parser, bool genDef) const
{
    if (genDef)
        parser.GenerateDefaults(*this);
    std::string hold;
    if (terms.size() == 0)
        return FloatToString(scalar);
    if (scalar != 0)
        hold += FloatToString(scalar) + " + ";

    hold += terms[0].ToString(parser, false);
    for (int i = 1; i < terms.size(); i++)
        hold += " + " + terms[i].ToString(parser, false);
    return hold;
}

SymExp SymExp::Eval(const SymExpTable& table) const
{
    SymExp exp(scalar);
    for (int i = 0; i < terms.size(); i++)
        exp += terms[i].Eval(table);

    exp.Simplify();
    return exp;
}

float SymExp::SclEval(const SymExpTable& table) const
{
    float sum(scalar);
    for (int i = 0; i < terms.size(); i++)
        sum += terms[i].SclEval(table);

    return sum;
}

SymExp SymExp::Eval(const std::vector<int> valIds, const std::vector<float> values) const
{
    SymExp exp(scalar);
    for (int i = 0; i < terms.size(); i++)
        exp += terms[i].Eval(valIds, values);

    exp.Simplify();
    return exp;
}

float SymExp::SclEval(const std::vector<int> valIds, const std::vector<float> values) const
{
    float sum(scalar);
    for (int i = 0; i < terms.size(); i++)
        sum += terms[i].SclEval(valIds, values);

    return sum;
}

std::vector<SymExp> SymExp::Gradient() const { return Gradient(GetIds()); }
void SymExp::Gradient(std::vector<SymExp>& grad) const { Gradient(GetIds(), grad); }

std::vector<SymExp> SymExp::Gradient(const std::vector<int>& ids) const
{
    std::vector<SymExp> grad; grad.resize(ids.size());
    for (int i = 0; i < terms.size(); i++)
    {
        std::vector<Product>& prodGrad = _prodGrad; terms[i].Gradient(ids, prodGrad);
        for (int i = 0; i < ids.size(); i++)
        {
            if (prodGrad[i].ids.size() == 0)
                grad[i].scalar += prodGrad[i].coeff;
            else
                grad[i].terms.push_back(mov(prodGrad[i]));
        }
    }

    return mov(grad);
}
void SymExp::Gradient(const std::vector<int>& ids, std::vector<SymExp>& grad) const
{
    grad.clear();
    grad.resize(ids.size()); //gonna need to implement some memory system so SymExp vectors don't go out of scope on dealloc
    for (int i = 0; i < terms.size(); i++)
    {
        std::vector<Product>& prodGrad = _prodGrad; terms[i].Gradient(ids, prodGrad);
        for (int i = 0; i < ids.size(); i++)
        {
            if (prodGrad[i].ids.size() == 0)
                grad[i].scalar += prodGrad[i].coeff;
            else
                grad[i].terms.push_back(mov(prodGrad[i]));
        }
    }
}

std::vector<int> GetIds(const std::vector<SymExp>& exps)
{
    std::vector<int> ids;
    for (int i = 0; i < exps.size(); i++)
    {
        const std::vector<Product>& terms = exps[i].terms;
        for (int i = 0; i < terms.size(); i++)
        {
            const std::vector<int>& termIds = terms[i].ids;
            for (int j = 0; j < termIds.size(); j++)
            {
                int ind = BinSearch(ids, termIds[j]);
                if (ids.size() == ind || ids[ind] != ids[j])
                    ids.insert(ids.begin()+ind, termIds[j]);
            }
        }
    }
    return mov(ids);
}

Vector2D<SymExp> Gradient(const std::vector<SymExp>& exps){ return Gradient(exps, GetIds(exps));}
void Gradient(const std::vector<SymExp>& exps, Vector2D<SymExp>& grad){ Gradient(exps, GetIds(exps), grad);}

Vector2D<SymExp> Gradient(const std::vector<SymExp>& exps, const std::vector<int>& gradIds) //vector of vecs, first component is the gradient for each exp
{
    Vector2D<SymExp> grad(gradIds.size(), exps.size());
    for (int i = 0; i < exps.size(); i++)
    {
        std::vector<SymExp>& expGrad = _expGrad; exps[i].Gradient(gradIds, expGrad);
        for (int j = 0; j < gradIds.size(); j++)
            grad.At(j,i) = mov(expGrad[j]);
    }

    return mov(grad);
}
void Gradient(const std::vector<SymExp>& exps, const std::vector<int>& gradIds, Vector2D<SymExp>& grad) //vector of vecs, first component is the gradient for each exp
{
    grad.clear();
    grad.SetDim(gradIds.size(), exps.size());
    for (int i = 0; i < exps.size(); i++)
    {
        std::vector<SymExp>& expGrad = _expGrad; exps[i].Gradient(gradIds, expGrad);
        for (int j = 0; j < gradIds.size(); j++)
            grad.At(j,i) = mov(expGrad[j]);
    }
}

std::vector<float> SclEvalVec(const std::vector<SymExp>& exps, const SymExpTable& table)
{
    std::vector<float> hold; hold.resize(exps.size());
    for (int i = 0; i < exps.size(); i++)
        hold[i] = exps[i].SclEval(table);
    return mov(hold);
}
std::vector<float> SclEvalVec(const std::vector<SymExp>& exps, const std::vector<int>& valIds, const std::vector<float>& values)
{
    std::vector<float> hold; hold.resize(exps.size());
    for (int i = 0; i < exps.size(); i++)
        hold[i] = exps[i].SclEval(valIds, values);
    return mov(hold);
}
Vector2D<float> SclEvalVec2D(const Vector2D<SymExp>& exps, const SymExpTable& table)
{
    Vector2D<float> hold(exps.width, exps.height);
    for (int i = 0; i < hold.size(); i++)
        hold[i] = exps[i].SclEval(table);
    return mov(hold);
}
Vector2D<float> SclEvalVec2D(const Vector2D<SymExp>& exps, const std::vector<int>& valIds, const std::vector<float>& values)
{
    Vector2D<float> hold(exps.width,exps.height);;
    for (int i = 0; i < hold.size(); i++)
        hold[i] = exps[i].SclEval(valIds, values);
    return mov(hold);
}

std::vector<float> NMnTo1(const SymExp& poly, const std::vector<int> valIds, const std::vector<float> initial, const float threshold)
{
    Vector<float>& testValue = _testValue; testValue = initial;
    Vector<float>& lastValue = _lastValue; lastValue = initial;
    Vector<float>& difVec = _difVec; difVec = initial;
    float dif = threshold+1;

    Vector<SymExp>& grad = _grad1D; grad = poly.Gradient(valIds);

    int count = 0;
    Vector<float>& gradEval = _gradEval1D; gradEval = initial;
    while (dif > threshold && dif < 250 && count < 40)
    {
        float eval = poly.SclEval(valIds, testValue);
        gradEval = SclEvalVec(grad, valIds, testValue);

        lastValue = testValue;
        float mult = eval/Dot(gradEval,gradEval);//P(x)/delP*delP
        SubEq(testValue, Mul(gradEval,mult));
        //testValue.R() -= gradEval.R()*mult;

        Sub(testValue, lastValue, difVec);
        //difVec = testValue.R()-lastValue;
        dif = SumAbs(difVec);

        count++;
    }

    return mov(testValue);
}

std::vector<float> NMnTom(const std::vector<SymExp>& polys, const std::vector<int> valIds, const std::vector<float> initial, const float threshold)
{
    Vector<float>& testValue = _testValue; testValue = initial;
    Vector<float>& lastValue = _lastValue; lastValue = initial;
    Vector<float>& difVec = _difVec; difVec = initial;
    float dif = threshold+1;

    Vector2D<SymExp>& grad = _grad2D; grad = Gradient(polys, valIds);

    int count = 0;
    Vector2D<float>& gradEval = _gradEval2D; gradEval.SetDim(grad.width, grad.height);
    Vector2D<float>& orthoVBasis = _orthoVBasis; orthoVBasis.SetDim(initial.size(), initial.size());
    Vector2D<float>& gradEvalOrthoBasis = _gradEvalOrthoBasis; gradEvalOrthoBasis.SetDim(grad.width, grad.height);
    Vector<float>& eval = _eval; eval.resize(polys.size());
    Vector<float>& change = _change; change = initial; //should staticify these variables and add descriptions
    while (dif > threshold && dif < 250 && count < 80)
    {
        eval = SclEvalVec(polys, valIds, testValue);
        gradEval = SclEvalVec2D(grad,valIds,testValue);
        orthoVBasis = CoeffFromGramSchmidt(gradEval,gradEvalOrthoBasis);

        lastValue = testValue;
        
        for (int i = 0; i < initial.size(); i++) //calculate change in orthobasis
            //change[i] = eval.Dot(gradEvalOrthoBasis.GetCol(i))/gradEvalOrthoBasis.GetCol(i).Dot(gradEvalOrthoBasis.GetCol(i));
            change[i] = Dot(eval,gradEvalOrthoBasis.GetCol(i))/Dot(gradEvalOrthoBasis.GetCol(i),gradEvalOrthoBasis.GetCol(i));

        for (int i = initial.size()-1; i >= 0; i--) //transform change into original basis
            for (int j = i+1; j < initial.size(); j++)
                change[i] += change[j]*orthoVBasis.At(i,j);
          
        SubEq(testValue, change);
        //testValue -= change;

        Sub(testValue, lastValue, difVec);
        //difVec = testValue-lastValue;
        dif = SumAbs(difVec);

        count++;
    }

    return mov(testValue);
}

//name technically describes how function works, but not it's intent lol
Vector2D<float> CoeffFromGramSchmidt(const Vector2D<float>& vecMap, Vector2D<float>& outNewMap) { return CoeffFromGramSchmidt(vecMap, &outNewMap); }

Vector2D<float> CoeffFromGramSchmidt(const Vector2D<float>& vecMap, Vector2D<float>* outNewMap)
{
    Vector2D<float>& hold = _hold; hold.SetDim(vecMap.width, vecMap.width); //number of coords equations, max equation involves all coords
    Vector2D<float>& newMap = _newMap; newMap = vecMap;
    for (int i = 0; i < hold.width; i++) //get equation for each coord (equivalently, set vecs in new map such that all are ortho to column i)
    {
        hold.At(i,i) = 1; //when going from x' to x, xi = 1x'i + m1x(i+1) + ...
        VectorRef<float> vi = newMap.GetCol(i); //shouldn't copy vecs, should use ptr to location
        float div = Dot(vi,vi);
        for (int j = i+1; j < hold.height; j++) //apply Gram-schidt to all vectors in terms of vi
        { //technically gram schmidt is typically done by iterating on each vector then past vectors instead of each vector then future vectors, but it is identical, just moving computation around
            float projCoeff = -Dot(vi,newMap.GetCol(j))/div;
            hold.At(i,j) = projCoeff;
            for (int k = 0; k < vecMap.height; k++)
                newMap.At(j,k) += vi[k]*projCoeff; //update newMap with partially GramSchmidted vectors
        }
    }
    
    if (outNewMap != nullptr)
        *outNewMap = newMap;
    return hold;
}

void ResetGlobals()
{
    _testValue = Vector<float>(); _testValue.shrink_to_fit();
    _lastValue = Vector<float>(); _lastValue.shrink_to_fit();
    _difVec = Vector<float>(); _difVec.shrink_to_fit();
    _grad1D = Vector<SymExp>(); _grad1D.shrink_to_fit();
    _grad2D = Vector<SymExp>(); _grad2D.shrink_to_fit();
    _gradEval1D = Vector<float>(); _gradEval1D.shrink_to_fit();
    _gradEval2D = Vector2D<float>(); _gradEval2D.shrink_to_fit();
    _orthoVBasis = Vector<float>(); _orthoVBasis.shrink_to_fit();
    _gradEvalOrthoBasis = Vector<float>(); _gradEvalOrthoBasis.shrink_to_fit();
    _eval = Vector<float>(); _eval.shrink_to_fit();
    _change = Vector<float>(); _change.shrink_to_fit();
    _hold = Vector2D<float>(); _hold.shrink_to_fit();
    _newMap = Vector2D<float>(); _newMap.shrink_to_fit();
    _prodGrad = std::vector<Product>(); _prodGrad.shrink_to_fit();
    _expGrad = std::vector<SymExp>(); _expGrad.shrink_to_fit();
}

/*std::vector<SymExp> ComplexPoly(const SymExp exp)
{
    std::vector<int> ids = exp.GetIds();
    if (ids.size() == 0)
    {
        std::vector<SymExp> hold; hold.push_back(exp); hold.push_back(SymExp(0));
        return hold;
    }
    std::vector<int> cIds; cIds.resize(ids.size());
    int idsInd = 0;
    int i = 0;
    for (; idsInd < ids.size(); i++)
        if (i < ids[idsInd])
            cIds.push_back(i);
        else
            idsInd++;
    while (cIds.size() < ids.size())
    {
        cIds.push_back(i);
        i++;
    }
    return ComplexPoly(exp, ids, cIds);
}

std::vector<SymExp> ComplexPoly(const SymExp exp, SymParser& parse)
{
    std::vector<int> ids = exp.GetIds();
    if (ids.size() == 0)
    {
        std::vector<SymExp> hold; hold.push_back(exp); hold.push_back(SymExp(0));
        return hold;
    }
    std::vector<int> cIds; cIds.resize(ids.size());
    int idsInd = 0;
    int i = 0;
    for (; idsInd < ids.size(); i++)
        if (i < ids[idsInd])
            cIds.push_back(i);
        else
            idsInd++;
    while (cIds.size() < ids.size())
    {
        cIds.push_back(i);
        i++;
    }
    return ComplexPoly(exp, parse, ids, cIds);
}

std::vector<SymExp> ComplexPoly(const SymExp exp, SymParser& parse, const std::vector<int> ids, const std::vector<int> cIds)
{
    for (int i = 0; i < ids.size(); i++)
    {
        std::string name = parse.SearchId(ids[i]);
        parse.Add(cIds[i], name+"c");
    }

    return ComplexPoly(exp, ids, cIds);
}

std::vector<SymExp> ComplexPoly(const SymExp exp, const std::vector<int> ids, const std::vector<int> cIds)
{
    std::vector<SymExp> hold; hold.resize(2);
    int holdTermInd = 0;
    for (int i = 0; exp.terms.size(); i++)
    {
        const Product& term = exp.terms[i];
        for (int i = 0; i < ids.size(); i++)
        {
            int idInd = BinSearch(term.ids, ids[i]);
            if (idInd = term.ids.size())
            {
                hold[0].terms.push_back(term);
            }
            int idPow = term.pows[idInd];
        }
    }
    unfinished //error
    return hold;
} */


//SymExpTable ===================================================
void SymExpTable::Add(int key, const SymExp& exp)
{
    int ind = BinSearch(lookup, key);
    if (lookup.size() == ind)
    {
        lookup.push_back(key);
        exps.push_back(exp);
        return;
    }
    if (lookup[ind] == key)
        exps[ind] = exp;
    else
    {
        lookup.insert(lookup.begin()+ind, key);
        exps.insert(exps.begin()+ind, exp);
    }
}

//SymParser ====================================================
std::string LettersFromNum(int num)
{
    std::string hold;
    for (; true; num /= 26)
    {
        hold.push_back('a'+(num%26));
        if (num < 26)
            break;
    }
    return hold;
}

void SymParser::Add(int key, std::string name)
{
    int ind = BinSearch(lookup, key);
    if (lookup.size() == ind)
    {
        lookup.push_back(key);
        names.push_back(name);
        return;
    }
    if (lookup[ind] == key)
        names[ind] = name;
    else
    {
        lookup.insert(lookup.begin()+ind, key);
        names.insert(names.begin()+ind, name);
    }
}

void SymParser::GenerateDefaults(const Product& prod)
{
    int name = 0;
    std::string nameS = LettersFromNum(name);
    while (true)
    {
        if (SearchName(nameS) == -1)
            break;
        else
        {
            name++;
            nameS = LettersFromNum(name);
        }
    }

    const std::vector<int>& ids = prod.ids;
    for (int j = 0; j < ids.size(); j++)
    {
        int ind = BinSearch(lookup, ids[j]);
        if (lookup.size() == ind || lookup[ind] != ids[j])
        {
            lookup.insert(lookup.begin()+ind, ids[j]);
            names.insert(names.begin()+ind, nameS);
            name++; nameS = LettersFromNum(name);
            while (true)
            {
                if (SearchName(nameS) == -1)
                    break;
                else
                {
                    name++;
                    nameS = LettersFromNum(name);
                }
            }
        }

    }
}

void SymParser::GenerateDefaults(const SymExp& symExp)
{
    //could somewhat simplify by using GetIds first, then iterating over ids

    const std::vector<Product>& ts = symExp.terms;
    int name = 0;
    std::string nameS = LettersFromNum(name);
    //determine first name
    while (true)
    {
        if (SearchName(nameS) == -1)
            break;
        else
        {
            name++;
            nameS = LettersFromNum(name);
        }
    }
    //find ids
    for (int i = 0; i < ts.size(); i++)
    {
        const std::vector<int>& ids = ts[i].ids;
        for (int j = 0; j < ids.size(); j++)
        {
            int ind = BinSearch(lookup, ids[j]);
            if (lookup.size() == ind || lookup[ind] != ids[j])
            {
                lookup.insert(lookup.begin()+ind, ids[j]);
                names.insert(names.begin()+ind, nameS);
                name++; nameS = LettersFromNum(name);
                while (true)
                {
                    if (SearchName(nameS) == -1)
                        break;
                    else
                    {
                        name++;
                        nameS = LettersFromNum(name);
                    }
                }
            }
        }
    }
}

std::string SymParser::SearchId(int id) const
{
    int index = BinSearch(lookup, id);
    if (lookup[index] == id)
        return names[index];
    else
        return std::to_string(id);
}
int SymParser::SearchName(const std::string& name) const
{
    for (int i = 0; i < names.size(); i++)
        if (names[i] == name)
            return i;
    return -1;
}


std::string FloatToString(const float& floatRef)
{
    int f = int(floatRef);
    int i = 0;
    while (f != 0)
    { i++; f /= 10; }
    std::stringstream hStream;
    hStream << std::setprecision(5+i) << floatRef;
    return hStream.str();
}