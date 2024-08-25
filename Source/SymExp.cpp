#include "SymExp.hpp"

#include <sstream>
#include <iomanip>
#include <Helpers.hpp>

#define float float

template<typename T>
inline T&& mov(T& t) { return reinterpret_cast<T&&>(t); }

static thread_local Vector<float> _testValue;
static thread_local Vector<float> _lastValue;
static thread_local Vector<float> _difVec;
static thread_local Vector<SymExp> _grad1D;
static thread_local Vector2D<SymExp> _grad2D;
static thread_local Vector<float> _gradEval1D;
static thread_local Vector2D<float> _gradEval2D;
static thread_local Vector2D<float> _vBasis;
static thread_local Vector2D<float> _wBasis;
static thread_local Vector<float> _wDivForV;
static thread_local Vector<float> _eval;
//static thread_local Vector2D<float> _hold;
//static thread_local Vector2D<float> _newMap;
static thread_local std::vector<Product> _prodGrad;
static thread_local std::vector<SymExp> _expGrad;


int BinSearch(std::vector<int>& vec, int val) //value at index is greater than or equal to val, can return index of vec.size()
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
int BinSearch(std::vector<int>& vec, int lower, int upper, int val)
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
    pows.resize(pows.size()-dif);
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

SymExp Product::operator+(Product& prod)
{
    SymExp hold;
    hold.terms.push_back(*this);
    hold.terms.push_back(prod);
    return hold;
}

SymExp Product::operator+(Product&& prod)
{
    return operator+(prod);
}

SymExp Product::operator-(Product& prod)
{
    SymExp hold;
    hold.terms.push_back(*this);
    hold.terms.push_back(prod);
    hold.terms.back().coeff *= -1;
    return hold;
}

SymExp Product::operator-(Product&& prod)
{
    return operator-(prod);
}

Product Product::operator*(Product& prod)
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

SymExp Product::operator*(Product&& prod)
{
    return operator*(prod);
}

/*Product Product::operator/(Product& prod)
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
}*/

bool Product::CheckType(Product& prod)
{
    if (ids.size() != prod.ids.size())
        return false;
    for (int i = 0; i < ids.size(); i++)
        if (ids[i] != prod.ids[i] || pows[i] != prod.pows[i])
            return false;
    return true;
}

SymExp Product::Eval(SymExpTable& table)
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

float Product::SclEval(SymExpTable& table)
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

SymExp Product::Eval(std::vector<int> valIds, std::vector<float> values)
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

float Product::SclEval(std::vector<int> valIds, std::vector<float> values)
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


std::string Product::ToString() { SymParser p; return ToString(p, true); }

std::string Product::ToString(SymParser& parser, bool genDef)
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

std::vector<int> Product::GetIds()
{
    return std::vector<int>(ids);
}

std::vector<Product> Product::Gradient() { std::vector<int> ids = GetIds(); return Gradient(ids); }
void Product::Gradient(std::vector<Product>& grad) { std::vector<int> ids = GetIds(); Gradient(ids, grad); }

std::vector<Product> Product::Gradient(std::vector<int>& gradIds)
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
void Product::Gradient(std::vector<int>& gradIds, std::vector<Product>& grad)
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
SymExp SymExp::operator+(SymExp& symExp)
{
    std::vector<Product> termsP;
    termsP.insert(termsP.end(),terms.begin(),terms.end());
    termsP.insert(termsP.end(),symExp.terms.begin(),symExp.terms.end());
    return SymExp(scalar+symExp.scalar, termsP);
}

SymExp SymExp::operator-(SymExp& symExp)
{
    std::vector<Product> termsP;
    termsP.insert(termsP.end(),terms.begin(),terms.end());
    int ind = termsP.size();
    termsP.insert(termsP.end(),symExp.terms.begin(),symExp.terms.end());
    for (; ind < termsP.size(); ind++)
        termsP[ind].coeff *= -1;
    return SymExp(scalar+symExp.scalar, termsP);
}

SymExp SymExp::operator*(SymExp& symExp)
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

/*SymExp SymExp::operator/(Product& prod)
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
}*/

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

std::vector<int> SymExp::GetIds()
{
    std::vector<int> ids;
    for (int i = 0; i < terms.size(); i++)
    {
        std::vector<int>& termIds = terms[i].ids;
        for (int j = 0; j < termIds.size(); j++)
        {
            int ind = BinSearch(ids, termIds[j]);
            if (ids.size() == ind || ids[ind] != termIds[j])
                ids.insert(ids.begin()+ind, termIds[j]);
        }
    }
    return mov(ids);
}

std::string SymExp::ToString() { SymParser p; return ToString(p, true); }

std::string SymExp::ToString(SymParser& parser, bool genDef)
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

SymExp SymExp::Eval(SymExpTable& table)
{
    SymExp exp(scalar);
    for (int i = 0; i < terms.size(); i++)
        exp += terms[i].Eval(table);

    exp.Simplify();
    return exp;
}

float SymExp::SclEval(SymExpTable& table)
{
    float sum(scalar);
    for (int i = 0; i < terms.size(); i++)
        sum += terms[i].SclEval(table);

    return sum;
}

SymExp SymExp::Eval(std::vector<int> valIds, std::vector<float> values)
{
    SymExp exp(scalar);
    for (int i = 0; i < terms.size(); i++)
        exp += terms[i].Eval(valIds, values);

    exp.Simplify();
    return exp;
}

float SymExp::SclEval(std::vector<int> valIds, std::vector<float> values)
{
    float sum(scalar);
    for (int i = 0; i < terms.size(); i++)
        sum += terms[i].SclEval(valIds, values);

    return sum;
}

std::vector<SymExp> SymExp::Gradient() { std::vector<int> ids = GetIds(); return Gradient(ids); }
void SymExp::Gradient(std::vector<SymExp>& grad) { std::vector<int> ids = GetIds(); Gradient(ids, grad); }

std::vector<SymExp> SymExp::Gradient(std::vector<int>& ids)
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
void SymExp::Gradient(std::vector<int>& ids, std::vector<SymExp>& grad)
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

std::vector<int> GetIds(std::vector<SymExp>& exps)
{
    std::vector<int> ids;
    for (int i = 0; i < exps.size(); i++)
    {
        std::vector<Product>& terms = exps[i].terms;
        for (int i = 0; i < terms.size(); i++)
        {
            std::vector<int>& termIds = terms[i].ids;
            for (int j = 0; j < termIds.size(); j++)
            {
                int ind = BinSearch(ids, termIds[j]);
                if (ids.size() == ind || ids[ind] != termIds[j])
                    ids.insert(ids.begin()+ind, termIds[j]);
            }
        }
    }
    return mov(ids);
}

Vector2D<SymExp> Gradient(std::vector<SymExp>& exps){ std::vector<int> ids = GetIds(exps); return Gradient(exps, ids);}
void Gradient(std::vector<SymExp>& exps, Vector2D<SymExp>& grad){ std::vector<int> ids = GetIds(exps); Gradient(exps, ids, grad);}

Vector2D<SymExp> Gradient(std::vector<SymExp>& exps, std::vector<int>& gradIds) //vector of vecs, first component is the partial derivative of all exps to the first id
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
void Gradient(std::vector<SymExp>& exps, std::vector<int>& gradIds, Vector2D<SymExp>& grad) //vector of vecs, first component is the partial derivative of all exps to the first id
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

std::vector<float> SclEvalVec(std::vector<SymExp>& exps, SymExpTable& table)
{
    std::vector<float> hold; hold.resize(exps.size());
    for (int i = 0; i < exps.size(); i++)
        hold[i] = exps[i].SclEval(table);
    return mov(hold);
}
std::vector<float> SclEvalVec(std::vector<SymExp>& exps, std::vector<int>& valIds, std::vector<float>& values)
{
    std::vector<float> hold; hold.resize(exps.size());
    for (int i = 0; i < exps.size(); i++)
        hold[i] = exps[i].SclEval(valIds, values);
    return mov(hold);
}
Vector2D<float> SclEvalVec2D(Vector2D<SymExp>& exps, SymExpTable& table)
{
    Vector2D<float> hold(exps.width, exps.height);
    for (int i = 0; i < hold.size(); i++)
        hold[i] = exps[i].SclEval(table);
    return mov(hold);
}
Vector2D<float> SclEvalVec2D(Vector2D<SymExp>& exps, std::vector<int>& valIds, std::vector<float>& values)
{
    Vector2D<float> hold(exps.width,exps.height);;
    for (int i = 0; i < hold.size(); i++)
        hold[i] = exps[i].SclEval(valIds, values);
    return mov(hold);
}

std::vector<float> NMnTo1(SymExp& poly, std::vector<int> valIds, std::vector<float> initial, float threshold)
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
        MulEq(gradEval, mult);
        SubEq(testValue, gradEval); //creates new vec here :P
        //testValue.R() -= gradEval.R()*mult;

        Sub(testValue, lastValue, difVec);
        //difVec = testValue.R()-lastValue;
        dif = SumAbs(difVec);

        count++;
    }

    return testValue;
}

std::vector<float> NMnTom(std::vector<SymExp>& polys, std::vector<int> valIds, std::vector<float> initial, float threshold)
{
    Vector<float>& testValue = _testValue; testValue = initial;
    Vector<float>& lastValue = _lastValue; lastValue = initial;
    Vector<float>& difVec = _difVec; difVec = initial;
    float dif = threshold+1;

    Vector2D<SymExp>& grad = _grad2D; grad = Gradient(polys, valIds);
    
    int count = 0;
    float evalSum = 2;
    Vector2D<float>& gradEval = _gradEval2D; gradEval.SetDim(initial.size(), polys.size());
    Vector2D<float>& vBasis = _vBasis; vBasis.SetDim(initial.size(),initial.size());
    Vector2D<float>& wBasis = _wBasis; wBasis.SetDim(initial.size(),polys.size());
    Vector<float>& wDivForV = _wDivForV; wDivForV.resize(initial.size());
    Vector<float>& eval = _eval; eval.resize(polys.size());
    while ((dif > threshold || evalSum > 0.0001f) && dif < 250 && count < 80)
    {
        eval = SclEvalVec(polys, valIds, testValue);
        gradEval = SclEvalVec2D(grad,valIds,testValue);

        evalSum = 0;
        for (int i = 0; i < polys.size(); i++)
            evalSum += glm::abs(eval[i]);

        //prevents floating point error when some ws are very close to 0, while other's aren't
        float maxWDiv = -1; //will be set on first w, a bad first w is possible

        lastValue = testValue;

        //transform gradEval into new basis with ortho wi
        int vCount = 0;
        for (int i = 0; i < polys.size() && vCount < initial.size(); i++) //i is the grad row used, vCount is the number of vis currently found
        {
            //calculate starting vi and wi
            VectorRef<float> initialV = vBasis.GetCol(vCount);
            for (int j = 0; j < initial.size(); j++)
                initialV[j] = gradEval.At(j,i);
            VectorRef<float> initialW = wBasis.GetCol(vCount);
            for (int i = 0; i < polys.size(); i++)
                initialW[i] = 0;

            //calculate initialW from vs and ws
            for (int i = 0; i < initial.size(); i++)
            {
                VectorRef<float> wTemp = gradEval.GetCol(i);
                for (int j = 0; j < polys.size(); j++)
                    initialW[j] += wTemp[j]*initialV[i];
            }

            //antiproject wi and vi
            for (int i = vCount-1; i >= 0; i--)
            {
                VectorRef<float> antiV = vBasis.GetCol(i);
                VectorRef<float> antiW = wBasis.GetCol(i);
                float coeff = Dot(initialW, antiW) / wDivForV[i];
                //project wi
                for (int i = 0; i < polys.size(); i++)
                    initialW[i] -= antiW[i]*coeff;
                //project vi
                for (int i = 0; i < initial.size(); i++)
                    initialV[i] -= antiV[i]*coeff;
            }
            
            wDivForV[vCount] = Dot(initialW, initialW);

            if (wDivForV[vCount] > maxWDiv)
                maxWDiv = wDivForV[vCount];

            //if wDivForV = 0, then vi and wi would be 0, so this vector should be skipped 
            if (wDivForV[vCount] > 0.000001f*maxWDiv || wDivForV[vCount] < -0.000001f*maxWDiv) //maybe not ferr, but approaching a 0 derivative causes problems
                vCount++;
        }

        //apply change
        for (int i = 0; i < vCount; i++)
        {
            VectorRef<float> wi = wBasis.GetCol(i);
            VectorRef<float> vi = vBasis.GetCol(i);
            float coeff = Dot(eval, wi)/wDivForV[i];
            for (int i = 0; i < initial.size(); i++)
                testValue[i] -= coeff*vi[i];
        }
        

        Sub(testValue, lastValue, difVec);
        //difVec = testValue-lastValue;
        dif = SumAbs(difVec);

        count++;
    }
    return testValue;
}

//original, incorrect implementation. Might be identical to complex NM on applicable functions (confirmed for z^3 - 1 = 0)
/*std::vector<float> NMnTom(const std::vector<SymExp>& polys, const std::vector<int> valIds, const std::vector<float> initial, const float threshold)
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
    Vector<float>& change = _change; change = initial;
    while (dif > threshold && dif < 250 && count < 80)
    {
        eval = SclEvalVec(polys, valIds, testValue);
        gradEval = SclEvalVec2D(grad,valIds,testValue);
        orthoVBasis = CoeffFromGramSchmidt(gradEval,gradEvalOrthoBasis);

        lastValue = testValue;
        
        err //algor is wrong, can't handle x^2 + y^2 - 1 = 0 , 0 = 0
            //smth to do with <P(x)*proj(delP)_vi / (proj(delP)_vi)^2>, specifically proj part 
        for (int i = 0; i < initial.size(); i++) //calculate change in orthobasis
        {
            change[i] = Dot(eval, gradEvalOrthoBasis.GetCol(i))/Dot(gradEvalOrthoBasis.GetCol(i), gradEvalOrthoBasis.GetCol(i));
            //if (isnan(change[i])) nans could happen in CoeffFromGramSchimdt
            //    change[i] = 0;
        }

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

    return testValue;
}*/

//name technically describes how function works, but not it's intent lol
/*
Vector2D<float> CoeffFromGramSchmidt(const Vector2D<float>& vecMap, Vector2D<float>& outNewMap) { return CoeffFromGramSchmidt(vecMap, &outNewMap); }

Vector2D<float> CoeffFromGramSchmidt(const Vector2D<float>& vecMap, Vector2D<float>* outNewMap)
{
    Vector2D<float>& hold = _hold; hold.SetDim(vecMap.width, vecMap.width); //number of coords equations, max equation involves all coords
    Vector2D<float>& newMap = _newMap; newMap = vecMap;
    for (int i = 0; i < hold.width; i++) //get equation for each coord (equivalently, set vecs in new map such that all are ortho to column i)
    {
        hold.At(i,i) = 1; //when going from x' to x, xi = 1x'i + m1x(i+1) + ...
        VectorRef<float> vi = newMap.GetCol(i);
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
*/

void SymExpResetGlobals()
{
    _testValue = Vector<float>(); _testValue.shrink_to_fit();
    _lastValue = Vector<float>(); _lastValue.shrink_to_fit();
    _difVec = Vector<float>(); _difVec.shrink_to_fit();
    _grad1D = Vector<SymExp>(); _grad1D.shrink_to_fit();
    _grad2D = Vector<SymExp>(); _grad2D.shrink_to_fit();
    _gradEval1D = Vector<float>(); _gradEval1D.shrink_to_fit();
    _gradEval2D = Vector2D<float>(); _gradEval2D.shrink_to_fit();
    _vBasis = Vector<float>(); _vBasis.shrink_to_fit();
    _wBasis = Vector2D<float>(); _wBasis.shrink_to_fit();
    _wDivForV = Vector<float>(); _wDivForV.shrink_to_fit();
    _eval = Vector<float>(); _eval.shrink_to_fit();
    //_hold = Vector2D<float>(); _hold.shrink_to_fit();
    //_newMap = Vector2D<float>(); _newMap.shrink_to_fit();
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
void SymExpTable::Add(int key, SymExp exp)
{
    int ind = BinSearch(lookup, key);
    if (lookup.size() == ind)
    {
        lookup.push_back(key);
        exps.push_back(mov(exp));
        return;
    }
    if (lookup[ind] == key)
        exps[ind] = mov(exp);
    else
    {
        lookup.insert(lookup.begin()+ind, key);
        exps.insert(exps.begin()+ind, mov(exp));
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

void SymParser::GenerateDefaults(Product& prod)
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

    std::vector<int>& ids = prod.ids;
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

void SymParser::GenerateDefaults(SymExp& symExp)
{
    //could somewhat simplify by using GetIds first, then iterating over ids

    std::vector<Product>& ts = symExp.terms;
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
        std::vector<int>& ids = ts[i].ids;
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

std::string SymParser::SearchId(int id)
{
    int index = BinSearch(lookup, id);
    if (lookup[index] == id)
        return names[index];
    else
        return std::to_string(id);
}
int SymParser::SearchName(std::string& name)
{
    for (int i = 0; i < names.size(); i++)
        if (names[i] == name)
            return i;
    return -1;
}


std::string FloatToString(float& floatRef)
{
    int f = int(floatRef);
    int i = 0;
    while (f != 0)
    { i++; f /= 10; }
    std::stringstream hStream;
    hStream << std::setprecision(5+i) << floatRef;
    return hStream.str();
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

std::vector<SymExp> ComplexPoly(SymExp& exp)
{
    std::vector<int> ids = exp.GetIds();
    return ComplexPoly(exp, ids);
}

std::vector<SymExp> ComplexPoly(SymExp& exp, std::vector<int> ids)
{
    std::vector<SymExp> hold; hold.resize(2);
    std::vector<int> newIds; newIds.resize(ids.size());
    int newId = 0;
    int newIndex = 0;
    int iId = 0;
    for (int i = 0; i < ids.size() && newIndex < ids.size()+1;)
    {
        if (ids[i] > newId)
        {
            if (newIndex == ids.size())
                iId = newId;
            else
                newIds[newIndex] = newId;
            newIndex++;
        }
        else
            i++;
        newId++;
    }
    for (; newIndex < ids.size()+1; newIndex++)
    {
        if (newIndex == ids.size())
            iId = newId;
        else
            newIds[newIndex] = newId;
        newId++;
    }
    SymExpTable table;
    SymExp sub; sub.terms.push_back(Product(1, 0, 1)); sub.terms.push_back(Product(1, 0, 1)); sub.terms[1].MultId(iId, 1); //sub = x_r + x_i*i
    for (int i = 0; i < ids.size(); i++)
    {
        sub.terms[0].ids[0] = ids[i];
        sub.terms[1].ids[0] = newIds[i];
        table.Add(ids[i], sub);
    }
    hold[0] = exp.Eval(table);
    for (int i = hold[0].terms.size()-1; i >= 0; i--)
    {
        Product& term = hold[0].terms[i];
        int iIndex = -1;
        for (int i = 0; i < term.ids.size(); i++)
        {
            if (term.ids[i] == iId)
            {
                iIndex = i;
                break;
            }
        }
        if (iIndex == -1)
            continue;
        if (term.pows[iIndex] % 4 == 0)
        {
            term.ids.erase(term.ids.begin()+iIndex);
            term.pows.erase(term.pows.begin()+iIndex);
        }
        else if (term.pows[iIndex] % 4 == 2)
        {
            term.coeff *= -1;
            term.ids.erase(term.ids.begin()+iIndex);
            term.pows.erase(term.pows.begin()+iIndex);
        }
        else if (term.pows[iIndex] % 4 == 1)
        {
            term.ids.erase(term.ids.begin()+iIndex);
            term.pows.erase(term.pows.begin()+iIndex);
            hold[1].terms.push_back(mov(term));
            hold[0].terms.erase(hold[0].terms.begin()+i);
        }
        else
        {
            term.coeff *= -1;
            term.ids.erase(term.ids.begin()+iIndex);
            term.pows.erase(term.pows.begin()+iIndex);
            hold[1].terms.push_back(mov(term));
            hold[0].terms.erase(hold[0].terms.begin()+i);
        }
    }
    return hold;
}

std::vector<SymExp> ComplexPoly(SymExp& exp, int iId)
{
    std::vector<int> ids = exp.GetIds();
    int ind = BinSearch(ids, iId);
    if (ind != ids.size())
        ids.erase(ids.begin()+ind);
    return ComplexPoly(exp, ids, iId);
}

std::vector<SymExp> ComplexPoly(SymExp& exp, std::vector<int> ids, int iId)
{
    std::vector<SymExp> hold; hold.resize(2);
    std::vector<int> newIds; newIds.resize(ids.size());
    int newId = 0;
    int newIndex = 0;
    for (int i = 0; i < ids.size() && newIndex < ids.size();)
    {
        if (ids[i] > newId)
        {
            if (newId == iId)
            {
                newId++;
                continue;
            }
            newIds[newIndex] = newId;
            newIndex++;
        }
        else
            i++;
        newId++;
    }
    for (; newIndex < ids.size(); newIndex++)
    {
        if (newId == iId)
            newId++;
        
        newIds[newIndex] = newId;
        newId++;
    }
    SymExpTable table;
    SymExp sub; sub.terms.push_back(Product(1, 0, 1)); sub.terms.push_back(Product(1, 0, 1)); sub.terms[1].MultId(iId, 1); //sub = x_r + x_i*i
    for (int i = 0; i < ids.size(); i++)
    {
        sub.terms[0].ids[0] = ids[i];
        if (newIds[i] < iId)
        {
            sub.terms[1].ids[0] = newIds[i];
            sub.terms[1].ids[1] = iId;
        }
        else
        {
            sub.terms[1].ids[0] = iId;
            sub.terms[1].ids[1] = newIds[i];
        }
        table.Add(ids[i], sub);
    }
    hold[0] = exp.Eval(table);
    for (int i = hold[0].terms.size()-1; i >= 0; i--)
    {
        Product& term = hold[0].terms[i];
        int iIndex = -1;
        for (int i = 0; i < term.ids.size(); i++)
        {
            if (term.ids[i] == iId)
            {
                iIndex = i;
                break;
            }
        }
        if (iIndex == -1)
            continue;
        if (term.pows[iIndex] % 4 == 0)
        {
            term.ids.erase(term.ids.begin()+iIndex);
            term.pows.erase(term.pows.begin()+iIndex);
        }
        else if (term.pows[iIndex] % 4 == 2)
        {
            term.coeff *= -1;
            term.ids.erase(term.ids.begin()+iIndex);
            term.pows.erase(term.pows.begin()+iIndex);
        }
        else if (term.pows[iIndex] % 4 == 1)
        {
            term.ids.erase(term.ids.begin()+iIndex);
            term.pows.erase(term.pows.begin()+iIndex);
            hold[1].terms.push_back(mov(term));
            hold[0].terms.erase(hold[0].terms.begin()+i);
        }
        else
        {
            term.coeff *= -1;
            term.ids.erase(term.ids.begin()+iIndex);
            term.pows.erase(term.pows.begin()+iIndex);
            hold[1].terms.push_back(mov(term));
            hold[0].terms.erase(hold[0].terms.begin()+i);
        }
    }
    return hold;
}

#undef float