
void SclEval(global uchar* flatSym,  local float* coords, local float* out, int expC, int idC);

kernel void NMnTom(global uchar* flatSym, global uchar* gradient, 
    global float* inCoords, global float* outCoords, global int* outIter,
    local float* testValueBuf, local float* lastValueBuf,
    local float* evalBuf, local float* gradEvalBuf,
    local float* vBasisBuf, local float* wBasisBuf, local float* wDivForVBuf,
    float threshold, int maxCount, int expC, int idC)
{
    //setup local vectors
    int localId = get_local_id(0);
    local float* testValue = testValueBuf+(idC*localId);
    local float* lastValue = lastValueBuf+(idC*localId);
    local float* eval = evalBuf+(expC*localId);
    local float* gradEval = gradEvalBuf+(idC*expC*localId);
    local float* vBasis = vBasisBuf+(idC*idC*localId);
    local float* wBasis = wBasisBuf+(expC*expC*localId);
    local float* wDivForV = wDivForVBuf+(idC*localId);

    for (int i = 0; i < idC; i++)
        testValue[i] = inCoords[(get_global_id(0)*idC)+i];

    float dif = threshold+1;
    float evalSum = 0;
    int count = 0;
    bool success = false;
    while (dif < 2500000 && count < maxCount)
    {
        SclEval(flatSym, testValue, eval, expC, idC);
        SclEval(gradient, testValue, gradEval, idC*expC, idC);

        evalSum = 0;
        for (int i = 0; i < expC; i++)
            evalSum += fabs(eval[i]);

        int normMult = 0; //prevents floating point overflow (for positive exponent, negative exponent is unlikely and can be caused from ferror anyway)
        for (int i = 0; i < idC; i++)
            for (int j = 0; j < expC; j++)
                normMult = max(normMult, ilogb(gradEval[(i*expC)+j]));

        for (int i = 0; i < expC; i++)
        {
            eval[i] = ldexp(eval[i],-normMult);
            for (int j = 0; j < idC; j++)
                gradEval[(j*expC)+i] = ldexp(gradEval[(j*expC)+i],-normMult);
        }

        for (int i = 0; i < idC; i++)
            lastValue[i] = testValue[i];

        //transform gradEval into new basis with ortho wi
        int vCount = 0;
        for (int i = 0; i < expC && vCount < idC; i++) //i is the grad row used, vCount is the number of vis currently found
        {
            //calculate starting vi and wi
            local float* initialV = vBasis+(vCount*idC);
            for (int j = 0; j < idC; j++)
                initialV[j] = gradEval[(j*expC)+i];
            local float* initialW = wBasis+(vCount*expC);
            for (int i = 0; i < expC; i++)
                initialW[i] = 0;

            //calculate initialW from vs and ws
            for (int i = 0; i < idC; i++)
            {
                local float* wTemp = gradEval+(i*expC);
                for (int j = 0; j < expC; j++)
                    initialW[j] += wTemp[j]*initialV[i];
            }

            //antiproject wi and vi
            for (int i = vCount-1; i >= 0; i--)
            {
                local float* antiV = vBasis+(i*idC);
                local float* antiW = wBasis+(i*expC);
                float coeff = 0;
                for (int i = 0; i < expC; i++)
                    coeff += initialW[i]*antiW[i];
                coeff /= wDivForV[i];
                //project wi
                for (int i = 0; i < expC; i++)
                    initialW[i] -= antiW[i]*coeff;
                //project vi
                for (int i = 0; i < idC; i++)
                    initialV[i] -= antiV[i]*coeff;
            }

            wDivForV[vCount] = 0;
            for (int i = 0; i < expC; i++)
                wDivForV[vCount] += initialW[i]*initialW[i];
            //if wDivForV = 0, then vi and wi would be 0, so this vector should be skipped 
            if (wDivForV[vCount] > 0 || wDivForV[vCount] < -0) //ferr is relevant
                vCount++;
        }

        //apply change
        for (int i = 0; i < vCount; i++)
        {
            local float* wi = wBasis+(i*expC);
            local float* vi = vBasis+(i*idC);
            float coeff = 0;
            for (int i = 0; i < expC; i++)
                coeff += eval[i]*wi[i];
            coeff /= wDivForV[i];
            for (int i = 0; i < idC; i++)
                testValue[i] -= coeff*vi[i];
        }

        dif = 0;
        for (int i = 0; i < idC; i++)
            dif += fabs(testValue[i]-lastValue[i]);
        
        if (dif <= threshold && evalSum <= 0.01)
        {
            success = true;
            break;
        }
        count++;
    }

    if (success == false)
        count = -1;

    outIter[get_global_id(0)] = count;
    for (int i = 0; i < idC; i++)
        outCoords[(get_global_id(0)*idC)+i] = testValue[i];
}

kernel void NMnTomRect(global uchar* flatSym, global uchar* gradient, 
    global float* initial, global float* base1, global float* base2, global float* outCoords, global int* outIter,
    local float* testValueBuf, local float* lastValueBuf,
    local float* evalBuf, local float* gradEvalBuf,
    local float* vBasisBuf, local float* wBasisBuf, local float* wDivForVBuf,
    float threshold, int maxCount, int expC, int idC, int width, int height)
{
    //setup local vectors
    int localId = get_local_id(0);
    local float* testValue = testValueBuf+(idC*localId);
    local float* lastValue = lastValueBuf+(idC*localId);
    local float* eval = evalBuf+(expC*localId);
    local float* gradEval = gradEvalBuf+(idC*expC*localId);
    local float* vBasis = vBasisBuf+(idC*idC*localId);
    local float* wBasis = wBasisBuf+(expC*expC*localId);
    local float* wDivForV = wDivForVBuf+(idC*localId);

    for (int i = 0; i < idC; i++)
        testValue[i] = initial[i]+(base1[i]*(get_global_id(0)%width))+(base2[i]*(height-(get_global_id(0)/width)));

    float dif = threshold+1;
    float evalSum = 0;
    int count = 0;
    bool success = false;
    while (dif < 2500000 && count < maxCount)
    {
        SclEval(flatSym, testValue, eval, expC, idC);
        SclEval(gradient, testValue, gradEval, idC*expC, idC);

        int normMult = 0; //prevents floating point overflow (for positive exponent, negative exponent is unlikely and can be caused from ferror anyway)
        for (int i = 0; i < idC; i++)
            for (int j = 0; j < expC; j++)
                normMult = max(normMult, ilogb(gradEval[(i*expC)+j]));

        for (int i = 0; i < expC; i++)
        {
            eval[i] = ldexp(eval[i],-normMult);
            for (int j = 0; j < idC; j++)
                gradEval[(j*expC)+i] = ldexp(gradEval[(j*expC)+i],-normMult);
        }

        //prevents floating point error when some ws are very close to 0, while other's aren't
        float maxWDiv = -1; //will be set on first w, a bad first w is possible
        
        
        for (int i = 0; i < idC; i++)
            lastValue[i] = testValue[i];

        //transform gradEval into new basis with ortho wi
        int vCount = 0;
        for (int i = 0; i < expC && vCount < idC; i++) //i is the grad row used, vCount is the number of vis currently found
        {
            //calculate starting vi and wi
            local float* initialV = vBasis+(vCount*idC);
            for (int j = 0; j < idC; j++)
                initialV[j] = gradEval[(j*expC)+i];
            local float* initialW = wBasis+(vCount*expC);
            for (int i = 0; i < expC; i++)
                initialW[i] = 0;

            //calculate initialW from vs and ws
            for (int i = 0; i < idC; i++)
            {
                local float* wTemp = gradEval+(i*expC);
                for (int j = 0; j < expC; j++)
                    initialW[j] += wTemp[j]*initialV[i];
            }

            //antiproject wi and vi
            for (int i = vCount-1; i >= 0; i--)
            {
                local float* antiV = vBasis+(i*idC);
                local float* antiW = wBasis+(i*expC);
                float coeff = 0;
                for (int i = 0; i < expC; i++)
                    coeff += initialW[i]*antiW[i];
                coeff /= wDivForV[i];
                //project wi
                for (int i = 0; i < expC; i++)
                    initialW[i] -= antiW[i]*coeff;
                //project vi
                for (int i = 0; i < idC; i++)
                    initialV[i] -= antiV[i]*coeff;
            }

            wDivForV[vCount] = 0;
            for (int i = 0; i < expC; i++)
                wDivForV[vCount] += initialW[i]*initialW[i];

            if (wDivForV[vCount] > maxWDiv)
                maxWDiv = wDivForV[vCount];
            maxWDiv = 0; //idk sometimes helps, sometimes doesn't

            //if wDivForV = 0, then vi and wi would be 0, so this vector should be skipped 
            if (wDivForV[vCount] > 0.0000001*maxWDiv || wDivForV[vCount] < -0.0000001*maxWDiv) //ferr is relevant
                vCount++;
        }

        //apply change
        for (int i = 0; i < vCount; i++)
        {
            local float* wi = wBasis+(i*expC);
            local float* vi = vBasis+(i*idC);
            float coeff = 0;
            for (int i = 0; i < expC; i++)
                coeff += eval[i]*wi[i];
            coeff /= wDivForV[i];
            for (int i = 0; i < idC; i++)
                testValue[i] -= coeff*vi[i];
        }

        dif = 0;
        for (int i = 0; i < idC; i++)
            dif += fabs(testValue[i]-lastValue[i]);
        evalSum = 0;
        for (int i = 0; i < expC; i++)
            evalSum += fabs(eval[i]);

        if (dif <= threshold && evalSum <= 0.01)
        {
            success = true;
            break;
        }

        count++;
    }

    if (success == false)
        count = -1;

    outIter[get_global_id(0)] = count;
    for (int i = 0; i < idC; i++)
        outCoords[(get_global_id(0)*idC)+i] = testValue[i];
}

const int soi = sizeof(int);
const int sof = sizeof(float);

int iAt(global uchar* ptr)
{
    return *((int*)(ptr));
}
float fAt(global uchar* ptr)
{
    return *((float*)(ptr));
}

void SclEval(global uchar* flatSym,  local float* coords, local float* out, int expC, int idC)
{
    //skip id info and exp indexes
    global uchar* bufIndex = flatSym+ ( soi*(1+idC+1+expC) );


    for (int i = 0; i < expC; i++)
    {
        out[i] = fAt(bufIndex); bufIndex += sof;
        int prodCount = iAt(bufIndex); bufIndex += soi;

        for (int j = 0; j < prodCount; j++)
        {
            float prodVal = fAt(bufIndex); bufIndex += sof;
            int facCount = iAt(bufIndex); bufIndex += soi;

            for (int k = 0; k < facCount; k++)
            {
                float facVal = coords[iAt(bufIndex)]; bufIndex += soi;
                int pow = iAt(bufIndex); bufIndex += soi;
                for (int p = 0; p < pow; p++)
                    prodVal *= facVal;
            }

            out[i] += prodVal;
        }
    }
}  