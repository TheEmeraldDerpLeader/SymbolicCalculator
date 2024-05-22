
void SclEval(global uchar* flatSym,  local float* coords, local float* out, int expC, int idC);

kernel void NMnTom(global uchar* flatSym, global uchar* gradient, 
    global float* inCoords, global float* outCoords,
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
    local float* gradEval = gradEvalBuf+(expC*idC*localId);
    local float* vBasis = vBasisBuf+(idC*idC*localId);
    local float* wBasis = wBasisBuf+(idC*expC*localId);
    local float* wDivForV = wDivForV+(idC*localId);

    for (int i = 0; i < idC; i++)
        testValue[i] = inCoords[(get_global_id(0)*idC)+i];

    //SclEval works
    SclEval(flatSym, testValue, eval, expC, idC);   

    for (int i = 0; i < idC; i++)
        outCoords[(get_global_id(0)*idC)+i] = eval[i];
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