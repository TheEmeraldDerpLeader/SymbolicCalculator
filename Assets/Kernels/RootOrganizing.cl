
kernel void AssociateCoords( 
    global float* inCoords, global float* inRoots, global int* inIters, global int* outIds,
    local float* coordBuf,
    int idC, int rootC)
{
    if (inIters[get_global_id(0)] == -1)
    {
        outIds[get_global_id(0)] = -2;
        return;
    }

    //setup local vectors
    int localId = get_local_id(0);
    local float* coord = coordBuf+(idC*localId);

    for (int i = 0; i < idC; i++)
        coord[i] = inCoords[(get_global_id(0)*idC)+i];

    
    for (int i = 0; i < rootC; i++)
    {
        float dif = 0;
        for (int j = 0; j < idC; j++)
            dif += fabs(coord[j]-inRoots[(i*idC)+j]);
        if (dif <= 0.1)
        {
            outIds[get_global_id(0)] = i;
            return;
        } 
    }

    outIds[get_global_id(0)] = -1;
}

kernel void FindNewRoots( 
    global int* inIds, global int* outRootIndexes,
    int checkSpan, int offset)
{

    global int* initial = inIds+(get_global_id(0)*offset);
    for (int i = 0; i < checkSpan; i++)
    {
        if (initial[i] == -1)
        {
            outRootIndexes[get_global_id(0)] = (get_global_id(0)*offset)+i;
            return;
        }
    }
    outRootIndexes[get_global_id(0)] = -1;
    
}

kernel void ColorTex( 
    global int* inId, global int* inIters, global uchar* outTex)
{
    int glob = get_global_id(0);
    int id = inId[glob];

    if (id >= 18)
        id %= 18;

    if (inIters[glob] >= 0 && inIters[glob] <= 2)
    {
        outTex[(glob*4)+0] = 255; outTex[(glob*4)+1] = 255; outTex[(glob*4)+2] = 255; outTex[(glob*4)+3] = 255; return;
    }

    float iterDif = 1.0-(inIters[glob]/60.0);

    switch (id)
    {
        default: outTex[(glob*4)+0] = 0; outTex[(glob*4)+1] = 0; outTex[(glob*4)+2] = 0; outTex[(glob*4)+3] = 255; break;
        case -1: outTex[(glob*4)+0] = 180; outTex[(glob*4)+1] = 180; outTex[(glob*4)+2] = 180; outTex[(glob*4)+3] = 255; break;
        case 0: outTex[(glob*4)+0] = 255; outTex[(glob*4)+1] = 0; outTex[(glob*4)+2] = 0; outTex[(glob*4)+3] = 255; break;
        case 1: outTex[(glob*4)+0] = 0; outTex[(glob*4)+1] = 255; outTex[(glob*4)+2] = 0; outTex[(glob*4)+3] = 255; break;
        case 2: outTex[(glob*4)+0] = 0; outTex[(glob*4)+1] = 0; outTex[(glob*4)+2] = 255; outTex[(glob*4)+3] = 255; break;
        case 3: outTex[(glob*4)+0] = 0; outTex[(glob*4)+1] = 255; outTex[(glob*4)+2] = 255; outTex[(glob*4)+3] = 255; break;
        case 4: outTex[(glob*4)+0] = 255; outTex[(glob*4)+1] = 0; outTex[(glob*4)+2] = 255; outTex[(glob*4)+3] = 255; break;
        case 5: outTex[(glob*4)+0] = 255; outTex[(glob*4)+1] = 255; outTex[(glob*4)+2] = 0; outTex[(glob*4)+3] = 255; break;
        case 6: outTex[(glob*4)+0] = 128; outTex[(glob*4)+1] = 0; outTex[(glob*4)+2] = 0; outTex[(glob*4)+3] = 255; break;
        case 7: outTex[(glob*4)+0] = 0; outTex[(glob*4)+1] = 128; outTex[(glob*4)+2] = 0; outTex[(glob*4)+3] = 255; break;
        case 8: outTex[(glob*4)+0] = 0; outTex[(glob*4)+1] = 0; outTex[(glob*4)+2] = 128; outTex[(glob*4)+3] = 255; break;
        case 9: outTex[(glob*4)+0] = 0; outTex[(glob*4)+1] = 128; outTex[(glob*4)+2] = 128; outTex[(glob*4)+3] = 255; break;
        case 10: outTex[(glob*4)+0] = 128; outTex[(glob*4)+1] = 0; outTex[(glob*4)+2] = 128; outTex[(glob*4)+3] = 255; break;
        case 11: outTex[(glob*4)+0] = 128; outTex[(glob*4)+1] = 128; outTex[(glob*4)+2] = 0; outTex[(glob*4)+3] = 255; break;
        case 12: outTex[(glob*4)+0] = 180; outTex[(glob*4)+1] = 60; outTex[(glob*4)+2] = 60; outTex[(glob*4)+3] = 255; break;
        case 13: outTex[(glob*4)+0] = 60; outTex[(glob*4)+1] = 180; outTex[(glob*4)+2] = 60; outTex[(glob*4)+3] = 255; break;
        case 14: outTex[(glob*4)+0] = 60; outTex[(glob*4)+1] = 60; outTex[(glob*4)+2] = 180; outTex[(glob*4)+3] = 255; break;
        case 15: outTex[(glob*4)+0] = 60; outTex[(glob*4)+1] = 180; outTex[(glob*4)+2] = 180; outTex[(glob*4)+3] = 255; break;
        case 16: outTex[(glob*4)+0] = 180; outTex[(glob*4)+1] = 60; outTex[(glob*4)+2] = 180; outTex[(glob*4)+3] = 255; break;
        case 17: outTex[(glob*4)+0] = 180; outTex[(glob*4)+1] = 180; outTex[(glob*4)+2] = 60; outTex[(glob*4)+3] = 255; break;
    }
    outTex[(glob*4)+0] *= iterDif; outTex[(glob*4)+1] *= iterDif; outTex[(glob*4)+2] *= iterDif; 
}