#include <stdint.h>
#include <stdio.h>
#include "tree_utils.h"

void int_to_bin_digit(int64_t in, int count, int* out)
{
    /* assert: count <= sizeof(int)*CHAR_BIT */
    unsigned int mask = 1U << (count-1);
    int i;
    for (i = 0; i < count; i++) {
        out[i] = (in & mask) ? 1 : 0;
        in <<= 1;
    }
}

void printRealMatrix(double *orthoMatrix, int rows, int cols){

    for(int i = 0; i < rows; i++){
        printf("\n");
        for(int j = 0; j < cols; j++){
            printf("%10.2f ",orthoMatrix[i*cols + j]);
        }
    }
    printf("\n");
}

#include <stdio.h>
#include <stdint.h>
#include <math.h>

inline double logbinom(double n, double k) {
    return lgamma(n+1)-lgamma(n-k+1)-lgamma(k+1);
}
inline double binom(double n, double k) {
    return exp(logbinom(n,k));
}

void getncsfs1(int *inpnsomo, int *inpms, int *outncsfs){
    int nsomo = *inpnsomo;
    int ms = *inpms;
    int nparcoupl = (nsomo + ms)/2;
    *outncsfs = binom(nsomo, nparcoupl);
}

void getncsfs(int NSOMO, int MS, int *outncsfs){
    int nparcoupl = (NSOMO + MS)/2;
    int nparcouplp1 = ((NSOMO + MS)/2)+1;
    double tmpndets=0.0;
    tmpndets = binom(NSOMO, nparcoupl);
    (*outncsfs) = round(tmpndets - binom(NSOMO, nparcouplp1));
}

#include <stdio.h>
#include <stdint.h>

void printCFGList(int32_t *inpNint, int32_t *inpNcfgs, int64_t *cfglist){
    int Ncfgs    = *inpNcfgs;
    int N_int    = *inpNint;
    int digit[MAX_SOMO];
    int64_t cfg1,cfg2;
    int nsomo=4;
    int ms=0;
    int ncsfs=0;
    printf("In 64 printcfglist\n");
    printf("Ncfgs = %d Nint=%d\n",Ncfgs, N_int);
    printf(" 1-- %ld \n -- %ld \n",cfglist[0*(2*Ncfgs) + 0*(Ncfgs) + 0], cfglist[0*(2*Ncfgs) + 1*(Ncfgs) + 0]);
    for(int i = 0; i < 15; i++){
        cfg1 = cfglist[1 + i*2];
        cfg2 = cfglist[0 + i*2];
        printf("%d> domo=%ld somo=%ld\n",i,cfg1,cfg2);
        int_to_bin_digit(cfg2,18,digit);
        for(int j=0;j<18;j++)
            printf("%d ",digit[j]);
        printf("\n");
    }
    getncsfs1(&nsomo,&ms,&ncsfs);
    printf("Nsomos = %d\n",ncsfs);
}

#include <stdio.h>

void getBFIndexList(int NSOMO, int *BF1, int *IdxListBF1){
    int Iidx;
    int Jidx;
    int BFcopy[NSOMO];

    int dictidx[2];
    dictidx[0] = -1;
    dictidx[1] =  1;

    for(int i = 0; i < NSOMO; i++)
        BFcopy[i] = BF1[i];

    for(int i = 0; i < NSOMO; i++){
        Iidx = i;
        if(BFcopy[i] == 0){
            int countN1=0;
            for(int j = i+1; j < NSOMO; j++){
                Jidx = j;
                countN1 = countN1 + dictidx[BFcopy[j]];
                if(countN1 > 0){
                    break;
                }
            }
            BFcopy[Iidx] = -1;
            BFcopy[Jidx] = -1;
            IdxListBF1[Jidx] = Iidx;
            IdxListBF1[Iidx] = Jidx;
        }
    }

}

#include <stdio.h>

void getIslands(int NSOMO, int *BF1, int *BF2, int *nislands, int *phasefactor){

    // Get BF ids
    int *IdxListBF1 = malloc(NSOMO * sizeof(int));
    int *IdxListBF2 = malloc(NSOMO * sizeof(int));

    getBFIndexList(NSOMO, BF1, IdxListBF1);
    getBFIndexList(NSOMO, BF2, IdxListBF2);

    int maxcount=0;
    *nislands = 0;
    *phasefactor = 1;

    int BF1copy[NSOMO];
    for(int i = 0; i < NSOMO; i++)
        BF1copy[i] = IdxListBF1[i];
    int BF2copy[NSOMO];
    for(int i = 0; i < NSOMO; i++)
        BF2copy[i] = IdxListBF2[i];

    for(int i = 0; i < NSOMO; i++){
        int thisId = i;
        int nextId = BF1copy[i];
        maxcount = 0;
        while(BF1copy[thisId] != -1 && maxcount < 20){
            if(maxcount==0) *nislands += 1;
            if(maxcount==19) *nislands -= 1;

            maxcount++;

            // First the bra
            nextId = BF1copy[thisId];
            BF1copy[thisId] = -1;
            BF1copy[nextId] = -1;
            //printf("\n(%d) %d> %d -> %d\n",i,maxcount,thisId,nextId);

            // Get the phase factor bra
            if(nextId < thisId) *phasefactor *= -1;

            // Then the ket
            thisId = BF2copy[nextId];
            BF2copy[thisId] = -1;
            BF2copy[nextId] = -1;
            //printf("\n(%d) %d> %d -> %d\n",i,maxcount,nextId,thisId);

            // Get the phase factor bra
            if(nextId > thisId) *phasefactor *= -1;

        }
        //printf("\nBF1\n");
        //for(int j = 0; j < NSOMO; j++)
        //    printf("%d ",BF1copy[j]);
        //printf("\nBF2\n");
        //for(int j = 0; j < NSOMO; j++)
        //    printf("%d ",BF2copy[j]);
    }

    // Garbage collection
    free(IdxListBF1);
    free(IdxListBF2);

}

void getOverlapMatrix(int64_t Isomo, int64_t MS, double **overlapMatrixptr, int *rows, int *cols, int *NSOMOout){

    int NBF = 0;
    int NSOMO = 0;

    Tree bftree = (Tree){  .rootNode = NULL, .NBF = -1 };
    bftree.rootNode = malloc(sizeof(Node));
    (*bftree.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    generateAllBFs(Isomo, MS, &bftree, &NBF, &NSOMO);

    *NSOMOout = NSOMO;

    //printTreeDriver(&bftree, NSOMO);

    // Initialize overlap matrix
    (*overlapMatrixptr) = malloc(NBF*NBF*sizeof(double));
    (*rows) = NBF;
    (*cols) = NBF;

    double *overlapMatrix = (*overlapMatrixptr);

    //// initialize Matrix
    //for(int i = 0; i < NBF; i++)
    //    for(int j = 0; j < NBF; j++)
    //        overlapMatrix[i*NBF + j] = 0.0;

    int addI = 0;
    int addJ = 0;
    int *BF1 = malloc(MAX_SOMO * sizeof(int));
    int *BF2 = malloc(MAX_SOMO * sizeof(int));
    int *IdxListBF1 = malloc(MAX_SOMO * sizeof(int));
    int *IdxListBF2 = malloc(MAX_SOMO * sizeof(int));

    int g = 0;
    g = (NSOMO - MS)/2;
    printf("NBFs = %d NSOMOs = %d MS = %ld g = %d\n",NBF,NSOMO,MS,g);

    int nislands; // Note that nislands < g always
    int phasefactor;

    int dictPhase[2];

    dictPhase[0] = 1;
    dictPhase[1] =-1;


    // Set block elements
    for(int i = 0; i < NBF; i++){
        addI = i;
        getIthBFDriver(&bftree, NSOMO, addI, BF1);
        getBFIndexList(NSOMO, BF1, IdxListBF1);

        //printf("addI : %d > ",addI);
        //for(int k=0;k<NSOMO;k++)
        //    printf("%d ",BF1[k]);
        //printf("\n");

        for(int j = 0; j < NBF; j++){
            addJ = j;
            getIthBFDriver(&bftree, NSOMO, addJ, BF2);
            getBFIndexList(NSOMO, BF2, IdxListBF2);
            //printf("addJ : %d > ",addI);
            //for(int k=0;k<NSOMO;k++)
            //    printf("%d ",BF2[k]);
            //printf("\n");

            // Get the i and r factors
            getIslands(NSOMO, BF1, BF2, &nislands, &phasefactor);

            //printf("(%d, %d) is=%d ph=%d fac=%10.15f\n",addI, addJ, nislands, phasefactor, phasefactor*1.0/(1 << (g-nislands)));

            overlapMatrix[i*NBF + j] = -1.0*phasefactor / (1 << (g - nislands));
        }
    }

    // Garbage collection
    free(BF1);
    free(IdxListBF1);
    free(BF2);
    free(IdxListBF2);

}

void getSetBits(int64_t n, int *nsetbits){
    int count = 0;
    while(n){
        count += n & 1;
        n >>= 1;
    }
    *nsetbits = count;
}

void generateAllBFs(int64_t Isomo, int64_t MS, Tree *bftree, int *NBF, int *NSOMO){
    getSetBits(Isomo, NSOMO);
    buildTreeDriver(bftree, *NSOMO, MS, NBF);
}

void gramSchmidt(double *overlapMatrix, int rows, int cols, double *orthoMatrix){

    // vector
    double norm = 0.0;
    orthoMatrix[(rows-1)*cols + cols-1] = 1.0;
    for(int i = cols-2; i > -1; i--){ orthoMatrix[(rows-1)*cols + i] = 0.0; }

    // Gram-Schmidt loop
    for(int i = rows-2; i > -1; i--){
        for(int k = cols-1; k > -1; k--){ orthoMatrix[(i)*cols + k] = 0.0; }
        orthoMatrix[i*cols + i] = 1.0;
        for(int j = rows-1; j > i; j--){
            for(int k = rows-1; k >= j; k--){
                orthoMatrix[i*cols + j] += -1.0*orthoMatrix[j*cols + k]*overlapMatrix[i*cols + k];
            }
        }

        // Normalization
        norm = 0.0;
        for(int j = rows-1; j >= i; j--){
            norm += orthoMatrix[i*cols + j]*orthoMatrix[i*cols + j];
        }
        norm = sqrt(norm);
        for(int j = rows-1; j >= i; j--){
            orthoMatrix[i*cols + j] /= norm;
        }

    }

}

void getbftodetfunction(Tree *dettree, int NSOMO, int MS, int *BF1, double *rowvec){
    int npairs = 1 << ((NSOMO - MS)/2);
    int idxp = 0;
    int idxq = 0;
    int *detslist = malloc(npairs*NSOMO*sizeof(int));
    double *phaselist = malloc(npairs*sizeof(double));
    for(int i=0;i<npairs;i++)
        phaselist[i] = 1.0;
    int shft = npairs;
    int donepq[NSOMO];
    double fac = 1.0;
    for(int i = 0; i < NSOMO; i++)
        donepq[i] = 0.0;
    //for(int i = 0; i < NSOMO; i++)
    //    printf("%d) %d\n",i,BF1[i]);

    for(int i = 0; i < NSOMO; i++){
        idxp = BF1[i];
        idxq = BF1[idxp];
        //printf("idxp=%d idxq=%d\n",idxp,idxq);
        // Do one pair only once
        if(donepq[idxp] > 0.0 || donepq[idxq] > 0.0) continue;
        fac *= 2.0;
        donepq[idxp] = 1.0;
        donepq[idxq] = 1.0;
        for(int j = 0; j < npairs; j = j + shft){
            for(int k = 0; k < shft/2; k++){
                detslist[(k+j)*NSOMO + idxp] = 1;
                detslist[(k+j)*NSOMO + idxq] = 0;
            }
            for(int k = shft/2; k < shft; k++){
                detslist[(k+j)*NSOMO + idxp] = 0;
                detslist[(k+j)*NSOMO + idxq] = 1;
                phaselist[k+j] *=-1;
            }
        }
        shft /= 2;
    }

    // Now get the addresses
    int inpdet[NSOMO];
    int addr = -1;
    for(int i = 0; i < npairs; i++){
        for(int j = 0; j < NSOMO; j++)
            inpdet[j] = detslist[i*NSOMO + j];
        findAddofDetDriver(dettree, NSOMO, inpdet, &addr);
        //rowvec[addr] = 1.0 * phaselist[i]/sqrt(fac);
        // Upon transformation from
        // SOMO to DET basis,
        // all dets have the same phase
        rowvec[addr] = 1.0/sqrt(fac);
    }

    free(detslist);
    free(phaselist);
}

void convertBFtoDetBasis(int64_t Isomo, int MS, double **bftodetmatrixptr, int *rows, int *cols){

    int NSOMO=0;
    getSetBits(Isomo, &NSOMO);
    int ndets = 0;
    int NBF = 0;
    double dNSOMO = NSOMO*1.0;
    double nalpha = (NSOMO + MS)/2.0;
    ndets = (int)binom(dNSOMO, nalpha);
    //printf("Ndets = %d\n",ndets);

    Tree dettree = (Tree){  .rootNode = NULL, .NBF = -1 };
    dettree.rootNode = malloc(sizeof(Node));
    (*dettree.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    genDetBasis(&dettree, Isomo, MS, &ndets);

    //printTreeDriver(&dettree, NSOMO);
    //printf("Ndets = %d\n",ndets);

    //int addr = -1;
    //int inpdet[NSOMO];
    //inpdet[0] = 1;
    //inpdet[1] = 1;
    //inpdet[2] = 1;
    //inpdet[3] = 0;
    //inpdet[4] = 0;
    //inpdet[5] = 0;

    //findAddofDetDriver(&dettree, NSOMO, inpdet, &addr);

    int detlist[ndets];
    getDetlistDriver(&dettree, NSOMO, detlist);

    //printf("\n");
    //for(int i=0;i<ndets;i++)
    //    printf("%d ",detlist[i]);
    //printf("\n");

    //printf("addr of det=%d\n",addr);

    // Prepare BFs
    Tree bftree = (Tree){  .rootNode = NULL, .NBF = -1 };
    bftree.rootNode = malloc(sizeof(Node));
    (*bftree.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    generateAllBFs(Isomo, MS, &bftree, &NBF, &NSOMO);

    //printf("in convert NBFs = %d ndets=%d\n",NBF,ndets);

    // Initialize transformation matrix
    (*bftodetmatrixptr) = malloc(NBF*ndets*sizeof(double));
    (*rows) = NBF;
    (*cols) = ndets;

    double *bftodetmatrix = (*bftodetmatrixptr);

    // Build BF to det matrix
    int addI = 0;
    int addJ = 0;
    double rowvec[ndets];
    for(int i=0;i<ndets;i++)
        rowvec[i]=0.0;
    int *BF1 = malloc(MAX_SOMO * sizeof(int));
    int *BF2 = malloc(MAX_SOMO * sizeof(int));
    int *IdxListBF1 = malloc(MAX_SOMO * sizeof(int));
    int *IdxListBF2 = malloc(MAX_SOMO * sizeof(int));

    for(int i = 0; i < NBF; i++){
        addI = i;
        getIthBFDriver(&bftree, NSOMO, addI, BF1);
        getBFIndexList(NSOMO, BF1, IdxListBF1);


        //printf("addI : %d > ",addI);
        //for(int k=0;k<NSOMO;k++)
        //    printf("%d ",BF1[k]);
        //printf("\n");

        // Get ith row
        getbftodetfunction(&dettree, NSOMO, MS, IdxListBF1, rowvec);

        //printf("---%d---\n",i);
        //for(int k=0;k<ndets;k++)
        //    printf("%10.4f ",rowvec[k]);
        //printf("\n");

        //printf("(%d, %d) is=%d ph=%d fac=%10.15f\n",addI, addJ, nislands, phasefactor, phasefactor*1.0/(1 << (g-nislands)));

        for(int j = 0; j < ndets; j++)
            bftodetmatrix[i*ndets + j] = rowvec[j];

        for(int k=0;k<ndets;k++)
            rowvec[k]=0.0;
    }

    // Garbage collection
    free(BF1);
    free(IdxListBF1);
    free(BF2);
    free(IdxListBF2);

}

unsigned int shftbit(int num, int p){
    unsigned int maskleft = ~(0 | ((1<<p)-1));
    unsigned int maskright = ((1<<(p-1))-1);
    int numleft = num & maskleft;
    int numright = num & maskright;
    numleft = numleft >> 1;
    return(numleft | numright);
};

int getphase(int num, int p, int q, int nmo){
    // CSF: 1 1 1 1 1 1 1 1 1 1
    // DET: 1 1 0 0 1 1 0 0 1 0
    //        |         |
    //        p         q
    //        |         |
    // CSF: 1 1 1 1 1 1 1 1 1 1
    // DET: 1 0 0 0 1 1 1 0 1 0
    unsigned int maskleft = ~(0 | ((1<<q)-1));
    unsigned int maskright = ((1<<(p-1))-1);
    unsigned int maskmo = ((1<<nmo)-1);
    int numleft = num & maskleft;
    int numleftright = numleft & maskright;
    int nalpha = __builtin_popcount(numleftright & maskmo);
    int nbeta = p-q-1 - nalpha;
    //printf("na=%d nb=%d ( %d, %d)\n",nalpha,nbeta,p,q);
    int maskatp = (1<<(p-1));
    int nelecalphaatp = __builtin_popcount(num & maskatp);
    int maskatq = (1<<(q-1));
    int nelecalphaatq = __builtin_popcount(num & maskatq);
    //printf("nelecalphaatp=%d\n",nelecalphaatp);
    int nfermions = nelecalphaatp == 0 ? nbeta : nalpha;
    int phase = nfermions % 2 == 0 ? 1 : -1;
    if(nelecalphaatp == nelecalphaatq) phase *= -1;
    return(phase);
};

void calcMEdetpair(int *detlistI, int *detlistJ, int orbI, int orbJ, int Isomo, int Jsomo, int ndetI, int ndetJ, int NMO, double *matelemdetbasis){


    int maskI;
    int nelecatI;
    unsigned int maskleft;
    unsigned int maskright;
    unsigned int psomo;
    unsigned int qsomo;

    // E(p,q) |I> = cpq |J>


    int p,q; // The two orbitals p is always > q.
    p = orbI >= orbJ ? orbI : orbJ;
    q = orbI >= orbJ ? orbJ : orbI;

    // Find the corresponding case
    // 1. NdetI > NdetJ  (SOMO -> SOMO)
    // 2. NdetI < NdetJ  (DOMO -> VMO)
    // 3. NdetI == NdetJ (SOMO -> VMO and DOMO -> SOMO)

    // Converting the above four cases into int:
    int case_type = abs(ndetI - ndetJ) == 0 ? 3 : (ndetI > ndetJ ? 1 : 2);

    switch (case_type){
        case 1:
            // SOMO -> SOMO
            printf("SOMO->SOMO\n");
            // Find the orbital ids in model space
            maskleft = (0 | ((1<<(p))-1));
            maskright =(0 | ((1<<(q-1))-1));
            psomo = __builtin_popcount(Isomo & maskleft);
            qsomo = q == 1 ? 1 : __builtin_popcount(Isomo & maskright);
            p = psomo >= qsomo ? psomo : qsomo;
            q = psomo >= qsomo ? qsomo : psomo;

            //printf("I=%d J=%d (%d %d)\n",Isomo,Jsomo,p,q);

            //printf("SOMO->SOMO\n");
            //printf("\np=%d q=%d  (%d %d)\n",q,p,psomo,qsomo);
            for(int i=0;i<ndetI;i++){
                int idet = detlistI[i];
                //printf("leading test "BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(idet));                // Calculate phase
                int phase = getphase(idet,p,q,NMO);
                // Shift bits for
                idet = shftbit(shftbit(detlistI[i],q),p-1);
                //printf(" -> "BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(idet));
                //printf(" %d\n",phase);
                for(int j=0;j<ndetJ;j++){
                    int jdet = (detlistJ[j]);
                    if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                }
            }
            break;
        case 2:
            // DOMO -> VMO
            printf("DOMO->VMO\n");
            // Find the orbital ids in model space
            maskleft = (0 | ((1<<(p))-1));
            maskright =(0 | ((1<<(q-1))-1));
            psomo = __builtin_popcount(Jsomo & maskleft);
            qsomo = q == 1 ? 1 : __builtin_popcount(Jsomo & maskright);
            p = psomo >= qsomo ? psomo : qsomo;
            q = psomo >= qsomo ? qsomo : psomo;

            //printf("I=%d J=%d (%d %d)\n",Isomo,Jsomo,p,q);

            for(int i=0;i<ndetI;i++){
                // Get phase
                int idet = detlistI[i];
                for(int j=0;j<ndetJ;j++){
                    int jdet = (detlistJ[j]);
                    // Calculate phase
                    // Note -1 phase difference wrt to SOMO->SOMO
                    int phase = 1*getphase(jdet,p,q,NMO);
                    // Shift bits for I
                    jdet = shftbit(shftbit(detlistJ[j],q),p-1);
                    if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                }
            }
            break;
        case 3:
            // (SOMO -> VMO or DOMO -> SOMO)
            printf("SOMO->VMO and DOMO->SOMO\n");
            // Find the orbital ids in model space
            maskleft = (0 | ((1<<(p))-1));
            maskright =(0 | ((1<<(q-1))-1));
            psomo = __builtin_popcount(Isomo & maskleft);
            qsomo = q == 1 ? 1 : __builtin_popcount(Isomo & maskright);
            p = psomo >= qsomo ? psomo : qsomo;
            q = psomo >= qsomo ? qsomo : psomo;

            //printf("I=%d J=%d (>%d %d)\n",Isomo,Jsomo,p,q);
            for(int i=0;i<ndetI;i++){
                // Get phase
                int idet = detlistI[i];
                //printf("leading test "BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(idet));
                // Calculate phase
                int phase = 1*getphase(idet,p,q,NMO);
                // Shift bits for I
                idet = shftbit(detlistI[i],p);
                //printf(" -> "BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(idet));
                //printf("\n");
                for(int j=0;j<ndetJ;j++){
                //printf("\tleading test "BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(detlistJ[j]));
                    int jdet = shftbit(detlistJ[j],q);
                //printf("\t -> "BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(jdet));
                //printf("\n");
                //printf("(%d  %d) -> %d\n",i,j,phase);
                    if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
                }
            }

            //maskI = (1 << (orbI-1));
            //nelecatI = __builtin_popcount(Isomo & maskI);
            ////printf("maskI=%d nelecatI=%d \n",maskI,Isomo);

            //switch(nelecatI){
            //    case 1:
            //        // SOMO -> VMO
            //        printf("SOMO->VMO\n");
            //        // Find the orbital ids in model space
            //        maskleft = (0 | ((1<<(p))-1));
            //        maskright =(0 | ((1<<(q-1))-1));
            //        psomo = __builtin_popcount(Isomo & maskleft);
            //        qsomo = q == 1 ? 1 : __builtin_popcount(Isomo & maskright);
            //        p = psomo >= qsomo ? psomo : qsomo;
            //        q = psomo >= qsomo ? qsomo : psomo;

            //        //printf("I=%d J=%d (>%d %d)\n",Isomo,Jsomo,p,q);
            //        for(int i=0;i<ndetI;i++){
            //            // Get phase
            //            int idet = detlistI[i];
            //            //printf("leading test "BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(idet));
            //            // Calculate phase
            //            int phase = 1*getphase(idet,p,q,NMO);
            //            // Shift bits for I
            //            idet = shftbit(detlistI[i],p);
            //            //printf(" -> "BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(idet));
            //            //printf("\n");
            //            for(int j=0;j<ndetJ;j++){
            //            //printf("\tleading test "BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(detlistJ[j]));
            //                int jdet = shftbit(detlistJ[j],q);
            //            //printf("\t -> "BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(jdet));
            //            //printf("\n");
            //            //printf("(%d  %d) -> %d\n",i,j,phase);
            //                if(idet == jdet) matelemdetbasis[i*ndetJ + j] = 1.0*phase;
            //            }
            //        }
            //        break;
            //    case 0:
            //        // DOMO -> SOMO
            //        printf("DOMO->SOMO\n");
            //        // Find the orbital ids in model space
            //        maskleft = (0 | ((1<<(p))-1));
            //        maskright =(0 | ((1<<(q-1))-1));
            //        psomo = __builtin_popcount(Isomo & maskleft);
            //        qsomo = q == 1 ? 1 : __builtin_popcount(Isomo & maskright);
            //        p = psomo >= qsomo ? psomo : qsomo;
            //        q = psomo >= qsomo ? qsomo : psomo;

            //        printf("I=%d J=%d (%d %d)\n",Isomo,Jsomo,p,q);
            //        break;
            //    default:
            //        printf("Something went wrong in calcME\n");
            //        break;
            //}

            break;
        default:
            printf("Something is wrong in calc ME\n");
            break;
    } // end select
    //printRealMatrix(matelemdetbasis,ndetI,ndetJ);

}

void callcalcMEij(int Isomo, int Jsomo, int orbI, int orbJ, int MS, int NMO, double **ApqIJptr, int *rowsA, int *colsA){
    // Get dets for I
    int ndetI;
    int ndetJ;

    // Get detlist
    int NSOMOI=0;
    int NSOMOJ=0;
    getSetBits(Isomo, &NSOMOI);
    getSetBits(Jsomo, &NSOMOJ);

    Tree dettreeI = (Tree){  .rootNode = NULL, .NBF = -1 };
    dettreeI.rootNode = malloc(sizeof(Node));
    (*dettreeI.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    genDetBasis(&dettreeI, Isomo, MS, &ndetI);


    Tree dettreeJ = (Tree){  .rootNode = NULL, .NBF = -1 };
    dettreeJ.rootNode = malloc(sizeof(Node));
    (*dettreeJ.rootNode) = (Node){ .C0 = NULL, .C1 = NULL, .PREV = NULL, .addr = 0, .cpl = -1, .iSOMO = -1};

    genDetBasis(&dettreeJ, Jsomo, MS, &ndetJ);

    int detlistI[ndetI];
    int detlistJ[ndetJ];

    // Get detlist
    getDetlistDriver(&dettreeI, NSOMOI, detlistI);
    getDetlistDriver(&dettreeJ, NSOMOJ, detlistJ);
    // printdets I
    //printf("Idets\n");
    //for(int i=0;i<ndetI;i++){
    //    printf("leading test "BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(detlistI[i]));
    //    printf("\n");
    //}

    (*ApqIJptr) = malloc(ndetI*ndetJ*sizeof(double));
    (*rowsA) = ndetI;
    (*colsA) = ndetJ;
    //printf("ndetI=%d ndetJ=%d\n",ndetI,ndetJ);

    double *matelemdetbasis = (*ApqIJptr);

    for(int i=0;i<ndetI;i++)
        for(int j=0;j<ndetJ;j++)
            matelemdetbasis[i*ndetJ + j]=0.0;

    // Garbage collection
    calcMEdetpair(detlistI, detlistJ, orbI, orbJ, Isomo, Jsomo, ndetI, ndetJ, NMO, matelemdetbasis);

    //printRealMatrix(matelemdetbasis, ndetI, ndetJ);

    // Garbage collection
}

void getApqIJMatrixDims(int64_t Isomo, int64_t Jsomo, int64_t MS, int32_t *rowsout, int32_t *colsout){
    int NSOMOI=0;
    int NSOMOJ=0;
    printf("Isomo=%ld Jsomo=%ld\n",Isomo,Jsomo);
    getSetBits(Isomo, &NSOMOI);
    getSetBits(Jsomo, &NSOMOJ);
    printf("NsomoI=%d NsomoJ=%d\n",NSOMOI,NSOMOJ);
    int NBFI=0;
    int NBFJ=0;
    getncsfs(NSOMOI, MS, &NBFI);
    getncsfs(NSOMOJ, MS, &NBFJ);
    (*rowsout) = NBFI;
    (*colsout) = NBFJ;
    printf("\t >> %d %d\n",NBFI,NBFJ);
}

void getApqIJMatrixDriver(int64_t Isomo, int64_t Jsomo, int orbp, int orbq, int64_t MS, int64_t NMO, double **CSFICSFJApqIJptr, int *rowsout, int *colsout){

    double *overlapMatrixI;
    double *overlapMatrixJ;
    double *orthoMatrixI;
    double *orthoMatrixJ;
    double *bftodetmatrixI;
    double *bftodetmatrixJ;
    double *ApqIJ;
    int NSOMO=0;

    /***********************************
                   Doing I
    ************************************/
    // Fill matrix
    int rowsI = 0;
    int colsI = 0;

    getOverlapMatrix(Isomo, MS, &overlapMatrixI, &rowsI, &colsI, &NSOMO);

    //printf("\nDone Overlap Matrix I\n");
    //printRealMatrix(overlapMatrixI, rowsI, colsI);
    //printf("\nDone Overlap Matrix I\n");

    orthoMatrixI = malloc(rowsI*colsI*sizeof(double));

    gramSchmidt(overlapMatrixI, rowsI, colsI, orthoMatrixI);

    //printf("\nDone Gram-Schmidt orthonormalization I\n");
    //printRealMatrix(orthoMatrixI, rowsI, colsI);
    //printf("\nGen det basis I \n");

    int rowsbftodetI, colsbftodetI;

    convertBFtoDetBasis(Isomo, MS, &bftodetmatrixI, &rowsbftodetI, &colsbftodetI);

    //printf("\nBF to det I\n");
    //printRealMatrix(bftodetmatrixI, rowsbftodetI, colsbftodetI);
    //printf("\nBF to det I\n");

    /***********************************
                   Doing J
    ************************************/

    int rowsJ = 0;
    int colsJ = 0;
    // Fill matrix
    getOverlapMatrix(Jsomo, MS, &overlapMatrixJ, &rowsJ, &colsJ, &NSOMO);

    //printf("\nDone overlap J\n");
    //printRealMatrix(overlapMatrixJ, rowsJ, colsJ);
    //printf("\nDone overlap J\n");

    orthoMatrixJ = malloc(rowsJ*colsJ*sizeof(double));

    gramSchmidt(overlapMatrixJ, rowsJ, colsJ, orthoMatrixJ);

    //printf("\nDone Gram-Schmidt orthonormalization\n");
    //printRealMatrix(orthoMatrixJ, rowsJ, colsJ);
    //printf("\nDone Gram-Schmidt orthonormalization\n");


    int rowsbftodetJ, colsbftodetJ;

    convertBFtoDetBasis(Jsomo, MS, &bftodetmatrixJ, &rowsbftodetJ, &colsbftodetJ);

    //printf("dims BFtoDetJ rowsbftodetJ=%d colsbftodetJ=%d\n",rowsbftodetJ,colsbftodetJ);

    //printf("\nGen det basis J \n");
    //printRealMatrix(bftodetmatrixJ, rowsbftodetJ, colsbftodetJ);
    //printf("\nGen det basis  J\n");

    int rowsA = 0;
    int colsA = 0;

    callcalcMEij(Isomo, Jsomo, orbp, orbq, MS, NMO, &ApqIJ, &rowsA, &colsA);

    //printf("Done MEij\n");
    //printRealMatrix(ApqIJ, rowsA, colsA);
    //printf("Done MEij\n");

    // Final ME in BF basis

    // First transform I in bf basis
    double *bfIApqIJ = malloc(rowsbftodetI*colsA*sizeof(double));

    int transA=false;
    int transB=false;
    callBlasMatxMat(bftodetmatrixI, rowsbftodetI, colsbftodetI, ApqIJ, rowsA, colsA, bfIApqIJ, transA, transB);

    //printf("Done blas BFI\n");
    //printRealMatrix(bfIApqIJ, colsI, colsA);

    // now transform I in csf basis
    double *CSFIApqIJ = malloc(rowsI*colsA*sizeof(double));
    transA = false;
    transB = false;
    callBlasMatxMat(orthoMatrixI, rowsI, colsI, bfIApqIJ, colsI, colsA, CSFIApqIJ, transA, transB);

    //printf("Done blas CSFI\n");
    //printRealMatrix(CSFIApqIJ, rowsI, colsA);
    //printf("Done blas CSFI\n");

    // now transform J in BF basis
    double *CSFIbfJApqIJ = malloc(rowsI*rowsbftodetJ*sizeof(double));
    //printf("rowsI = %d colsA=%d | rowsbftodetJ=%d colsbftodetJ=%d\n",rowsI,colsA,rowsbftodetJ,colsbftodetJ);
    transA = false;
    transB = true;
    callBlasMatxMat(CSFIApqIJ, rowsI, colsA, bftodetmatrixJ, rowsbftodetJ, colsbftodetJ, CSFIbfJApqIJ, transA, transB);

    //printf("Done blas BFJ\n");
    //printRealMatrix(CSFIbfJApqIJ, rowsI, rowsbftodetJ);
    //printf("Done blas BFJ\n");

    // now transform J in CSF basis
    (*CSFICSFJApqIJptr) = malloc(rowsI*rowsJ*sizeof(double));
    (*rowsout) = rowsI;
    (*colsout) = rowsJ;

    double *CSFICSFJApqIJ = (*CSFICSFJApqIJptr);
    transA = false;
    transB = true;
    callBlasMatxMat(CSFIbfJApqIJ, rowsI, rowsbftodetJ, orthoMatrixJ, rowsJ, colsJ, CSFICSFJApqIJ, transA, transB);

    //printf("ME CSF basis\n");
    //printRealMatrix(CSFICSFJApqIJ, rowsI, rowsJ);


    // Garbage collection
    free(overlapMatrixI);
    free(overlapMatrixJ);
    free(orthoMatrixI);
    free(orthoMatrixJ);
    free(bftodetmatrixI);
    free(bftodetmatrixJ);
    free(ApqIJ);
    free(bfIApqIJ);
    free(CSFIApqIJ);
    free(CSFIbfJApqIJ);
}

void getApqIJMatrixDriverArrayInp(int64_t Isomo, int64_t Jsomo, int32_t orbp, int32_t orbq, int64_t MS, int64_t NMO, double *CSFICSFJApqIJ, int32_t rowsmax, int32_t colsmax){

    double *overlapMatrixI;
    double *overlapMatrixJ;
    double *orthoMatrixI;
    double *orthoMatrixJ;
    double *bftodetmatrixI;
    double *bftodetmatrixJ;
    double *ApqIJ;
    int NSOMO=0;

    /***********************************
                   Doing I
    ************************************/
    // Fill matrix
    int rowsI = 0;
    int colsI = 0;

    getOverlapMatrix(Isomo, MS, &overlapMatrixI, &rowsI, &colsI, &NSOMO);

    //printf("\nDone Overlap Matrix I\n");
    //printRealMatrix(overlapMatrixI, rowsI, colsI);
    //printf("\nDone Overlap Matrix I\n");

    orthoMatrixI = malloc(rowsI*colsI*sizeof(double));

    gramSchmidt(overlapMatrixI, rowsI, colsI, orthoMatrixI);

    //printf("\nDone Gram-Schmidt orthonormalization I\n");
    //printRealMatrix(orthoMatrixI, rowsI, colsI);
    //printf("\nGen det basis I \n");

    int rowsbftodetI, colsbftodetI;

    convertBFtoDetBasis(Isomo, MS, &bftodetmatrixI, &rowsbftodetI, &colsbftodetI);

    //printf("\nBF to det I\n");
    //printRealMatrix(bftodetmatrixI, rowsbftodetI, colsbftodetI);
    //printf("\nBF to det I\n");

    /***********************************
                   Doing J
    ************************************/

    int rowsJ = 0;
    int colsJ = 0;
    // Fill matrix
    getOverlapMatrix(Jsomo, MS, &overlapMatrixJ, &rowsJ, &colsJ, &NSOMO);

    //printf("\nDone overlap J\n");
    //printRealMatrix(overlapMatrixJ, rowsJ, colsJ);
    //printf("\nDone overlap J\n");

    orthoMatrixJ = malloc(rowsJ*colsJ*sizeof(double));

    gramSchmidt(overlapMatrixJ, rowsJ, colsJ, orthoMatrixJ);

    //printf("\nDone Gram-Schmidt orthonormalization\n");
    //printRealMatrix(orthoMatrixJ, rowsJ, colsJ);
    //printf("\nDone Gram-Schmidt orthonormalization\n");


    int rowsbftodetJ, colsbftodetJ;

    convertBFtoDetBasis(Jsomo, MS, &bftodetmatrixJ, &rowsbftodetJ, &colsbftodetJ);

    //printf("dims BFtoDetJ rowsbftodetJ=%d colsbftodetJ=%d\n",rowsbftodetJ,colsbftodetJ);

    //printf("\nGen det basis J \n");
    //printRealMatrix(bftodetmatrixJ, rowsbftodetJ, colsbftodetJ);
    //printf("\nGen det basis  J\n");

    int rowsA = 0;
    int colsA = 0;

    callcalcMEij(Isomo, Jsomo, orbp, orbq, MS, NMO, &ApqIJ, &rowsA, &colsA);

    //printf("Done MEij\n");
    //printRealMatrix(ApqIJ, rowsA, colsA);
    //printf("Done MEij\n");

    // Final ME in BF basis

    // First transform I in bf basis
    double *bfIApqIJ = malloc(rowsbftodetI*colsA*sizeof(double));

    int transA=false;
    int transB=false;
    callBlasMatxMat(bftodetmatrixI, rowsbftodetI, colsbftodetI, ApqIJ, rowsA, colsA, bfIApqIJ, transA, transB);

    //printf("Done blas BFI\n");
    //printRealMatrix(bfIApqIJ, colsI, colsA);

    // now transform I in csf basis
    double *CSFIApqIJ = malloc(rowsI*colsA*sizeof(double));
    transA = false;
    transB = false;
    callBlasMatxMat(orthoMatrixI, rowsI, colsI, bfIApqIJ, colsI, colsA, CSFIApqIJ, transA, transB);

    //printf("Done blas CSFI\n");
    //printRealMatrix(CSFIApqIJ, rowsI, colsA);
    //printf("Done blas CSFI\n");

    // now transform J in BF basis
    double *CSFIbfJApqIJ = malloc(rowsI*rowsbftodetJ*sizeof(double));
    //printf("rowsI = %d colsA=%d | rowsbftodetJ=%d colsbftodetJ=%d\n",rowsI,colsA,rowsbftodetJ,colsbftodetJ);
    transA = false;
    transB = true;
    callBlasMatxMat(CSFIApqIJ, rowsI, colsA, bftodetmatrixJ, rowsbftodetJ, colsbftodetJ, CSFIbfJApqIJ, transA, transB);

    //printf("Done blas BFJ\n");
    //printRealMatrix(CSFIbfJApqIJ, rowsI, rowsbftodetJ);
    //printf("Done blas BFJ\n");

    // now transform J in CSF basis
    //(*CSFICSFJApqIJptr) = malloc(rowsI*rowsJ*sizeof(double));
    //(*rowsout) = rowsI;
    //(*colsout) = rowsJ;

    //double *CSFICSFJApqIJ = (*CSFICSFJApqIJptr);
    transA = false;
    transB = true;
    callBlasMatxMat(CSFIbfJApqIJ, rowsI, rowsbftodetJ, orthoMatrixJ, rowsJ, colsJ, CSFICSFJApqIJ, transA, transB);

    //printf("ME CSF basis\n");
    //printRealMatrix(CSFICSFJApqIJ, rowsI, rowsJ);


    // Garbage collection
    free(overlapMatrixI);
    free(overlapMatrixJ);
    free(orthoMatrixI);
    free(orthoMatrixJ);
    free(bftodetmatrixI);
    free(bftodetmatrixJ);
    free(ApqIJ);
    free(bfIApqIJ);
    free(CSFIApqIJ);
    free(CSFIbfJApqIJ);
}

void calculateMETypeSOMOSOMO(int *BF1, int *BF2, int moi, int moj, double *factor, int *phasefactor){

    // Calculate the factor following rules in the table
    // find the type

    if(BF1[moi] == moj){
        // Type I
        (*factor) = sqrt(2.0);
        (*phasefactor) = 1;
    }
    else{
        if(BF1[moi] != moi || BF1[moj] != moj){
            // Type II, III and IV
            (*factor) = 1.0/sqrt(2.0);
            (*phasefactor) =-1;
        }
        else{
            // Type V
            (*factor) = 0.0;
            (*phasefactor) = 1;
        }
    }


}
