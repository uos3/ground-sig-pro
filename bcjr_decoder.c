/**
 * @file bcjr_decoder.c
 * @author your name (you@domain.com)
 * @brief BCJR algorithm is implemented to estimate the original bit sequence sent. 
 * @version 0.1
 * @date 2022-02-15
 * 
 * Refer to GDP 17-18
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* note, the second half of the table must only be the second transition to each state */

static int translte8[16][4] =  {
    {1-1,        1-1,        0,          0},
    {2-1,        5-1,        0,          0},
    {3-1,        6-1,        0,          1},
    {4-1,        2-1,        0,          1},
    {5-1,        3-1,        0,          1},
    {6-1,        7-1,        0,          1},
    {7-1,        8-1,        0,          0},
    {8-1,        4-1,        0,          0},
    {1-1,        5-1,        1,          1},
    {2-1,        1-1,        1,          1},
    {3-1,        2-1,        1,          0},
    {4-1,        6-1,        1,          0},
    {5-1,        7-1,        1,          0},
    {6-1,        3-1,        1,          0},
    {7-1,        4-1,        1,          1},
    {8-1,        8-1,        1,          1}
};

/**
 * @brief The statement chooses the biggest number out of 2 different variables. Note __typeof__ 
 * used if writing a header file that must work when included in ISO C programs.
 */
#ifndef max
 #define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#endif


/**
 * @brief Variation of complex BCJR algorithm where all multiplications are turned into sums to reduce
 * computation complexity, where log-MAP algorithm is being used. 
 */
__inline  double maxstar(double a, double b)
{
    /* fabs() takes a single argument (in double) and returns absolute value of that number */
    return  max(a,b) + log(1.0+exp(-fabs(a-b)));
}

/**
 * @brief BCJR is an algorithm for maximum a posteriori (MAP) decoding of error correcting codes (ECC) 
 * defined on trellises (Covolution Code), where instead of just going forward as Viterbi decoding it 
 * goes forward and backwards getting a posteriori log likelihood ratio for each message width.
 * 
 * While Viterbi decoding used ML (finds the codeword which is closest to distance), while MAP gives 
 * LLR for each bit given all the received information (soft output) produces probabilities on the bits, 
 * not just the likely bits.
 * 
 * More info:
 * https://www.semanticscholar.org/paper/From-BCJR-to-turbo-decoding%3A-MAP-algorithms-made-Abrantes/6e7f04a3c3de4957238e2fadfdabdd0d272b362e
 * 
 * @param uncoded_in 
 * @param coded_in1 
 * @param uncoded_out 
 * @param len 
 * @param trans 
 * @param state_count 
 * @param trans_count 
 * @param last_state 
 */
void bcjr_decoder(double *uncoded_in, double *coded_in1, double *uncoded_out,
        int len, int trans[][4], int state_count, int trans_count, int last_state)
{

    int i,j;
    int transc = trans_count;
    int set0,set1;
    double p1,p0;
    int a,b;
    double **gammas;
    double **alphas;
    double **betas;
    double temp[16];
    int states = state_count;

    /* calculate gammas 
    gamma is a probability representing present (received symbol is y_k at time k and current state is S_k=s)
    */
    gammas = (double**) malloc (transc*sizeof(double *));
    for (i = 0; i < transc; i++)
        gammas[i] = (double*) malloc(len*sizeof(double));

    for (i = 0; i < transc; i++)
    {
        if (trans[i][2] && trans[i][3]){
            for (j = 0; j < len; j++){
                gammas[i][j] = uncoded_in[j] + coded_in1[j];
            }
        }
        else if (trans[i][2]){
            for (j = 0; j < len; j++)
                gammas[i][j] = uncoded_in[j];
        }
        else if (trans[i][3]){
            for (j = 0; j < len; j++)
                gammas[i][j] = coded_in1[j];
        }
        else{
            for (j = 0; j < len; j++)
                gammas[i][j] = 0;
        }
    }



    /* set and initialise memory */
    alphas = (double**) malloc (states*sizeof(double *));
    for (i = 0; i < states; i++)
        alphas[i] = (double*) malloc(len*sizeof(double));
    betas = (double**) malloc (states*sizeof(double *));
    for (i = 0; i < states; i++)
        betas[i] = (double*) malloc(len*sizeof(double));


    /* forward recursion by computing forward probabilities (alphas), computed as sequence y is received, calculated
    * by going forward from the beginning to the end of the trellis
    */
    alphas[0][0] = 0;        /* first state */
    for (i = 1; i < states; i++)
        alphas[i][0] = -9000;
     for (i = 1; i < len; i++)
     {
         for (j = 0; j < states; j++)
             temp[trans[j][1]] = alphas[trans[j][0]][i-1] + gammas[j][i-1];
         for (j = states; j < transc; j++)
             alphas[trans[j][1]][i] = maxstar( temp[trans[j][1]],  alphas[trans[j][0]][i-1] + gammas[j][i-1]  );
     }



    /* backwards recursion by computing backwward probabilities(betas), only computed after whole sequence y is received  */
    /* double betas[8][len]; */
    if (last_state < 1 || last_state > states){
        for (i = 0; i < states; i++)
            betas[i][len-1] = 0;     /* end state unknown */
    }
    else
    {
        for (i = 0; i < states; i++)
            betas[i][len-1] = -9000;
        betas[last_state-1][len-1] = 0;
    }
    for (i = len-2; i >= 0; i--)
    {
        for (j = 0; j < states; j++)
            temp[trans[j][0]] = betas[trans[j][1]][i+1] + gammas[j][i+1];
        for (j = states; j < transc; j++)
            betas[trans[j][0]][i] = maxstar( temp[trans[j][0]],  betas[trans[j][1]][i+1] + gammas[j][i+1]  );
    }


    /*
     * if (trans_prob != (double*)0)
     * {
     * mexPrintf("trans probs\n");
     * for (i=0;i<transc;i++)
     * {
     * mexPrintf("%f ",trans_prob[i]);
     *
     *
     * }
     * }
     * mexPrintf("GAMMAS\n");
     * for (i=0;i<transc;i++)
     * {
     * for (j=0;j<len;j++)
     * mexPrintf("%f ",gammas[i][j]);
     * mexPrintf("\n");
     *
     * }
     * mexPrintf("ALPHAS\n");
     * for (i=0;i<states;i++)
     * {
     * for (j=0;j<len;j++)
     * mexPrintf("%f ",alphas[i][j]);
     * mexPrintf("\n");
     *
     * }
     * mexPrintf("\n");
     * mexPrintf("\n");
     * mexPrintf("BETAS\n");
     * for (i=0;i<states;i++)
     * {
     * for (j=0;j<len;j++)
     * mexPrintf("%f ",betas[i][j]);
     * mexPrintf("\n");
     *
     * }
     * mexPrintf("\n");
     */

    /* deltas */
    /* reuse gammas memory */
    for (i = 0; i < transc; i++)
    {
        a = trans[i][0];
        b = trans[i][1];
        for (j = 0; j < len; j++)
            gammas[i][j] += alphas[a][j] + betas[b][j];
    }


    /* extrinisic uncoded llr 
    log-likelyhood Ratio (LLR) 
    */
    for (j = 0; j < len; j++)
    {
        set0 = 0;
        set1 = 0;

        for (i = 0; i < transc; i++)
        {
            if (trans[i][2]){
                if (set1)
                    p1 = maxstar(p1,gammas[i][j]);
                else
                {
                    set1 = 1;
                    p1 = gammas[i][j];
                }
            }
            else{
                if (set0)
                    p0 = maxstar(p0,gammas[i][j]);
                else
                {
                    set0 = 1;
                    p0 = gammas[i][j];
                }
            }
        }
        uncoded_out[j] = p1 - p0 - uncoded_in[j];
    }


    /* Gammas are erased from memory*/
    for (i = 0; i < transc; i++){
        free(gammas[i]);
    }
    free(gammas);

    /* Alphas are erased from memory*/
    for (i = 0; i < states; i++){
        free(alphas[i]);
    }
    free(alphas);

    /* Betas are erased from memory*/
    for (i = 0; i < states; i++){
        free(betas[i]);
    }
    free(betas);


}
