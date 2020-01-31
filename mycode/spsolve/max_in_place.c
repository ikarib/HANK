/*====================================================================
 *
 * spsolve.c  updates a data matrix in-place with the max value
 *                 of the matrix and a 2nd matrix of the same size
 *
 * The calling syntax is:
 *
 *		max_in_place(matrix1, matrix2)
 *
 * matrix1 will be updated with the maximal values from corresponding
 * indices of the 2 matrices
 *
 * Both inputs must be double 2D real non-sparse matrices of same size
 *
 * Yair Altman 2018-07-18
 * http://UndocumentedMatlab.com/blog/multi-threaded-mex
 *
 * Adapted from original work by Dirk-Jan Kroon
 * http://mathworks.com/matlabcentral/profile/authors/1097878-dirk-jan-kroon
 *
 *==================================================================*/

#include <math.h>
#include "mex.h"

/* undef needed for LCC compiler */
#undef EXTERN_C
#ifdef _WIN32
    #include <windows.h>
    #include <process.h>
#else
    #include <pthread.h>
#endif

/* Input Arguments */
#define	hMatrix1	prhs[0]
#define	hMatrix2	prhs[1]

/* Macros */
#if !defined(MAX)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

/* Main processing loop function */
void main_loop(const mxArray *prhs[], int startIdx, int endIdx)
{
    /* Assign pointers to the various parameters */
    double *p1 = mxGetPr(hMatrix1);
    double *p2 = mxGetPr(hMatrix2);
    /*mexPrintf("comparing matrix idx = %d : %d\n", startIdx, endIdx); */

    /* Loop through all matrix coordinates */
    for (int idx=startIdx; idx<=endIdx; idx++)
    {
        /* Update hMatrix1 with the maximal value of hMatrix1,hMatrix2 */
        if (p1[idx] < p2[idx]) {
            p1[idx] = p2[idx];
        }
    }
}

/* Computation function in threads */
#ifdef _WIN32
  unsigned __stdcall thread_func(void *ThreadArgs_) {
#else
  void thread_func(void *ThreadArgs_) {
#endif
    double **ThreadArgs = ThreadArgs_;  /* void* => double** */
    const mxArray** prhs = (const mxArray**) ThreadArgs[0];

    int ThreadID = (int) ThreadArgs[1][0];
    int startIdx = (int) ThreadArgs[2][0];
    int endIdx   = (int) ThreadArgs[3][0];
    /*mexPrintf("Starting thread #%d: idx=%d:%d\n", ThreadID, startIdx, endIdx); */

    /* Run the main processing function */
    main_loop(prhs, startIdx, endIdx);

    /* Explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
    #ifdef _WIN32
        _endthreadex( 0 );
        return 0;
    #else
        pthread_exit(NULL);
    #endif
}

/* Validate inputs */
void validateInputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments */
    if (nrhs != 2) { 
        mexErrMsgIdAndTxt( "MATLAB:max_in_place:invalidNumInputs",
                "2 input arguments are required, but %d were provided.", nrhs);
    }

    if (nlhs > 0) {
        mexErrMsgIdAndTxt( "MATLAB:max_in_place:maxlhs",
                "Too many output arguments (no outputs are allowed)."); 
    }

    /* Ensure the inputs are of the correct data type */
    if (!mxIsDouble(hMatrix1) || mxIsComplex(hMatrix1) || mxIsSparse(hMatrix1)) {
        mexErrMsgIdAndTxt( "MATLAB:max_in_place:invalidInputType",
                "1st input must be a 2D double real non-sparse input matrix.");
    }

    if (!mxIsDouble(hMatrix2) || mxIsComplex(hMatrix2) || mxIsSparse(hMatrix2)) {
        mexErrMsgIdAndTxt( "MATLAB:max_in_place:invalidInputType",
                "2nd input must be a 2D double real non-sparse input matrix.");
    }

    /* Ensure the inputs have the same size (not necessarily same dimensions) */ 
    size_t n1 = mxGetNumberOfElements(hMatrix1);
    size_t n2 = mxGetNumberOfElements(hMatrix2);
    if (n1 != n2) {
        mexErrMsgIdAndTxt( "MATLAB:max_in_place:unequalInputsSize",
                "Both input matrices must have the same number of elements.");
    }
}

/* Main entry function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{
    /* Validate the inputs */
    validateInputs(nlhs, plhs, nrhs, prhs);

    /* Quick bail-out in the trivial case of empty inputs */
    if (mxIsEmpty(hMatrix1))  return;

    /* Get the number of threads from the Matlab engine (maxNumCompThreads) */
    mxArray *matlabCallOut[1] = {0};
    mxArray *matlabCallIn[1]  = {0};
    mexCallMATLAB(1, matlabCallOut, 0, matlabCallIn, "maxNumCompThreads");
    double *Nthreadsd = mxGetPr(matlabCallOut[0]);
    int Nthreads = (int) Nthreadsd[0];

    /* Get the number of elements to process */
    size_t n1 = mxGetNumberOfElements(hMatrix1);

    if (Nthreads == 1) {

        /* Process the inputs directly (not via a thread) */
        main_loop(prhs, 0, n1-1);

    } else {  /* multi-threaded */

        /* Allocate memory for handles of worker threads */
        #ifdef _WIN32
            HANDLE    *ThreadList = (HANDLE*)   malloc(Nthreads*sizeof(HANDLE));
        #else
            pthread_t *ThreadList = (pthread_t*)malloc(Nthreads*sizeof(pthread_t));
        #endif

        /* Allocate memory for the thread arguments (attributes) */
        double **ThreadID, **ThreadStartIdx, **ThreadEndIdx, ***ThreadArgs;
        double *ThreadID1, *ThreadStartIdx1, *ThreadEndIdx1, **ThreadArgs1;

        ThreadID       = (double **) malloc( Nthreads* sizeof(double *) );
        ThreadStartIdx = (double **) malloc( Nthreads* sizeof(double *) );
        ThreadEndIdx   = (double **) malloc( Nthreads* sizeof(double *) );
        ThreadArgs     = (double ***)malloc( Nthreads* sizeof(double **) );

        /* Launch the requested number of threads */
        int i;
        int threadBlockSize = ceil( ((double)n1) / Nthreads );
        for (i=0; i<Nthreads; i++)
        {
            /* Create thread ID */
            ThreadID1 = (double *)malloc( 1* sizeof(double) );
            ThreadID1[0] = i;
            ThreadID[i] = ThreadID1;

            /* Compute start/end indexes for this thread */
            ThreadStartIdx1 = (double *) malloc( sizeof(double) );
            ThreadStartIdx1[0] = i * threadBlockSize;
            ThreadStartIdx[i] = ThreadStartIdx1;

            ThreadEndIdx1 = (double *) malloc( sizeof(double) );
            ThreadEndIdx1[0] = MIN((i+1)*threadBlockSize, n1) - 1;
            ThreadEndIdx[i] = ThreadEndIdx1;

            /* Prepare thread input args */
            ThreadArgs1 = (double **) malloc( 4* sizeof(double*) );
            ThreadArgs1[0] = (double *) prhs;
            ThreadArgs1[1] = ThreadID[i];
            ThreadArgs1[2] = ThreadStartIdx[i];
            ThreadArgs1[3] = ThreadEndIdx[i];

            ThreadArgs[i] = ThreadArgs1;

            /* Launch the thread with its associated args */
            #ifdef _WIN32
                ThreadList[i] = (HANDLE)_beginthreadex(NULL, 0, &thread_func, ThreadArgs[i], 0, NULL);
            #else
                pthread_create ((pthread_t*)&ThreadList[i], NULL, (void *) &thread_func, ThreadArgs[i]);
            #endif
        }

        /* Wait for all the treads to finish working */
        #ifdef _WIN32
            for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
            for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }
        #else
            for (i=0; i<Nthreads; i++) { pthread_join(ThreadList[i],NULL); }
        #endif

        /* Free the memory resources allocated for the threads */
        for (i=0; i<Nthreads; i++)
        {
            free(ThreadArgs[i]);
            free(ThreadID[i]);
            free(ThreadStartIdx[i]);
            free(ThreadEndIdx[i]);
        }

        free(ThreadArgs);
        free(ThreadID );
        free(ThreadStartIdx);
        free(ThreadEndIdx);
        free(ThreadList);
    }

    return;
}
