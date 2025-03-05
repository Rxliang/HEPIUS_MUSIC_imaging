#if WIN
#include <windows.h>
#include <process.h>  
#endif
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include "matrix.h"
#include "mex.h"

/* Input Arguments */

#define	I prhs[0]
#define	PRFIn prhs[1]
#define	dopPRFIn prhs[2]
#define	ttnaIn prhs[3]
#define	nDopPRIs prhs[4]
#define	fileSpecIn prhs[5]
#define datenumIn prhs[6]
#define originIn prhs[7]
#define PDeltaIn prhs[8]
#define freqIn prhs[9]
#define closeFileAfterWriteIn prhs[10]
#define instanceNoIn prhs[11]
#define modeIn prhs[12]
#define freqTWIn prhs[13]
#define fileVersionIn prhs[14]

//#define DEBUG 1

/* Output Arguments */

typedef struct {
  float *iMat;
  size_t len;
  FILE *fid;
  int threadID;
  int closeFileAfterWriting;
  unsigned int instanceNo;
} saveIQStructType;

static volatile int WaitForThread1;
static volatile int WaitForThread2;

static float *im;

#if WIN
unsigned __stdcall save_iq(saveIQStructType *saveIQStruct) {
#else
unsigned int save_iq(saveIQStructType *saveIQStruct) {
#endif

#if DEBUG
  printf("len =  %d, fid = %x ", saveIQStruct->len,  saveIQStruct->fid);
  printf("part =  %x\n", saveIQStruct->iqMat);    
#endif
  size_t ret;
  
  if (saveIQStruct->iMat != NULL) {
    ret = fwrite(saveIQStruct->iMat, sizeof(float),
		 saveIQStruct->len,
	   saveIQStruct->fid);
  }

  if (ret==0)
    printf("*** error writing file \n");
  
  fflush(saveIQStruct->fid);
  
  if (saveIQStruct->closeFileAfterWriting) {
    fclose(saveIQStruct->fid);
#if DEBUG
    printf("instance %u files closed\n", saveIQStruct->instanceNo);
#endif
  }

#if DEBUG
  printf("thread %d done \n",saveIQStruct->threadID );
#endif
  
  switch (saveIQStruct->threadID) {
  case 1:
    WaitForThread1=0;
    break;
  case 2:
    WaitForThread2=0;
  }
#if WIN
  _endthreadex( 0 );
#endif
  return 0;
}

// Note, if closeFileAfterWrite is 0, newer filenames will be ignored,
// and data will be appended to file specified at last time
// firstRun==1. consult stresssaveiq.m for usage and testing

// iq local storage is used because Matlab will sometimes
// throw access violation exceptions when reading the input
// copy inside the detached threads.

void mexFunction(
	int		nlhs,
	mxArray	*plhs[],
	int		nrhs,
	const mxArray *prhs[]
	)
{

  mwSize fileVersion;
#if WIN
#else
  typedef int HANDLE;
#endif

  static HANDLE *ThreadList; // Handles to the worker threads
  
  static int firstRun=1;
  static FILE *fid[1];
  static mwSize lenAlloc;
  static int cnt=0;

  static saveIQStructType saveIQStruct1;
  static char outFilename1[1024];

  char *mode_ptr;
  double *freqTW_ptr;
  
  if ((nrhs!=15)) {
    printf("Number of input arguments = %u\n", nrhs);
    mexErrMsgTxt("saveim.c: Wrong number of arguments input\n");
    return;
  }

  fileVersion = (mwSize)(*((double*)mxGetPr(fileVersionIn)));
  mode_ptr = (char*)mxGetChars(modeIn);
  freqTW_ptr = (double*)mxGetPr(freqTWIn);

  if (mxIsComplex(I)) {
    mexErrMsgTxt("saveim.c: Image data must be real.\n");
    return;
  }
  
  if (  WaitForThread1 ) {
    printf("*** thread (%d) not done yet. frame block skipped. *** \n",
	   WaitForThread1);
    return;
  }
  
  double *origin_ptr = (double*)mxGetPr(originIn);
  double *PDelta_ptr = (double*)mxGetPr(PDeltaIn);
  double *freq_ptr = (double*)mxGetPr(freqIn);
  double *PRF_ptr = (double*)mxGetPr(PRFIn);
  double *dopPRF_ptr = (double*)mxGetPr(dopPRFIn);
  double *ttna_ptr = (double*)mxGetPr(ttnaIn);
  double *nDopPRIs_ptr = (double*)mxGetPr(nDopPRIs);    
  double PRF = PRF_ptr[0]; 
  double dopPRF = dopPRF_ptr[0];
  double *datenum_ptr = (double*)mxGetPr(datenumIn);
  double *closeFileAfterWrite_ptr = (double*)mxGetPr(closeFileAfterWriteIn);
  double *instanceNo_ptr = (double*)mxGetPr(instanceNoIn);
  unsigned int instanceNo = (unsigned int)(*instanceNo_ptr);
  unsigned int closeFileAfterWrite = (*closeFileAfterWrite_ptr!=0);

#if DEBUG
  printf("instance %u started with firstRun = %d\n", instanceNo, firstRun);
#endif      
 
  if ( (dopPRF == 0) & (!firstRun) ) {
#if DEBUG
      printf("Closing files because dopPRF == 0 and !firstRun\n");
#endif      
      fclose(fid[0]);
      free(ThreadList);
      free(im);
      firstRun=1;
      return;
  }
   
  // avoid case where this function persists after VSX stopped and close all
  if (dopPRF < 0) {
#if DEBUG
    printf("saveim.c: dopPRF < 0, force init\n");
    printf("WaitForThread1 = %d\n",  WaitForThread1);
    printf("closeFileAfterWrite = %d\n", closeFileAfterWrite);
#endif
    firstRun=1;
    dopPRF=-dopPRF;
  }
    
  size_t nDim = mxGetNumberOfDimensions(I);

  if ( ! ((nDim != 2) | (nDim != 3)) ) {
    mexErrMsgTxt(
      "Image matrix needs to be 2D (one frame) or 3D (multiple frame) \n");
  }
    
  mwSize* dimIQPre =  mxGetDimensions(I);
  mwSize* dimdateNum =  mxGetDimensions(datenumIn);    
				
  mwSize dimIQ[3];
  dimIQ[2] = 1;      

  for (int i = 0; i < nDim; i++) {
    dimIQ[i] = dimIQPre[i];
  }

  mwSize len = dimIQ[0]*dimIQ[1]*dimIQ[2];

  if (!firstRun & (len > lenAlloc)) {
    free(im);
    im = (float*)malloc(len*sizeof(float));
#if DEBUG    
    printf("Image data size increase from length: %d to %d, reallocated.\n",
	   lenAlloc, len);
#endif
    lenAlloc = len;
  }

#if DEBUG
  printf("DopPRF: %lf ", dopPRF);
  printf("freq: %lf ", *freq_ptr);  
  printf("IQ dimension: ");
  for (int i = 0; i < 3; i++) {
    printf("[%d]", dimIQ[i]);
  }
  printf("\n");

  printf("dateNum dimension: %d x %d\n", dimdateNum[0], dimdateNum[1]);
  
#endif
    
  if ((dimdateNum[0] != 1) | (dimdateNum[1] != dimIQ[2])) {
    printf("dateNum has %d columns, IQMat has %d frames.\n", 
	   dimdateNum[1], dimIQ[2]);
    mexErrMsgTxt("dateNum needs to be same length as number of frames.\n");
  }

  // if we get here, all threads are done from any previous calls
  // if we are closing file after writing, we need to open new files now
  if (firstRun | closeFileAfterWrite) {
    // local copies of IQ, to avoid access violation exceptions from Matlab
    im = (float*)malloc(len*sizeof(float));

    if (im==NULL)
      mexErrMsgTxt("saveim: Error allocating im\n");
    
    mwSize strln;
    strln = mxGetN(fileSpecIn)+1;
    mxGetString(fileSpecIn, outFilename1, strln);
    
    printf("Image output file: %s\n", outFilename1);
				
    fid[0] = fopen(outFilename1, "wb");

    if ((fid[0] == NULL) ){
      printf("fid[0]: %x \n", fid[0]);
      mexErrMsgTxt("saveim: Could not open file for writing!\n");
    }
	
#if DEBUG
    printf("instance %u files opened: fid[0] = %x \n",
	   instanceNo, fid[0]);
#endif
    
    //printf("cnt = %d, len =  %d, fid1 = %u, fid2 = %u \n", cnt, len, fid[0], fid[1]);
    
    for (int i = 0; i < 1; i++) {
      fwrite(&fileVersion, sizeof(fileVersion), 1, fid[i]);
      fwrite(dimIQ, sizeof(mwSize), 3, fid[i]);
      fwrite(&PRF, sizeof(double), 1, fid[i]);
      fwrite(&dopPRF, sizeof(double), 1, fid[i]);
      fwrite(ttna_ptr, sizeof(double), 1, fid[i]);
      fwrite(nDopPRIs_ptr, sizeof(double), 1, fid[i]);    
      fwrite(datenum_ptr, sizeof(double), dimIQ[2], fid[i]);
      fwrite(origin_ptr, sizeof(double), 3, fid[i]);
      fwrite(PDelta_ptr, sizeof(double), 3, fid[i]);
      fwrite(freq_ptr, sizeof(double), 1, fid[i]);
      fwrite(mode_ptr, sizeof(char), 1, fid[i]);
      fwrite(freqTW_ptr, sizeof(double), 1, fid[i]);
      }
    }
    
    saveIQStruct1.fid = fid[0];
    saveIQStruct1.threadID = 1;    
    saveIQStruct1.len = len;
    saveIQStruct1.instanceNo = instanceNo;

  if (firstRun) {
    ThreadList = (HANDLE*)malloc(1*sizeof( HANDLE ));
    firstRun=0;
  }

  #if DEBUG
  printf("b4 copy imIn: %x ", imIn);
  #endif

  float *imIn = (float*)mxGetPr(I);
    
#if DEBUG  
  printf("Will now copy %d floats\n", len);
#endif
  
  memcpy(im, imIn, len*sizeof(float));
  
  saveIQStruct1.iMat = im;

#if DEBUG
  printf("iMat: %x ", saveIQStruct1.iMat);
#endif
  
  // closing files is done in the thread
  saveIQStruct1.closeFileAfterWriting = closeFileAfterWrite;
 
  /*  printf("WaitForThread1 = %d, WaitForThread2 = %d \n",
      WaitForThread1, WaitForThread2);*/

  WaitForThread1 = 1;
 
#if WIN
  // Reserve room for handles of threads in ThreadList
  ThreadList[0] = (HANDLE)_beginthreadex( NULL, 0, &save_iq,
					  &saveIQStruct1 , 0, NULL );
  // this does not kill the thread, only removes the handle
  CloseHandle( ThreadList[0] );
#else
  save_iq(&saveIQStruct1);
#endif    
  return;
}
