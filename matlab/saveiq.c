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

#define	IQ prhs[0]
#define	PRFIn prhs[1]
#define	dopPRFIn prhs[2]
#define	ttnaIn prhs[3]
#define	nDopPRIs prhs[4]
#define	fileSpecInRe prhs[5]
#define	fileSpecInIm prhs[6]
#define datenumIn prhs[7]
#define originIn prhs[8]
#define PDeltaIn prhs[9]
#define freqIn prhs[10]
#define closeFileAfterWriteIn prhs[11]
#define instanceNoIn prhs[12]
#define modeIn prhs[13]
#define freqTWIn prhs[14]
#define lateralSubSampleIndexIn prhs[15]
#define rangeSubSampleIndexIn prhs[16]
#define fileVersionIn prhs[17]

//#define DEBUG 1

/* Output Arguments */

typedef struct {
  float *iqMat;
  size_t len;
  FILE *fid;
  int threadID;
  int closeFileAfterWriting;
  unsigned int instanceNo;
} saveIQStructType;

static volatile int WaitForThread1;
static volatile int WaitForThread2;

static float *iqRe, *iqIm;

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
  
  if (saveIQStruct->iqMat != NULL) {
    ret = fwrite(saveIQStruct->iqMat, sizeof(float),
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
  static FILE *fid[2];
  static mwSize lenAlloc;
  static int cnt=0;

  static saveIQStructType saveIQStruct1, saveIQStruct2;
  static char outFilename1[1024];
  static char outFilename2[1024];    

  char *mode_ptr;
  double *freqTW_ptr;
  uint16_t *lateralSubSampleIndex;
  uint16_t *rangeSubSampleIndex;
  mwSize numLateral, numRange;
  
  if ( (nrhs!=13) & (nrhs!=18)) {
    printf("Number of input arguments = %u vs 13 or 18\n", nrhs);
    mexErrMsgTxt("saveiq.c: Wrong number of arguments input\n");
    return;
  }

  if (nrhs!=18) {
    fileVersion = 2;
  }
  else {
    fileVersion = (mwSize)(*((double*)mxGetPr(fileVersionIn)));
    mode_ptr = (char*)mxGetChars(modeIn);
    freqTW_ptr = (double*)mxGetPr(freqTWIn);
    lateralSubSampleIndex = (uint16_t*)mxGetData(lateralSubSampleIndexIn);
    numLateral = mxGetN(lateralSubSampleIndexIn);
    rangeSubSampleIndex = (uint16_t*)mxGetData(rangeSubSampleIndexIn);
    numRange = mxGetN(rangeSubSampleIndexIn);
  }

  if (!mxIsComplex(IQ)) {
    mexErrMsgTxt("saveiq.c: IQ data must be complex\n");
    return;
  }
  
  if (  WaitForThread1 | WaitForThread2 ) {
    printf("*** threads (%d, %d) not done yet. frame block skipped. *** \n",
	   WaitForThread1, WaitForThread2);
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
      fclose(fid[1]);
      free(ThreadList);
      free(iqRe); free(iqIm);
      firstRun=1;
      return;
  }
   
  // avoid case where this function persists after VSX stopped and close all
  if (dopPRF < 0) {
#if DEBUG
    printf("saveiq.c: dopPRF < 0, force init\n");
    printf("WaitForThread1 = %d\n",  WaitForThread1);
    printf("WaitForThread2 = %d\n",  WaitForThread2);
    printf("closeFileAfterWrite = %d\n", closeFileAfterWrite);
#endif
    firstRun=1;
    dopPRF=-dopPRF;
  }
    
  size_t nDim = mxGetNumberOfDimensions(IQ);

  if ( ! ((nDim != 3) | (nDim != 4)) ) {
    mexErrMsgTxt(
      "IQ matrix needs to be 3D (one frame) or 4D (multiple frame) \n");
  }
    
  mwSize* dimIQPre =  mxGetDimensions(IQ);
  mwSize* dimdateNum =  mxGetDimensions(datenumIn);    
				
  mwSize dimIQ[4];
  dimIQ[3] = 1;      

  for (int i = 0; i < nDim; i++) {
    dimIQ[i] = dimIQPre[i];
  }

  mwSize len = dimIQ[0]*dimIQ[1]*dimIQ[2]*dimIQ[3];

  if (!firstRun & (len > lenAlloc)) {
    free(iqRe); free(iqIm);
    iqRe = (float*)malloc(len*sizeof(float));
    iqIm = (float*)malloc(len*sizeof(float));
#if DEBUG    
    printf("IQ data size increase from length: %d to %d, reallocated.\n",
	   lenAlloc, len);
#endif
    lenAlloc = len;
  }

#if DEBUG
  printf("DopPRF: %lf ", dopPRF);
  printf("freq: %lf ", *freq_ptr);  
  printf("IQ dimension: ");
  for (int i = 0; i < 4; i++) {
    printf("[%d]", dimIQ[i]);
  }
  printf("\n");

  printf("dateNum dimension: %d x %d\n", dimdateNum[0], dimdateNum[1]);
  
#endif
    
  if ((dimdateNum[0] != 1) | (dimdateNum[1] != dimIQ[nDim-1])) {
    printf("dateNum has %d columns, IQMat has %d frames.\n", 
	   dimdateNum[1], dimIQ[nDim-1]);
    mexErrMsgTxt("dateNum needs to be same length as number of frames.\n");
  }

  // if we get here, all threads are done from any previous calls
  // if we are closing file after writing, we need to open new files now
  if (firstRun | closeFileAfterWrite) {
    // local copies of IQ, to avoid access violation exceptions from Matlab
    iqRe = (float*)malloc(len*sizeof(float));
    iqIm = (float*)malloc(len*sizeof(float));

    if (iqRe==NULL)
      mexErrMsgTxt("saveiq: Error allocating iqRe\n");
    if (iqIm==NULL)
      mexErrMsgTxt("saveiq: Error allocating iqIm\n");
    
    mwSize strln;
    strln = mxGetN(fileSpecInRe)+1;
    mxGetString(fileSpecInRe, outFilename1, strln);
    strln = mxGetN(fileSpecInIm)+1;
    mxGetString(fileSpecInIm, outFilename2, strln);    
    
    printf("IQ output files: %s %s\n", outFilename1, outFilename2);
				
    fid[0] = fopen(outFilename1, "wb");
    fid[1] = fopen(outFilename2, "wb");    

    if ((fid[0] == NULL) | (fid[1] == NULL) ){
      printf("fid[0]: %x fid[1]: %x \n", fid[0], fid[1]);
      mexErrMsgTxt("saveiq: Could not open file for writing!\n");
    }
	
#if DEBUG
    printf("instance %u files opened: fid[0] = %x, fid[1] = %x \n",
	   instanceNo, fid[0], fid[1]);
#endif
    
    //printf("cnt = %d, len =  %d, fid1 = %u, fid2 = %u \n", cnt, len, fid[0], fid[1]);
    
    for (int i = 0; i < 2; i++) {
      fwrite(&fileVersion, sizeof(fileVersion), 1, fid[i]);
      fwrite(dimIQ, sizeof(mwSize), 4, fid[i]);
      fwrite(&PRF, sizeof(double), 1, fid[i]);
      fwrite(&dopPRF, sizeof(double), 1, fid[i]);
      fwrite(ttna_ptr, sizeof(double), 1, fid[i]);
      fwrite(nDopPRIs_ptr, sizeof(double), 1, fid[i]);    
      fwrite(datenum_ptr, sizeof(double), dimIQ[3], fid[i]);
      fwrite(origin_ptr, sizeof(double), 3, fid[i]);
      fwrite(PDelta_ptr, sizeof(double), 3, fid[i]);
      fwrite(freq_ptr, sizeof(double), 1, fid[i]);
      if (fileVersion==3) {
	fwrite(mode_ptr, sizeof(char), 1, fid[i]);
	fwrite(freqTW_ptr, sizeof(double), 1, fid[i]);
	fwrite(&numLateral, sizeof(mwSize), 1, fid[i]);
	fwrite(lateralSubSampleIndex, sizeof(uint16_t), numLateral, fid[i]);
	fwrite(&numRange, sizeof(mwSize), 1, fid[i]);
	fwrite(rangeSubSampleIndex, sizeof(uint16_t), numRange, fid[i]);
      }
    }
    
    saveIQStruct1.fid = fid[0];
    saveIQStruct1.threadID = 1;    
    saveIQStruct2.fid = fid[1];
    saveIQStruct2.threadID = 2;        
    saveIQStruct1.len = len;
    saveIQStruct2.len = len;
    saveIQStruct1.instanceNo = instanceNo;
    saveIQStruct2.instanceNo = instanceNo;    
      
  } // close file after write

  if (firstRun) {
    ThreadList = (HANDLE*)malloc(2*sizeof( HANDLE ));
    firstRun=0;
  }

  #if DEBUG
  printf("b4 copy iqRe: %x ", iqRe);
  printf("b4 copy iqIm: %x\n ", iqIm);  
  #endif

  float *realIQ = (float*)mxGetPr(IQ);
  float *imagIQ = (float*)mxGetPi(IQ);  
    
  /*  printf("tmpa: %x ", tmpa);
  printf("tmpb: %x\n ", tmpb);
  printf("tmpa[0]: %f ", tmpa[0]);
  printf("tmpa[%d]: %f\n ", len-1, tmpa[len-1]);
  printf("tmpb[0]: %f\n", tmpb[0]);
  printf("tmpb[%d]: %f\n ", len-1, tmpb[len-1]);

  printf("iqRe[0]: %f\n", iqRe[0]);
  printf("iqRe[%d]: %f\n ", len-1, iqRe[len-1]);
  printf("iqIm[0]: %f\n", iqIm[0]);
  printf("iqIm[%d]: %f\n ", len-1, iqIm[len-1]);
  */

#if DEBUG  
  printf("Will now copy %d floats\n", len);
#endif
  
  memcpy(iqRe, realIQ, len*sizeof(float));
  memcpy(iqIm, imagIQ, len*sizeof(float));  

  // why this causes access violation error during file
  // write (another thread) is a mystery, since it is
  // okay when not doing it with a var under VSX
  //saveIQStruct1.iqMat = (float*)mxGetPr(IQ);
  //saveIQStruct2.iqMat = (float*)mxGetPi(IQ);
  
  saveIQStruct1.iqMat = iqRe;
  saveIQStruct2.iqMat = iqIm;
  
  
#if DEBUG
  printf("iqMat real: %x ", saveIQStruct1.iqMat);
  printf("iqMat imag: %x\n ", saveIQStruct2.iqMat);  
#endif
  
  // closing files is done in the thread
  saveIQStruct1.closeFileAfterWriting = closeFileAfterWrite;
  saveIQStruct2.closeFileAfterWriting = closeFileAfterWrite;
 
  /*  printf("WaitForThread1 = %d, WaitForThread2 = %d \n",
      WaitForThread1, WaitForThread2);*/

  WaitForThread1 = 1;
  WaitForThread2 = 1;  
 
#if WIN
  // Reserve room for handles of threads in ThreadList
  ThreadList[0] = (HANDLE)_beginthreadex( NULL, 0, &save_iq,
					  &saveIQStruct1 , 0, NULL );
  ThreadList[1] = (HANDLE)_beginthreadex( NULL, 0, &save_iq,
					  &saveIQStruct2 , 0, NULL );
  // this does not kill the thread, only removes the handle
  CloseHandle( ThreadList[0] );
  CloseHandle( ThreadList[1] );    
#else
  save_iq(&saveIQStruct1);
  save_iq(&saveIQStruct2);
#endif    
  return;
}
