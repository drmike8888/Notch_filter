/****************************************************************************
 *                                                                          *
 *  Parse binary data files collected by heart monitor.  Output for gnuplot *
 *                                                                          *
 *                        Author = Mike Rosing                              *
 *                        Date = July 6, 2014                               *
 *                                                                          *
 *    Modified July 21, 2014 to output binary data after just 60 Hz low     *
 *  pass filter.  Attempts at removing baseline all failed, so use an       *
 *  algorithm that does not need baseline removed.                          *
 *    Modified Aug 1, 2014 to output binary data after 60 Hz notch filter   *
 *  applied to full data.  Keeps full 2kHz sample rate rather than down     *
 *  sampling.                                                               *
 *                                                                          *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

/*  These coefficients computed using IIR notch filter design algorithm
    of Pei and Tseng (1996).  See notch_filter.c 
 
1 coefficients are -1.9162329361  0.9507867324
3 coefficients are -1.6490444725  0.9530853152
5 coefficients are -1.1482741905  0.9535607368
7 coefficients are -0.4858941845  0.9538156135
9 coefficients are  0.2449035617  0.9540193348
11 coefficients are  0.9414624597  0.9542403314
13 coefficients are  1.5060265873  0.9545758641
15 coefficients are  1.8597673176  0.9554750804

*/
    
#define A11 (-1.9162329361)
#define A12 (0.9507867324)
#define A31 (-1.6490444725)
#define A32 (0.9530853152)
#define A51 (-1.1482741905)
#define A52 (0.9535607368)
#define A71 (-0.4858941845)
#define A72 (0.9538156135)
#define A91 (0.2449035617)
#define A92 (0.9540193348)
#define A111 (0.9414624597)
#define A112 (0.9542403314)
#define A131 (1.5060265873)
#define A132 (0.9545758641)
#define A151 (1.8597673176)
#define A152 (0.9554750804)

#define FFTSIZE 16384

#include "heart_low_pass.c"

int main(int argc, char*argv[])
{  
  FILE *in, *out;
  char filename[128], *dot;
  int i, j, k, value[9], index, chan[3];
  int numsamples, numgroups;
  double *rawdata, *rawin, *endraw;
  double *pass1, *passptr, *pass2;
  double bcoef[8][3], acoef[8][3];
  double stage[3][8][3];  // curve, harmonic, history (dimension meaning)

  fftw_plan plan, invrs;
  double *fftin, *mag, phase, *fftrvout, maxmag[32];
  fftw_complex *fftout, *fftrvin;

  in = fopen(argv[1], "r");
  if(!in)
  {
    printf("can't find file %s\n", argv[1]);
    exit(0);
  }
  strcpy(filename, argv[1]);
  dot = strstr( filename, ".data");
  strcpy( dot, ".raw");
  printf("output filename: %s\n", filename);
  out = fopen(filename, "w");
  if(!out)
  {
    printf("can't create %s\n", filename);
    exit(0);
  }

/*  assume we don't take more than 30 minutes worth of data  */

  rawdata = (double*) malloc(sizeof(double)*3*2000*30*60);
  rawin = &rawdata[3*3];  // zero out first sample

/*  Each channel uses 3 bytes.  There are 3 channels.  */

  while(!feof(in))
  {
    for(i=0; i<9; i++) value[i] = fgetc(in);
    chan[0] = (value[0] << 16) | (value[1] << 8) | value[2];
    if(chan[0] & 0x800000) chan[0] |= 0xff000000;
    chan[1] = (value[3] << 16) | (value[4] << 8) | value[5];
    if(chan[1] & 0x800000) chan[1] |= 0xff000000;
    chan[2] = (value[6] << 16) | (value[7] << 8) | value[8];
    if(chan[2] & 0x800000) chan[2] |= 0xff000000;
    *rawin = chan[0];
    rawin++;
    *rawin = chan[1];
    rawin++;
    *rawin = chan[2];
    rawin++;
  }
  fclose(in);
  for(i=0; i<9; i++) // zero fill initial data
    rawdata[i] = 0.0;
  endraw = rawin - 9; // ignore last values collected
  for(i=0; i<9; i++)
    endraw[i] = 0.0;
  endraw += 9;
  printf("data read in\n");
  rawin = rawdata;
  while (rawin < endraw)
  {
    fprintf(out, "%lf %lf %lf\n", rawin[0], rawin[1], rawin[2]);
    rawin += 3;
  }
  fclose(out);
  numsamples = (endraw - rawdata)/3;
  printf("saving %d samples\n", numsamples);


/*  compute IIR notch filter
    one notch for each odd harmonic of 60 Hz out to 900 Hz */

  dot = strstr( filename, ".raw");
  strcpy( dot, ".plt_ntch");
  printf("output filename: %s\n", filename);
  out = fopen(filename, "w");
  if(!out)
  {
    printf("can't create %s\n", filename);
    exit(0);
  }
  pass1 = (double*)malloc(sizeof(double)*(3*numsamples + 6*FILTER_TAP_NUM));
  rawin = rawdata + 3*3;
  passptr = pass1 + 3*(FILTER_TAP_NUM - 3);

// notch filter

  bcoef[0][0] = (1.0 + A12)/2.0;
  bcoef[0][1] = A11;
  bcoef[0][2] = (1.0 + A12)/2.0;
  acoef[0][0] = 1.0;
  acoef[0][1] = A11;
  acoef[0][2] = A12;

  bcoef[1][0] = (1.0 + A32)/2.0;
  bcoef[1][1] = A31;
  bcoef[1][2] = (1.0 + A32)/2.0;
  acoef[1][0] = 1.0;
  acoef[1][1] = A31;
  acoef[1][2] = A32;

  bcoef[2][0] = (1.0 + A52)/2.0;
  bcoef[2][1] = A51;
  bcoef[2][2] = (1.0 + A52)/2.0;
  acoef[2][0] = 1.0;
  acoef[2][1] = A51;
  acoef[2][2] = A52;

  bcoef[3][0] = (1.0 + A72)/2.0;
  bcoef[3][1] = A71;
  bcoef[3][2] = (1.0 + A72)/2.0;
  acoef[3][0] = 1.0;
  acoef[3][1] = A71;
  acoef[3][2] = A72;

  bcoef[4][0] = (1.0 + A92)/2.0;
  bcoef[4][1] = A91;
  bcoef[4][2] = (1.0 + A92)/2.0;
  acoef[4][0] = 1.0;
  acoef[4][1] = A91;
  acoef[4][2] = A92;

  bcoef[5][0] = (1.0 + A112)/2.0;
  bcoef[5][1] = A111;
  bcoef[5][2] = (1.0 + A112)/2.0;
  acoef[5][0] = 1.0;
  acoef[5][1] = A111;
  acoef[5][2] = A12;

  bcoef[6][0] = (1.0 + A132)/2.0;
  bcoef[6][1] = A131;
  bcoef[6][2] = (1.0 + A132)/2.0;
  acoef[6][0] = 1.0;
  acoef[6][1] = A131;
  acoef[6][2] = A132;

  bcoef[7][0] = (1.0 + A152)/2.0;
  bcoef[7][1] = A151;
  bcoef[7][2] = (1.0 + A152)/2.0;
  acoef[7][0] = 1.0;
  acoef[7][1] = A151;
  acoef[7][2] = A152;

  for(k=0; k<3; k++)
    for(i=0; i<8; i++)
      for(j=0; j<3; j++)
	stage[k][i][j] = 0.0;

  while(rawin < endraw)
  {
    for(j=0; j<3; j++)  // loop over each curve
    {

      stage[j][0][2] = 0.0;
      for(i=0; i<3; i++)
	stage[j][0][2] += rawin[-3*i + j]*bcoef[0][i];
      stage[j][0][2] -= acoef[0][1]*stage[j][0][1] + acoef[0][2]*stage[j][0][0];
      for(k=1; k<7; k++) // loop over each harmonic
      {
	stage[j][k][2] = 0.0;
	for(i=0; i<3; i++)
	  stage[j][k][2] += stage[j][k-1][2-i]*bcoef[k][i];
	stage[j][k][2] -= acoef[k][1]*stage[j][k][1] + acoef[k][2]*stage[j][k][0];
      }
      passptr[j] = 0.0;
      for(i=0; i<3; i++)
	passptr[j] += stage[j][6][2-i]*bcoef[7][i];
      passptr[j] -= acoef[7][1]*passptr[j-3] + acoef[7][2]*passptr[j-6];
      for(k=0; k<7; k++)
      {
	stage[j][k][0] = stage[j][k][1];
	stage[j][k][1] = stage[j][k][2];
      }
/*
      passptr[j] = 0.0;
      for(i=0; i<3; i++)
	passptr[j] += rawin[-3*i + j]*bcoef[0][i];
      passptr[j] -= acoef[0][1]*passptr[j-3] + acoef[0][2]*passptr[j-6];
*/      
    }
    passptr += 3;
    rawin += 3;
  }
  printf("done with pass 1\n");

/*  perform FIR low pass on remaining data  */

  endraw = passptr;
  pass2 = (double*)malloc(sizeof(double)*(3*numsamples + 3*FILTER_TAP_NUM));
  rawin = pass1;
  passptr = pass2;
  while(rawin < endraw + 3*FILTER_TAP_NUM)
  {
    for(j=0; j<3; j++)
    {
      passptr[j] = 0.0;
      for(i=0; i<FILTER_TAP_NUM; i++)
	passptr[j] += rawin[3*i + j]*filter_taps[i];
    }
    passptr += 3;
    rawin += 3;
  }
  printf("done with pass 2\n");

/*  save data for gnuplot  */

  rawin = pass2 + 6*FILTER_TAP_NUM;
  for(i=0; i<numsamples - 3*FILTER_TAP_NUM; i++)
  {
    fprintf(out, "%lf %lf %lf\n", rawin[0] , rawin[1], rawin[2]);
    rawin += 3;
  }
  fclose(out);

/*  save data for further processing  */

  dot = strstr( filename, ".plt_ntch");
  strcpy( dot, ".bin_ntch");
  printf("output filename: %s\n", filename);
  out = fopen(filename, "w");
  if(!out)
  {
    printf("can't create %s\n", filename);
    exit(0);
  }
  fwrite(pass2 + 6*FILTER_TAP_NUM, sizeof(double), (numsamples - FILTER_TAP_NUM)*3, out);
  fclose(out);
  exit(0);    // skip FFT's for now

/*  compute 1D FFT's on each wave form  */

  fftin = fftw_alloc_real(FFTSIZE);
  fftout = fftw_alloc_complex(FFTSIZE/2);
  mag = (double*)malloc(sizeof(double)*FFTSIZE/2);
  plan = fftw_plan_dft_r2c_1d(FFTSIZE, fftin, fftout, FFTW_MEASURE);
  rawin = pass1;
  k = 0;
  while(rawin < endraw - 3*FFTSIZE)
  {
    for(i=0; i<FFTSIZE; i++) fftin[i] = rawin[3*i + 2];
    fftw_execute(plan);
    strcpy(filename, argv[1]);
    dot = strstr( filename, ".data");
    sprintf( dot, "_%03d_fft.plt", k);
    printf("output filename: %s\n", filename);
    out = fopen(filename, "w");
    if(!out)
    {
      printf("can't create %s\n", filename);
      exit(0);
    }
    for(i=0; i<FFTSIZE/2; i++)
    {
      mag[i] = sqrt(fftout[i][0]*fftout[i][0] + fftout[i][1]*fftout[i][1]);
      phase = atan2(fftout[i][1], fftout[i][0]);
      fprintf(out, "%lg %lg\n", mag[i], phase);
    }
    fclose(out);
    k++;
    rawin += 3*FFTSIZE;
  }
  fftw_free(fftin);
  fftw_free(fftout);
  free(mag);
  free(pass1);
  free(pass2);
  free(rawdata);
}
