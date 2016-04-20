/***************************************************************
 *                                                             *
 *   test to see if iir notch filters might work on 60Hz at    *
 *   full sample rate.                                         *
 *                                                             *
 *             Author = Mike Rosing                            *
 *             Date = July 31, 2014                            *
 *                                                             *
 **************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main()
{
  int i, j;
  double omega, acoef[3], bcoef[3], hnr, hni, hdr, hdi, mag;
  double c1, c2, f, s1, s2;
  double omegan1, omegabw1;
  double omega1, omega2, theta1, theta2, beta1, beta2;
  double p[2], q[2][2], det, qinv[2][2];
  double am[2], tmp;

  for(j=1; j<16; j+=2)
  {
    omegan1 = 2.0*M_PI*60.0*j/2000.0;
    omegabw1 = 2.0*M_PI*15.0/2000.0;
    omega1 = omegan1 - omegabw1/2.0;
    omega2 = omegan1;
    theta1 = -M_PI/2.0;
    theta2 = -M_PI;
    beta1 = theta1/2.0 + omega1;
    beta2 = theta2/2.0 + omega2;
    p[0] = tan(beta1);
    p[1] = tan(beta2);
    q[0][0] = sin(omega1) - p[0]*cos(omega1);
    q[0][1] = sin(2.0*omega1) - p[0]*cos(2.0*omega1);
    q[1][0] = sin(omega2) - p[1]*cos(omega2);
    q[1][1] = sin(2.0*omega2) - p[1]*cos(2.0*omega2);
    det = q[0][0]*q[1][1] - q[0][1]*q[1][0];
    qinv[0][0] = q[1][1]/det;
    qinv[0][1] = -q[0][1]/det;
    qinv[1][0] = -q[1][0]/det;
    qinv[1][1] = q[0][0]/det;
    am[0] = qinv[0][0]*p[0] + qinv[0][1]*p[1];
    am[1] = qinv[1][0]*p[0] + qinv[1][1]*p[1];

    printf("%d coefficients are %13.10lf %13.10lf\n", j, am[0], am[1]);
  }
  exit(0);

  acoef[0] = 1.0 + am[1];
  acoef[1] = 2.0*am[0];
  acoef[2] = 1.0 + am[1];
  bcoef[0] = 1.0;
  bcoef[1] = am[0];
  bcoef[2] = am[1];
  for(i=1; i<250; i++)
  {
    f = i;
    omega = 2.0*M_PI*f/2000.0;
    c1 = cos(omega);
    c2 = cos(2.0*omega);
    s1 = sin(omega);
    s2 = sin(2.0*omega);
    hnr = acoef[0] + acoef[1]*c1 + acoef[2]*c2;
    hni = acoef[1]*s1 + acoef[2]*s2;
    hdr = bcoef[0] + bcoef[1]*c1 + bcoef[2]*c2;
    hdi = bcoef[1]*s1 + bcoef[2]*s2;
    mag = sqrt((hnr*hnr + hni*hni)/(hdr*hdr + hdi*hdi))/2.0;
    printf("%lf %lf\n", f, mag);
  }
}
