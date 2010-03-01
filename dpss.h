#ifndef _DPSS_H
#define _DPSS_H 1
/**
 * @file   dpss.h
 * @author Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:41:14 2010
 * 
 * @brief  C function to calculate discrete prolate spherical sequences.
 *
 * Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
 * Commons Attribution-Noncommercial-Share Alike 3.0 United States
 * License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/) * 
 */

/**
 * Computes the discrete prolate spherical sequences used in the
 * multitaper method power spectrum calculations.
 *
 * Inputs:
 * @param  npoints   the number of points in the window
 * @param  nw        the time-bandwidth product. Must be an integer or half-integer
 *                   (typical choices are 2, 5/2, 3, 7/2, or 4)
 * @param  k         how many DPSS vectors to return (up to npoints but k>nw*2-1 are not stable)
 *
 * Outputs:
 * @param  tapers  - k DPSS sequences in order of decreasing eigenvalue (size npoints*k)
 * @param  lambdas - k eigenvalues associated with each taper
 * @param  [outputs need to be allocated]
 *
 * @return
 *    0 for success
 *   -1 for invalid parameters (NW >= npoints/2, npoints < 0, k < 0, etc)
 *   -2 for failed eigenvalue solver (shouldn't ever happen)
 */
int dpss(double *tapers, double *lambda, int npoints, double NW, int k);


#endif /* _DPSS_H */

