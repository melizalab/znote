#ifndef _DPSS_H
#define _DPSS_H 1


/* taper generating functions */

/**
 * Computes the discrete prolate spherical sequences used in the
 * multitaper method power spectrum calculations.
 *
 * Inputs:
 *   npoints   the number of points in the window
 *   nw        the time-bandwidth product. Must be an integer or half-integer
 *             (typical choices are 2, 5/2, 3, 7/2, or 4)
 *   k         how many DPSS vectors to return (up to npoints but k>nw*2-1 are not stable)
 *
 * Outputs:
 *   tapers  - k DPSS sequences in order of decreasing eigenvalue (size npoints*k)
 *   lambdas - k eigenvalues associated with each taper
 *   [outputs need to be allocated]
 *
 * Returns:
 *    0 for success
 *   -1 for invalid parameters (NW >= npoints/2, npoints < 0, k < 0, etc)
 *   -2 for failed eigenvalue solver (shouldn't ever happen)
 */
int dpss(double *tapers, double *lambda, int npoints, double NW, int k);


#endif /* _DPSS_H */

