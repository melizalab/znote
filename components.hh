#ifndef COMPONENTS_H
#define COMPONENTS_H
/**
 * @file   components.hh
 * @author Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:33:47 2010
 * 
 * @brief  Functions to find and manipulated connected components
 * 
 * Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
 * Commons Attribution-Noncommercial-Share Alike 3.0 United States
 * License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).
 */
#include "common.hh"
#include <algorithm>

const int UNMARKED = -1;	/**< value of unmarked points in the label array */

/** 
 * Generate a list of coordinates in an ellipse around the center of a
 * rectangular array.
 * 
 * @param indices output vector of 2-D coordinates
 * @param n_row Number of rows in the array
 * @param n_col Number of columns in the array
 */
void ellipse_neighborhood(coord_list &indices, int n_row, int n_col);

/** 
 * Convert a rectangular array of search points into a list of
 * coordinates.
 * 
 * @param indices output vector of 2-D coordinates where <nmat> is
 * positive, relative to the center of the array
 * @param nmat input arrayl
 */
void neighborhood(coord_list &indices, const blitz::Array<bool,2> &nmat);

/** 
 * Calculate the connected components in a thresholded spectrogram.
 * Regarding a matrix I as an image, find all the connected components
 * in its non-zero region, mark all the non-zero pixels by the
 * components they belong to and all the zero pixels by -1, and return
 * the marks in another matrix M which has the same size as I.  A
 * connected component C of a set V of pixels is a maximal subset of
 * the latter, such that for any u, v in C, there are w(0), .... w(n)
 * in C, with u=w(0), v=w(n), and w(k) a neighbor of w(k-1) for k=1,
 * ...n.
 *
 * The neighborhood structure is specified by a list of
 * coordinates. Two pixels u=(u(0), u(1)) and v=(v(0), v(1)) are
 * neighbors if u-v is in <nbhd>.  The image is NOT considered a
 * torus, so there is no wrapping at the boundary.
 * 
 * @param spec_thresh  2-D boolean array indicating which points should be included in the search
 * @param neighborhood  List of points considered to be in the neighborhood of the current point
 * @param component_matrix  Output integer array indicating which points are in which components.
 * 
 * @return a vector of coordinate lists; each element of the vector is
 * a list of the points in the corresponding component.
 */
clist_vector matrix_cbranches(const imatrix &spec_thresh, const coord_list &neighborhood, 
			     imatrix &component_matrix);

/** 
 * Convert a label matrix into a collection of coordinate
 * lists. Letting <labels> be a 2D integer array L, each point L(i,j)
 * >= 0 belongs to the component L(i,j)
 * 
 * @param labels  Matrix of labels
 * 
 * @return a vector of coordinate lists; each element of the vector is
 * a list of the points in the corresponding elements.  If there are
 * empty elements (i.e. if the list of unique values in L is not
 * contiguous), there will be empty elements in this vector.
 */
clist_vector get_components(const imatrix &labels);

/** 
 * Remove a component from a label matrix by setting all the
 * coordinates to unmarked.
 * 
 * @param component_matrix  Matrix of labels
 * @param indices  List of coordinates in the element to be deleted.
 */
void delete_component(imatrix &component_matrix, const coord_list &indices);

/** 
 * Renumber the components of a matrix so that the unique values of
 * the matrix are contiguous.
 * 
 * @param component_matrix The matrix to be renumbered.
 */
void renumber_components(imatrix &component_matrix);

/** 
 * Calculate a mask from a component.  The mask is a matrix (W) with
 * reference point (ci,cj). Each M(i,j) set as follows.  For each
 * (k,l) in C, the weight it assigns to (i,j) is W(i-k+ci, j-l+cj) if
 * (i-k+ci, j-l+cj) is in the range allowed by the dimension of W.
 * Then M(i, j) is the maximum among the weights over all (k,l) with
 * I(k,l)==ID.  For (i,j) that do not get any assigned weights,
 * M(i,j)=0.
 * 
 * @param component C, list of component coordinates
 * @param mask W, mask weights
 * @param mask_center (ci,cj) the reference point of the mask
 * @param out M, output mask.
 */
void make_mask(const coord_list &component, const dmatrix &mask, const coord &mask_center, dmatrix &out);

/** 
 * Calculate sum of the masks for all the components.
 * 
 * @param components vector of coordinate lists, defining all the components
 * @param mask W, mask weights
 * @param mask_center (ci,cj)
 * @param out M, output mask
 * @see make_mask
 */
void mask_sum(const clist_vector &components, const dmatrix &mask, const coord &mask_center, dmatrix &out);

/** 
 * Calculate the bounds of a component. Given a list of coordinates,
 * returns a rectangular domain defined by the minimum values in each dimension.
 * 
 * @param coord_list list of coordinates, C_i(n_1,n_2,...)
 * 
 * @return RectDomain where lbound = min(C(n_1)), min(C(n_2),... and ubound = max(C(n_1))...
 */
template <int N_rank> inline
blitz::RectDomain<N_rank> component_bounds(const blitz::Array<blitz::TinyVector<int,N_rank>,1> &coord_list) 
{
	blitz::TinyVector<int,N_rank> lbound, ubound;
	for (int i = 0; i < N_rank; i++) {
		lbound(i) = blitz::min(coord_list.extractComponent(int(),i,N_rank));
		ubound(i) = blitz::max(coord_list.extractComponent(int(),i,N_rank));
	}
	return blitz::RectDomain<N_rank>(lbound, ubound);
}


/** 
 * Compute quantile of cumulative power in an array of any rank, such
 * that all the points below the output value constitute <frac> of the
 * total power.
 * 
 * @param data power values, in an N-D array
 * @param frac the quantile to compute
 * 
 * @return the power at the given quantile <frac>
 */
template <typename T, int N_rank>
T frac_thresh(const blitz::Array<T,N_rank> &data, double frac)
{
	unsigned int i;
	std::vector<T> val(data.size());
	std::copy(data.begin(), data.end(), val.begin());
 	std::sort(val.begin(), val.end());
	double sum = 0, max = (double)blitz::sum(data);
	for (i = 0; (i < val.size()) && (sum < frac); i++)
		sum += (double)val[i] / max;
	return val[i];
}

#endif
