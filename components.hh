#ifndef COMPONENTS_H
#define COMPONENTS_H

#include "common.hh"

const int UNMARKED = -1;

void ellipse_neighborhood(coord_list &indices, int n_row, int n_col);
void neighborhood(coord_list &indices, const blitz::Array<bool,2> &nmat);

clist_vector matrix_cbranches(const imatrix &spec_thresh, const coord_list &neighborhood, 
			     imatrix &component_matrix);

clist_vector get_components(const imatrix &labels);

void delete_component(imatrix &component_matrix, const coord_list &indices);
void renumber_components(imatrix &component_matrix);

void make_mask(const coord_list &component, const dmatrix &mask, const coord &mask_center, dmatrix &out);
void mask_sum(const clist_vector &components, const dmatrix &mask, const coord &mask_center, dmatrix &out);

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

#endif
