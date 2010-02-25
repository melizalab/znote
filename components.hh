#ifndef COMPONENTS_H
#define COMPONENTS_H

#include "math.hh"

const int UNMARKED = -1;

void ellipse_neighborhood(coord_list &indices, int n_row, int n_col);
void neighborhood(coord_list &indices, const blitz::Array<bool,2> &nmat);

clist_vector matrix_cbranches(const imatrix &spec_thresh, const coord_list &neighborhood, 
			     imatrix &component_matrix);

void delete_component(imatrix &component_matrix, const coord_list &indices);

#endif
