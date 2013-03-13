/**
 * @file   components.cc
 * @author Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:38:31 2010
 *
 * Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
 * Commons Attribution-Noncommercial-Share Alike 3.0 United States
 * License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).
 *
 */

#include "components.hh"
#include <set>

// try to determine if we need to include a 0.9 header
#if BZ_MAJOR_VERSION == 0 && BZ_MINOR_VERSION < 10
#include <blitz/tinyvec-et.h>
#endif

void
ellipse_neighborhood(coord_list &indices, int n_row, int n_col)
{
	blitz::firstIndex i;
	blitz::secondIndex j;
	double x(n_row), y(n_col);
	blitz::Array<bool,2> nmat(2 * n_row + 1, 2 * n_col + 1);
	nmat = (blitz::sqr((i - x)/x) + blitz::sqr((j - y)/y)) <= 1.0;

	neighborhood(indices, nmat);
}

void
neighborhood(coord_list &indices, const blitz::Array<bool,2> &nmat)
{
	blitz::find(indices, nmat);
	indices -= blitz::TinyVector<int,2>(nmat.rows() / 2, nmat.cols() / 2);
}

clist_vector
matrix_cbranches(const imatrix &I, const coord_list &N, imatrix &M)
{
	int nbranches = 0;
	coord_vector active_stack, node_stack;
	clist_vector branch_list;

	M.resize(I.shape());
	M = UNMARKED;
	for (int c = 0; c < I.cols(); c++) {
		for (int r = 0; r < I.rows(); r++) {
			coord start(r,c);
			if (I(start) && M(start) == UNMARKED) {
				M(start) = nbranches;
				active_stack.push_back(start);
				node_stack.clear();
				node_stack.push_back(start);
				// search for connected neighbors
				while (!active_stack.empty()) {
					coord curr(active_stack.back());
					active_stack.pop_back();
					int n_added = 0;
					for (size_t ni = 0; ni < N.size(); ni++) {
						coord search = curr + N(ni);
						if (inrange(I,search) && I(search) &&
						    M(search)==UNMARKED && blitz::all(curr!=search)) {
							M(search) = M(curr);
							active_stack.push_back(search);
							node_stack.push_back(search);
							n_added += 1;
						}
					}
				}
				// done searching for this component
				nbranches += 1;
				coord_list branch(node_stack.size());
				std::copy(node_stack.begin(), node_stack.end(), branch.begin());
				branch_list.push_back(branch);
			}
		}
	}
	return branch_list;
}

clist_vector
get_components(const imatrix &labels)
{
	int max_index = blitz::max(labels);
	std::vector<coord_vector> branchv(max_index+1);
	for (int r = 0; r < labels.rows(); r++) {
		for (int c = 0; c < labels.cols(); c++) {
			coord p(r,c);
			if (labels(p) != UNMARKED)
				branchv[labels(p)].push_back(p);
		}
	}

	clist_vector branches;
	for (int i = 0; i <= max_index; i++) {
		coord_list branch(branchv[i].size());
		std::copy(branchv[i].begin(), branchv[i].end(), branch.begin());
		branches.push_back(branch);
	}
	return branches;
}

void
delete_component(imatrix &M, const coord_list &branch)
{
	coord_list::const_iterator it;
	for (it = branch.begin(); it != branch.end(); it++)
		M(*it) = UNMARKED;
}

void
renumber_components(imatrix &M)
{
	int i;
	std::set<int> unique_vals(M.begin(),M.end());
	unique_vals.erase(-1);
	std::set<int>::const_reverse_iterator rit = unique_vals.rbegin();

 	ivector index(*rit+1);
 	for (i=unique_vals.size()-1; rit != unique_vals.rend(); --i, ++rit)
 		index(*rit) = i;

	imatrix::iterator it = M.begin();
	for (; it != M.end(); it++) {
		if (*it > -1)
			*it = index(*it);
	}
}

void
make_mask(const coord_list &component, const dmatrix &mask, const coord &mask_center, dmatrix &out)
{
	int j,k;
	coord_list::const_iterator it = component.begin();
	for (; it !=component.end(); it++) {
		for (j = 0; j < mask.rows(); j++) {
			for (k = 0; k < mask.cols(); k++) {
				coord tgt(*it + coord(j,k));
				if (inrange(out, tgt))
					out(tgt) = fmax(out(tgt),mask(j,k));
			}
		}
	}
}

void
mask_sum(const clist_vector &components, const dmatrix &mask, const coord &mask_center, dmatrix &out)
{
	out = 0;
	dmatrix temp(out.shape());
	for (unsigned int i = 0; i < components.size(); i++) {
		temp = 0;
		make_mask(components[i], mask, mask_center, temp);
		out += temp;
	}
}
