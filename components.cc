#include "components.hh"
#include "blitz_io.hh"
#include <blitz/tinyvec-et.h>

using namespace std;

void
ellipse_neighborhood(coord_list &indices, int n_row, int n_col)
{
	double x(n_row), y(n_col);
	blitz::Array<bool,2> nmat(2 * n_row + 1, 2 * n_col + 1);
	nmat = (blitz::sqr((blitz::tensor::i - x)/x) + blitz::sqr((blitz::tensor::j - y)/y)) <= 1.0;

	write_bin("nhd.bin",nmat);
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
	for (int r = 0; r < I.rows(); r++) {
		for (int c = 0; c < I.cols(); c++) {
			coord start(r,c);
			if (I(start) && M(start) == UNMARKED) {
				M(start) = nbranches;
				active_stack.push_back(start);
				node_stack.push_back(start);
				// search for connected neighbors
				while (!active_stack.empty()) {
					coord curr(active_stack.back());
					active_stack.pop_back();
					int n_added = 0;
					for (int ni = 0; ni < N.size(); ni++) {
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
