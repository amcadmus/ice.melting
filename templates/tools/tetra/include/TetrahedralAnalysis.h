#pragma once

#include <vector>
#include <cassert>

using namespace std;

class TetrahedralAnalysis 
{
public:
  TetrahedralAnalysis (const double rc);
  void compute (vector<double > & qq,
		const CellList & clist,
		const vector<double> & box,
		const vector<vector<double > > & com) const;
private:
  double rc, rc2;
}
    ;

TetrahedralAnalysis::
TetrahedralAnalysis (const double rc_) 
    : rc (rc_) , rc2(rc_*rc_)
{
}

void
TetrahedralAnalysis::
compute (vector<double > & qq,
	 const CellList & clist,
	 const vector<double> & box,
	 const vector<vector<double > > & com) const
{
  // for (unsigned ii = 0; ii < com.size(); ++ii){
  //   cout << com[ii][0] << " " 
  // 	 << com[ii][1] << " " 
  // 	 << com[ii][2] << " "  << endl;
  // }

  double rup = rc;
  int xiter = rup / clist.getCellSize().x;
  if (xiter * clist.getCellSize().x < rup) xiter ++;
  int yiter = rup / clist.getCellSize().y;
  if (yiter * clist.getCellSize().y < rup) yiter ++;
  int ziter = rup / clist.getCellSize().z;
  if (ziter * clist.getCellSize().z < rup) ziter ++;
  assert (xiter * clist.getCellSize().x >= rup);
  assert (yiter * clist.getCellSize().y >= rup);
  assert (ziter * clist.getCellSize().z >= rup);

  vector<vector<int > > nlist4 (com.size());
  vector<vector<pair<double,int> > > i_nei_info(com.size());

  IntVectorType nCell = clist.getNumCell();
  unsigned cellIndexUpper = unsigned(nCell.x * nCell.y * nCell.z);
  for (unsigned iCellIndex = 0;
       iCellIndex < cellIndexUpper;
       iCellIndex += 1){
    const vector<unsigned> & iCellList (clist.getList()[iCellIndex]);
    vector<unsigned > neighborCellIndex =
	clist.neighboringCellIndex (iCellIndex, IntVectorType (xiter, yiter, ziter));
    // loop of all neighboring cells of i
    for (unsigned iNeighborCellIndex = 0;
	 iNeighborCellIndex < neighborCellIndex.size();
	 ++iNeighborCellIndex){
      unsigned jCellIndex = neighborCellIndex[iNeighborCellIndex];
      const vector<unsigned> & jCellList (clist.getList()[jCellIndex]);
      bool sameCell (iCellIndex == jCellIndex);
      for (unsigned ii = 0; ii < iCellList.size(); ++ii){
	int i_index = iCellList[ii];
	for (unsigned jj = 0; jj < jCellList.size(); ++jj){
	  if (sameCell && ii == jj) continue;	    
	  int j_index = jCellList[jj];	    
	  vector<double > diff(3) ;
	  for (int dd = 0; dd < 3; ++dd) diff[dd] = com[i_index][dd] - com[j_index][dd];
	  for (int dd = 0; dd < 3; ++dd){
	    if      (diff[dd] < -.5 * box[dd]) diff[dd] += box[dd];
	    else if (diff[dd] >= .5 * box[dd]) diff[dd] -= box[dd];
	  }
	  double r2(0);
	  for (int dd = 0; dd < 3; ++dd) r2 += diff[dd] * diff[dd];
	  if (r2 < rc2) {
	    i_nei_info[i_index].push_back (pair<double, int> (r2, j_index));
	  }
	}
      }
    }
  }

  for (unsigned ii = 0; ii < i_nei_info.size(); ++ii){
    if (i_nei_info[ii].size() < 4){
      cout << i_nei_info[ii].size() << endl;
      cerr << ii << " has " << i_nei_info[ii].size() 
	   << " neighbors, less than 4, the rc is too small, exit" << endl;
      exit(1);
    }
    sort (i_nei_info[ii].begin(), i_nei_info[ii].end());
    // cout << ii << " " ;
    for (int dd = 0; dd < 4; ++dd){
      nlist4[ii].push_back (i_nei_info[ii][dd].second);
      // if (ii == 0) cout << sqrt(i_nei_info[ii][dd].first) << " " << i_nei_info[ii][dd].second << "    " ;
    }
    // if (ii == 0) cout << endl;
  }

  qq.resize (com.size());
  assert (nlist4.size() == com.size());
  for (unsigned ii = 0; ii < nlist4.size(); ++ii){
    assert (nlist4[ii].size() == 4);
    double sum = 0;
    for (unsigned jj = 0; jj < nlist4[ii].size() - 1; ++jj){
      for (unsigned kk = jj+1; kk < nlist4[ii].size(); ++kk){
	vector<double > r0(3), r1(3);
	for (int dd = 0; dd < 3; ++dd) r0[dd] = com[nlist4[ii][jj]][dd] - com[ii][dd];
	for (int dd = 0; dd < 3; ++dd){
	  if      (r0[dd] < -.5 * box[dd]) r0[dd] += box[dd];
	  else if (r0[dd] >= .5 * box[dd]) r0[dd] -= box[dd];
	}
	for (int dd = 0; dd < 3; ++dd) r1[dd] = com[nlist4[ii][kk]][dd] - com[ii][dd];
	for (int dd = 0; dd < 3; ++dd){
	  if      (r1[dd] < -.5 * box[dd]) r1[dd] += box[dd];
	  else if (r1[dd] >= .5 * box[dd]) r1[dd] -= box[dd];
	}
	double dr0(0), dr1(0), r0r1(0);
	for (unsigned dd = 0; dd < 3; ++dd){
	  dr0 += r0[dd] * r0[dd];
	  dr1 += r1[dd] * r1[dd];
	  r0r1 += r0[dd] * r1[dd];
	}
	double cosv = r0r1 / (sqrt(dr0) * sqrt(dr1));
	sum += (cosv + 1./3.) * (cosv + 1./3.);	
	// sum += (acos (cosv) - acos (-1./3.)) * (acos (cosv) - acos (-1./3.));
      }
    }
    qq[ii] = 1. - 3./8. * sum;
    // qq[ii] = 1. - 1./1.8165 * sqrt(sum);
  }
}


