
/*******************************************************************************************************
*  Copyright 2018 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  ("Alliance") under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
*  The Government retains for itself and others acting on its behalf a nonexclusive, paid-up,
*  irrevocable worldwide license in the software to reproduce, prepare derivative works, distribute
*  copies to the public, perform publicly and display publicly, and to permit others to do so.
*
*  Redistribution and use in source and binary forms, with or without modification, are permitted
*  provided that the following conditions are met:
*
*  1. Redistributions of source code must retain the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer.
*
*  2. Redistributions in binary form must reproduce the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer in the documentation and/or
*  other materials provided with the distribution.
*
*  3. The entire corresponding source code of any redistribution, with or without modification, by a
*  research entity, including but not limited to any contracting manager/operator of a United States
*  National Laboratory, any institution of higher learning, and any non-profit organization, must be
*  made publicly available under this license for as long as the redistribution is made available by
*  the research entity.
*
*  4. Redistribution of this software, without modification, must refer to the software by the same
*  designation. Redistribution of a modified version of this software (i) may not refer to the modified
*  version by the same designation, or by any confusingly similar designation, and (ii) must refer to
*  the underlying software originally provided by Alliance as "SolTrace". Except to comply with the 
*  foregoing, the term "SolTrace", or any confusingly similar designation may not be used to refer to 
*  any modified version of this software or any modified version of the underlying software originally 
*  provided by Alliance without the prior written consent of Alliance.
*
*  5. The name of the copyright holder, contributors, the United States Government, the United States
*  Department of Energy, or any of their employees may not be used to endorse or promote products
*  derived from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
*  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
*  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER,
*  CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
*  EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
*  IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
*  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************************************/


#ifndef _TREE_MESH_
#define _TREE_MESH_ 1

#include <vector>
#include <string>
using namespace std;


struct KDLayoutData
{
    double xlim[2];     //[min, max]
    double ylim[2];     //[min, max]

    double min_unit_dx;
    double min_unit_dy;

};

//-------------------------------------------------------------------------------------------------

class st_tree_node
{
	st_tree_node *m0, *m1;
	vector<void*> data;
	bool terminal;
    string address;
	
protected:
	st_tree_node *m_proc(string &key, int index);
	vector<st_tree_node*> m_get_children();
	
public:
    st_tree_node();
	void setup(st_tree_node *child0);
	void setup(st_tree_node *child0, st_tree_node *child1);
	void setup(void* data);
	bool is_terminal();
	vector<void*> *get_array();
	vector<void*> get_child_data();
    st_tree_node *get_child_node(bool upper);
    string get_address();
    void set_address( const string &addr);
};

//-------------------------------------------------------------------------------------------------

class st_opt_element : public st_tree_node
{
	double
		xr[2],
		yr[2];
    vector<void*> neighbors;

public:
	st_opt_element(){};
	void set_range(double xrlo, double xrhi, double yrlo, double yrhi);
	void set_range(double xr[2], double yr[2]);
	double *get_yr();
	double *get_xr();
	st_opt_element *process(string &key, int index);
	vector<st_opt_element*> get_children();
    void add_neighbor(void* ptr);
    vector<void*> *get_neighbor_data();
};

//-------------------------------------------------------------------------------------------------

class st_hash_tree
{
protected:
	KDLayoutData Data;
	vector<st_opt_element> nodes;
	st_opt_element head_node;
	int nx_req, ny_req; //The number of divisions required to achieve the desired resolution
	double log2inv;
	void create_node(st_opt_element &node, int index, string &binary, void *object, double *objsize=0);

public:
	st_hash_tree();
	void reset();
	bool create_mesh(KDLayoutData &Data);
	
    virtual void add_object(void *object, double locx, double locy, double *objsize=0);
    virtual string pos_to_binary_base(double x, double y);
	virtual string pos_to_binary(double x, double y, double res);
    virtual void binary_to_pos(const string &binary, double *x, double *y);
	
    void get_terminal_data(vector<vector<void*>*> &tdata);
	void get_terminal_nodes(vector<st_opt_element*> &tnodes);
	vector<st_opt_element>* get_all_nodes();
    void add_neighborhood_data(); //vector<void*> *local_dat, string &binary);
    bool get_all_data_at_loc(vector<void*> &data, double locx, double locy);
};


#endif
