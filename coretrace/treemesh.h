#ifndef _TREE_MESH_
#define _TREE_MESH_ 1

#include <vector>
#include <string>
using namespace std;


struct LayoutData
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
	LayoutData *Data;
	vector<st_opt_element> nodes;
	st_opt_element head_node;
	int nx_req, ny_req; //The number of divisions required to achieve the desired resolution
	double log2inv;
	void create_node(st_opt_element &node, int index, string &binary, void *object, double *objsize=0);

public:
	st_hash_tree();
	void reset();
	bool create_mesh(LayoutData *Data);
	
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
