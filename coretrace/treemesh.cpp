
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


#include "treemesh.h"
#include <math.h>
#include <algorithm>

#include <unordered_map>
#include <set>

#include <stdio.h>



//-------------------------------------------------------------------------------------------------

inline bool binary_add(const string &key, string &modkey, int val[2])
{
    /* 
    key is:
    "x0-y0-x1-y1-x2-y2..."
    without the dashes

    val is:
    <adder in x (-1,0,1)>,<adder in y>
    
    modkey is the sum of key and val

    Binary addition. The summed number will necessarily be of the same length as the value provided because of
    the resolution requirements for the binary tree.

    If this is adapted for some other purpose, the summed number should probably be shifted so that there's
    room for extra leading digits.
    */

    int lx=0;
    int ly=0;
    for(int i=0; i<key.size(); i+=2)
    {
        if(key.at(i) == 'x')
            break;
        lx++;
    }
    for(int i=1; i<key.size(); i+=2)
    {
        if(key.at(i) == 'x')
            break;
        ly++;
    }

    modkey = key;  //copy

    //deal with x
    int xadd = val[0];
    int c=0;
    bool first=true;
    for(int i=2*(lx-1); i>-1; i-=2)
    {
        int n = (int)(key.at(i) == '1');
        if(first)
            n += xadd;
        else
            n += c;
        first = false;

        if( n > 1 )
        {
            modkey.at(i) = '0';
            c = 1;
        }
        else if(n < 0)
        {
            modkey.at(i) = '1';
            c = -1;
        }
        else 
        {
            modkey.at(i) = n == 0 ? '0' : '1';   //0 or 1    
            c=0;
            break;
        }

        if( i==0 && c != 0 )
        {
            //out of bounds
            return false;
        }
    }

    //deal with y
    int yadd = val[1];
    c=0;
    first=true;
    for(int i=2*(ly-1)+1; i>0; i-=2)
    {
        int n = (int)(key.at(i) == '1');
        if(first)
            n += yadd;
        else
            n += c;
        first = false;

        if( n > 1 )
        {
            modkey.at(i) = '0';
            c = 1;
        }
        else if(n < 0)
        {
            modkey.at(i) = '1';
            c = -1;
        }
        else 
        {
            modkey[i] = n == 0 ? '0' : '1';   //0 or 1    
            c=0;
            break;
        }

        if( i==0 && c != 0 )
        {
            //out of bounds
            return false;
        }
    }

    return true;
}

st_tree_node::st_tree_node()
{
    //Initialize
    m0 = m1 = 0;
    terminal = false;
}

void st_tree_node::setup(st_tree_node *child0){
    //Set both children equal to specified node. Used for 'x' keys.
	terminal = false;
	m0 = child0;
	m1 = m0;
}
void st_tree_node::setup(st_tree_node *child0, st_tree_node *child1){
    //Set children nodes equal to node pointer values (as applicable). For ptrs to 0, don't set. 
	terminal = false;
	if(child0 != 0)
        m0 = child0;
	if(child1 != 0)
        m1 = child1;
}
void st_tree_node::setup(void* Data)
{
    //Set node as terminal and add data
	terminal = true;
    if(Data!=0)     //Don't try to add a null value
	    data.push_back(Data);
}
bool st_tree_node::is_terminal(){
	return terminal;
}
vector<void*> *st_tree_node::get_array(){
	return &data;
}
st_tree_node *st_tree_node::m_proc(string &key, int index){
    /* 
    ** Recursive **

    Process 'key' string at current index. Traverses tree structure. If no node is found at the specified index,
    return 0.

    */

	char c;
	try
	{
		c = key.at(index);
	}
	catch(...){
		return this;
	}
    if(terminal)
        return this;
    switch (c)
    {
    case 't':
        return this;
    case 'x':
    case '0':
        if(m0 != 0)
            return m0->m_proc(key, index+1);
        else
            return 0;
    case '1':
        if(m1 != 0)
            return m1->m_proc(key, index+1);
        else
            return 0;
    }

	return 0;
}
vector<st_tree_node*> st_tree_node::m_get_children(){
	vector<st_tree_node*> kids;
	if(! terminal){
		if( m0 == m1 ){
			kids.push_back(m0);
			vector<st_tree_node*> m0kids = m0->m_get_children();
			for(unsigned int i=0; i<m0kids.size(); i++)
				kids.push_back(m0kids.at(i));
		}
		else{
			kids.push_back(m0);
			kids.push_back(m1);
			vector<st_tree_node*> m0kids = m0->m_get_children();
			vector<st_tree_node*> m1kids = m1->m_get_children();
			for(unsigned int i=0; i<m0kids.size(); i++)
				kids.push_back(m0kids.at(i));
			for(unsigned int i=0; i<m1kids.size(); i++)
				kids.push_back(m1kids.at(i));
		}
	}
	return kids;
}
vector<void*> st_tree_node::get_child_data(){
	if(terminal){
		return data;
	}
	else{
		if(m0 == m1){
			return m0->get_child_data();
		}
		else{
			vector<void*> m0dat, m1dat, alldat;
			m0dat = m0->get_child_data();
			m1dat = m1->get_child_data();
			for(unsigned int i=0; i<m0dat.size(); i++)
				alldat.push_back(m0dat.at(i));
			for(unsigned int i=0; i<m1dat.size(); i++)
				alldat.push_back(m1dat.at(i));
			return alldat;
		}
	}

}

st_tree_node* st_tree_node::get_child_node(bool upper)
{
    return upper ? m1 : m0;
}

string st_tree_node::get_address()
{
    return address;
}

void st_tree_node::set_address(const string &addr)
{
    address = addr;
}


//-------------------------------------------------------------------------------------------------
void st_opt_element::set_range(double xrlo, double xrhi, double yrlo, double yrhi){
	xr[0] = xrlo;
	xr[1] = xrhi;
	yr[0] = yrlo;
	yr[1] = yrhi;
}

void st_opt_element::set_range(double *xri, double *yri){	
	for(int i=0; i<2; i++){
		xr[i] = xri[i];
		yr[i] = yri[i];
	}
}

st_opt_element *st_opt_element::process(string &key, int index){
	return (st_opt_element*)m_proc(key, index);
}

vector<st_opt_element*> st_opt_element::get_children(){
	vector<st_opt_element*> children;
	vector<st_tree_node*> m_children = m_get_children();
	for( vector<st_tree_node*>::iterator it = m_children.begin(); it != m_children.end(); it++)
		children.push_back( (st_opt_element*) *it );
	return children;
}

double *st_opt_element::get_yr(){return yr;}
double *st_opt_element::get_xr(){return xr;}

void st_opt_element::add_neighbor(void *ptr)
{
    neighbors.push_back(ptr);
}

vector<void*>* st_opt_element::get_neighbor_data()
{
    return &neighbors;
}

//-------------------------------------------------------------------------------------------------
st_hash_tree::st_hash_tree()
{
	log2inv = 1./log(2.);
}

bool st_hash_tree::create_mesh(KDLayoutData &data){
	/*
	Create a mesh of the heliostat field according to the performance surface
	provided by the 'integrals' class.
	*/
	Data = data;
        
	//Calculate min and max recursion levels based on user zone size limitations
	double dextx = (Data.xlim[1] - Data.xlim[0]);
	nx_req = (int)floor( log(dextx/Data.min_unit_dx)*log2inv );
    nx_req = nx_req < 1 ? 1 : nx_req;
	double dexty = (Data.ylim[1] - Data.ylim[0]);
	ny_req = (int)floor( log(dexty/Data.min_unit_dy)*log2inv );
    ny_req = ny_req < 1 ? 1 : ny_req;

	//estimate the maximum number of nodes and reserve memory
	int nmaxdivx = (int)pow(2., nx_req);		//maximum number of zones x
	int nmaxdivy = (int)pow(2., ny_req);	//max y
	int nmaxterm =  nmaxdivx * nmaxdivy;	//total max number of zones
	int maxreclevel = max(nmaxdivx, nmaxdivy); //worst case max recursion level
	int nmaxnodes = 0;
	for(int i=0; i<maxreclevel; i++)
		nmaxnodes += nmaxterm/(int)pow(2.,i);		//Add each level in the node tree
	
	//Try reserving the number of required nodes, catch any memory error
	try
	{
		nodes.reserve(nmaxnodes*2); //include a 100% buffer
	}
	catch(...)
	{
        //Memory error
        return false;
	}
	
    //set up the head node's range. This doesn't actually create the tree yet.
	head_node.set_range(Data.xlim[0], Data.xlim[1], Data.ylim[0], Data.ylim[1]);
    
    return true;
}

string st_hash_tree::pos_to_binary_base(double x, double y){
    double res = fmin(Data.min_unit_dx, Data.min_unit_dy);
	return pos_to_binary(x, y, res);
}

void st_hash_tree::reset(){
    Data.min_unit_dx = 
        Data.min_unit_dy = 
        Data.xlim[0] = Data.xlim[1] = Data.ylim[0] = Data.ylim[1] =
            std::numeric_limits<double>::quiet_NaN();
    
	head_node = st_opt_element();
	nodes.clear();
	nx_req = -1;
	ny_req = -1;
}

void st_hash_tree::get_terminal_data(vector<vector<void*>*> &retdata){

    for( vector<st_opt_element>::iterator it = nodes.begin(); it != nodes.end(); it++){
		if(! it->is_terminal() ) continue;
		retdata.push_back(it->get_array());
	}
}
void st_hash_tree::get_terminal_nodes(vector<st_opt_element*> &tnodes){

	tnodes.clear(); 

    for(int i=0; i<(int)nodes.size(); i++){
		if( nodes.at(i).is_terminal() )
			tnodes.push_back(&nodes.at(i));
	}

	return;
}

vector<st_opt_element>* st_hash_tree::get_all_nodes()
{
    return &nodes;
}

void st_hash_tree::add_neighborhood_data() 
{
    /* 
    Collect the data elements that surround addresses that already contain elements. If the neighboring elements
    of the same resolution are empty (no data), then create empty zones and flag as needing to be processed.
    */

    vector<st_opt_element*> tnodes;
    unordered_map< string, set<st_opt_element*> > empty_nb_map;

    get_terminal_nodes(tnodes);
    for(int k=0; k<tnodes.size(); k++)
    {
        st_opt_element* node = tnodes.at(k);
        int add[2];

        for(int i=-1; i<2; i++)
        {
            for(int j=-1;j<2; j++)
            {
                if( i == 0 && j == 0 )
                    continue;   //omit self zone

                //do binary addition to find the specified neighboring zone
                add[0] = i;
                add[1] = j;

                string nzone;

                if( binary_add(node->get_address(), nzone, add) )
                {
                    st_opt_element *el = head_node.process(nzone, 0);
                    if(el != 0)
                    {
                        vector<void*>* dat = el->get_array();
                        if( dat != 0 )
                        {
                            for(int k=0; k<dat->size(); k++)
                                node->get_neighbor_data()->push_back(dat->at(k));
                        }
                    }
                    else
                    {
                        empty_nb_map[ nzone ].insert( node );
                    }
                }
            }
        }
    }

    //Go back and add empty neighbors that were identified 
    for(unordered_map< string, set<st_opt_element*> >::iterator mp=empty_nb_map.begin(); mp != empty_nb_map.end(); mp++)
    {
        double x,y;
        //get center of zone
        binary_to_pos(mp->first, &x, &y);
        //get an appropriate element size to add
        double *xr = (*mp->second.begin())->get_xr();
        double *yr = (*mp->second.begin())->get_yr();
        double rmin[2];
        rmin[0] = (xr[1] - xr[0])*0.55; //make the "object" smaller than the previous zone size but larger than half the zone size
        rmin[1] = (yr[1] - yr[0])*0.55;
        //Add the element (actually just create the empty terminal zone)
        add_object(0, x, y, rmin);

        //now add neighbors directly to the new terminal zone
        st_opt_element* newnode = &nodes.back();
        for( set<st_opt_element*>::iterator nit = mp->second.begin(); nit != mp->second.end(); nit++)
        {
            for(int j=0; j<(*nit)->get_array()->size(); j++)
                newnode->get_neighbor_data()->push_back( (*nit)->get_array()->at(j) );
        }
    }
}


bool st_hash_tree::get_all_data_at_loc(vector<void*> &data, double locx, double locy)
{
    string bin = pos_to_binary_base(locx, locy);

    data.clear();

    st_opt_element* z = head_node.process(bin, 0);
    if( z != 0 )
    {
        vector<void*>* zd = z->get_array();
        if(zd != 0)
            for(int i=0; i<zd->size(); i++)
                data.push_back(zd->at(i));
        
        zd = z->get_neighbor_data();
        if(zd != 0)
            for(int i=0; i<zd->size(); i++)
                data.push_back(zd->at(i));

        if( !data.empty() )
            return true;
    }

    return false;
}

void st_hash_tree::create_node(st_opt_element &node, int index, string &binary, void *object, double *dprojected)
{
    /* 
    node        |   pointer to parent node
    index       |   current character location in the string 'binary'
    binary      |   string of characters indicating the path to follow to create a new node
    object      |   pointer to the object that is to be added to the new node
    dprojected  |   2-value array: 1st value is span in azimuthal direction, 2nd value spans zenith direction
    */
    
    
    //evaluate the derivatives at the center of the element.
    double 
		xr0 = node.get_xr()[0],
		xr1 = node.get_xr()[1],
		yr0 = node.get_yr()[0],
		yr1 = node.get_yr()[1],
		C0 = (xr0 + xr1)*0.5,
		C1 = (yr0 + yr1)*0.5;
    double dpx, dpy;
    if( dprojected == 0)
    {
        dpx = 0.;
        dpy = 0.;
    }
    else
    {
        dpx = dprojected[0];
        dpy = dprojected[1];
    }
    
    bool x_direction = (double)(index / 2) == ((double)index / 2.);

    int x_rec_level, y_rec_level;

    if(x_direction)
    {
        x_rec_level = index / 2;
        y_rec_level = x_rec_level - 1;
    }
    else
    {
        x_rec_level = y_rec_level = (index-1)/2;
    }


    //check to see if we're at the end of the string
    if(index > binary.size()-1)
    {
        //Add the object to the data vector and return
        node.setup(object);
        return;
    }

    //check to see if this node is terminal. If so, add the data here rather than continuing branching.
    if( node.is_terminal() )
    {
        node.setup(object);
        return;
    }

    //which direction split is indicated?
    bool split_pos = binary.at(index) == '1';


    if( x_direction )
    {
        if( x_rec_level < nx_req && dpx < (C0-xr0))
        {

            //does a node already exist in the proposed split direction?
            if( node.get_child_node( split_pos ) == 0 )
            {
                //we can split in the x direction
                nodes.push_back(st_opt_element());
                st_opt_element *m = &nodes.back();
                m->set_address( node.get_address()+binary.at(index) );

                if(split_pos)   //upper
                {
                    double xr[] = {C0, xr1};
                    double yr[] = {yr0, yr1};
                    m->set_range(xr, yr);
                    node.setup(0, m);

                    create_node(*m, index + 1, binary, object, dprojected);
                    return;
                }
                else        //lower
                {
                    double xr[] = {xr0, C0};
                    double yr[] = {yr0, yr1};
                    m->set_range(xr,yr);
                    node.setup(m, 0);

                    create_node(*m, index + 1, binary, object, dprojected);
                    return;
                }
            }
            else
            {
                //a node already exists, so just follow along without adding a new node in its place
                create_node(*(st_opt_element*)(node.get_child_node(split_pos)), index + 1, binary, object, dprojected);
                return;
            }
        }
        else
        {
            //no more splits allowed in the x direction
            //check if more y splits are allowed
            if( y_rec_level < ny_req && dpy < (C1 - yr0) )
            {
                if( node.get_child_node( split_pos ) == 0 )
                {
                    nodes.push_back(st_opt_element());
                    st_opt_element *m = &nodes.back();

                    double xr[] = {xr0, xr1};
                    double yr[] = {yr0, yr1};

                    m->set_range(xr, yr);
                    m->set_address( node.get_address() + "x" );
                    node.setup(m);
                
                    create_node(*m, index + 1, binary, object, dprojected);
                    return;
                }
                else
                {
                    create_node(*(st_opt_element*)(node.get_child_node(split_pos)), index + 1, binary, object, dprojected);
                    return;
                }
            }
            else
            {
                //no more y splits are allowed either, so add the object and return
                node.setup(object);
                return;
            }
        }
    }
    else    //Y direction
    {
        if( y_rec_level < ny_req && dpy < (C1 - yr0) )
        {

            //does a node already exist in the proposed split direction?
            if( node.get_child_node( split_pos ) == 0 )
            {
                //we can split in the y direction
                nodes.push_back(st_opt_element());
                st_opt_element *m = &nodes.back();
                m->set_address( node.get_address()+binary.at(index) );

                if(split_pos)   //upper
                {
                    double xr[] = {xr0, xr1};
                    double yr[] = {C1, yr1};
                    m->set_range(xr, yr);
                    node.setup(0, m);

                    create_node(*m, index + 1, binary, object, dprojected);
                    return;
                }
                else        //lower
                {
                    double xr[] = {xr0, xr1};
                    double yr[] = {yr0, C1};
                    m->set_range(xr,yr);
                    node.setup(m, 0);

                    create_node(*m, index + 1, binary, object, dprojected);
                    return;
                }
            }
            else
            {
                //a node already exists, so just follow along without adding a new node in its place
                create_node(*(st_opt_element*)(node.get_child_node(split_pos)), index + 1, binary, object, dprojected);
                return;
            }
        }
        else
        {
            //no more splits allowed in the y direction
            //check if more y splits are allowed
            if( x_rec_level < nx_req && dpx < (C0 - xr0) )
            {
                if( node.get_child_node( split_pos ) == 0 )
                {
                    nodes.push_back(st_opt_element());
                    st_opt_element *m = &nodes.back();

                    double xr[] = {xr0, xr1};
                    double yr[] = {yr0, yr1};

                    m->set_range(xr, yr);
                    m->set_address( node.get_address() + "x" );
                    node.setup(m);
                
                    create_node(*m, index + 1, binary, object, dprojected);
                    return;
                }
                else
                {
                    create_node(*(st_opt_element*)(node.get_child_node(split_pos)), index + 1, binary, object, dprojected);
                    return;
                }
            }
            else
            {
                //no more x splits are allowed either, so add the object and return
                node.setup(object);
                return;
            }
        }
    }

}

string st_hash_tree::pos_to_binary(double x, double y, double res){
	/*
	Convert an x-y position into a binary tag
	*/
	    
	string tag;
        
	bool x_mode = true; //start with radius
        
	double 
		y0 = Data.ylim[0],
		y1 = Data.ylim[1],
		x0 = Data.xlim[0],
		x1 = Data.xlim[1];
	    
	int nc = max(nx_req, ny_req)*2;
        
	for(int i=0; i<nc; i++){
		if( x_mode){
			double cx = (x0 + x1)*0.5;
			if(x > cx){
				x0 = cx;
				tag.append("1");
			}
			else{
				x1 = cx;
				tag.append("0");
			}
		}
		else{
			double cy = (y0 + y1)*0.5;
			if(y > cy){
				y0 = cy;
				tag.append("1");
			}
			else{
				y1 = cy;
				tag.append("0");
			}
		}
		x_mode = ! x_mode;
	}
	return tag;
}

void st_hash_tree::binary_to_pos(const string &binary, double *x, double *y)
{
    /* 
    Returns double[2] of x,y position of center of zone described by binary.
    */
    bool x_mode = true; //start with radius
        
	double 
		y0 = Data.ylim[0],
		y1 = Data.ylim[1],
		x0 = Data.xlim[0],
		x1 = Data.xlim[1];

    *x = (x0 + x1)/2.;
    *y = (y0 + y1)/2.;

    for(int i=0; i<binary.size(); i++)
    {
        if(x_mode)
        {
            switch(binary.at(i))
            {
            case '0':
                x1 = *x;
                *x = (x1 + x0)/2.;
                break;
            case '1':
                x0 = *x;
                *x = (x1 + x0)/2.;
                break;
            case 'x':
                break;
            }
        }
        else
        {
            switch(binary.at(i))
            {
            case '0':
                y1 = *y;
                *y = (y1 + y0)/2.;
                break;
            case '1':
                y0 = *y;
                *y = (y1 + y0)/2.;
                break;
            case 'x':
                break;
            }
        }

        x_mode = !x_mode;
    }

    return;
}

void st_hash_tree::add_object(void *object, double locx, double locy, double *objsize)
{
    string tag = pos_to_binary_base(locx, locy);
    create_node(head_node, 0, tag, object, objsize);
}
