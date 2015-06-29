#ifdef _WIN32
	#include <windows.h>
	#include <GL/gl.h>
	#include <GL/freeglut.h>
	#include <GL/glu.h>
#else	
	#include <GL/gl.h>
	#include <GL/glut.h>
	#include <GL/glu.h>

#endif

#include <pthread.h>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <ctype.h>
#include <string>
#include <cassert>
#include <vec3.h>
#include <mesh.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

struct camera{
  vec3f    pos;
  GLfloat elevation;
  GLfloat azimut;
  GLfloat zangle;
}cam;

mesh model;

GLfloat *vertices ;
GLuint *indices ;

int keyb[256]; //keyboard list

struct bvh_node{
	size_t idx;
	bvh_node *parent;

	aabb bb;
	bvh_node(const aabb &bb_):idx(0),bb(bb_){}
	
	virtual bool is_leaf() const  = 0;	
	//virtual size_t get_ntris() const = 0;	
	//virtual size_t count() = 0;
	
};

struct inner_node: public bvh_node{
	bvh_node * nodes_ptr[2];
	size_t   nodes_idx[2];
	inner_node(const aabb &bb_,
			   size_t node1_idx,
	           size_t node2_idx):
	           bvh_node(bb_){
		nodes_idx[0] = node1_idx;
		nodes_idx[1] = node2_idx;
	};
	
	virtual bool is_leaf() const{
		return false;
	}
	
	//virtual size_t get_ntris() const;
	
	//virtual size_t count();
	
};

struct leaf_node: public bvh_node{
	size_t start_idx, end_idx;

	
	leaf_node(const aabb &bb_,
			  size_t start_idx_, 
			  size_t end_idx_):bvh_node(bb_),
				start_idx(start_idx_),
				end_idx(end_idx_){}
					
	virtual bool is_leaf() const{
		return true;
	}
	
	//virtual size_t get_ntris() const;
	
	//virtual size_t count();
	
};


std::vector<bvh_node*>  nodes;

bvh_node *root;

bvh_node *current;


void draw_reference_frame();

void draw_plane(int Nx,int Nz);

void init(void);

bool load_model_ssv(const std::string filename);

bool load_model_nbin(const std::string filename);


void draw_model(void);

void draw_bounding_box(const aabb &bb);


void display(void);

void reshape (int w, int h);

void keyboardUp(unsigned char key, int x, int y);

void keyboardDown(unsigned char key, int x, int y);

void loop(int value);

std::string trim(const std::string& str, const std::string& whitespace = " \t");

bool get_valid_line(std::istream &ifile,std::string &line);

bool compare_extension(const std::string &filename,const std::string &extension);


int main(int argc, char** argv){	
	pthread_getconcurrency();
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB| GLUT_DEPTH);
	glutInitWindowSize (1024, 768); 
	glutInitWindowPosition (100, 100);
	glutCreateWindow (argv[0]);
	init ();
	if(argc >1){ 
		std::string filename = std::string(argv[1]);		
		if(compare_extension(filename,"ssv")){
			load_model_ssv(filename);
		}else if(compare_extension(filename,"bin")){
				assert(false);
		}else if(compare_extension(filename,"nbin")){
			load_model_nbin(filename);

		}else{
			std::cerr << "file format not recognized:  " <<filename << std::endl;
		}
	}
	
	glutDisplayFunc(display); 
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboardDown);
	glutKeyboardUpFunc(keyboardUp);
	glutTimerFunc(50,loop, 0);
	glutMainLoop();
	return 0;
}



static aabb get_tri_aabb(const vec3d &va,const vec3d &vb,const vec3d &vc){
	const vec3d min_vec = cwisemin(cwisemin(va,vb),vc);
	const vec3d max_vec = cwisemax(cwisemax(va,vb),vc);    
	return aabb(min_vec,max_vec);
}


aabb compute_aabb(const std::vector<vec3d> &vlist,
				  const std::vector<vec4i> &flist){

	
    aabb bb;
    for(size_t i = 0; i< flist.size(); i++){
		const vec3d va = vlist[flist[i].a];
		const vec3d vb = vlist[flist[i].b];
		const vec3d vc = vlist[flist[i].c];	
	
        bb.grow(get_tri_aabb(va,vb,vc));
    }

	return bb;
}


struct __attribute__((__packed__)) vertex3f{
	float x,y,z;
	
	void operator =(const vec3d &v){
		x = v.x;
		y = v.y;
		z = v.z;
	}
};

struct __attribute__((__packed__)) face4i{
	int a,b,c,d;
	
	void operator =(const vec4i &f){
		a = f.a;
		b = f.b;
		c = f.c;
		d = f.d;		
	}
};

namespace nvcuda_types{

	struct __attribute__((__packed__)) nv_int2{
		int x,y;
	};

	struct __attribute__((__packed__)) nv_float4{
		float x,y,z,w;
	};

}

bool load_model_nbin(const std::string filename){
	std::string line;
    size_t num_verts = 0;
    size_t num_faces = 0;
    size_t remapping_table_size =0;
    size_t node_list_size =0;
    
    std::ifstream ifile(filename.c_str());
		
	if (!ifile.is_open()) {		
		return false;
	}
	
	{	//read the number of vertex and faces
		get_valid_line(ifile,line);
	    std::stringstream ss(line);
	    ss >> num_verts >> num_faces;	    
	}
	
	vertex3f * packed_verts = new vertex3f[num_verts];
	face4i * packed_faces = new face4i[num_faces];
	
	ifile.read(reinterpret_cast<char*>(packed_verts),sizeof(vertex3f)*num_verts);
	ifile.read(reinterpret_cast<char*>(packed_faces),sizeof(face4i)*num_faces);

	//skip the new line character
	int  newline = ifile.get();
	assert(newline == '\n');
	
	
	
	
	
	model.vlist.reserve(num_verts);
	model.flist.reserve(num_faces);
	
	for(size_t i = 0 ; i < num_verts; i++){		
		model.vlist.push_back(
			make_vec3d(packed_verts[i].x,
			           packed_verts[i].y,
					   packed_verts[i].z
			)
		);
	}

	for(size_t i = 0 ; i < num_faces; i++){		
		model.flist.push_back(
			make_vec4i(packed_faces[i].a,
			           packed_faces[i].b,
					   packed_faces[i].c,
					   packed_faces[i].d
			)
		);
	}



	aabb bb = compute_aabb(model.vlist, model.flist);
	
	vec3d half_size = bb.get_half_size();
	vec3d center    = bb.get_center();
	double scale_factor = 10/min(half_size);
	
	//rescaling the vertices	
	for(size_t i = 0 ; i < num_verts; i++){
		model.vlist[i] = (model.vlist[i]-center)*scale_factor;
	}	
	//recompute the aabb
	model.bb = compute_aabb(model.vlist, model.flist);
	delete [] packed_verts;
	delete [] packed_faces;
	

	
	
	//copy the array of vertices
	vertices = new GLfloat[3*num_verts];
	indices  = new GLuint[3*num_faces];
	
	for(size_t i =0; i < num_verts;i++){
		vertices[3*i]   = model.vlist[i].x; 
		vertices[3*i+1] = model.vlist[i].y;
		vertices[3*i+2] = model.vlist[i].z;		
	}
	 
	for(size_t i =0; i < num_faces;i++){
		indices[3*i]   = model.flist[i].a; 
		indices[3*i+1] = model.flist[i].b; 
		indices[3*i+2] = model.flist[i].c; 
	
	} 
	
	
	{//read the size of remapping table and nonde list
		get_valid_line(ifile,line);
	    std::stringstream ss(line);
	    ss >> remapping_table_size >> node_list_size;	    
	}
	
	int32_t * remapping_table = new int32_t[remapping_table_size];
	ifile.read(reinterpret_cast<char*>(remapping_table),sizeof(int32_t)*remapping_table_size);

	delete []remapping_table;
	
	
	using namespace nvcuda_types;
	nv_int2   * N0 = new nv_int2[node_list_size];
	nv_float4 * N1 = new nv_float4[node_list_size];	
	nv_float4 * N2 = new nv_float4[node_list_size];	
	nv_float4 * N3 = new nv_float4[node_list_size];


	ifile.read(reinterpret_cast<char*>(N0),sizeof(nv_int2)*node_list_size);
	ifile.read(reinterpret_cast<char*>(N1),sizeof(nv_float4)*node_list_size);
	ifile.read(reinterpret_cast<char*>(N2),sizeof(nv_float4)*node_list_size);
	ifile.read(reinterpret_cast<char*>(N3),sizeof(nv_float4)*node_list_size);

	
	
	nodes.resize(node_list_size);
	for(size_t i =0; i < node_list_size; i++){
		nodes[i] = NULL;
	}
	
	nodes[0] = new inner_node(model.bb,0,0);
	//read the node list
	for(size_t i =0; i < node_list_size; i++){
		if(!nodes[i]->is_leaf()){
			inner_node * inode = reinterpret_cast<inner_node*>(nodes[i]);
			
			aabb bb1(
				make_vec3d(
					N1[i].x,N1[i].z,N3[i].x
				),
				make_vec3d(
					N1[i].y,N1[i].w,N3[i].y
				)
			);
			
			bb1.min_vec = (bb1.min_vec - center)*scale_factor;
			bb1.max_vec = (bb1.max_vec - center)*scale_factor;
			
			aabb bb2(
				make_vec3d(
					N2[i].x,N2[i].z,N3[i].z
				),
				make_vec3d(
					N2[i].y,N2[i].w,N3[i].w
				)
			);

			bb2.min_vec = (bb2.min_vec - center)*scale_factor;
			bb2.max_vec = (bb2.max_vec - center)*scale_factor;
						
			if(N0[i].x > 0){
				inode->nodes_idx[0] = N0[i].x;
				nodes[N0[i].x] = new inner_node(bb1,0,0);
			}else{
				inode->nodes_idx[0]= ~N0[i].x;
				nodes[inode->nodes_idx[0]] = new leaf_node(bb1,0,0);				
			}

			if(N0[i].y > 0){
				inode->nodes_idx[1] = N0[i].y;
				nodes[N0[i].y] = new inner_node(bb2,0,0);
			}else{
				inode->nodes_idx[1]= ~N0[i].y;
				nodes[inode->nodes_idx[1]] = new leaf_node(bb2,0,0);				
			}

			
		}else{
			leaf_node * lnode = reinterpret_cast<leaf_node*>(nodes[i]);			
			lnode->start_idx = N0[i].x;
			lnode->end_idx = N0[i].y;
		}
		
	}

	//link the pointers of the nodes 
	for(size_t i =0; i < node_list_size; i++){
		if(!nodes[i]->is_leaf()){
			inner_node *inode = reinterpret_cast<inner_node*>(nodes[i]);
			inode->nodes_ptr[0] = nodes[inode->nodes_idx[0]];
			inode->nodes_ptr[1] = nodes[inode->nodes_idx[1]];
			nodes[inode->nodes_idx[0]]->parent = inode;
			nodes[inode->nodes_idx[1]]->parent = inode;
		}
	}



	root = nodes[0];
	root->parent = NULL;
	
	
	current = root;
	


	return true;
}




bool load_model_ssv(const std::string filename){
	
	std::string line;
    size_t num_verts = 0;
    size_t num_faces = 0;
    size_t remapping_table_size =0;
    size_t node_list_size =0;
    
    std::ifstream ifile(filename.c_str());
		
	if (!ifile.is_open()) {		
		return false;
	}
	
	{	//read the number of vertex and faces
		get_valid_line(ifile,line);
	    std::stringstream ss(line);
	    ss >> num_verts >> num_faces;	    
	}
	
	model.vlist.reserve(num_verts);
	model.flist.reserve(num_faces);
	
	for(size_t i = 0 ; i < num_verts; i++){
		double x,y,z;
		get_valid_line(ifile,line);
	    std::stringstream ss(line);
		ss >> x>>y>>z;
		model.vlist.push_back(make_vec3d(x,y,z));
	}

	for(size_t i = 0 ; i < num_faces; i++){
		int a,b,c,d;
		get_valid_line(ifile,line);
	    std::stringstream ss(line);
		ss >> a>>b>>c>>d;
		model.flist.push_back(make_vec4i(a,b,c,d));
	}
	
	aabb bb = compute_aabb(model.vlist, model.flist);
	
	vec3d half_size = bb.get_half_size();
	vec3d center    = bb.get_center();
    
	double scale_factor = (min(half_size) !=0 )?10/min(half_size):1;
	
	//rescaling the vertices	
	for(size_t i = 0 ; i < num_verts; i++){
		model.vlist[i] = (model.vlist[i]-center)*scale_factor;
	}	
	//recompute the aabb
	model.bb = compute_aabb(model.vlist, model.flist);
	
	{//read the size of remapping table and nonde list
		get_valid_line(ifile,line);
	    std::stringstream ss(line);
	    ss >> remapping_table_size >> node_list_size;	    
	}
	
	//skip the remapping table
	get_valid_line(ifile,line);
	
	
	nodes.reserve(node_list_size);
	//read the node list
	for(size_t i =0; i < node_list_size; i++){
		get_valid_line(ifile,line);
	    std::stringstream ss(line);
	    int node_type ;
	    ss >>node_type;
	    switch(node_type){
			case 0:{ //inner node 
				// 0 node1_idx node2_idx aabb min (x y z) aabb max (x y z)
				size_t node1_idx,node2_idx;
				aabb bb;
				
				ss >> node1_idx >> node2_idx 
				   >> bb.min_vec.x >> bb.min_vec.y >> bb.min_vec.z
				   >> bb.max_vec.x >> bb.max_vec.y >> bb.max_vec.z;
				
				bb.min_vec = (bb.min_vec - center)*scale_factor;
				bb.max_vec = (bb.max_vec - center)*scale_factor;

				//bb.scale(scale_factor);

				nodes.push_back(new inner_node(bb,node1_idx,node2_idx));

			}break;
			case 1:{ //leaf node
				// 1 start_idx end_idx aabb min (x y z) aabb max (x y z)
				size_t start_index,end_index;
				aabb bb;
				ss >> start_index >> end_index 
				   >> bb.min_vec.x >> bb.min_vec.y >> bb.min_vec.z
				   >> bb.max_vec.x >> bb.max_vec.y >> bb.max_vec.z;

				bb.min_vec = (bb.min_vec - center)*scale_factor;
				bb.max_vec = (bb.max_vec - center)*scale_factor;


  				//bb.scale(scale_factor);

				nodes.push_back(new leaf_node(bb,start_index,end_index));				   
			}break;
		}
	}

	//link the pointers of the nodes 
	for(size_t i =0; i < node_list_size; i++){
		if(!nodes[i]->is_leaf()){
			inner_node *inode = reinterpret_cast<inner_node*>(nodes[i]);
			inode->nodes_ptr[0] = nodes[inode->nodes_idx[0]];
			inode->nodes_ptr[1] = nodes[inode->nodes_idx[1]];
			nodes[inode->nodes_idx[0]]->parent = inode;
			nodes[inode->nodes_idx[1]]->parent = inode;
		}
	}



	root = nodes[0];
	root->parent = NULL;
	
	
	current = root;
	
	
	//copy the array of vertices
	vertices = new GLfloat[3*num_verts];
	indices  = new GLuint[3*num_faces];
	
	for(size_t i =0; i < num_verts;i++){
		vertices[3*i]   = model.vlist[i].x; 
		vertices[3*i+1] = model.vlist[i].y;
		vertices[3*i+2] = model.vlist[i].z;		
	}
	 
	for(size_t i =0; i < num_faces;i++){
		indices[3*i]   = model.flist[i].a; 
		indices[3*i+1] = model.flist[i].b; 
		indices[3*i+2] = model.flist[i].c; 
	} 
	
	return true; 
}


void init(void){
    glClearColor (0.0f, 0.0f, 0.0f, 0.0f);
    glClearDepth(1.0f);				
    glDepthFunc(GL_LESS);			        
    glEnable(GL_DEPTH_TEST);		       
    glShadeModel (GL_SMOOTH);
    srand((unsigned int) time(NULL));
    cam.pos = make_vec3f(0,3,15);

    cam.elevation = 0;
    cam.azimut = 0;
}


vec3f get_colour(GLfloat v,GLfloat vmin,GLfloat vmax)
{
   vec3f c = {1.0f,1.0f,1.0f}; // white
   GLfloat dv;

   if (v < vmin)
      v = vmin;
   if (v > vmax)
      v = vmax;
   dv = vmax - vmin;

   if (v < (vmin + 0.25f * dv)) {
      c.x = 0;
      c.y = 4 * (v - vmin) / dv;
   } else if (v < (vmin + 0.5f * dv)) {
      c.x = 0;
      c.z = 1 + 4 * (vmin + 0.25f * dv - v) / dv;
   } else if (v < (vmin + 0.75f * dv)) {
      c.x = 4 * (v - vmin - 0.5f * dv) / dv;
      c.z = 0;
   } else {
      c.y = 1 + 4 * (vmin + 0.75f * dv - v) / dv;
      c.z = 0;
   }

   return c;
}


void draw_model(void){
    
	if(model.vlist.size() != 0){	
		// activate and specify pointer to vertex array
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, vertices);

		// draw a cube
		glDrawElements(GL_TRIANGLES, 3*model.flist.size(), GL_UNSIGNED_INT, indices);

		// deactivate vertex arrays after drawing
		glDisableClientState(GL_VERTEX_ARRAY);

		
		
// 		glBegin(GL_TRIANGLES);
// 
// 		for(size_t i = 0; i < model.flist.size(); i++){
// 			vec3d  a = model.vlist[model.flist[i].a];
// 			vec3d  b = model.vlist[model.flist[i].b];
// 			vec3d  c = model.vlist[model.flist[i].c];
//             std::cout << a << " " << b << " " << c << " " <<std::endl;
// 			////GLfloat ka = 1/(rho1[faces[i].x]*rho2[faces[i].x]);
// 			////GLfloat kb = 1/(rho1[faces[i].y]*rho2[faces[i].y]);
// 			////GLfloat kc = 1/(rho1[faces[i].z]*rho2[faces[i].z]);
// 			////vec3f ca = get_colour(ka,min_k,max_k);
// 			////vec3f cb = get_colour(kb,min_k,max_k);
// 			////vec3f cc = get_colour(kc,min_k,max_k);
// 			//////printf("%f %f %f \n",ca.x, ca.y, ca.z);
// 
// 			////glColor3f (ca.x, ca.y, ca.z);
// 			glVertex3f(a.x, a.y, a.z);
// 					
//      		////glColor3f (cb.x, cb.y, cb.z);
// 			glVertex3f(b.x, b.y, b.z);
// 			
// 			////glColor3f (cc.x, cc.y, cc.z);
// 			glVertex3f(c.x, c.y, c.z);
// 
// 		}	
// 		glEnd();
        
	}
}



void display(void){
    glClear (GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

    glPushMatrix();

    //posizionamento telecamera su montatura azimutale
    glRotatef(cam.elevation,1.0,0.0,0.0);
    glRotatef(cam.azimut,0.0,1.0,0.0);   
    glTranslatef(-cam.pos.x,-cam.pos.y,-cam.pos.z);
    glRotatef(-90,1.0,0.0,0.0);
    // assi
    glPushMatrix();
    glScalef(1.0,1.0,1.0);
    draw_reference_frame();
    glPopMatrix();  
   
    //disegna piano
    glColor3f (1.0, 1.0, 0.0);  
    //draw_plane(100,100);
	
	draw_model();
	
    glColor3f (1.0, 0.0, 0.0);  
	draw_bounding_box(current->bb);
	
	if(!current->is_leaf()){
		inner_node *inode = reinterpret_cast<inner_node*>(current);

		glColor3f (0.0, 1.0, 0.0);  
		draw_bounding_box(inode->nodes_ptr[0]->bb);
		glColor3f (0.0, 0.0, 1.0);  
		draw_bounding_box(inode->nodes_ptr[1]->bb);
	}
	
    glPopMatrix();
    
	
	
    glutSwapBuffers();
}

void reshape (int w, int h){
    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(65.0f,(GLfloat) w/((GLfloat)h),0.01f,1000.0f);
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();
}

void keyboardUp(unsigned char key, int x, int y){ 
	
	bool show_index = false; 	
	keyb[key] = 0;
	
	if(key == 'v'){
		if(!current->is_leaf()){			
			current = reinterpret_cast<inner_node*>(current)->nodes_ptr[0];
			show_index = true;
		}
    }
    if(key ==  'b'){
   		if(!current->is_leaf()){			
			current = reinterpret_cast<inner_node*>(current)->nodes_ptr[1];
			show_index = true;
		}
    }
    if(key ==  'p'){
		if(current->parent != NULL){
			current = current->parent;
			show_index = true;
		}
    }	
    
    if( show_index && !current->is_leaf() ){
		std::cout << "node1 idx:" <<  reinterpret_cast<inner_node*>(current)->nodes_idx[0] 
				  << "node2 idx:" <<  reinterpret_cast<inner_node*>(current)->nodes_idx[1]
				  << std::endl;
	}
    
}
void keyboardDown(unsigned char key, int x, int y){

  keyb[key] = 1;   
}

void loop(int value){
    GLfloat dx,dz;
    if(keyb['w']){
        dx = sinf(cam.azimut*M_PI/180.0f)/10.0f;
        dz = -cosf(cam.azimut*M_PI/180.0f)/10.0f;
        cam.pos.x +=dx;
        cam.pos.z +=dz;
    }
    if(keyb['s']){
        dx = sinf(cam.azimut*M_PI/180.0f)/10.0f;
        dz = -cosf(cam.azimut*M_PI/180.0f)/10.0f;
        cam.pos.x -= dx;
        cam.pos.z -= dz;
    }
    
    if(keyb['a']){
		dx = -cosf(cam.azimut*M_PI/180.0f)/10.0f;
        dz = -sinf(cam.azimut*M_PI/180.0f)/10.0f;
        cam.pos.x +=dx;
        cam.pos.z +=dz;
    }
    if(keyb['d']){
        dx = -cosf(cam.azimut*M_PI/180.0f)/10.0f;
        dz = -sinf(cam.azimut*M_PI/180.0f)/10.0f;
        
        cam.pos.x -= dx;
        cam.pos.z -= dz;
    }

    
    if(keyb['m']){
        cam.pos.y+=0.1;        
    }
    if(keyb['n']){
        cam.pos.y-=0.1;        
    }
    if(keyb['i']){
        cam.elevation+=0.5f;
    }
    if(keyb['k']){
        cam.elevation-=0.5f;
    }
    if(keyb['j']){
        cam.azimut-=0.5f;
    }
    if(keyb['l']){
        cam.azimut+=0.5f;
    }
   
    
    if(keyb[27]==1){
    	exit(0);
    }
    glutPostRedisplay();
    glutTimerFunc(30,loop, 0);
}

void draw_reference_frame(){
  glLineWidth(1.5f);
  glBegin(GL_LINES);
  //asse x;
  glColor3f (1.0f, .0f, 0.0f);  
  glVertex3f(0.0f,0.0f,0.0f);
  glVertex3f(1.0f,0.0f,0.0f);
  //asse y
  glColor3f (0.0f, 1.0f, 0.0f); 
  glVertex3f(0.0f,0.0f,0.0f);
  glVertex3f(0.0f,1.0,0.0f);
  //asse z
  glColor3f (0.0f, 0.0f, 1.0f); 
  glVertex3f(0.0f,0.0f,0.0f);
  glVertex3f(0.0f,0.0f,1.0f);

  glEnd();
  glLineWidth(1.0f);
}


void draw_segment(const vec3d A, const vec3d B){
    glVertex3f(A.x,A.y,A.z);
    glVertex3f(B.x,B.y,B.z);
}

void draw_bounding_box(const aabb &bb){
	glLineWidth(1.5f);
        
	glBegin(GL_LINES);    
	const vec3d  A = bb.min_vec;
	const vec3d  B = make_vec3d(bb.max_vec.x,bb.min_vec.y,bb.min_vec.z);
	const vec3d  C = make_vec3d(bb.max_vec.x,bb.min_vec.y,bb.max_vec.z);
	const vec3d  D = make_vec3d(bb.min_vec.x,bb.min_vec.y,bb.max_vec.z);
	const vec3d  E = make_vec3d(bb.min_vec.x,bb.max_vec.y,bb.min_vec.z);	
	const vec3d  F = make_vec3d(bb.max_vec.x,bb.max_vec.y,bb.min_vec.z);	
	const vec3d  G = bb.max_vec;
	const vec3d  H = make_vec3d(bb.min_vec.x,bb.max_vec.y,bb.max_vec.z);	
	//glColor3f (0.5f, 1.0f, 1.0f);  
	
	draw_segment(A,B);
	draw_segment(B,C);
	draw_segment(C,D);
	draw_segment(D,A);
	draw_segment(A,E);
    draw_segment(D,H);
    draw_segment(B,F);
    draw_segment(C,G);
    
    draw_segment(E,F);
    draw_segment(F,G);
    draw_segment(G,H);
    draw_segment(H,E);
    
    glEnd();
	glLineWidth(1.0f);        
}

void draw_plane(int Nx,int Nz){
    int i,j;
    GLfloat x,z;
  
    for(i = 0;i<2*Nx;i++){
        for(j = 0;j<2*Nz;j++){
            x = i-Nx+0.5f;
            z = j-Nz+0.5f;
            if((i%2 ==0 && j %2==0) || (i%2 ==1 && j %2 ==1)){
                glColor3f (1.0f, 1.0f, 1.0f);
            }else{
	            glColor3f (0.0f,0.8f,0.9f);
	        }
            glBegin(GL_TRIANGLES);
            glVertex3f(x-0.5f,0.0f,z+0.5f); 
            glVertex3f(x+0.5f,0.0f,z-0.5f);
            glVertex3f(x-0.5f,0.0f,z-0.5f);            
            glVertex3f(x-0.5f,0.0f,z+0.5f);
            glVertex3f(x+0.5f,0.0f,z+0.5f); 
            glVertex3f(x+0.5f,0.0f,z-0.5f);                
            glEnd();
        }
    }
}


/***** TRIM ****/
std::string trim(const std::string& str, const std::string& whitespace){
    const size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const size_t strEnd = str.find_last_not_of(whitespace);
    const size_t strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


bool get_valid_line(std::istream &ifile,std::string &line){
    do{
        std::getline(ifile, line);
        line = trim(line); // trim leading and trailing white spaces
    }while(ifile.good() && (line[0]=='#'  || line.size()==0));
    
    return ifile.good();
}


//TODO: this doesn't handle correctly the filepath
bool compare_extension(const std::string &filename, 
					   const std::string &extension){
						   
	return filename.substr(filename.find_last_of(".") + 1) == extension;
}

