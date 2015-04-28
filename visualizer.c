

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
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#define PI 3.14159265358979f
#define MAX_LINE 1024


struct vec3f_tag{ GLfloat x,y,z; };
struct vec3i_tag{ int x,y,z; };
typedef struct vec3f_tag vec3f;
typedef struct vec3i_tag vec3i;

vec3f make_vec3f(GLfloat x,GLfloat y,GLfloat z);
vec3i make_vec3i(int x,int y,int z);
vec3f addv(const vec3f *v1, const vec3f *v2);
GLfloat sq_dist(const vec3f *v1, const vec3f *v2);
GLfloat my_fminf(GLfloat a,GLfloat b);

vec3f *verts   = NULL; 
vec3f *normals = NULL; 
vec3f *X1      = NULL; 
vec3f *X2      = NULL; 
GLfloat *rho1  = NULL;
GLfloat *rho2  = NULL;
vec3i *faces   = NULL;
int num_verts = 0;
int num_faces = 0;
GLfloat max_k;
GLfloat min_k;
GLfloat min_len;

struct camera{
  vec3f    pos;
  GLfloat elevation;
  GLfloat azimut;
  GLfloat zangle;
}cam;

int keyb[256]; //keyboard list

void draw_reference_frame();

void draw_plane(int Nx,int Nz);

void init(void);

void load_model(const char *filename);

void draw_model(void);

void display(void);

void reshape (int w, int h);

void keyboardUp(unsigned char key, int x, int y);

void keyboardDown(unsigned char key, int x, int y);

void loop(int value);

int get_line(FILE *ifile, char *line, int max_line);

void trim(char *line);

int main(int argc, char** argv){	
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB| GLUT_DEPTH);
	glutInitWindowSize (1024, 768); 
	glutInitWindowPosition (100, 100);
	glutCreateWindow (argv[0]);
	init ();
	if(argc >1) load_model(argv[1]);

	glutDisplayFunc(display); 
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboardDown);
	glutKeyboardUpFunc(keyboardUp);
	glutTimerFunc(50,loop, 0);
	glutMainLoop();
	return 0;
}

int get_valid_line(FILE *ifile,char *line,int max_line){
	int len;
	while((len = get_line(ifile,line,max_line)) > 0){
		trim(line);
		
		if(line[0] != '#' && strlen(line) > 0)
			break; //skip comment
	}
	return len;
}

void load_model(const char *filename){
	FILE *ifile;
	char *line;
	int i ,len;
	line = malloc(sizeof(char)*MAX_LINE);
	ifile = fopen(filename,"r");
	len = get_valid_line(ifile,line,MAX_LINE);
	sscanf(line,"%d %d",&num_verts,&num_faces);
	printf("num_verts: %d num_faces %d\n",num_verts,num_faces);
	verts = malloc(num_verts*sizeof(vec3f));
	normals = malloc(num_verts*sizeof(vec3f));
	X1 = malloc(num_verts*sizeof(vec3f));
	X2 = malloc(num_verts*sizeof(vec3f));
	rho1 = malloc(num_verts*sizeof(GLfloat));
	rho2 = malloc(num_verts*sizeof(GLfloat));
	faces = malloc(num_faces*sizeof(vec3i));
	for(i =0; i< num_verts; i++){
		GLfloat x,y,z, nx,ny,nz, X1_x,X1_y,X1_z, rrho1,rrho2;
		
		len = get_valid_line(ifile,line,MAX_LINE);
		assert(len >0);
		//printf("%s\n",line);
		sscanf(line,"%f %f %f    %f %f %f    %f %f %f  %f %f",&x,&y,&z,&nx,&ny,&nz,&X1_x,&X1_y,&X1_z,&rrho1,&rrho2);
		verts[i]   = make_vec3f(x,y,z);
		normals[i] = make_vec3f(nx,ny,nz);
		X1[i]      = make_vec3f(X1_x,X1_y,X1_z);
		rho1[i]    = 1.0f/rrho1;
		rho2[i]    = 1.0f/rrho2;
		
		//printf("%f %f %f\n",verts[i].x,verts[i].y,verts[i].z);
	}
	//compute the min and the max curvature 
	min_k = max_k = 1/(rho1[0]*rho2[0]);
	for(i =0; i< num_verts; i++){
		GLfloat  k = 1/(rho1[i]*rho2[i]);
		if(k < min_k ) min_k = k;
		if(k > max_k ) max_k = k;
	}

	printf("min k: %f max k %f\n",min_k,max_k);
	
	for(i =0; i< num_faces; i++){
		int ia,ib,ic;
		len = get_valid_line(ifile,line,MAX_LINE);
		//printf("%s\n",line);

		assert(len >0);

		sscanf(line,"%d %d %d   ",&ia,&ib,&ic);
		faces[i] = make_vec3i(ia,ib,ic);
		//printf("%d %d %d\n",faces[i].x,faces[i].y,faces[i].z);
	}
	
    min_len = sq_dist(&verts[faces[0].x], &verts[faces[1].y]);
    //compute the min edge_size
    for(i =0; i< num_faces; i++){        
        
        GLfloat sq_len_a = sq_dist(&verts[faces[i].x], &verts[faces[i].y]);
        GLfloat sq_len_b = sq_dist(&verts[faces[i].y], &verts[faces[i].z]);
        GLfloat sq_len_c = sq_dist(&verts[faces[i].z], &verts[faces[i].x]);
        
        min_len = my_fminf(min_len,my_fminf(sq_len_a,my_fminf(sq_len_b,sq_len_c)));
        
    }
    min_len = sqrtf(min_len);
    printf("min edge len: %f \n",min_len);
	fclose(ifile);
	free(line);
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
	int i;
    GLfloat linelen = min_len*0.7f;
	if(verts != NULL){
		glBegin(GL_TRIANGLES);

		for(i = 0; i < num_faces; i++){
			vec3f  a = verts[faces[i].x];
			vec3f  b = verts[faces[i].y];
			vec3f  c = verts[faces[i].z];
			GLfloat ka = 1/(rho1[faces[i].x]*rho2[faces[i].x]);
			GLfloat kb = 1/(rho1[faces[i].y]*rho2[faces[i].y]);
			GLfloat kc = 1/(rho1[faces[i].z]*rho2[faces[i].z]);
			vec3f ca = get_colour(ka,min_k,max_k);
			vec3f cb = get_colour(kb,min_k,max_k);
			vec3f cc = get_colour(kc,min_k,max_k);
			//printf("%f %f %f \n",ca.x, ca.y, ca.z);

			glColor3f (ca.x, ca.y, ca.z);
			glVertex3f(a.x, a.y, a.z);
					
     		glColor3f (cb.x, cb.y, cb.z);
			glVertex3f(b.x, b.y, b.z);
			
			glColor3f (cc.x, cc.y, cc.z);
			glVertex3f(c.x, c.y, c.z);

		}	
		glEnd();
        //draw normals
        glLineWidth(1.5f);
        
        glBegin(GL_LINES);
        for(i =0; i< num_verts; i++){ 
			vec3f  a = verts[i];
   			vec3f  n = normals[i];        
            glColor3f (0.5f, 1.0f, 1.0f);  
            glVertex3f(a.x,a.y,a.z);
            glVertex3f(a.x+linelen*n.x,a.y+linelen*n.y,a.z+linelen*n.z);
        }
        glEnd();
        //draw_X1
        glBegin(GL_LINES);
        for(i =0; i< num_verts; i++){ 
			vec3f  a = verts[i];
   			vec3f  x1 = X1[i];        
            glColor3f (0.1f, 0.1f, 0.1f);  
            glVertex3f(a.x,a.y,a.z);
            glVertex3f(a.x+linelen*x1.x,a.y+linelen*x1.y,a.z+linelen*x1.z);
        }
        glEnd();
        
        glLineWidth(1.0f);
        
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
	keyb[key] = 0;
}
void keyboardDown(unsigned char key, int x, int y){

  keyb[key] = 1;   
}

void loop(int value){
    GLfloat dx,dz;
    if(keyb['w']){
        dx = sinf(cam.azimut*PI/180.0f)/10.0f;
        dz = -cosf(cam.azimut*PI/180.0f)/10.0f;
        cam.pos.x +=dx;
        cam.pos.z +=dz;
    }
    if(keyb['s']){
        dx = sinf(cam.azimut*PI/180.0f)/10.0f;
        dz = -cosf(cam.azimut*PI/180.0f)/10.0f;
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


int get_line(FILE *ifile, char *line, int max_line){
	int c,i;
	i = 0;
	while( --max_line > 0 && (c = getc(ifile))!= EOF && c != '\n'){
		line[i++] = c;
	}
	if(c == '\n')
		line[i++] = c;
	
	line[i] = '\0';
	return i;
}

void trim(char *line){
	int p,i;
	p = i = 0;
	while(line[p] && isspace(line[p])) p++;	
	while(line[i] = line[p]){ i++; p++;};
	do{i--;}while(isspace(line[i]));
	line[++i] = '\0';
}



/*** VEC3f library */
vec3f make_vec3f(GLfloat x,GLfloat y,GLfloat z){
	vec3f v = {x,y,z}; 
	return v; 
}

vec3i make_vec3i(int x,int y,int z){
	vec3i v = {x,y,z}; 
	return v; 
}

vec3f addv(const vec3f *v1, const vec3f *v2){	
	return  make_vec3f(v1->x+v2->x,v1->y+v2->y,v1->z+v2->z); 
}

GLfloat sq_dist(const vec3f *v1, const vec3f *v2){
    return (v1->x - v2->x)*(v1->x - v2->x) + 
           (v1->y - v2->y)*(v1->y - v2->y) +
           (v1->z - v2->z)*(v1->z - v2->z);
}

GLfloat my_fminf(GLfloat a,GLfloat b){
    return (a <b)?a:b;
}

