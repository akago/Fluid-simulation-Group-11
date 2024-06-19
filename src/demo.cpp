/*
  ======================================================================
   demo.c --- protoype to show off the simple solver
  ----------------------------------------------------------------------
   Author : Jos Stam (jstam@aw.sgi.com)
   Creation Date : Jan 9 2003

   Description:

	This code is a simple prototype that demonstrates how to use the
	code provided in my GDC2003 paper entitles "Real-Time Fluid Dynamics
	for Games". This code uses OpenGL and GLUT for graphics and interface

  =======================================================================
*/
#include <gfx/vec2.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#include <iostream>
#include <limits>

#include "ParticleSystem.h"
#include "SampleSystems.h"
#include "Particle.h"
#include "Force.h"
#include "imageio.h"
#include "Constraint.h"
#include "LinearSolver.h"
#include "RigidBody.h"


/* macros */

#define IX(i,j) ((i)+(N+2)*(j))


/* external definitions (from solver.c) */

//extern void dens_step ( int N, float * x, float * x0, float * u, float * v, int * mask, float diff, float dt );
//extern void vel_step ( int N, float * u, float * v, float * u0, float * v0, int * mask, float visc, float dt );

/* global variables */

static int N;
static float dt, diff, visc;
static float force, source;
static int dvel, drb;

static float * u, * v, * u_prev, * v_prev;
static float * dens, * dens_prev;
static int * internal_bd;
static float* bnd_vel_u, * bnd_vel_v;

static Particle *mouseParticle;
static SpringForce *mouseSpringForce;
static double mouse_kd;
static double mouse_ks;
static ParticleSystem* particleSystem;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int omx, omy, mx, my;


/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/


static void free_data ( void )
{
	if ( u ) free ( u );
	if ( v ) free ( v );
	if ( u_prev ) free ( u_prev );
	if ( v_prev ) free ( v_prev );
	if ( dens ) free ( dens );
	if ( dens_prev ) free ( dens_prev );
	if ( internal_bd ) free ( internal_bd );
	if (bnd_vel_u) free(bnd_vel_u);
	if (bnd_vel_v) free(bnd_vel_v);

	delete particleSystem;
	delete mouseParticle;
	delete mouseSpringForce;
	delete Constraint::GlobalJ;
	delete Constraint::GlobalJdot;
}

static void init_internal_bnd( void ){
	int i, size=(N+2)*(N+2);

	memset(internal_bd, 0, size * sizeof(int));
	for ( i=0; i<=N+1 ; i++ ) {
		internal_bd[IX(0  ,i)] = 1;
		internal_bd[IX(N+1,i)] = 1;
		internal_bd[IX(i,0  )] = 1;
		internal_bd[IX(i,N+1)] = 1;
	}
	
}

static void init_particle_system(void)
{
	particleSystem = rigid1(dt);
	mouse_kd = 0.1;
	mouse_ks = 0.001;
	particleSystem->reset();
}

static void reset_bndAndvel(void) {
	int i, size = (N + 2)*(N + 2);

	for (i = 0; i < size; i++) {
		bnd_vel_u[i] = 0.0f;
		bnd_vel_v[i] = 0.0f;
		if (internal_bd[i] > 1) {
			internal_bd[i] = 0;
		}
	}
}

static void clear_data ( void )
{
	int i, size=(N+2)*(N+2);
	
	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
	}
	init_internal_bnd();
	particleSystem->reset();
}

static int allocate_data ( void )
{
	int size = (N+2)*(N+2);

	u			= (float *) malloc ( size*sizeof(float) );
	v			= (float *) malloc ( size*sizeof(float) );
	u_prev		= (float *) malloc ( size*sizeof(float) );
	v_prev		= (float *) malloc ( size*sizeof(float) );
	dens		= (float *) malloc ( size*sizeof(float) );	
	dens_prev	= (float *) malloc ( size*sizeof(float) );
	internal_bd = (int *) malloc ( size*sizeof(int) );
	bnd_vel_u = (float *)malloc(size * sizeof(float));
	bnd_vel_v = (float *)malloc(size * sizeof(float));
	mouseParticle = new Particle(Vec2f(0, 0));

	if ( !u || !v || !u_prev || !v_prev || !dens || !dens_prev || !internal_bd || !bnd_vel_u || !bnd_vel_v) {
		fprintf ( stderr, "cannot allocate data\n" );
		return ( 0 );
	}

	init_internal_bnd();
	init_particle_system();
	reset_bndAndvel();

	return ( 1 );
}


/*
  ----------------------------------------------------------------------
   OpenGL specific drawing routines
  ----------------------------------------------------------------------
*/

static void pre_display ( void )
{
	glViewport ( 0, 0, win_x, win_y );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ();
	gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
}

static void post_display ( void )
{
	glutSwapBuffers ();
}

static void draw_velocity ( void )
{
	int i, j;
	float x, y, h;

	h = 1.0f/N;

	glColor3f ( 1.0f, 1.0f, 1.0f );
	glLineWidth ( 1.0f );

	glBegin ( GL_LINES );

		for ( i=1 ; i<=N ; i++ ) {
			x = (i-0.5f)*h;
			for ( j=1 ; j<=N ; j++ ) {
				y = (j-0.5f)*h;

				glVertex2f ( x, y );
				glVertex2f ( x+u[IX(i,j)], y+v[IX(i,j)] );
			}
		}

	glEnd ();
}

static void draw_density ( void )
{
	int i, j;
	float x, y, h, d00, d01, d10, d11;

	h = 1.0f/N;

	glBegin ( GL_QUADS );

		for ( i=0 ; i<=N ; i++ ) {
			x = (i-0.5f)*h;
			for ( j=0 ; j<=N ; j++ ) {
				y = (j-0.5f)*h;
				if (internal_bd[IX(i,j)] == 0){

					d00 = dens[IX(i,j)];
					d01 = dens[IX(i,j+1)]; 
					d10 = dens[IX(i+1,j)];
					d11 = dens[IX(i+1,j+1)];

					glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
					glColor3f ( d10, d10, d10 ); glVertex2f ( x+h, y );
					glColor3f ( d11, d11, d11 ); glVertex2f ( x+h, y+h );
					glColor3f ( d01, d01, d01 ); glVertex2f ( x, y+h );
				}
				else if (internal_bd[IX(i, j)] == 2) {
					glColor3f(1.0f, 0.0f, 0.0f);
					glVertex2f(x, y);
					glVertex2f(x + h, y);
					glVertex2f(x + h, y + h);
					glVertex2f(x, y + h);
				}
				else if (internal_bd[IX(i, j)] == 1) {
					glColor3f(0.5f, 1.0f, 0.5f); 
					glVertex2f ( x, y );
					glVertex2f ( x+h, y );
					glVertex2f ( x+h, y+h );
					glVertex2f ( x, y+h );
				}
			}
		}

	glEnd ();
}

/*
  ----------------------------------------------------------------------
   relates mouse movements to force sources
  ----------------------------------------------------------------------
*/

static void get_from_UI ( float * d, float * u, float * v, int* int_bnd_mask)
{
	int i, j, size = (N+2)*(N+2);

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = d[i] = 0.0f;
	}

	if ( !mouse_down[0] && !mouse_down[2] && !mouse_down[1] ) return;

	i = (int)((       mx /(float)win_x)*N+1);
	j = (int)(((win_y-my)/(float)win_y)*N+1);

	if ( i<1 || i>N || j<1 || j>N ) return;

	if ( mouse_down[0] && !drb ) {
		u[IX(i,j)] = force * (mx-omx);
		v[IX(i,j)] = force * (omy-my);
	}

	if ( mouse_down[1] ) {
		int_bnd_mask[IX(i,j)] = 1;
		int_bnd_mask[IX(i+1,j)] = 1;
		int_bnd_mask[IX(i,j+1)] = 1;
		int_bnd_mask[IX(i+1,j+1)] = 1;
	}

	if ( mouse_down[2] ) {
		d[IX(i,j)] = source;
	}

	omx = mx;
	omy = my;

	return;
}

/*
  ----------------------------------------------------------------------
   GLUT callback routines
  ----------------------------------------------------------------------
*/

static void key_func ( unsigned char key, int x, int y )
{
	switch ( key )
	{
		case 'c':
		case 'C':
			clear_data ();
			break;

		case 'q':
		case 'Q':
			free_data ();
			exit ( 0 );
			break;

		case 'v':
		case 'V':
			dvel = !dvel;
			break;
		case 'x':
		case 'X':
			drb = !drb;
			break;
	}
}

static void mouse_func ( int button, int state, int x, int y )
{
	omx = mx = x;
	omx = my = y;

	mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func ( int x, int y )
{
	mx = x;
	my = y;
}

static void reshape_func ( int width, int height )
{
	glutSetWindow ( win_id );
	glutReshapeWindow ( width, height );

	win_x = width;
	win_y = height;
}

bool isbdvelallzero() {
	int size = (N + 2)*(N + 2);
	for (int i = 0; i < size; i++) {
		if (bnd_vel_u[i] != 0.0 || bnd_vel_v[i] != 0.0) {
			return false;
		}
	}
	return true;
}
bool isvelallzero() {
	int size = (N + 2)*(N + 2);
	for (int i = 0; i < size; i++) {
		printf("%d,%d\n", u[i], v[i]);
		if (u[i] != 0.0 || v[i] != 0.0) {
			return false;
		}
	}
	return true;
}

static void idle_func ( void )
{
	
	get_from_UI ( dens_prev, u_prev, v_prev, internal_bd);
	reset_bndAndvel();
	if (mouse_down[0] && drb) {
		mouseParticle->m_Position[0] = (2.0*mx / win_x) - 1;
		mouseParticle->m_Position[1] = -(2.0*my / win_y) + 1;

		Particle *closestParticle;
		float closestDistanceSquared = std::numeric_limits<float>::max();

		for (auto p : particleSystem->getParticles()) {
			float dx = (mouseParticle->m_Position[0] - p->m_Position[0]);
			float dy = (mouseParticle->m_Position[1] - p->m_Position[1]);
			float distanceSquared = dx * dx + dy * dy;
			if (distanceSquared < closestDistanceSquared) {
				closestParticle = p;
				closestDistanceSquared = distanceSquared;
			}
		}
		// Also consider rb vertices
		for (auto rb : particleSystem->getRigids()) {
			for (auto p : rb->m_Vertices) {
				float dx = (mouseParticle->m_Position[0] - (p->m_Position[0] + rb->m_Position[0]));
				float dy = (mouseParticle->m_Position[1] - (p->m_Position[1] + rb->m_Position[1]));
				float distanceSquared = dx * dx + dy * dy;
				if (distanceSquared < closestDistanceSquared) {
					closestParticle = p;
					closestDistanceSquared = distanceSquared;
				}
			}
		}

		mouseSpringForce = new SpringForce(mouseParticle, closestParticle, 0, mouse_kd, mouse_ks);
		
		particleSystem->addForce(mouseSpringForce);
		particleSystem->simulationStep();
		particleSystem->removeLastForce();
		delete mouseSpringForce;
	}
	else {
		particleSystem->simulationStep();
	}
	
	particleSystem->projectRigidBodies(N, internal_bd, bnd_vel_u, bnd_vel_v);
	
	vel_step ( N, u, v, u_prev, v_prev, internal_bd, bnd_vel_u, bnd_vel_v, visc, dt );
	dens_step ( N, dens, dens_prev, u, v, internal_bd, bnd_vel_u, bnd_vel_v, diff, dt );

	glutSetWindow ( win_id );
	glutPostRedisplay ();
}

static void display_func ( void )
{
	pre_display ();

		if ( dvel ) draw_velocity ();
		else		draw_density ();

		particleSystem->drawWalls();
		particleSystem->drawConstraints();
		particleSystem->drawForces();
		particleSystem->drawParticles();
		particleSystem->drawRigids();

	post_display ();
}


/*
  ----------------------------------------------------------------------
   open_glut_window --- open a glut compatible window and set callbacks
  ----------------------------------------------------------------------
*/

static void open_glut_window ( void )
{
	glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

	glutInitWindowPosition ( 0, 0 );
	glutInitWindowSize ( win_x, win_y );
	win_id = glutCreateWindow ( "Alias | wavefront" );

	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();

	pre_display ();

	glutKeyboardFunc ( key_func );
	glutMouseFunc ( mouse_func );
	glutMotionFunc ( motion_func );
	glutReshapeFunc ( reshape_func );
	glutIdleFunc ( idle_func );
	glutDisplayFunc ( display_func );
}


/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main ( int argc, char ** argv )
{
	glutInit ( &argc, argv );

	if ( argc != 1 && argc != 6 ) {
		fprintf ( stderr, "usage : %s N dt diff visc force source\n", argv[0] );
		fprintf ( stderr, "where:\n" );\
		fprintf ( stderr, "\t N      : grid resolution\n" );
		fprintf ( stderr, "\t dt     : time step\n" );
		fprintf ( stderr, "\t diff   : diffusion rate of the density\n" );
		fprintf ( stderr, "\t visc   : viscosity of the fluid\n" );
		fprintf ( stderr, "\t force  : scales the mouse movement that generate a force\n" );
		fprintf ( stderr, "\t source : amount of density that will be deposited\n" );
		exit ( 1 );
	}

	if ( argc == 1 ) {
		N = 160;
		dt = 0.1f;
		diff = 0.0f;
		visc = 0.0f;
		force = 5.0f;
		source = 100.0f;
		fprintf ( stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n",
			N, dt, diff, visc, force, source );
	} else {
		N = atoi(argv[1]);
		dt = atof(argv[2]);
		diff = atof(argv[3]);
		visc = atof(argv[4]);
		force = atof(argv[5]);
		source = atof(argv[6]);
	}

	printf ( "\n\nHow to use this demo:\n\n" );
	printf ( "\t Add densities with the right mouse button\n" );
	printf ( "\t Add boundaries with the middle mouse button\n" );
	printf("\t Use the 'x' key to switch between using the left mouse button to control the rigid body or add velocity. ");
	printf ( "\t Toggle density/velocity display with the 'v' key\n" );
	printf ( "\t Clear the simulation by pressing the 'c' key\n" );
	printf ( "\t Quit by pressing the 'q' key\n" );

	dvel = 0;
	drb = 0;

	if ( !allocate_data () ) exit ( 1 );
	clear_data ();

	win_x = 1024;
	win_y = 1024;
	open_glut_window ();

	glutMainLoop ();

	exit ( 0 );
}