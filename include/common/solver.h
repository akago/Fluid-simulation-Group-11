#include <stdlib.h>

void add_source ( int N, float * x, float * s, float dt );
void set_bnd ( int N, int b, float * x );
void lin_solve ( int N, int b, float * x, float * x0, float a, float c );
void diffuse ( int N, int b, float * x, float * x0, float diff, float dt );
void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt );
void project ( int N, float * u, float * v, float * p, float * div );
void vorticity_confinement(int N, float * u, float * v, float * u0, float * v0, float epsilon);


void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt );
void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt );
