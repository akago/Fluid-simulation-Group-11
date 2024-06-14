#include <stdlib.h>
#include <stdbool.h>

void add_source ( int N, float * x, float * s, float dt );
void set_bnd (  int N, int b, float * x, bool* mask );
void lin_solve ( int N, int b, float * x, float * x0, bool*mask, float a, float c  );
void diffuse ( int N, int b, float * x, float * x0, bool* mask, float diff, float dt );
void advect ( int N, int b, float * d, float * d0, float * u, float * v, bool* mask, float dt );
void project ( int N, float * u, float * v, float * p, float * div, bool* mask );
void vorticity_confinement(int N, float * u, float * v, float * u0, float * v0, float epsilon);


void dens_step ( int N, float * x, float * x0, float * u, float * v, bool * mask, float diff, float dt );
void vel_step ( int N, float * u, float * v, float * u0, float * v0, bool * mask, float visc, float dt );
