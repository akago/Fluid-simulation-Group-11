#include "Constraint.h"

GlobalMatrix* Constraint::GlobalJ = nullptr;
GlobalMatrix* Constraint::GlobalJdot = nullptr;
int Constraint::global_cons_num = 0;
double Constraint::ks = 0.01;
double Constraint::kd = 0.02;

