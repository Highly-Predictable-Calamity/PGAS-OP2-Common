//
// auto-generated by op2.py
//

// global constants
extern float gam;
extern float gm1;
extern float cfl;
extern float eps;
extern float mach;
extern float alpha;
extern float qinf[4];

// header
#include "op_lib_c.h"

void op_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){}
// user kernel files
#include "save_soln_acckernel.c"
#include "adt_calc_acckernel.c"
#include "res_calc_acckernel.c"
#include "bres_calc_acckernel.c"
#include "update_acckernel.c"
