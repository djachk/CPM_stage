#ifndef ALLOCATE_DEF
#define ALLOCATE_DEF

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "operation.h"

typedef struct cell {
	double *xcoord;
	double *ycoord;
	double *temp;
	int *targetarea;
	int *area;
	double *area_constraint;
	int *celltype;
	int *perimeter;
	int *surf_energy; //cells.sur_energy=J*cells.perimeter pour une mousse, mais pas pour un tissu.
	//int *oldperimeter; Unused
	int *nneighbours;
	int **neighbours;
	int *t_debut_division;
	int *interphase;
	int *vient_de_diviser;
	int *vient_de_diviser_pour_affichage;
	int *vient_de_naitre;
	int *duree_de_vie;
	int *morte;
	int *vient_de_mourir;
} Cell;

int* AllocateLineInterface(int n);

Cell AllocateCells(int n, int maxneighbour);
void FreeCells(Cell cells, int n);
void Duplicate(Cell copie, Cell original, int maxcells);

int PutCell(TYPE **plane, int y, int x, TYPE m, int ncol, int nrow, int side1, int side2);
void InitBubblePlane(int init_config, float fillfactor,int nrow,int ncol, int target_area, double a1, double a2, TYPE **state, Cell cells, int sliding, double area_constraint1, 
	int interphase1, int duree_de_vie1, int* nb_cellules, int* nb_cellules_vivantes,int* nb_cellules1, int* nb_cellules2, int maxcells,  int division_cellulaire, int apoptose);
int AssignNormalTargetarea(int mean, double mu2adim, int minimum);
int GeneratePolydispersity(int polydispersity, int blob, int maxcells, double fillfactor, int nrow, int ncol, int target_area, double targetareamu2, int target_area2, double alpha, Cell cells, double area_constraint2, 
	int interphase2, int duree_de_vie2, int* nb_cellules, int* nb_cellules1, int* nb_cellules2);

#endif
