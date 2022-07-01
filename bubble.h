#ifndef BUBBLE_DEF

#define BUBBLE_DEF
#define M_PI 3.14159265358979323846

#include <math.h>
#include <limits.h>
#include <omp.h>
#include <string.h>

#include "allocate.h"

typedef struct point point;

struct point {
    int x;
    int y;
};


double BubbleHamiltonian(int time, int dispersetime, int num, int neighbour_energy, int neighbour_copy, int neighbour_connected, int ncol, int nrow, TYPE **state, int mediumcell, int heat_bath, int** Jarray, Cell cells, double area_constraint , double temperature,
    int* nb_cellules_vivantes, int* nb_cellules_mortes, int* nb_cellules_condamnees, 
    int* nb_cellules1_condamnees, int* nb_cellules2_condamnees, int* nb_cellules1, int* nb_cellules2, int maxcells, int* pile_labels_libres, int* taille_labels_libres, int* commencer_division, int* nb_cellules_tuees);

bool test_only_two(TYPE **state, bool condwrap, int neighbour_copy, int pixel, int icandidate, int ncol, int nrow, int x, int y);

bool connected_4(TYPE **state, int nb_nei_id, bool condwrap, int pixel, int ncol, int nrow, int x, int y);

bool connected_8(TYPE **state, int nb_nei_id, bool condwrap, int pixel, int ncol, int nrow, int x, int y);

void ComputeLineInterface(int ttime, FILE* interfacefp, Cell cells, int ncol, int nrow, TYPE **state, int* line_interface);

void dessiner_tableau(int taille, int tab[taille][taille], int x0, int y0, int num_cell);
point trouver_point_dedans_boucle(int taille, int tab[taille][taille], int num_cell);
int tester_connexe_boucle(int taille_tab, int tab[taille_tab][taille_tab], int num_cell, int x_deb_tache, int y_deb_tache, int area_cell );
int tester_cellule_connexe_boucle(Cell cells,TYPE** state, int num_cell, int cote_carre, int nrow, int ncol);


int surface_par_type_cellule(Cell cells, int type_cellule, int nb_cellules);

void Diviser(Cell cells, int num_cell, int ttime, int nrow, int ncol, TYPE** state, int* nb_cellules, int* nb_cellules_vivantes, 
    int* nb_cellules_mortes, int* nb_cellules_par_division, int* nb_cellules1_par_division, int* nb_cellules2_par_division, int* nb_cellules1, int* nb_cellules2, 
    int* nb_cellules_tuees, int* nb_cellules_condamnees, int* nb_cellules1_condamnees, int* nb_cellules2_condamnees, int maxcells, int duree_de_vie1, 
    int duree_de_vie2, int interphase1, int interphase2, int* pile_labels_libres, int* taille_labels_libres, int cote_carre, int* nb_divisions_acceptees, int* nb_divisions_refusees);

void Condamner_cellule(Cell cells, int num_cell, int* nb_cellules_condamnees, int *nb_cellules1_condamnees, int *nb_cellules2_condamnees, int* nb_cellules_vivantes);

void Empiler_cellule_vide(Cell cells, int num_cell, int* nb_cellules_vivantes, 
	int* nb_cellules_mortes, int* nb_cellules1, int* nb_cellules2, 
    int* nb_cellules_tuees, int * nb_cellules_condamnees, int * nb_cellules1_condamnees, int * nb_cellules2_condamnees, int maxcells, int* pile_labels_libres, int* taille_labels_libres, int origine);

void ComputePerimeter(int maxcells, Cell cells, int ncol, int nrow, TYPE **state, int mediumcell, int neighbour_energy);
void FindNeighbours(int maxcells, Cell cells, int ncol, int nrow, TYPE **state, int mediumcell, int neighbour_connnected, int maxneighbours, int* side_interf12, int* side_interf10, int* side_interf20);
double ComputeEnergy(int maxcells, Cell cells, int ncol, int nrow, TYPE **state, int neighbour_energy, int** Jarray, double area_constraint);
void ComputeCenterCoords(Cell cells, int ncol, int nrow, TYPE **state, int nbrcells, int mediumcell);
int ComputeBoundary(Cell cells, int ncol, int nrow, TYPE **state, int neighbour_energy);
void AffichageCouleurs(int affichage, Cell cells, int ncol, int nrow, char *subdirname, char *subdirnameRAW, TYPE **state, TYPE **nstate, int montrer_division);
#endif
