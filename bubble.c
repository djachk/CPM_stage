#include "bubble.h"

//             reseau         4     |  6  |  8  |          16           |    20     |                 36
static const int nx[36]={ 0, 0, 1,-1, 1,-1, 1,-1, 1, 1,-1,-1, 2, 2,-2,-2, 0, 0, 2,-2, 2,-2, 0,-1,-2,-3,-3,-3,-3, 0, 1, 2, 3, 3, 3, 3};
static const int ny[36]={ 1,-1, 0, 0, 1,-1,-1, 1, 2,-2, 2,-2, 1,-1, 1,-1, 2,-2, 0, 0, 2,-2,-3,-3,-3,-3, 0,-1,-2, 3, 3, 3, 3, 0, 1, 2};

//distance between central pixel and neigborhood (neighbour_energy=8 only for now)
//static const double distance[8]={1,1,1,1,sqrt(2),sqrt(2),sqrt(2),sqrt(2)};


double BubbleHamiltonian(int time, int dispersetime, int num, int neighbour_energy, int neighbour_copy, int neighbour_connected, int ncol, int nrow, TYPE **state, int mediumcell, int heat_bath, int** Jarray, Cell cells, double area_constraint , double temperature,
			int* nb_cellules_vivantes, int* nb_cellules_mortes, int* nb_cellules_condamnees, 
			int* nb_cellules1_condamnees, int* nb_cellules2_condamnees, int* nb_cellules1, int* nb_cellules2, int maxcells, int* pile_labels_libres, int* taille_labels_libres, int* commencer_division, int* nb_cellules_tuees)
{
	bool sub_flag, nei_flag;
	double sum_E=0;
	int pixel, icandidate, ineighbour, y, x, pos_neg;
	int state_copy[neighbour_copy];
	
	for(int i=0;i<ncol*nrow;i++) {
		pixel=state[y=(int)(aleatoire(num)*nrow) +1][x=(int)(aleatoire(num)*ncol) +1]; //CASH PLANE are from 1 to ncol or nrow.
		if (pixel==-1) printf("ALERTE PIXEL APRES STATE %d\n",pixel);

		//flag to pass PeriodicWrap funtion if far of boundary.
		bool condwrap = true;
		if (x>neighbour_energy/4+2 && x<ncol-2-neighbour_energy/4 && y>neighbour_energy/4+2 && y<nrow-2-neighbour_energy/4)
			condwrap = false;
		
		pos_neg = 0 ; 
		
		
		// State of neighbours.
		int nb_nei_id = 0; // nombre de voisins identiques au pixel central
		if (condwrap) {
			for (int k=0 ; k < neighbour_copy ; k++) {
				if (state[PeriodicWrap(y+ny[k],nrow)][PeriodicWrap(x+nx[k],ncol)] == pixel) {
					nb_nei_id += 1;
				}
			}
		}
		else {
			for (int k=0 ; k < neighbour_copy ; k++) {
				if (state[y+ny[k]][x+nx[k]] == pixel) {
					nb_nei_id += 1;
				}
			}
		}
		if (nb_nei_id < neighbour_copy) { //on ne perd pas de temps si tous les voisins sont identiques au central (pixel non situé à une interface)
			for (int k=0 ; k<neighbour_copy ; k++) {
//				possible_state_copy[k][0] = mediumcell-1;
				state_copy[k]=-1;
			}
			int s;
			//pos_neg contient le nombre de pixels identiques au pixel central dans le voisinage de copie
			if (condwrap) {
				//state_copy[0]=state[PeriodicWrap(y+ny[0], nrow)][PeriodicWrap(x+nx[0], ncol)];
				for (int iposs=0 ; iposs<neighbour_copy ; iposs++) { 
					s = state[PeriodicWrap(y+ny[iposs], nrow)][PeriodicWrap(x+nx[iposs], ncol)];
					//if (cells.condamnee[s]) continue;
					for( int jposs=0; jposs<pos_neg+1; jposs++) {
						if (s == state_copy[jposs]) //valeur déjà répertoriée
							break;
						if (state_copy[jposs] == -1){ //valeur à ajouter
							state_copy[jposs] = s;
							pos_neg++;
							break;
						}
					}
					//AddIfNotIn(possible_state_copy, s, &pos_neg);
				}
			}
			else {
				//state_copy[0]=state[y+ny[0]][x+nx[0]];
				for (int iposs=0 ; iposs<neighbour_copy ; iposs++) {
					s = state[y+ny[iposs]][x+nx[iposs]];
					//if (cells.condamnee[s]) continue;
					for( int jposs=0; jposs<pos_neg+1; jposs++) {
						if (s == state_copy[jposs]) //valeur déjà répertoriée
							break;
						if (state_copy[jposs] == -1){ //valeur à ajouter
							state_copy[jposs] = s;
							pos_neg++;
							break;
						}
					}
					//AddIfNotIn(possible_state_copy, s, &pos_neg);
				}
			}
			//flag with value true if local connectivity (so, can try a metropolis step).
			nei_flag = true;
			if (pixel != mediumcell) {
			//pixel != mediumcell => mediumcell peut se fragmenter. Nécessaire pour la résorption aux premiers instants. 
			//Egalement nécessaire pour un système non confluent avec cellules pouvant s'agréger. 
				//on teste la connexité locale de la bulle candidate:
				nei_flag = connected_4(state, nb_nei_id, condwrap, pixel, ncol, nrow, x, y);
			}
			if (nei_flag) {						
				double delta_E = 0;
				int delta_surf_E=0;
				int delta_perimeter=0;
				icandidate=state_copy[(int)(aleatoire(num)*pos_neg)];
				if (icandidate != pixel) { //else, delta_E = 0 and we don't make calcul.
					if (condwrap) {
						for(int j=0;j<neighbour_energy;j++){
							ineighbour = state[PeriodicWrap(y+ny[j],nrow)][PeriodicWrap(x+nx[j],ncol)];
							if(pixel != ineighbour) {
								delta_E -= Jarray[cells.celltype[pixel]][cells.celltype[ineighbour]];//pixel never saved.
								delta_perimeter--;
							}
							if(icandidate != ineighbour) {
								delta_E += Jarray[cells.celltype[icandidate]][cells.celltype[ineighbour]];
								delta_perimeter++;
							}
						}
					}
					else {
						for(int j=0;j<neighbour_energy;j++){
							ineighbour = state[y+ny[j]][x+nx[j]];
							if(pixel != ineighbour) {
								delta_E -= Jarray[cells.celltype[pixel]][cells.celltype[ineighbour]];//pixel never saved.
								delta_perimeter--;
							}
							if(icandidate != ineighbour) {
								delta_E += Jarray[cells.celltype[icandidate]][cells.celltype[ineighbour]];
								delta_perimeter++;
							}
						}						
					}
					delta_surf_E=(int)(delta_E); //on stocke la variation d'énergie de surface		
					if(cells.targetarea[pixel] > 0){
						delta_E+=(cells.area_constraint[pixel]/cells.targetarea[pixel])*(-2*(cells.area[pixel]-cells.targetarea[pixel])+1);
					}
					if(cells.targetarea[icandidate] > 0){
						delta_E+=(cells.area_constraint[icandidate]/cells.targetarea[icandidate])*(2*(cells.area[icandidate]-cells.targetarea[icandidate])+1);
					}
					sub_flag = false; // flag qui devient true si le changement de valeur du pixel est autorisé par le metropolis step.	
					if ( delta_E < 0. || aleatoire(num) < exp(-delta_E/temperature) ){
						sub_flag = true;
					}
					if (cells.condamnee[pixel]==1) sub_flag=true;  //cellule à détruire
					if(sub_flag){
						//on teste maintenant la connexité locale de la bulle "cible". Si celle-ci ne l'était pas, la modification de la valeur du pixel entrainerait
						//la création d'une arche, et donc celle-ci deviendrait multiplement connexe, ce qu'on veut empêcher.
						//sub_flag=false; // flag qui devient true si la valeur du pixel est effectivement modifiée.
						bool update_flag=false;
						int voisins_identiques = 0;
						if (condwrap) {
							for(int j=0;j<neighbour_connected;j++) {
		                        if (state[PeriodicWrap(y+ny[j], nrow)][PeriodicWrap(x+nx[j], ncol)] == icandidate)
		                            voisins_identiques += 1;
		                    }
		                }	                
		                else {
							for(int j=0;j<neighbour_connected;j++) {
		                        if (state[y+ny[j]][x+nx[j]] == icandidate)
									voisins_identiques += 1;
								}
		                }
	                    if (pixel != mediumcell){
							update_flag=connected_4(state, voisins_identiques, condwrap, icandidate, ncol, nrow, x, y);
						}
						else { //si pixel appartient au médium, on veut qu'il puisse se résorber. Inutile pour un système confluent.
		                    if (connected_4(state, nb_nei_id, condwrap, pixel, ncol, nrow, x, y)){
		                    	update_flag=connected_4(state, voisins_identiques, condwrap, icandidate, ncol, nrow, x, y);	
		                    }
							else { //workaround pour limiter l'apparition de fragments de médium dans les cellules
								update_flag = connected_8(state, voisins_identiques, condwrap, icandidate, ncol, nrow, x, y);
								//workaround alternatif:
								//update_flag = !(test_only_two(state,condwrap, 36, pixel, icandidate, ncol, nrow, x, y));
							}
						}	
						if (cells.condamnee[icandidate]==1) update_flag=false;  //cellule condamnee ne doit pas gagner de terrain
						if (update_flag==true){
							cells.area[pixel]--;
							if (cells.area[0]==0 && !(*commencer_division))  (*commencer_division)=1;
							if (cells.area[pixel]==0) {
								//DC: cette cellule n'a plus aucun pixel...
								Empiler_cellule_vide(cells, pixel, nb_cellules_vivantes, 
									nb_cellules_mortes, nb_cellules1, nb_cellules2, nb_cellules_tuees, nb_cellules_condamnees, 
									nb_cellules1_condamnees, nb_cellules2_condamnees, maxcells, pile_labels_libres, taille_labels_libres, 1);
							}
							cells.area[icandidate]++;							
							cells.surf_energy[pixel] += delta_surf_E;
							cells.surf_energy[icandidate] += delta_surf_E;
							cells.perimeter[pixel] += delta_perimeter;
							cells.perimeter[icandidate] += delta_perimeter;
							//cells.area[mediumcell]=0;							
							sum_E += delta_E;
							state[y][x]=icandidate;
						}
					}
				}
			}
		}
	} 	
	return sum_E;
}


bool test_only_two(TYPE **state, bool condwrap, int neighbour_copy, int pixel, int icandidate, int ncol, int nrow, int x, int y) {
	//check that there are only two values in the neighborhood
	bool only_two = true;
	if (condwrap) {
		for (int k=0 ; k < neighbour_copy ; k++) {
			if ((state[PeriodicWrap(y+ny[k],nrow)][PeriodicWrap(x+nx[k],ncol)] != pixel) || (state[PeriodicWrap(y+ny[k],nrow)][PeriodicWrap(x+nx[k],ncol)] != icandidate)) {
				only_two = false;
				break;
			}
		}
	}
	else {
		for (int k=0 ; k < neighbour_copy ; k++) {
			if ((state[y+ny[k]][x+nx[k]] != pixel) || (state[y+ny[k]][x+nx[k]] != icandidate)) {
				only_two = false;
				break;
			}
		}
	}
	return only_two;	
}

bool connected_4(TYPE **state, int nb_nei_id, bool condwrap, int pixel, int ncol, int nrow, int x, int y) {
	//return true; //-> enlève la connexité.
	bool nei_flag = true;
	if (nb_nei_id > 1) {
		bool N, E, S, O, NE, NO, SO, SE;
		if (condwrap) {
			E=(state[y][PeriodicWrap(x+1, ncol)] == pixel);
			O=(state[y][PeriodicWrap(x-1, ncol)] == pixel);
			S=(state[PeriodicWrap(y+1, nrow)][x] == pixel);
			N=(state[PeriodicWrap(y-1, nrow)][x] == pixel);
			NE=(state[PeriodicWrap(y-1, nrow)][PeriodicWrap(x+1, ncol)] == pixel);
			NO=(state[PeriodicWrap(y-1, nrow)][PeriodicWrap(x-1, ncol)] == pixel);
			SE=(state[PeriodicWrap(y+1, nrow)][PeriodicWrap(x+1, ncol)] == pixel);
			SO=(state[PeriodicWrap(y+1, nrow)][PeriodicWrap(x-1, ncol)] == pixel);
		} else {
			E=(state[y][x+1] == pixel);
			O=(state[y][x-1] == pixel);
			S=(state[y+1][x] == pixel);
			N=(state[y-1][x] == pixel);
			NE=(state[y-1][x+1] == pixel);
			NO=(state[y-1][x-1] == pixel);
			SE=(state[y+1][x+1] == pixel);
			SO=(state[y+1][x-1] == pixel);
		}
		/* Ici, avoir le même état sous-entend « que le nœuds central considéré ».
		 * Cas à 1 voisins : forcement connexe ;
		 * Cas à 4 voisins : non considéré : pas de changement possible ;
		 * Cas à 2 voisins : soit ils sont opposés, donc ce n'est pas localement connexe,
		 * 					 soit ils sont proches, alors si le nœud diagonal entre eux a le même état, c'est connexe. Sinon, non.
		 * Cas à 3 voisins : on regarde quel voisin n'a pas le même état. Si les deux nœuds diagonaux opposés ont le même état,
		 * 					 c'est connexe. Sinon, non.
		 * On cherche les cas **non** connexes pour leur attribuer nei_flag = false.
		 * Comme il est plus facile de regarder les cas connexes, on s'assure que l'on est **pas** dans un cas connexe.
		 * Retour sur le cas à 3 voisins : l'un (et un seul) des quatres cas ne respecte pas la condition par le voisinage direct (E, O, N ou S).
		 * 		  Alors, s'il la respecte par les diagnonaux, c'est connexe, on ne rentre **pas** dans le « if » (d'où le « ! » avant la condition).
		 * 		  S'il ne la respecte pas, c'est bien **non** connexe, on rentre dans le « if ».
		 * */
		if ( (nb_nei_id == 2 && !( (N&&E&&NE) || (N&&O&&NO) || (S&&E&&SE) || (S&&O&&SO) ) ) || (nb_nei_id == 3 && !( (S||(NE&&NO)) && (E||(NO&&SO)) && (N||(SE&&SO)) && (O||(SE&&NE)) ) )) {
			nei_flag = false;
		}
	}
	return nei_flag;
}

bool connected_8(TYPE **state, int nb_nei_id, bool condwrap, int pixel, int ncol, int nrow, int x, int y) {
	//Check connectivity on Moore neighborhood
	bool nei_flag = true;
	if (nb_nei_id > 1) {
		bool N, E, S, O, NE, NO, SO, SE;
		if (condwrap) {
			E=(state[y][PeriodicWrap(x+1, ncol)] == pixel);
			O=(state[y][PeriodicWrap(x-1, ncol)] == pixel);
			S=(state[PeriodicWrap(y+1, nrow)][x] == pixel);
			N=(state[PeriodicWrap(y-1, nrow)][x] == pixel);
			NE=(state[PeriodicWrap(y-1, nrow)][PeriodicWrap(x+1, ncol)] == pixel);
			NO=(state[PeriodicWrap(y-1, nrow)][PeriodicWrap(x-1, ncol)] == pixel);
			SE=(state[PeriodicWrap(y+1, nrow)][PeriodicWrap(x+1, ncol)] == pixel);
			SO=(state[PeriodicWrap(y+1, nrow)][PeriodicWrap(x-1, ncol)] == pixel);
		} else {
			E=(state[y][x+1] == pixel);
			O=(state[y][x-1] == pixel);
			S=(state[y+1][x] == pixel);
			N=(state[y-1][x] == pixel);
			NE=(state[y-1][x+1] == pixel);
			NO=(state[y-1][x-1] == pixel);
			SE=(state[y+1][x+1] == pixel);
			SO=(state[y+1][x-1] == pixel);
		}
		//if ( (nb_nei_id == 1 && !( (N && !(SO)) || (N && !(SE)) || (S && !(NO)) || (S && !(NE)) || (O && !(NE)) || (O && !(SE)) || (E && !(NO)) || (E && !(SO)) )  ) || (nb_nei_id == 2 && !( (N&&E&&NE && !(SO)) || (N&&O&&NO && !(SE)) || (S&&E&&SE && !(NO)) || (S&&O&&SO && !(NE)) ||  (N&&S && !(SO)) || (N&&S && !(SE)) || (N&&S && !(NO)) || (N&&S && !(NE)) || (O&&E && !(NO)) || (O&&E && !(SO)) || (O&&E && !(NE)) || (O&&E && !(SE)) ) ) || (nb_nei_id == 3 && !( (S||(NE&&NO)) && (E||(NO&&SO)) && (N||(SE&&SO)) && (O||(SE&&NE)) ) )) {
		//if ( (nb_nei_id == 1 && ( (N && (SO)) || (N && (SE)) || (S && (NO)) || (S && (NE)) || (O && (NE)) || (O && (SE)) || (E && (NO)) || (E && (SO)) )  ) || (nb_nei_id == 2 && !( (N&&E&&NE && !(SO)) || (N&&O&&NO && !(SE)) || (S&&E&&SE && !(NO)) || (S&&O&&SO && !(NE)) ||  (N&&S && !(SO) && !(SE) && !(NO) && !(NE)) || (O&&E && !(NO) && !(SO) && !(NE) && !(SE)) ) ) || (nb_nei_id == 3 && !( (S||(NE&&NO)) && (E||(NO&&SO)) && (N||(SE&&SO)) && (O||(SE&&NE)) ) )) {
		if ( (nb_nei_id == 1 && !( (N && !(SO) && !(SE) ) || (S && !(NO) && !(NE) ) || (O && !(NE) && !(SE) ) || (E && !(NO) && !(SO) ) )  ) || (nb_nei_id == 2 && !( (N&&E&&NE && !(SO)) || (N&&O&&NO && !(SE)) || (S&&E&&SE && !(NO)) || (S&&O&&SO && !(NE)) ||  (N&&S && !(SO) && !(SE) && !(NO) && !(NE)) || (O&&E && !(NO) && !(SO) && !(NE) && !(SE)) ) ) || (nb_nei_id == 3 && !( (S||(NE&&NO)) && (E||(NO&&SO)) && (N||(SE&&SO)) && (O||(SE&&NE)) ) )) {
			nei_flag = false;

			
		}
	}
	return nei_flag;
}


void ComputePerimeter(int maxcells, Cell cells, int ncol, int nrow, TYPE **state, int mediumcell, int neighbour_energy){
	//Calcule le périmètre de chaque bulle.
	//Pour avoir le périmètre réel en pixel, il faut diviser par le paramètre line_to_area correspondant au 
	//neighbour_energy choisi (11.3 si neighbour_enery=20)
	int k,ii;

	for(k=0;k<maxcells;k++){
		//cells.oldperimeter[k]=cells.perimeter[k];
		cells.perimeter[k]=0;
	}

	PLANE(	//The PLANE() macro defines a loop over the whole field, i.e., i=1, 2, .., nrow and j=1, 2, .., ncol.
	if(cells.celltype[(k=state[i][j])]>mediumcell){
		for(ii=0;ii<neighbour_energy;ii++){
			if(state[PeriodicWrap(i+ny[ii],nrow)][PeriodicWrap(j+nx[ii],ncol)]!=k)
				//cells.perimeter[k]=cells.perimeter[k];
				cells.perimeter[k]++;
		}
	}
	);
}

double ComputeEnergy(int maxcells, Cell cells, int ncol, int nrow, TYPE **state, int neighbour_energy, int** Jarray, double area_constraint){
	//Calcule l'énergie du système.
//	
	double boundary_energy = 0;
	int k, cible, ii;
	double total_energy=0;
	PLANE( //The PLANE() macro defines a loop over the whole field, i.e., i=1, 2, .., nrow and j=1, 2, .., ncol.
		k=state[i][j];
		for(ii=0;ii<neighbour_energy;ii++){
			//Si on trouve un voisins différents avec un type différent, c'est une frontière.
			cible = state[PeriodicWrap(i+ny[ii],nrow)][PeriodicWrap(j+nx[ii],ncol)];
			if(cible!=k){
				boundary_energy += Jarray[cells.celltype[cible]][cells.celltype[k]];
				//fprintf(stderr,"%d\t",Jarray[cells.celltype[cible]][cells.celltype[k]]%1600);
			}
		}
	);
	if((int)(boundary_energy)%2 != 0){
		fprintf(stderr,"Error, Check of Parity of twice Boundary Energy failed = %f\n", boundary_energy);
		exit(EXIT_FAILURE);
	}
	
	boundary_energy=boundary_energy/2; //on divise par deux car chaque couplage est compté deux fois
	total_energy=boundary_energy;
	//fprintf(stderr,"%lf\t%lf\n",boundary_energy,total_energy);
	
	for(k=1;k<maxcells;k++) {//on ajoute l'énergie de compression
		if(cells.targetarea[k]!=0) {
			total_energy+=cells.area_constraint[k]*pow(cells.targetarea[k]-cells.area[k],2)/cells.targetarea[k];
		}
	}
	
	return total_energy;

}

void ComputeCenterCoords(Cell cells, int ncol, int nrow, TYPE **state, int nbrcells, int mediumcell) {
//Calcule les centres géométriques des cellules. Contrairement à la méthode itérative à chaque MCS, cette fonction ne
//peut donner que les coordonnées modulo ncol et nrow, et n'est donc pas adaptée à l'étude de la diffusion.

//Il vaudrait mieux définir Area*xcoord: car ce sont toujours des entiers, contrairement aux coords qui sont des doubles.	
	int airehautgauche[nbrcells], airehautdroite[nbrcells], airebasgauche[nbrcells], airebasdroite[nbrcells];
	int xgauche[nbrcells], xdroite[nbrcells], ybas[nbrcells], yhaut[nbrcells];
	memset(airehautgauche, 0, sizeof(airehautgauche));
	memset(airehautdroite, 0, sizeof(airehautdroite));
	memset(airebasgauche, 0, sizeof(airebasgauche));
	memset(airebasdroite, 0, sizeof(airebasdroite));
	memset(xgauche, 0, sizeof(xgauche));
	memset(xdroite, 0, sizeof(xdroite));
	memset(ybas, 0, sizeof(ybas));
	memset(yhaut, 0, sizeof(yhaut));
	int i, j, k, airegauche, airedroite, airehaut, airebas;
	for(i=1;i<=nrow/2;i++){
		for(j=1;j<=ncol/2;j++){
			k=state[i][j];
			xgauche[k] += j;
			yhaut[k] += i;
			airehautgauche[k]++;			
		}
		for(j=ncol/2+1;j<=ncol;j++){
			k=state[i][j];
			xdroite[k] += j;
			yhaut[k] += i;
			airehautdroite[k]++;			
		}		
	}
	for(i=nrow/2+1;i<=nrow;i++){
		for(j=1;j<=ncol/2;j++){
			k=state[i][j];
			xgauche[k] += j;
			ybas[k] += i;
			airebasgauche[k]++;			
		}
		for(j=ncol/2+1;j<=ncol;j++){
			k=state[i][j];
			xdroite[k] += j;
			ybas[k] += i;
			airebasdroite[k]++;			
		}		
	}
	for(k=mediumcell+1;k<nbrcells;k++){ //A cause des conditions aux limites périodiques, il faut regrouper les domaines cellulaires partagés.
		if (cells.celltype[k] != 0) {
			airehaut=airehautgauche[k]+airehautdroite[k];
			airebas=airebasgauche[k]+airebasdroite[k];
			airegauche=airehautgauche[k]+airebasgauche[k];
			airedroite=airehautdroite[k]+airebasdroite[k];
			if(airehaut + airebas != cells.area[k] || cells.area[k]==0){
				fprintf(stderr,"error: il y a une fuite d'aire de la cellule n°%d ! dont le type est %d et l'aire est %d\n", \
					k, cells.celltype[k], cells.area[k]);
				//exit(EXIT_FAILURE);
			}
			if (2*xdroite[k]*airegauche - 2*xgauche[k]*airedroite > airegauche*airedroite*ncol){
			//if (xdroite[k]/airedroite - xgauche[k]/airegauche > ncol/2){ //critère pour déterminer si la cellule est à cheval sur les bords.
				if(airegauche > airedroite)//on déplace la plus petite des deux moitiés
					xdroite[k]=xdroite[k]-ncol*airedroite;
				else
					xgauche[k]=xgauche[k]+ncol*airegauche;
			}
			if (2*ybas[k]*airehaut - 2*yhaut[k]*airebas > airehaut*airebas*nrow){
			//if (ybas[k]/airebas - yhaut[k]/airehaut > nrow/2){
				if(airehaut > airebas)
					ybas[k]=ybas[k]-nrow*airebas;
				else
					yhaut[k]=yhaut[k]+nrow*airehaut;
			}
			cells.xcoord[k]=(double)(xgauche[k]+xdroite[k])/cells.area[k]; //il faudra reitrer la division par cells.area pour ne travailler qu'avec des entiers
			cells.ycoord[k]=(double)(ybas[k]+yhaut[k])/cells.area[k];
		}
	}
}

int ComputeBoundary(Cell cells, int ncol, int nrow, TYPE **state, int neighbour_energy) {
	//Calcule la longueur de l'interface entre les deux types cellulaires. Le renvoie.
	//Pour avoir la longueur réelle en pixel, il faut diviser par le paramètre line_to_area correspondant au 
	//neighbour_energy choisi (11.3 si neighbour_energy=20).
	int interface = 0, k, cible, ii;
	PLANE(
		k=state[i][j];
		for(ii=0;ii<neighbour_energy;ii++){
			//Si on trouve un voisins différents avec un type différent, c'est une frontière.
			cible = state[PeriodicWrap(i+ny[ii],nrow)][PeriodicWrap(j+nx[ii],ncol)];
			if(cible!=k){
				if(cells.celltype[cible] != cells.celltype[k] ){
					interface ++;
					//break;
				}
			}
		}
	);
	return interface;
}


// DC: calcule l'interface et l'écrit dans un fichier
void ComputeLineInterface(int ttime, FILE* interfacefp, Cell cells, int ncol, int nrow, TYPE **state, int* line_interface) {
	//calcule et ecrit dans un fichier l'interface entre les deux tissus cellulaires
	int k, k_bas;
	line_interface[0]=ttime;
	for(int j=1;j<=ncol;j++){
		for(int i=1;i<=nrow;i++){
			if (i>nrow/10 && i<(9*nrow)/10) {
				k = cells.celltype[state[i][j]];
				k_bas = cells.celltype[state[i+1][j]];
				if (k==1 && k_bas==2 ){
					if (line_interface[j] != 0) printf("Attention interface double\n");
					line_interface[j]=nrow-i;  //lecture y=f(x)
					break;
				}
			}
		}
	}
	//ecriture dans fichier output
	for (int j=0; j<ncol+1; j++){
		fprintf(interfacefp, "%d ", line_interface[j]);
	}
	fprintf(interfacefp, "\n");
	//remise à zéro
	for (int j=0; j<ncol+1; j++){
		line_interface[j]=0;
	}		
}


int surface_tache(int taille, int deja_compte[taille][taille]) {
    int res=0;
    for(int i=0;i<taille;i++) {
        for(int j=0;j<taille;j++)   {
            if (deja_compte[i][j]) res++;
        }
    } 
    return res;    
}


#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

void dessiner_tableau(int taille, int tab[taille][taille], int x0, int y0, int num_cell) {
    for(int i=0;i<taille;i++) {
         for(int j=0;j<taille;j++)   {
			 if (i==x0 && j==y0){
				 printf("%sX ",KRED);
				 //printf("X ");
				 printf("%s", KNRM);
			 }else {
				 if (tab[i][j]==num_cell) {
					printf("%s%d ",KRED,1);	
					printf("%s", KNRM);				
					//printf("%d ",1);

				 }else{

					 printf("%d ",0);
				 }
			 }
         }
         printf("\n");
    }  
	printf("\n\n");    
}

void dessiner_tableau_deux_cellules(int taille, int tab[taille][taille], int x0, int y0, int num_cell, int num_new_cell) {
    for(int i=0;i<taille;i++) {
         for(int j=0;j<taille;j++)   {
			//  if (i==x0 && j==y0){
			// 	 printf("%sX ",KRED);
			// 	 printf("%s", KNRM);
			//  }else {
				 if (tab[i][j]==num_cell) {
					printf("%s%d ",KRED,1);	
					printf("%s", KNRM);				
				 }else if (tab[i][j]==num_new_cell){
					printf("%s%d ",KGRN,2);	
					printf("%s", KNRM);						 				 
				 } else{ 
					 printf("%d ",0);
				 }
			//} //
         }
         printf("\n");
    }  
	printf("\n\n");    
}

// DC: teste la connexité d'une cellule
int tester_cellule_connexe_boucle(Cell cells,TYPE** state, int num_cell, int cote_carre, int nrow, int ncol){
	int x0 = (int)cells.xcoord[num_cell];
	if (x0==0) x0=1;
	int y0 = (int)cells.ycoord[num_cell];
	if (y0==0) y0=1;
	int x=y0,y=x0;  //attention à l'inversion	
	int R = (int) sqrt(cells.area[num_cell]);
	//printf("R=%d\n",R);
	if(R==0) printf("ATTENTION R=0\n");
	int D =cote_carre*R;
	int taille_tab=2*D+1;
	int tab[taille_tab][taille_tab];
	//on recopie state dans tab
	for(int i=0; i<taille_tab; i++){
		int ip=x+i-D; if (ip>=(nrow+1)) ip=ip-nrow; if (ip<=0) ip=nrow+ip; 
		for(int j=0; j<taille_tab; j++){
			int jp=y+j-D; if (jp>=(ncol+1)) jp=jp-ncol; if (jp<=0) jp=ncol+jp;
			tab[i][j]=state[ip][jp];
		}
	}
	//printf("dessin de la cellule en cours de test de connexité, cellule: %d, area: %d\n", num_cell, cells.area[num_cell]);
	//dessiner_tableau(taille_tab, tab, D, D, num_cell );

	point p_num=trouver_point_dedans_boucle(taille_tab, tab, num_cell);
	if (p_num.x<=0 || p_num.y<=0 ) {
		//printf("pas reussi à trouver des points dedans dans connexité cellule, p.x=%d, p.y=%d\n", p_num.x, p_num.y);
		return 0;
	}
	int x_deb_tache_num=p_num.x;
	int y_deb_tache_num=p_num.y;
	int area_cell_num=cells.area[num_cell];	
	int connexe_num=tester_connexe_boucle(taille_tab, tab, num_cell, x_deb_tache_num, y_deb_tache_num, area_cell_num );
	if (!connexe_num) {
		//printf("ALERTE: la cellule %d n'est pas connexe!!\n",num_cell);
		//dessiner_tableau(taille_tab, tab, x_deb_tache_num, y_deb_tache_num, num_cell );
		;
	}
	return connexe_num;
	
}

void empiler_pour_connexite(int taille_pile_connexite, point pile_connexite[taille_pile_connexite], int* dans_pile_connexite, point pixel) {
    (*dans_pile_connexite)++;
    if ((*dans_pile_connexite>taille_pile_connexite-1)) printf("attention, pile connexite saturée!!\n");
    pile_connexite[(*dans_pile_connexite)-1]=pixel;
}

point depiler_pour_connexite(int taille_pile_connexite, point pile_connexite[taille_pile_connexite], int* dans_pile_connexite) {
    point res=pile_connexite[(*dans_pile_connexite)-1];
    (*dans_pile_connexite)--;
    if ((*dans_pile_connexite<0)) printf("attention, dans_pile_connexite négatif!!\n");
    return res;
}


// DC: teste la connexité d'une cellule à oartir d'un pixel. On fait l'extension connexe par une boucle à partie d'un pixel et on vérifie qu'on obtient bien sa surface totale
int tester_connexe_boucle(int taille_tab, int tab[taille_tab][taille_tab], int num_cell, int x_deb_tache, int y_deb_tache, int area_cell) {
	if (x_deb_tache<0 || y_deb_tache<0) {
		printf("x ou y du début de tache négatif!!\n");
	}

    int tab_copie[taille_tab][taille_tab];
	for(int i=0;i<taille_tab;i++) {
		for(int j=0;j<taille_tab;j++) {
			if (tab[i][j]==num_cell) {
				tab_copie[i][j]=1;   
			}else{
				tab_copie[i][j]=0;
			}
		} 
	}

    int inclus[taille_tab][taille_tab];
	for(int i=0;i<taille_tab;i++) {
		for(int j=0;j<taille_tab;j++) {
			inclus[i][j]=0;
		} 
	}

	int taille_pile_connexite=500;
	int dans_pile_connexite=0;
	point pile_connexite[taille_pile_connexite];
	for(int i=0; i<taille_pile_connexite;i++) {pile_connexite[i].x=0; pile_connexite[i].y=0;}
	point p_depart={x_deb_tache,y_deb_tache};
	empiler_pour_connexite(taille_pile_connexite, pile_connexite, &dans_pile_connexite,p_depart);
	inclus[x_deb_tache][y_deb_tache]=1;
	while(dans_pile_connexite>0){
		point p=depiler_pour_connexite(taille_pile_connexite, pile_connexite, &dans_pile_connexite);
		int i=p.x, j=p.y;
		int ip=i+1, jp=j;					
		if(tab_copie[ip][jp] && !inclus[ip][jp]){
			point pp={ip,jp};
			empiler_pour_connexite(taille_pile_connexite, pile_connexite, &dans_pile_connexite,pp);
			inclus[ip][jp]=1;
		}
		ip=i-1; jp=j;					
		if(tab_copie[ip][jp] && !inclus[ip][jp]){
			point pp={ip,jp};
			empiler_pour_connexite(taille_pile_connexite, pile_connexite, &dans_pile_connexite,pp);
			inclus[ip][jp]=1;
		}
		ip=i; jp=j+1;					
		if(tab_copie[ip][jp] && !inclus[ip][jp]){
			point pp={ip,jp};
			empiler_pour_connexite(taille_pile_connexite, pile_connexite, &dans_pile_connexite,pp);
			inclus[ip][jp]=1;
		}
		ip=i; jp=j-1;					
		if(tab_copie[ip][jp] && !inclus[ip][jp]){
			point pp={ip,jp};
			empiler_pour_connexite(taille_pile_connexite, pile_connexite, &dans_pile_connexite,pp);
			inclus[ip][jp]=1;
		}	
	}

	int surface=surface_tache(taille_tab, inclus);
	if (surface==area_cell) {
		//printf("connexite ok pour la cellule %d, area et surface = %d\n",num_cell,area_cell);
		return 1;
	}else{
		//printf("connexite incorrecte, area=%d et surface = %d\n",area_cell, surface);
		//afficher_allume(taille_tab,tab_copie,x_deb_tache,y_deb_tache);
		return 0;
	}
}


int dedans(int taille, int tab[taille][taille], int x0, int y0, int num_cell) {
	//if (tab[x0-1][y0-1]!=num_cell) return 0;
	if (tab[x0-1][y0]!=num_cell) return 0;
	//if (tab[x0-1][y0+1]!=num_cell) return 0;
	if (tab[x0][y0-1]!=num_cell) return 0;
	if (tab[x0][y0+1]!=num_cell) return 0;
	//if (tab[x0+1][y0-1]!=num_cell) return 0;
	if (tab[x0+1][y0]!=num_cell) return 0;
	//if (tab[x0+1][y0+1]!=num_cell) return 0;	
	return 1;
}


// DC: trouver un pixel appartenant à une cellule
point trouver_point_dedans_boucle(int taille, int tab[taille][taille], int num_cell) {
	point p={0,0};
	for (int i=1; i<taille-1; i++) {
		for (int j=1; j<taille-1; j++) {
			if (tab[i][j]==num_cell) {
				p.x=i; p.y=j;
				return p;
			}
		}
	}
	return p;
}

// teste si la division d'une cellule par une droite d'angle theta passant par son centre donne deux cellules filles connexes
int tester_division(int taille_tab, int tab[taille_tab][taille_tab],int D, int num_cell,double theta, int* nb_cellules, int area_cell, int impress) {
	
	int debut_tache_trouve_num=0, x_deb_tache_num=0, y_deb_tache_num=0, area_cell_num=0;
	int debut_tache_trouve_num_new=0, x_deb_tache_num_new=0, y_deb_tache_num_new=0, area_cell_num_new=0;
	area_cell_num=area_cell;
	area_cell_num_new=0;

    int tab_copie[taille_tab][taille_tab];
	for(int i=0;i<taille_tab;i++) {
		for(int j=0;j<taille_tab;j++) {
			tab_copie[i][j]=tab[i][j];
		} 
	}

	int num_new_cell=(*nb_cellules) + 1;
	for(int i=0; i<taille_tab; i++){
		for(int j=0; j<taille_tab; j++){
			if (tab_copie[i][j]==num_cell) {
				if (((i-D)*cos(theta) + (j-D)*sin(theta)) > 0) {
					tab_copie[i][j]=num_new_cell;
					area_cell_num_new++;
					area_cell_num--;
				}
			}
		}
	}
	point p_num=trouver_point_dedans_boucle(taille_tab, tab_copie, num_cell);
	point p_num_new=trouver_point_dedans_boucle(taille_tab, tab_copie, num_new_cell);
	if (p_num.x<=0 || p_num.y<=0 || p_num_new.x<=0 || p_num_new.y<=0) {
		//printf("pas reussi à trouver des points dedans, p.x=%d, p.y=%d, p_new.x=%d, p_new.y=%d\n", p_num.x, p_num.y, p_num_new.x, p_num_new.y);
		return 0;
	}
	x_deb_tache_num=p_num.x;
	y_deb_tache_num=p_num.y;
	x_deb_tache_num_new=p_num_new.x;
	y_deb_tache_num_new=p_num_new.y;	
	int connexe_num=tester_connexe_boucle(taille_tab, tab_copie, num_cell, x_deb_tache_num, y_deb_tache_num, area_cell_num );
	int connexe_num_new=tester_connexe_boucle(taille_tab, tab_copie, num_new_cell, x_deb_tache_num_new, y_deb_tache_num_new, area_cell_num_new );
	// if (connexe_num || connexe_num_new) printf("#############\n");
	// if (connexe_num) printf("CONNEXITE 1: OK\n") ;
	// if (connexe_num_new) printf("CONNEXITE 2: OK\n\n") ;

	if (impress && (!connexe_num || !connexe_num_new)){
		printf("Une des deux cellules_filles non connexe, voici la cellule mère %d avec theta=%.2lf degrés\n", num_cell, (theta/M_PI)*180.0);
		dessiner_tableau_deux_cellules(taille_tab, tab_copie, D, D, num_cell, num_new_cell);
	}
	
	return (connexe_num && connexe_num_new);
}


void Diviser(Cell cells, int num_cell, int ttime, int nrow, int ncol, TYPE** state, int* nb_cellules, int* nb_cellules_vivantes, 
	int* nb_cellules_mortes, int* nb_cellules_par_division, int* nb_cellules1_par_division, int* nb_cellules2_par_division, int* nb_cellules1, int* nb_cellules2, 
	int* nb_cellules_tuees, int* nb_cellules_condamnees, int* nb_cellules1_condamnees, int* nb_cellules2_condamnees, int maxcells, int duree_de_vie1, 
	int duree_de_vie2, int interphase1, int interphase2, int* pile_labels_libres, int* taille_labels_libres, int cote_carre, int* nb_divisions_acceptees, int* nb_divisions_refusees) {

	if (cells.celltype[num_cell]==0) return;
	if (cells.condamnee[num_cell]==1) return;
	if ((*nb_cellules)>=(maxcells-1) && ((*taille_labels_libres)==0)) return;
	int x0 = (int)cells.xcoord[num_cell];
	if (x0==0) x0=1;
	int y0 = (int)cells.ycoord[num_cell];
	if (y0==0) y0=1;
	int x=y0,y=x0;  //attention à l'inversion
	if (state[y0][x0] != num_cell) {
		//printf("attention dans Diviser, centre cellule %d mal positionné, je ne divise pas cette cellule\n", num_cell);
		// printf("attention dans Diviser, centre cellule %d dont type est %d mal positionné, je divise quand même cette cellule\n", 
		// 	num_cell, cells.celltype[num_cell]);
		;		
		//return;
	}
	//printf("je divise la cellule %d dont l'age est %d\n",num_cell, ttime - cells.t_debut_division[num_cell]);

	double theta = aleatoire(0)*2* M_PI;
	int R = (int) sqrt(cells.area[num_cell]);
	//printf("R=%d\n",R);
	if(R==0) printf("ATTENTION R=0\n");
	int D =cote_carre*R;
	int taille_tab=2*D+1;
	int tab[taille_tab][taille_tab];
	//on recopie state dans tab
	for(int i=0; i<taille_tab; i++){
		int ip=x+i-D; if (ip>=(nrow+1)) ip=ip-nrow; if (ip<=0) ip=nrow+ip; 
		for(int j=0; j<taille_tab; j++){
			int jp=y+j-D; if (jp>=(ncol+1)) jp=jp-ncol; if (jp<=0) jp=ncol+jp;
			tab[i][j]=state[ip][jp];
		}
	}
	int division_ok=tester_division(taille_tab, tab, D, num_cell,theta, nb_cellules, cells.area[num_cell],0);
	//DC: on essaie d'autres angles, e, faisant tourner la droite de pi/5
	double thetap=theta;
	if (!division_ok){ //DC: on va faire bouger theta
		//printf("Division refusée pour la cellule %d\n", num_cell);
		for(int i=1;i<=4;i++) {
			thetap=theta + i*M_PI/5;
			//printf("Nouvelle tentative: %d pour la cellule %d\n", i, num_cell);
			division_ok=tester_division(taille_tab, tab, D, num_cell,thetap, nb_cellules, cells.area[num_cell],0);
			if (division_ok) {
				//printf("J'arrive à diviser la cellule %d à la tentative %d\n\n",num_cell,i);
				theta=thetap;
				break;
			}else{
				//printf("Ca n'a pas marché pour la cellule %d à la tentative %d\n", num_cell, i);
				;
			}
		}
	}

	if (!division_ok) {
		//printf("					division refusee pour la cellule %d, area: %d\n",num_cell, cells.area[num_cell]);
		(*nb_divisions_refusees)++;
		cells.division_refusee[num_cell]=1;
		//dessiner_tableau(taille_tab, tab, D, D, num_cell);
		//printf("test de connexite: %d, cellule: %d, area: %d\n", tester_cellule_connexe_boucle(cells,state,num_cell,cote_carre,nrow,ncol), num_cell, cells.area[num_cell]);
		//DC: on rejoue la scène...
		// for(int i=0;i<=4;i++) {
		// 	thetap=theta + i*M_PI/5;
		// 	printf("Nouvelle tentative: %d pour la cellule %d\n", i, num_cell);
		// 	division_ok=tester_division(taille_tab, tab, D, num_cell,thetap, nb_cellules, cells.area[num_cell],1);
		// }
		return;
	}else{
		//printf("division acceptée pour la cellule %d, area: %d\n",num_cell, cells.area[num_cell]);
		(*nb_divisions_acceptees)++;
		cells.division_refusee[num_cell]=0;
		//dessiner_tableau(taille_tab, tab, D, D, num_cell);
		//printf("test de connexite: %d, cellule: %d, area: %d\n", tester_cellule_connexe_boucle(cells,state,num_cell,cote_carre,nrow,ncol), num_cell, cells.area[num_cell]);		
		;
	}

	int num_new_cell=0;
	if ((*taille_labels_libres)==0) {
		num_new_cell=(*nb_cellules) + 1;
		(*nb_cellules)++;
	} else {
		num_new_cell=pile_labels_libres[(*taille_labels_libres)-1];
		(*taille_labels_libres)--;
	}
	
	
	(*nb_cellules_vivantes)++;
	(*nb_cellules_par_division)++;
	if (cells.celltype[num_cell] == 1) {
		(*nb_cellules1)++;
		(*nb_cellules1_par_division)++;
	}
	else if (cells.celltype[num_cell] == 2) {
		(*nb_cellules2)++;
		(*nb_cellules2_par_division)++;
	}

	cells.targetarea[num_new_cell] = cells.targetarea[num_cell];
	cells.targetarea_original[num_new_cell] = cells.targetarea_original[num_cell];
	cells.area[num_new_cell] = 0;
	cells.area_constraint[num_new_cell] = cells.area_constraint[num_cell];
	cells.celltype[num_new_cell] = cells.celltype[num_cell];
	cells.t_debut_division[num_new_cell] = ttime; //cells.t_debut_division[num_cell]; //ttime;
	//cells.t_debut_division[num_cell]=ttime;
	cells.t_derniere_division[num_new_cell] = ttime;
	cells.t_derniere_division[num_cell] = ttime;
	if (cells.celltype[num_new_cell]==1)   cells.interphase[num_new_cell] =  (-100 + (int)(aleatoire(0)*200)) + interphase1;   //+ cells.interphase[num_cell];
	if (cells.celltype[num_new_cell]==2)   cells.interphase[num_new_cell] =  (-100 + (int)(aleatoire(0)*200)) + interphase2;   //+ cells.interphase[num_cell];
	//if (cells.interphase[num_new_cell] < 100) cells.interphase[num_new_cell]=100;
	cells.condamnee[num_new_cell]=0;
	cells.division_refusee[num_new_cell]=0;
	if (cells.celltype[num_new_cell]==1) {
		cells.duree_de_vie[num_new_cell] = duree_de_vie1; //(int) (aleatoire(0)*2*duree_de_vie1) + 1;
	}else if (cells.celltype[num_new_cell]==2) {
		cells.duree_de_vie[num_new_cell] = duree_de_vie2; //(int) (aleatoire(0)*2*duree_de_vie2) + 1;
	}

	for(int i=-D; i<=D; i++){
		int ip=x+i; if (ip>=(nrow+1)) ip=ip-nrow; if (ip<=0) ip=nrow+ip; 
		for(int j=-D; j<=D; j++){
			int jp=y+j; if (jp>=(ncol+1)) jp=jp-ncol; if (jp<=0) jp=ncol+jp;
			if (state[ip][jp]==num_cell) {
				if ((i*cos(theta) + j*sin(theta)) > 0) {
					state[ip][jp]=num_new_cell;
					cells.area[num_cell]--;
					cells.area[num_new_cell]++;
				}
			}
		}
	}

	//cells.vient_de_diviser[num_cell]=1;
	cells.vient_de_diviser[num_new_cell]=1;
	cells.vient_de_diviser_pour_affichage[num_new_cell]=1;
	cells.vient_de_diviser_pour_affichage[num_cell]=1;
	cells.vient_de_naitre[num_new_cell]=1;
	if (cells.area[num_new_cell]==0) {
		Empiler_cellule_vide(cells, num_new_cell, nb_cellules_vivantes, 
			nb_cellules_mortes, nb_cellules1, nb_cellules2, nb_cellules_tuees, nb_cellules_condamnees, nb_cellules1_condamnees,nb_cellules2_condamnees, maxcells, pile_labels_libres, taille_labels_libres, 0);
	}
	if (cells.area[num_cell]==0) {
		Empiler_cellule_vide(cells, num_cell, nb_cellules_vivantes, 
			nb_cellules_mortes, nb_cellules1, nb_cellules2, nb_cellules_tuees, nb_cellules_condamnees, nb_cellules1_condamnees,nb_cellules2_condamnees, maxcells, pile_labels_libres, taille_labels_libres, 0);
	}
	//printf("j'ai fait une division, cellule= %d, new_cellule= %d, nb_cellules=%d\n", num_cell, num_new_cell, *nb_cellules);
}

// DC:calcule la surface totale occupée par un type de cellules
int surface_par_type_cellule(Cell cells, int type_cellule, int nb_cellules) {
	int res=0;
	for (int i=0; i<nb_cellules;i++) {
		if (cells.celltype[i]==type_cellule && cells.condamnee[i]==0) {
			res+=cells.area[i];
		}
	}
	return res;
}


void Condamner_cellule(Cell cells, int num_cell, int *nb_cellules_condamnees, 
					int *nb_cellules1_condamnees, int *nb_cellules2_condamnees, int* nb_cellules_vivantes) {
	//printf("je condamne la cellule: %d, de type: %d, de surface: %d\n", num_cell, cells.celltype[num_cell], cells.area[num_cell]);
	if (cells.celltype[num_cell]==0 || num_cell==0) return;
	if (cells.condamnee[num_cell]==1) return;
	cells.condamnee[num_cell]=1;
	cells.targetarea[num_cell]=0; //utile??
	(*nb_cellules_condamnees)++;

	if (cells.celltype[num_cell]==1){
		(*nb_cellules1_condamnees)++;
	} else if (cells.celltype[num_cell]==2) {
		(*nb_cellules2_condamnees)++;
	}
	if ((*nb_cellules_condamnees)==(*nb_cellules_vivantes)) {
		printf("ALERTE: nb_cellules_condammees=nb_cellules_vivantes soit %d\n", (*nb_cellules_condamnees));
	}
}


void Empiler_cellule_vide(Cell cells, int num_cell, int* nb_cellules_vivantes, 
	int* nb_cellules_mortes, int* nb_cellules1, int* nb_cellules2, 
	int* nb_cellules_tuees, int* nb_cellules_condamnees, int* nb_cellules1_condamnees, int* nb_cellules2_condamnees, int maxcells, int* pile_labels_libres, int* taille_labels_libres, int origine) {
		
		//printf("j'empile la cellule: %d, area: %d, type: %d, origine: %d\n",num_cell, cells.area[num_cell], cells.celltype[num_cell], origine);	
		// printf("nb_cellules_vivantes: %d, nb_cellules1: %d, nb_cellules2 %d\n", (*nb_cellules_vivantes), (*nb_cellules1), (*nb_cellules2) );
		// if ((*nb_cellules_vivantes) != (*nb_cellules1) + (*nb_cellules2)) printf("ALERTE!\n");
		if (cells.celltype[num_cell]==0 || num_cell==0) return;
		(*taille_labels_libres)++;
		if ((*taille_labels_libres) >= 2*maxcells-1) { //sizeof(pile_labels_libres)/sizeof(int)-1) {
			printf("attention, pile saturéé dans Empiler_cellule_vide!!!\n");
			return;
		}
		pile_labels_libres[(*taille_labels_libres)-1] = num_cell;
				
		(*nb_cellules_vivantes)--;
		(*nb_cellules_mortes)++;
		if (cells.celltype[num_cell] == 1) {
			(*nb_cellules1)--;
		}
		else if (cells.celltype[num_cell] == 2) {
			(*nb_cellules2)--;
		}

		if (cells.condamnee[num_cell]) {
			(*nb_cellules_tuees)++;
			(*nb_cellules_condamnees)--;
			if (cells.celltype[num_cell]==1) {
				(*nb_cellules1_condamnees)--;
			} else if (cells.celltype[num_cell]==2) {
				(*nb_cellules2_condamnees)--;
			} 
		}


		cells.celltype[num_cell]=0;
		cells.condamnee[num_cell]=0;
		cells.vient_de_mourir[num_cell]=1;		
	}



void FindNeighbours(int maxcells, Cell cells, int ncol, int nrow, TYPE **state, int mediumcell, int neighbour_connected, int maxneighbours, int* side_interf12, int* side_interf10, int* side_interf20) {
	
	//En toute rigueur neighbour_connected devrait correpondre au voisinage d'adjacence, z=4. Mais du fait de la discrétisation, deux noeuds trivalents sont parfois fusionnés
	// en un seul neud quadrivalent. On a alors un problème: certains cotés on une longueur nulle et ne sont plus comptabilisés. On trouve alors qu'il y a plus de cellules
	// à 5 côtés qu'à 7 côtés ! Pour éviter ce pb, on étend la def de voisins aux 8 pixels (z=8) les plus proches.
	int k,ii,iii,nouv,neigh;
	for(k=0;k<maxcells;k++) {
		cells.nneighbours[k] = 0;
		cells.temp[k] = 0;
	}
	//*side_interf = 0;

	PLANE(	//The PLANE() macro defines a loop over the whole field, i.e., i=1, 2, .., nrow and j=1, 2, .., ncol.
	if(cells.celltype[(k=state[i][j])]>mediumcell){
		for(ii=0;ii<neighbour_connected;ii++){
			if((neigh=state[PeriodicWrap(i+ny[ii],nrow)][PeriodicWrap(j+nx[ii],ncol)])!=k){
				nouv=1;
				for(iii=0;iii<cells.nneighbours[k];iii++){//on vérifie que ce voisin n'est pas déjà compté
					if(neigh==cells.neighbours[k][iii]){
						nouv=0;
						break;
					}
				}
				if(nouv==1){
					if (neigh>mediumcell)
						cells.temp[k] += 1./(3./3.72*(sqrt(cells.area[k])-sqrt(cells.area[neigh]))*(sqrt(cells.area[k])-sqrt(cells.area[neigh]))/(sqrt(cells.area[k])*sqrt(cells.area[neigh])*(sqrt(cells.area[k])+sqrt(cells.area[neigh]))));
					
					cells.nneighbours[k]++;
					if(cells.nneighbours[k]>=maxneighbours){
						printf("trop de voisins! =%d pour la cellule %d\n",cells.nneighbours[k], k);
						fprintf(stderr,"Error, FindNeighbours: increase MAXNEIGHBOURS\n");
						exit(EXIT_FAILURE);
					}
					cells.neighbours[k][iii]=neigh;
					//Calcul de la longueur d'interface entre différents types cellulaires mesurée en nombres de côtés. Il faut diviser par 2 le résultat (qui doit donc être pair) pour avoir le nbr de  bord de Plateua.
					/*//On ne tient pas compte de l'interface avec le mediumcell.
					if((cells.celltype[k] != cells.celltype[neigh]) && (cells.celltype[neigh]>mediumcell)) {
						*side_interf12=*side_interf12 + 1;
					}*/
					if(cells.celltype[k] ==1) {
						if(cells.celltype[neigh]==2) {
							*side_interf12=*side_interf12 + 1;
						}
						else if (cells.celltype[neigh]==mediumcell) {
							*side_interf10=*side_interf10 + 1;
						}
					}
					else if ((cells.celltype[k] ==2) && (cells.celltype[neigh]==mediumcell)) {
						*side_interf20=*side_interf20 + 1;
					}
										
				}
			}
		}
	}
	);
}


void AffichageCouleurs(int affichage, Cell cells, int ncol, int nrow, char *subdirname, char *subdirnameRAW ,TYPE **state, TYPE **nstate, int montrer_division) {
	//couleurs 
	//affichage=0 for white monochrome, 1 for color with number of neighbours, 2 for color with celltypes, 3 for five different colors (large with 6, 7 neighbours, small with 5, 6 neighbours + default).
		 // Default : try (but easiliy fail!) to make a different color for each bubble.
	int n = 0;
	if (affichage == 1) {
		//Use of Cash plane, need i and j are no define variable.
		PLANE(
			n = cells.nneighbours[state[i][j]];
			if(n<3 || n>9)
				nstate[i][j]=0;
			else
				nstate[i][j]=n-2;
		);
	}
	else if (affichage == 2) {
		PLANE(
			if(state[i][j]!=0){
				if(cells.celltype[state[i][j]]==1)
					nstate[i][j]=4;  //4  10=rouge
				else
					nstate[i][j]=8;  //8 
			}
			else 
				nstate[i][j]=1;
			
			if (montrer_division){
				if (cells.vient_de_diviser_pour_affichage[state[i][j]]==1) {
					//printf("cette cellule s'est divisée\n");
					nstate[i][j]=10;
					//cells.vient_de_diviser[state[i][j]] = 0;
				} 
			}
		);
		PLANE(
			if (nstate[i][j]==10) cells.vient_de_diviser_pour_affichage[state[i][j]] = 0;
		);
	}
	else if (affichage == 3) {
		PLANE(
			n = cells.nneighbours[state[i][j]];
			if(state[i][j]!=0){
				if (cells.celltype[state[i][j]]==1 && n == 6)
					nstate[i][j]=4;
				else if (cells.targetarea[state[i][j]]==2 && n == 6)
					nstate[i][j]=7;
				else if (cells.targetarea[state[i][j]]==1 && (n == 5 || n == 7))
					nstate[i][j]=9;
				else if (cells.targetarea[state[i][j]]==2 && (n == 5 || n == 7))
					nstate[i][j]=5;
				else
					nstate[i][j]=2;
			}
			else
				nstate[i][j]=1;
		);
	}
	else if (affichage == 0) { //monochrome
				PLANE(
					nstate[i][j]=WHITE;
				);
			}
	else {
		PLANE(
			nstate[i][j]=state[i][j]+1;
		);
	}
	//plan affichage, plan cellules, yoffset, xoffset, coloroffset
	CellPlaneDisplay(nstate,state,0,0,0);
	//CellPlaneRAW2(subdirnameRAW,state,state,0);

	PlaneRAW2(subdirnameRAW,state,state,0);
	CellPlanePNG2(subdirname,nstate,state,0); //nstate définit dans quelle couleur une bulle apparaitra, tandis que state définit la position des bulles
}