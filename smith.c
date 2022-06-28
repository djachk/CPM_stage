#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct point point;

struct point {
    int x;
    int y;
};



int compter(int taille, int tab[taille][taille]) {
    int res=0;
    for(int i=0;i<taille;i++){
        for(int j=0;j<taille;j++){
            if (tab[i][j]==1) res++;
        }
    }
    return res;
}

int compter_allumes(int taille, int allume[taille][taille]){
    int res=0;
    for(int i=0;i<taille;i++) {
         for(int j=0;j<taille;j++)   {
             if (allume[i][j]==1) res++;
         }
    }  
    return res;  
}

void allumer_contours(int taille, int tab[taille][taille], int allume[taille][taille], int num_cell) {
    //remettre à zéro
    for(int i=0;i<taille;i++) {
         for(int j=0;j<taille;j++)   {
             allume[i][j]=0;
         }
    }
    for(int i=1;i<taille-1;i++) {
         for(int j=1;j<taille-1;j++)
            if (tab[i][j]==num_cell) {
                if (tab[i+1][j]!=num_cell || tab[i-1][j]!=num_cell || tab[i][j+1]!=num_cell || tab[i][j-1]!=num_cell) {
                    allume[i][j]=1;
                }
            } 
    }
    for(int i=0;i<taille;i++) {
        if (tab[i][0]==num_cell) allume[i][0]=1;
    }
    for(int i=0;i<taille;i++) {
        if (tab[i][taille-1]==num_cell) allume[i][taille-1]=1;
    }
    for(int j=0;j<taille;j++) {
        if (tab[0][j]==num_cell) allume[0][j]=1;
    }    
    for(int j=0;j<taille;j++) {
        if (tab[taille-1][j]==num_cell) allume[taille-1][j]=1;
    }     
}

void empiler(int taille_pile, point pile[taille_pile], int* haut_pile, point val) {
    (*haut_pile)++;
    if ((*haut_pile>taille_pile-1)) printf("attention, pile saturée!!\n");
    pile[(*haut_pile)-1]=val;
}

point depiler(int taille_pile, point pile[taille_pile], int* haut_pile) {
    point res=pile[(*haut_pile)-1];
    (*haut_pile)--;
    if ((*haut_pile<0)) printf("attention, haut_pile négatif!!\n");
    return res;
}

int interieur(int x, int y, int taille, int tab[taille][taille], int allume[taille][taille], int deja_compte[taille][taille]){
    if (x<0 || x>=taille || y<0 || y>=taille) return 0;
    if (allume[x][y] && !deja_compte[x][y]) deja_compte[x][y]=1;
    int res=!allume[x][y];
    return res;
}

int interieur_gauche(int x, int y, int taille, int tab[taille][taille], int allume[taille][taille], int deja_compte[taille][taille]){
    int res= interieur(x,y,taille,tab, allume, deja_compte) && !interieur(x-1,y,taille,tab,allume, deja_compte);
    return res;
}

int interieur_droit(int x, int y, int taille, int tab[taille][taille], int allume[taille][taille], int deja_compte[taille][taille]){
    int res= interieur(x,y,taille,tab, allume, deja_compte) && !interieur(x+1,y,taille,tab,allume, deja_compte);
    return res;
}

int Ajouts(int xg, int xd, int y, int taille, int tab[taille][taille], int allume[taille][taille], int taille_pile, point pile[taille_pile], int* haut_pile, int deja_compte[taille][taille] ){

    int xdd=xd;
    while(interieur(xdd,y,taille,tab,allume, deja_compte)) { 
        xdd++;
    }
    while(xg<=xdd) {
        while(xg<=xdd && !interieur(xdd,y,taille,tab,allume, deja_compte)) {
            xdd--;
        }
        if(xg<=xdd) {
            point p={xdd,y};
            empiler(taille_pile, pile, haut_pile, p);
            while(xg<=xdd && interieur(xdd,y,taille,tab,allume, deja_compte)) {
                xdd--;
            }
        }
    }
}

void afficher_allume(int taille, int allume[taille][taille]) {
    for(int i=0;i<taille;i++) {
         for(int j=0;j<taille;j++)   {
             printf("%d ",allume[i][j]);
         }
         printf("\n");
    }       
}

void completer_deja_compte(int taille, int tab[taille][taille], int deja_compte[taille][taille]) {
    int deja_compte_copie[taille][taille];
    for(int i=0;i<taille;i++) {
        for(int j=0;j<taille;j++)   {
            deja_compte_copie[i][j]=deja_compte[i][j] ;
        }
    }     
    for(int i=0;i<taille;i++) {
        for(int j=0;j<taille;j++)   {
            if (deja_compte[i][j]) {
                if (i>0 && i<taille-1 && j>0 && j<taille-1) {
                    if (tab[i+1][j]) deja_compte_copie[i+1][j]=1;
                    if (tab[i-1][j]) deja_compte_copie[i-1][j]=1; 
                    if (tab[i][j+1]) deja_compte_copie[i][j+1]=1;
                    if (tab[i][j-1]) deja_compte_copie[i][j-1]=1;
                }
                if (i==0 && j>0 && j<taille-1) {
                    if (tab[i+1][j]) deja_compte_copie[i+1][j]=1;
                    if (tab[i][j+1]) deja_compte_copie[i][j+1]=1;
                    if (tab[i][j-1]) deja_compte_copie[i][j-1]=1;
                }
                if (i==taille-1 && j>0 && j<taille-1) {
                    if (tab[i-1][j]) deja_compte_copie[i-1][j]=1; 
                    if (tab[i][j+1]) deja_compte_copie[i][j+1]=1;
                    if (tab[i][j-1]) deja_compte_copie[i][j-1]=1;
                } 
                if (i>0 && i<taille-1 && j==0) {
                    if (tab[i+1][j]) deja_compte_copie[i+1][j]=1;
                    if (tab[i-1][j]) deja_compte_copie[i-1][j]=1; 
                    if (tab[i][j+1]) deja_compte_copie[i][j+1]=1;
                }   
                if (i>0 && i<taille-1 && j==taille-1) {
                    if (tab[i+1][j]) deja_compte_copie[i+1][j]=1;
                    if (tab[i-1][j]) deja_compte_copie[i-1][j]=1; 
                    if (tab[i][j-1]) deja_compte_copie[i][j-1]=1;
                }                                         
            }
        }
    } 
    for(int i=0;i<taille;i++) {
        for(int j=0;j<taille;j++)   {
            deja_compte[i][j]=deja_compte_copie[i][j] ;
        }
    }      
    return ;    
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

int main(){
    int taille=20;
    int tab[20][20]={0};
    int allume[20][20]={0};
    int taille_pile=1000;
    point pile[1000]={0,0};
    int haut_pile=0;
    int deja_compte[20][20]={0};
    printf("%d\n", tab[5][8]);
    
    for (int i=0; i<20;i++) {
        // if (i==3) {
        //     for (int j=4; j<11; j++) tab[i][j]=1;
        // }
        if (i>=3 && i<=5) {
            for (int j=3; j<12; j++) tab[i][j]=1;
        }
        if (i>5 && i<15){
            for (int j=3; j<6; j++) tab[i][j]=1;
            for (int j=8; j<12; j++) tab[i][j]=1;
        }
        if (i>14 && i<19) {
            for (int j=15; j<17; j++) tab[i][j]=1;
        }
    }
    int res=compter(20,tab);
    printf("compte=%d\n",res);

    point point1={10,11};
    point point2={13,14};
    empiler(1000, pile, &haut_pile, point1 );
    empiler(1000, pile, &haut_pile, point2);
    printf("depil: %d\n", depiler(1000,pile,&haut_pile).x);
    printf("depil: %d\n", depiler(1000,pile,&haut_pile).x);

    allumer_contours(20,tab,allume,1);
    printf("contours=%d\n", compter_allumes(20,allume));
    allumer_contours(20,tab,allume,0);
    printf("contours=%d\n", compter_allumes(20,allume));

    //smith
    printf("debut-smith\n");
    printf("config initiale\n");
    afficher_allume(taille, tab);
    allumer_contours(20,tab,allume,1);
    printf("haut_pile=%d\n",haut_pile);
    printf("allume avant\n");
    afficher_allume(taille, allume);
    int x0=15, y0=16;
    //int x0=4, y0=6;
    if (!interieur(x0,y0,taille,tab,allume, deja_compte)) {
        printf("Point départ incorrect alors que tab[x0][y0]=%d et allume[x0][y0]=%d\n", tab[x0][y0], allume[x0][y0] );
        return 0;
    }
    int xd=x0, y, xg;
    while(!interieur_droit(xd,y0,taille,tab,allume, deja_compte)) {
        xd++;
    }
    point p={xd,y0};
    empiler(taille_pile, pile, &haut_pile, p);
    while(haut_pile>0) {
        point pc=depiler(taille_pile, pile, &haut_pile);
        xd=pc.x; y=pc.y; xg=xd;
        if (allume[xd][y]) continue;
        while (!interieur_gauche(xg,y,taille,tab,allume, deja_compte)) {
            xg--;
        }
        for(int i=xg;i<=xd;i++) {
            allume[i][y]=1;
            if (!deja_compte[i][y]) deja_compte[i][y]=1;
        }
        if (y>0) Ajouts(xg,xd,y-1,taille,tab,allume,taille_pile,pile,&haut_pile, deja_compte);
        if (y<taille-1) Ajouts(xg,xd,y+1,taille,tab,allume,taille_pile,pile,&haut_pile, deja_compte );
        
    }

    printf("allume apres\n");
    afficher_allume(taille, allume);

    completer_deja_compte(taille,tab,deja_compte);    

    printf("deja_compte\n");   
    afficher_allume(taille, deja_compte);
    printf("surface de la tache: %d\n", surface_tache(taille, deja_compte));     
    return 0;
}