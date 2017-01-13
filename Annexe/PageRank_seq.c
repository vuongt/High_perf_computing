
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// calculer le produit d'une matrice carrée n*n avec un vecteur de taille n
void produitMV(double **m, double *x, double *res, int n){
    int i,j;
    for (i=0;i<n;i++){
        res[i]=0;
        for (j=0;j<n;j++){
            res[i] += m[i][j]*x[j];
        }
    }
}

/*Calculer la puissant d'une matrice carrée A*/
void puissantMatrice(double **A, double ** res, int taille, int k){
    int i;
    for (i=0;i<k;i++){
        produitMM(A,A,res,taille);
    }
}

/*Calculer la norme d'un vectuer de taille n*/
double normeVecteur(double * x, int n){
    int j;
    double norme = 0;
    for(j=0;j<n;j++){
        norme += x[j]*x[j];
    }
    norme = (double) sqrt(norme);
    return norme;
}

/*Allocation dynamique d'une matrice*/
double ** allocation_m(int nbLignes,int nbColonnes){
    double ** ptr;

    ptr = malloc(nbLignes * sizeof(double));       //On alloue 'taille1' pointeurs.
    if(ptr == NULL){
        exit(EXIT_FAILURE); //Notifier l'erreur et de quitter le programme.
    }

    int i;

    for(i=0;i<nbLignes;i++){
        ptr[i] = malloc (nbColonnes * sizeof(double));
        if(ptr[i] == NULL){                              //En cas d'erreur d'allocation
         //Il faut libérer la mémoire déjà allouée
         for(i = i-1 ; i >= 0 ; i--)    //On parcourt la boucle dans l'ordre inverse pour libérer ce qui a déjà été alloué
              free(ptr[i]);
         free(ptr);  //Libérer le premier pointeur
         exit(EXIT_FAILURE);//notifier l'erreur et de quitter le programme.
     }
    }
    return ptr;

}

/*Libérer la mémoire*/
void * my_free_m(double **ptr, int nbLignes){
    int i;
    for(i=0 ; i < nbLignes ; i++){
        free(ptr[i]);
    }

    free(ptr);
    ptr = NULL;
    return NULL;
}

int main(int argc, char *argv[]){

    int const n = 4; // la dimension de la matrice
    double tolerance = 0.0001;

    //Création de matrice A
    double **A = allocation_m(n,n);

    //Allouer la mémoire pour le vecteur b
    double *b = malloc(n*sizeof(double));
    if (b==NULL){
        exit(EXIT_FAILURE);
    }
    double *tmp = malloc(n*sizeof(double));
    if (tmp == NULL){
        exit(EXIT_FAILURE);
    }
    int i,j;
    for (i=0;i<n;i++){
        b[i]= (double)1/n;
        tmp[i] = 0;
        for (j=0;j<n;j++){
            A[i][j] = 0;
        }
    }

    A[0][0] = 0;A[0][1] = (double)1/2 ;A[0][2] = (double)1/4 ;A[0][3] =0 ;
    A[1][0] = (double)1/3 ;A[1][1] = 0;A[1][2] = (double)1/4 ;A[1][3] = 0;
    A[2][0] = (double)1/3;A[2][1] =0 ;A[2][2] = (double)1/4;A[2][3] = 1;
    A[3][0] = (double)1/3;A[3][1] = (double)1/2;A[3][2] = (double)1/4;A[3][3] = 0;
    //**A = {{0,1/2,1/4,0},{1/3,0,1/4,0},{1/3,0,1/4,1},{1/3,1/2,1/4,0}};
    printf("Matrice A \n");
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
        printf("%lf ",A[i][j]);
        }
        printf("\n");
    }
    int satisfait = 0; // =1 si norme(bk+1) - norme(bk) soit inférieur la la valeur de tolérance définie ci-dessus
    while (satisfait==0){
    // Calculer la tmp = A*b
        produitMV(A,b,tmp,n);
        // Calculer la norme de A*b
        /*printf("Vecteur tmp: \n");
        for (i=0;i<n;i++){
            printf("%lf \n", tmp[i]);
        }*/
        double normeTmp = normeVecteur(tmp,n);
        double normeB = normeVecteur(b,n);
        printf("normeTmp = %lf \n",normeTmp);
        printf("normeB = %lf \n",normeB);
        if (fabs(normeTmp-normeB)<tolerance){
            satisfait = 1;
        }
        printf("Vecteur b: \n");
        for (i=0;i<n;i++){
            b[i] =(double) tmp[i]/normeTmp;
            printf("%lf \n", b[i]);
        }

    }

    my_free_m(A,n);
    free(b);
    free(tmp);
    return 0;
}
