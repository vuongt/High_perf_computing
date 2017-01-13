#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void produitPartielMM(double **m1, double **m2, double **res, int n);
void produitMV(double **m, double *x, double *res, int n);
void allocation(double *** tableau, int nbLignes,int nbColonnes);
void my_free(double ***tableau);
double normeVecteur(double * x, int n);
void Produit_2Matrices_Parallele(double **m1, double **m2, double **res, int n,int p, int rank);
void QRfactorization(double **A, double ** Q, double ** R, int n,int p,int rank);
void clone(double **A, double **B, int n);


int main(int argc, char *argv[]){
    int const n = 12; /*La taille de la matrice carrée*/
    int nbIteration = 100;
    int rank;       /* Mon rang*/
    int p;          /*Le nombre de processus*/
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    double ** A=NULL;
    double ** Q=NULL;
    double ** R=NULL;
    if (rank ==0){

        allocation(&A,n,n);
        allocation(&Q,n,n);
        allocation(&R,n,n);
        //Initialisation
        A[0][0] = 50.000000; A[0][1] = 0.000000   ;A[0][2] = 2.000000    ; A[0][3] = 0.000000   ;A[0][4]= 0.000000  ;A[0][5] =10.000000  ; A[0][6] = 0.000000  ; A[0][7]=0.000000   ; A[0][8] =  26.000000 ; A[0][9]= 0.000000 ;A[0][10]=  0.000000  ;A[0][11] = 1.000000;
        A[1][0] = 0.000000 ; A[1][1] =  17.000000 ;A[1][2] =   0.000000  ; A[1][3] =  0.000000  ;A[1][4]=  0.000000 ;A[1][5] = 0.000000  ; A[1][6] = 4.000000  ; A[1][7]=0.000000   ; A[1][8] =  29.000000 ; A[1][9]= 0.000000 ;A[1][10]=  41.000000 ;A[1][11] =  0.000000;
        A[2][0] = 2.000000 ; A[2][1] =  0.000000  ;A[2][2] =  199.000000 ; A[2][3] =  0.000000  ;A[2][4]=  0.000000 ;A[2][5] = 14.000000 ; A[2][6] =  16.000000; A[2][7]=  70.000000; A[2][8] =  14.000000 ; A[2][9]= 8.000000 ;A[2][10]=  28.000000 ;A[2][11] =  19.000000;
        A[3][0] = 0.000000 ; A[3][1] =  0.000000  ;A[3][2] =  0.000000   ; A[3][3] =  26.000000 ;A[3][4]=   0.000000;A[3][5] =  0.000000 ; A[3][6] =  0.000000 ; A[3][7]= 40.000000 ; A[3][8] =  0.000000  ; A[3][9]=3.000000  ;A[3][10]= 5.000000   ;A[3][11] =1.000000;
        A[4][0] = 0.000000 ; A[4][1] =  0.000000  ;A[4][2] =  0.000000   ; A[4][3] =  0.000000  ;A[4][4]=  1.000000 ;A[4][5] = 0.000000  ; A[4][6] = 0.000000  ; A[4][7]=0.000000   ; A[4][8] =  7.000000  ; A[4][9]=0.000000  ;A[4][10]= 0.000000   ;A[4][11] =0.000000;
        A[5][0] = 10.000000; A[5][1] =   0.000000 ;A[5][2] =   14.000000 ; A[5][3] =  0.000000  ;A[5][4]=  0.000000 ;A[5][5] = 74.000000 ; A[5][6] =  0.000000 ; A[5][7]= 0.000000  ; A[5][8] =  50.000000 ; A[5][9]= 0.000000 ;A[5][10]=  0.000000  ;A[5][11] = 3.000000;
        A[6][0] = 0.000000 ; A[6][1] =  4.000000  ;A[6][2] =  16.000000  ; A[6][3] =  0.000000  ;A[6][4]=  0.000000 ;A[6][5] = 0.000000  ; A[6][6] = 5.000000  ; A[6][7]=0.000000   ; A[6][8] =  6.000000  ; A[6][9]=2.000000  ;A[6][10]= 9.000000   ;A[6][11] =2.000000;
        A[7][0] = 0.000000 ; A[7][1] =  0.000000  ;A[7][2] =  70.000000  ; A[7][3] =  40.000000 ;A[7][4]=   0.000000;A[7][5] =  0.000000 ; A[7][6] =  0.000000 ; A[7][7]= 114.000000; A[7][8] =  0.000000  ; A[7][9]=0.000000  ;A[7][10]= 12.000000  ;A[7][11] = 7.000000;
        A[8][0] = 26.000000; A[8][1] =   29.000000;A[8][2] =    14.000000; A[8][3] =  0.000000  ;A[8][4]=  7.000000 ;A[8][5] = 50.000000 ; A[8][6] =  6.000000 ; A[8][7]= 0.000000  ; A[8][8] =  161.000000; A[8][9]=  0.000000;A[8][10]=   79.000000;A[8][11] =   13.000000;
        A[9][0] = 0.000000 ; A[9][1] =  0.000000  ;A[9][2] =  8.000000   ; A[9][3] =  3.000000  ;A[9][4]=  0.000000 ;A[9][5] = 0.000000  ; A[9][6] = 2.000000  ; A[9][7]=0.000000   ; A[9][8] =  0.000000  ; A[9][9]=10.000000 ;A[9][10]=  0.000000  ;A[9][11] = 4.000000;
        A[10][0] = 0.000000 ; A[10][1] =  41.000000 ;A[10][2] =   28.000000 ; A[10][3] =  5.000000  ;A[10][4]=  0.000000 ;A[10][5] = 0.000000  ; A[10][6] = 9.000000  ; A[10][7]=12.000000  ; A[10][8] =  79.000000 ; A[10][9]= 0.000000 ;A[10][10]=  123.000000;A[10][11] =   0.000000;
        A[11][0] = 1.000000 ; A[11][1] =  0.000000  ;A[11][2] =  19.000000  ; A[11][3] =  1.000000  ;A[11][4]=  0.000000 ;A[11][5] = 3.000000  ; A[11][6] = 2.000000  ; A[11][7]=7.000000   ; A[11][8] =  13.000000 ; A[11][9]= 4.000000 ;A[11][10]=  0.000000  ;A[11][11] = 68.00000;
        int i,j;
        //Afficher la matrice A
        printf("Matrice A \n");
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }

    }

    int m;
    for (m=1;m<nbIteration;m++){
        QRfactorization(A,Q,R,n,p,rank);
        Produit_2Matrices_Parallele(R,Q,A,n,p,rank);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    //Afficher les résultats
    if (rank ==0){
        printf("Matrice A diagonalisée \n");
        int i,j;
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
        printf("Les valeurs propres approchées de la matrice A sont: \n");
        for (i=0;i<n;i++){
            printf("%lf \n ", A[i][i]);
        }
        my_free(&R);
        my_free(&A);
        my_free(&Q);
    }
    MPI_Finalize();
    return 0;
}


/**
Cette fonction calcule le produit de deux matrices m1 et m2 en parallèle puis stocker le résultat dans la matrice res
n est la taille des matrices
p est le nombre de processus
rank est le rang du processus qui appelle la fonction
*/
void Produit_2Matrices_Parallele(double **m1, double **m2, double **res, int n,int p, int rank){
    int n_local;    /*La taille de la matrice locale pour chaque processeur*/
    int tag1= rand();
    int tag2 = rand();
    int tag3 = rand();

    MPI_Status status;  /*valeur de retour pour le récepteur*/
    /*On doit donc diviser la matrice globale m1,m2 en p matrices carrées extraites.
    Supposons que p = q*q
    On a donc n_local = n/q
    Supposons dans un premier temps que sqrt(p) est un nombre entier et que n est divisable par q*/
    int q = (int)sqrt(p); // le nombre de block par dimension
    if (q*q !=p) {
        printf("La taille de la matrice de connectivité est %d, entrer un nombre de processus de la forme q*q tel que q divise n \n",n);
        exit(EXIT_FAILURE);
    }
    if (n%q!=0) {
        printf("La taille de la matrice de connectivité est %d, entrer un nombre de processus de la forme q*q tel que q divise n \n",n);
        exit(EXIT_FAILURE);
    }

    n_local = (int)n/q;    // La taille de chaque block
	int I = 1 + rank%q;
	int J = 1 + rank/q;
	//printf("Processus de rang %d, le couple d'indice associé est: I = %d,J = %d \n",rank, I, J);

	/*Créer les commutateur ligne et colonne*/
    MPI_Comm commutateur_ligne;
    MPI_Comm commutateur_colonne;
	MPI_Comm_split(MPI_COMM_WORLD,I,J,&commutateur_ligne);
    MPI_Comm_split(MPI_COMM_WORLD,J,I,&commutateur_colonne);

    /*Allouer la mémoire et initialiser les matrices locales pour chaque processeur*/
    double ** m1_local;
    allocation(&m1_local,n_local,n_local);
    double ** m2_local;
    allocation(&m2_local,n_local,n_local);
    double ** res_local;
    allocation(&res_local,n_local,n_local);
    int i,j;
    for (i=0;i<n_local;i++){
            for (j=0;j<n_local;j++){
                m1_local[i][j] = 0;
                m2_local[i][j] = 0;
                res_local[i][j] = 0;
            }
        }

    if (rank == 0){
    /**La mission du processus 0 est de diviser les matrice m1 et m2 en blocs puis
    * envoyer les blocs aux autres processus en fonction de leur rang.
    * A la fin, il récupère les matrice res_local
    * et reconstruit la matrice du résultat global */

        /*Diviser la matrice gobale et envoyer vers les processus*/
		int k; // le rang du destinataire
		/*copy les blocs dans les matrices m1_local et m2_local
		* m1_local stock le bloc (I,J) de m1
		* m2_local stock le bloc (I,J) de m2
		* I = 1 + k%q;
        * J = 1 + k/q;
		*/
		for (k = 1; k<p; k++){
            for (i=0;i<n_local;i++){
                for (j=0;j<n_local;j++){
                    m1_local[i][j] = m1[(k%q)*n/q +i][(k/q)*n/q + j];
                    m2_local[i][j] = m2[(k%q)*n/q +i][(k/q)*n/q + j];
                }
            }
            //Envoyer vers le processus k
            MPI_Send(&(m1_local[0][0]),n_local*n_local, MPI_DOUBLE,k,tag1,MPI_COMM_WORLD);
            MPI_Send(&(m2_local[0][0]),n_local*n_local, MPI_DOUBLE,k,tag2,MPI_COMM_WORLD);

            
		}
		/*En fin, lors que tous les blocs sont envoyés,
		* le processus 0 copy les blocs (0,0) de m1 et m2
		* pour effectuer un calcul local lui-même
		*/
        for (i=0;i<n_local;i++){
            for (j=0;j<n_local;j++){
                m1_local[i][j] = m1[i][j];
                m2_local[i][j] = m2[i][j];
            }
        }
		
    }
    else {
        
        //recevoir les matrices extraites stocké dans m1_local et m2_local
        MPI_Recv(&(m1_local[0][0]), n_local*n_local, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &status);
        MPI_Recv(&(m2_local[0][0]), n_local*n_local, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD, &status);
       
    }
    //uncomment les lignes suivantes pour vérifier les ranks ligne et colonne
    /*
    int rank_ligne;
    MPI_Comm_rank(commutateur_ligne, &rank_ligne);
    int rank_colonne;
    MPI_Comm_rank(commutateur_colonne, &rank_colonne);
    printf("This is processus %d, my rank in commutateur ligne is: %d, my rank in commutateur colonne is: %d \n", rank, rank_ligne, rank_colonne);
	*/
    double ** m1_temp;
    allocation(&m1_temp,n_local,n_local); // matrice temporaire pour stocker les messages provenant des autres processus
    double ** m2_temp;
    allocation(&m2_temp,n_local,n_local); //matrice temporaire pour stocker les messages provenant des autres processus
        for (i=0;i<n_local;i++){
            for (j=0;j<n_local;j++){
                m1_temp[i][j] = 0;
                m2_temp[i][j] = 0;
            }
        }

    /**Processus de l'envoie et de la reception des blocs dans la même ligne ou même colonne*/
    int K; // Numéro du bloc à envoyer/recevoir bloc != rang !
    /*le rang est est entre 0 et q-1
        *rank = numéro du bloc -1*/
    // a chaque itération K, c'est le processus K-1 qui envoie/recevoir les données
    for (K=1;K<=q;K++){
        //Si la matrice à computer est la matrice locale, m1_temp = m1_local
        if (K == J){
            int i,j;
            for (i=0;i<n_local;i++){
                for (j=0;j<n_local;j++){
                    m1_temp[i][j] = m1_local[i][j];
                }
            }
        }
        // si rank = K-1 alors m1_temp = m1_local est envoyé aux autres dans la même ligne
        // si rank != K-1 alors la matrice reçue est stocké dans m1_temp
        MPI_Bcast(&(m1_temp[0][0]), n_local*n_local, MPI_DOUBLE, K-1, commutateur_ligne);
        if (K == I){
            int i,j;
            for (i=0;i<n_local;i++){
                for (j=0;j<n_local;j++){
                    m2_temp[i][j] = m2_local[i][j];
                }
            }
        }
        // si rank = K-1 alors m2_temp = m2_local est envoyé aux autres dans la même colonne
        // si rank != K-1 alors la matrice reçue est stocké dans m2_temp
        MPI_Bcast(&(m2_temp[0][0]), n_local*n_local, MPI_DOUBLE, K-1, commutateur_colonne);


        produitPartielMM(m1_temp,m2_temp,res_local,n_local);
    }

    MPI_Barrier(MPI_COMM_WORLD);


    if (rank!=0){
        //envoyer les résultats locaux vers le processus 0
        MPI_Send(&(res_local[0][0]),n_local*n_local, MPI_DOUBLE,0,tag3,MPI_COMM_WORLD);
    }
    else{
        //double **res;
        //allocation(&res,n,n);
        int i,j;
        //Cloner le résultat calculer par le processus 0 à la matrice durésultat global
        for (i=0;i<n_local;i++){
            for (j=0;j<n_local;j++){
                res[i][j] = res_local[i][j];
            }
        }

        int k; // rang du processus qui envoie le résultat
        for (k=1;k<p;k++){
            MPI_Recv(&(res_local[0][0]), n_local*n_local, MPI_DOUBLE, k, tag3, MPI_COMM_WORLD, &status);

            //Cloner le résultat local au résultat global
            for (i=0;i<n_local;i++){
                for (j=0;j<n_local;j++){
                    res[(k%q)*n/q +i][(k/q)*n/q + j] = res_local[i][j];
                }
            }
        }
        //Afficher le résultat final
        /*
        printf("Matrice res = m1*m2 :\n");
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
            printf("%lf ",res[i][j]);
            }
            printf("\n");
        }
        */
        //my_free(&res);
    }
    my_free(&m1_temp);
    my_free(&m2_temp);

    my_free(&m1_local);
    my_free(&m2_local);
    my_free(&res_local);
    MPI_Comm_free(&commutateur_colonne);
    MPI_Comm_free(&commutateur_ligne);

}

/** Cette fonction calcule le produit de 2 matrices carrées de taille n (m1 et m2)
 * et AJOUTER le résultat dans la matrice res
 * res = res + m1*m2
 */
void produitPartielMM(double **m1, double **m2,double **res, int n){
    int i,j,k;
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            for (k=0;k<n;k++){
                res[i][j] += m1[i][k]*m2[k][j];
            }
        }
    }
}

/**calculer le produit d'une matrice carrée n*n avec un vecteur de taille n
*/
void produitMV(double **m, double *x, double *res, int n){
    int i,j;
    for (i=0;i<n;i++){
        res[i]=0;
        for (j=0;j<n;j++){
            res[i] += m[i][j]*x[j];
        }
    }
}

/**Cette fonction retourne un tableau de taille nbLinges*nbColonnes
 * alloué en mémoire CONTIGUE
 */
void allocation(double *** tableau, int nbLignes,int nbColonnes){

    /* allouer nbColonnes*nbLignes cases continues dans la mémoire */
    double *p = (double *)malloc(nbColonnes*nbLignes*sizeof(double));
    if (p==NULL) {exit(EXIT_FAILURE);}

    /* Allouer le pointeur des lignes */
    (*tableau) = (double **)malloc(nbLignes*sizeof(double*));
    if (*tableau==NULL) {
       free(p);
       exit(EXIT_FAILURE);
    }
    /* faire corespondant entre le pointeur et
    la mémoire allouée */
    int i;
    for (i=0; i<nbLignes; i++){
       (*tableau)[i] = &(p[i*nbColonnes]);
    }
}

/** Cette fonction libère la mémoire associée à un tableau */
void my_free(double ***tableau){
    /*Libérer la mémoire par l'adresse
    du premier élément du tableau*/
    free(&((*tableau)[0][0]));

    /*Libérer le pointeur*/
    free(*tableau);
}


/**Calculer la norme d'un vecteur de taille n*/
double normeVecteur(double * x, int n){
    int j;
    double norme = 0;
    for(j=0;j<n;j++){
        norme += x[j]*x[j];
    }
    norme = (double) sqrt(norme);
    return norme;
}
//Copier de dest = source
void clone(double **dest, double **source, int n){
    int i,j;
    for (i=0;i<n;i++){
            for (j=0;j<n;j++){
                dest[i][j] = source[i][j];
            }
        }

}

/**
Cette fonction calcule la décomposition QR d'une matrice A et écrire les matrices de la décomposition dans Q et R
n est la taille des matrices
p est le nombre de processus
rank est le rang du processus qui appelle la fonction
(c.f rapport section "La décomposition QR")
*/
void QRfactorization(double **A, double ** Q, double ** R, int n,int p, int rank){
    int k;
    double *x; 
    double *u;
    double **H = NULL; // La matrice de réflexion de Householder
    double **Q_temp = NULL; //Matrice temporaire pour stocker résultat après chaque opération
    if (rank == 0){
        allocation(&H,n,n);
        allocation(&Q_temp,n,n);

        //Initialiser Q = I
        int i,j;
        for (i=0;i<n;i++){
            for(j=0;j<n;j++){
                if (i==j){
                    Q[i][j] = 1;
                }
                else {
                    Q[i][j] = 0;
                }
            }
        }
    }
    for (k=0; k<n-1;k++){
        if (rank ==0){
            //Extraire le vecteur x à partir d'une colonne de A. x est de taille n-k
            x = malloc((n-k)*sizeof(double));
            u = malloc((n-k)*sizeof(double));
            int i,j;
            for (j=0;j<(n-k);j++){
                x[j]=A[k+j][k];
            }

            //Calculer u normé associé à la matrice de Householder de x
            int signx1 = x[0] < 0 ? (-1) : 1;
            for (i=0;i<(n-k);i++){u[i] = 0; }
            u[0] = signx1 * normeVecteur(x,n-k);
            for (j=0;j<(n-k);j++){
                u[j]+= x[j];
            }
            double normeU = normeVecteur(u,n-k);
            for (j=0;j<(n-k);j++){
                u[j] = u[j]/normeU ;
            }

            //Calculer la matrice Hk
            // réinitialiser H =I
            for (i=0;i<n;i++){
                for(j=0;j<n;j++){
                    if (i==j){
                        H[i][j] = 1;
                    }
                    else {
                        H[i][j] = 0;
                    }
                }
            }
            //écrire sur le bloc en ba à droite de H
            for (i=0;i<n-k;i++){
                for (j=0;j<n-k;j++){
                    if (i==j){
                        H[i+k][j+k] = 1-(2*u[i]*u[j]);
                    }
                    else{
                        H[i+k][j+k] = -2*u[i]*u[j];
                    }
                }
            }
        }
        // R = H*A
        Produit_2Matrices_Parallele(H,A,R,n,p,rank);

        //A = R
        if (rank ==0){clone(A,R,n);}

        //Q_temp = Q*H
        Produit_2Matrices_Parallele(Q,H,Q_temp,n,p,rank);

        //Q = Q_temp
        if (rank ==0){
            clone(Q,Q_temp,n);
            free(x);
            free(u);
        }
    }
    if (rank ==0){
        my_free(&H);
        my_free(&Q_temp);
    }

}
