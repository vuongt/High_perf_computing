#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

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

/**Calculer la norme d'un vectuer de taille n*/
double normeVecteur(double * x, int n){
    int j;
    double norme = 0;
    for(j=0;j<n;j++){
        norme += x[j]*x[j];
    }
    norme = (double) sqrt(norme);
    return norme;
}

/**Cette fonction retourne un tableau de taille nbLinges*nbColonnes
 * alloué en mémoire CONTIGUE
 */
void allocation(double *** tableau, int nbLignes,int nbColonnes){

    /* allouer nbColonnes*nbLignes cases continues dans la mémoire */
    double *p = (double *)malloc(nbColonnes*nbLignes*sizeof(double));
    if (p==NULL) {exit(EXIT_FAILURE);}

    /* allouer les pointers de lignes dans la mémoire */
    (*tableau) = (double **)malloc(nbLignes*sizeof(double*));
    if (*tableau==NULL) {
       free(p);
       exit(EXIT_FAILURE);
    }
    /* faire correspondant entre les pointeur du tableau et la mémoire allouée avant */
    int i;
    for (i=0; i<nbLignes; i++){
       (*tableau)[i] = &(p[i*nbColonnes]);
    }

}
/** Cette fonction libère la mémoire associée à un tableau */
void my_free(double ***tableau){
    free(&((*tableau)[0][0]));
    free(*tableau);
}

int main(int argc, char *argv[]){
    int const n = 12; /*La taille de la matrice carrée de connectivité du net*/
    int nbIteration = 50;
    int rank;
    int p;          /*Le nombre de processus*/
    int n_local;    /*La taille de la matrice locale pour chaque processeur*/
    int tag1= 1;
    int tag2 = 2;
    int tag3 = 3;

    MPI_Status status; /*valeur de retour pour le récepteur*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

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
    printf("Le nombre de blocs par dimension est: %d \n", q);
    printf("La taille de chaque bloc est: %d \n", n_local);
	int I = 1 + rank%q;
	int J = 1 + rank/q;
	printf("Processus de rang %d, le couple d'indice associé est: I = %d,J = %d \n",rank, I, J);

    MPI_Comm commutateur_ligne;
    MPI_Comm commutateur_colonne;
	MPI_Comm_split(MPI_COMM_WORLD,I,J,&commutateur_ligne);
    MPI_Comm_split(MPI_COMM_WORLD,J,I,&commutateur_colonne);

    /*Allouer la mémoire et initialiser les matrices locales pour chaque processeur*/
    /*Il est nécessaire d'allouer 2 matrices locales même s'elles sont identiques au début
    * parce que ma matrice m1_locale va etre remplacé en A², A³,..
    * dans la boucle d'itération
    */
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
    /**La mission du processus 0 est de créer la matrice A,
    * la diviser en blocs puis
    * envoyer les blocs aux autres processus en fonction de leur rang.
    * A la fin, il récupère les matrice res_local
    * et reconstruit la matrice du résultat global */
    
        double ** A;
        allocation(&A,n,n);


        int i,j;
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
                A[i][j] = 0;
            }
        }
        //Initialisation
        A[0][0] = (double)0  ;  A[0][1] = (double)0  ;  A[0][2] = (double)1/8;  A[0][3] = (double)0;  A[0][4] = (double)1/8;  A[0][5] = (double)1/3;  A[0][6] = (double) 1/2;  A[0][7] = (double) 0  ;  A[0][8] = (double) 0  ;  A[0][9] = (double) 0  ; A[0][10] = (double) 0  ;  A[0][11] = (double) 0  ;
        A[1][0] = (double)1/6;  A[1][1] = (double)0  ;  A[1][2] = (double)0  ;  A[1][3] = (double)0;  A[1][4] = (double)1/8;  A[1][5] = (double)0;  A[1][6] = (double) 1/2;  A[1][7] = (double) 0  ;  A[1][8] = (double) 1/5  ;  A[1][9] = (double) 1/4  ; A[1][10] = (double) 1/3 ;  A[1][11] = (double) 0  ;
        A[2][0] = (double)0  ;  A[2][1] = (double)0  ;  A[2][2] = (double)1/8;  A[2][3] = (double)1;  A[2][4] = (double)1/8;  A[2][5] = (double)0  ;  A[2][6] = (double) 0  ;  A[2][7] = (double) 1/5;  A[2][8] = (double) 0  ;  A[2][9] = (double) 0  ; A[2][10] = (double) 0  ;  A[2][11] = (double) 0  ;
        A[3][0] = (double)1/6;  A[3][1] = (double)1/4;  A[3][2] = (double)1/8;  A[3][3] = (double)0;  A[3][4] = (double)0  ;  A[3][5] = (double)0  ;  A[3][6] = (double) 0  ;  A[3][7] = (double) 1/5;  A[3][8] = (double) 1/5  ;  A[3][9] = (double) 1/4  ; A[3][10] = (double) 0  ;  A[3][11] = (double) 0  ;
        A[4][0] = (double)0  ;  A[4][1] = (double)1/4;  A[4][2] = (double)0  ;  A[4][3] = (double)0;  A[4][4] = (double)1/8;  A[4][5] = (double)0  ;  A[4][6] = (double) 0  ;  A[4][7] = (double) 0;  A[4][8] = (double) 0  ;  A[4][9] = (double) 0  ; A[4][10] = (double) 1/3 ;  A[4][11] = (double) 1/2 ;
        A[5][0] = (double)0  ;  A[5][1] = (double)1/4;  A[5][2] = (double)1/8;  A[5][3] = (double)0;  A[5][4] = (double)1/8;  A[5][5] = (double)0  ;  A[5][6] = (double) 0  ;  A[5][7] = (double) 1/5;  A[5][8] = (double) 0  ;  A[5][9] = (double) 0  ; A[5][10] = (double) 0  ;  A[5][11] = (double) 0  ;
        A[6][0] = (double)1/6;  A[6][1] = (double)1/4;  A[6][2] = (double)1/8;  A[6][3] = (double)0;  A[6][4] = (double)0  ;  A[6][5] = (double)1/3;  A[6][6] = (double) 0  ;  A[6][7] = (double) 0  ;  A[6][8] = (double) 1/5  ;  A[6][9] = (double) 1/4  ; A[6][10] = (double) 0  ;  A[6][11] = (double) 0  ;
        A[7][0] = (double)0  ;  A[7][1] = (double)0  ;  A[7][2] = (double)1/8;  A[7][3] = (double)0;  A[7][4] = (double)0  ;  A[7][5] = (double)0  ;  A[7][6] = (double) 0  ;  A[7][7] = (double) 0;  A[7][8] = (double) 0  ;  A[7][9] = (double) 0  ; A[7][10] = (double) 0  ;  A[7][11] = (double) 0  ;
        A[8][0] = (double)1/6;  A[8][1] = (double)0  ;  A[8][2] = (double)0  ;  A[8][3] = (double)0;  A[8][4] = (double)1/8;  A[8][5] = (double)0  ;  A[8][6] = (double) 0  ;  A[8][7] = (double) 0;  A[8][8] = (double) 1/5  ;  A[8][9] = (double) 1/4  ; A[8][10] = (double) 1/3 ;  A[8][11] = (double) 0  ;
        A[9][0] = (double)0  ;  A[9][1] = (double)0  ;  A[9][2] = (double)1/8;  A[9][3] = (double)0;  A[9][4] = (double)1/8;  A[9][5] = (double)0  ;  A[9][6] = (double) 0  ;  A[9][7] = (double) 0;  A[9][8] = (double) 0  ;  A[9][9] = (double) 0  ; A[9][10] = (double) 0  ;  A[9][11] = (double) 1/2 ;
        A[10][0] = (double)1/6; A[10][1] = (double)0  ; A[10][2] = (double)1/8; A[10][3] = (double)0; A[10][4] = (double)0  ; A[10][5] = (double)0  ; A[10][6] = (double) 0  ; A[10][7] = (double) 1/5; A[10][8] = (double) 0  ; A[10][9] = (double) 0  ;A[10][10] = (double) 0  ; A[10][11] = (double) 0  ;
        A[11][0] = (double)1/6; A[11][1] = (double)0  ; A[11][2] = (double)0  ; A[11][3] = (double)0; A[11][4] = (double)1/8; A[11][5] = (double)1/3  ; A[11][6] = (double) 0  ; A[11][7] = (double) 1/5; A[11][8] = (double) 1/5  ; A[11][9] = (double) 0  ;A[11][10] = (double) 0  ; A[11][11] = (double) 0  ;

        //**A = {{0,1/2,1/4,0},{1/3,0,1/4,0},{1/3,0,1/4,1},{1/3,1/2,1/4,0}};
        //Afficher la matrice A
        printf("Matrice A \n");
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
            printf("%lf ",A[i][j]);
            }
            printf("\n");
        }

        /*Diviser la matrice A et envoyer vers les processus*/
		int k; // le rang du destinataire
		for (k = 1; k<p; k++){
            /*copy les blocs dans les matrices m1_local et m2_local
            * m1_local et m2_local stock le bloc (I,J) de A
            * I = 1 + k%q;
            * J = 1 + k/q;
            */
            for (i=0;i<n_local;i++){
                for (j=0;j<n_local;j++){
                    m1_local[i][j] = A[(k%q)*n/q +i][(k/q)*n/q + j];
                    m2_local[i][j] = A[(k%q)*n/q +i][(k/q)*n/q + j];
                }
            }
            //Envoyer vers le processus k
            MPI_Send(&(m1_local[0][0]),n_local*n_local, MPI_DOUBLE,k,tag1,MPI_COMM_WORLD);
            MPI_Send(&(m2_local[0][0]),n_local*n_local, MPI_DOUBLE,k,tag2,MPI_COMM_WORLD);
            //printf("Matrice extraite envoyée ver le processus %d \n", k);

            /******VERIFICATION: UNCOMMENT POUR AFFICHER m1_local ET m_2 local ENVOYÉS******
            printf("Envoie vers processus %d. Matrice m1 local \n", k);
            for (i=0;i<n_local;i++){
                for (j=0;j<n_local;j++){
                printf("%lf ",m1_local[i][j]);
                }
                printf("\n");
            }
            printf("Envoie vers processus %d. Matrice m2 local \n", k);
            for (i=0;i<n_local;i++){
                for (j=0;j<n_local;j++){
                printf("%lf ",m2_local[i][j]);
                }
                printf("\n");
            }
            /********************************/
		}
		/*En fin, lors que tous les blocs sont envoyés,
		* le processus 0 copy les blocs (0,0) de m1 et m2
		* pour effectuer un calcul local lui-même
		*/
        for (i=0;i<n_local;i++){
            for (j=0;j<n_local;j++){
                m1_local[i][j] = A[i][j];
                m2_local[i][j] = A[i][j];
            }
        }
        my_free(&A);
    }
    else {
        //printf("coucou from processus %d \n", rank);
        // recevoir les matrices extraites
        //stocker les matrice extraite dans m1_local et m2_local
        MPI_Recv(&(m1_local[0][0]), n_local*n_local, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &status);
        MPI_Recv(&(m2_local[0][0]), n_local*n_local, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD, &status);

        /******VERIFICATION : UNCOMMENT POUR AFFICHER m1_local ET m2_local REÇUS********
        int i,j;
        printf("This is processus %d. Matrice m1 local \n", rank);
        for (i=0;i<n_local;i++){
            for (j=0;j<n_local;j++){
            printf("%lf ",m1_local[i][j]);
            }
            printf("\n");
        }
        printf("This is processus %d. Matrice m2 local \n", rank);
        for (i=0;i<n_local;i++){
            for (j=0;j<n_local;j++){
            printf("%lf ",m2_local[i][j]);
            }
            printf("\n");
        }
        /************************************/
    }

    /*****************BROARDCASTING ET CALCUL DU RESULTAT LOCAL*******************/
    int rank_ligne;
    MPI_Comm_rank(commutateur_ligne, &rank_ligne);
    int rank_colonne;
    MPI_Comm_rank(commutateur_colonne, &rank_colonne);
    //uncomment la ligne suivante pour vérifier les ranks ligne et colonne
    //printf("This is processus %d, my rank in commutateur ligne is: %d, my rank in commutateur colonne is: %d \n", rank, rank_ligne, rank_colonne);
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
    int m;
    //On répète le calcul nbIteration fois
    for (m=1;m<=nbIteration;m++){
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

        //remplacer m1_local par res_local puis continuer la boucle sur nbIteration
        for (i=0;i<n_local;i++){
            for (j=0;j<n_local;j++){
                m1_local[i][j] = res_local[i][j];
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /***************FINI CALCUL DU RESULTAT LOCAL****************/


    /*****************Recupérer les résulats locaux*****************/
    if (rank!=0){
    //Envoyer les résultats locaux vers le processus 0
        MPI_Send(&(res_local[0][0]),n_local*n_local, MPI_DOUBLE,0,tag3,MPI_COMM_WORLD);
    }
    else{
        double **res;
        allocation(&res,n,n);
        int i,j;
        //Cloner le résultat
        for (i=0;i<n_local;i++){
            for (j=0;j<n_local;j++){
                res[i][j] = res_local[i][j];
            }
        }
        int k; // rang du processus qui envoie le résultat
        for (k=1;k<p;k++){
            MPI_Recv(&(res_local[0][0]), n_local*n_local, MPI_DOUBLE, k, tag3, MPI_COMM_WORLD, &status);
            /*******************VERIFICATION Résultat local reçu**************/
            printf("This is processus %d receive from %d. res_local:\n", rank,k);
            for (i=0;i<n_local;i++){
                for (j=0;j<n_local;j++){
                printf("%lf ",res_local[i][j]);
                }
                printf("\n");
            }
            /************************************/

            //Cloner le résultat local au résultat global
            for (i=0;i<n_local;i++){
                for (j=0;j<n_local;j++){
                    res[(k%q)*n/q +i][(k/q)*n/q + j] = res_local[i][j];
                }
            }
        }
        printf("Matrice résultat A puissance %d:\n",nbIteration);
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
            printf("%lf ",res[i][j]);
            }
            printf("\n");
        }

        //Vecteur b est le vecteur propre à chercher de A, contient les score d'importance
        double *b = malloc(n*sizeof(double));
        if (b==NULL){
            exit(EXIT_FAILURE);
        }
        // tmp est un vecteur temporaire de taille n
        double *tmp = malloc(n*sizeof(double));
        if (tmp == NULL){
            exit(EXIT_FAILURE);
        }

        //Initiation les coefficients de b à 1/n et tmp à 0
        for (i=0;i<n;i++){
            b[i]= (double)1/n;
            tmp[i] = 0;
        }

        //temp = A^N*b
        produitMV(res,b,tmp,n);

        double normeTmp = normeVecteur(tmp,n);
        // b = tmp/norme(tmp)
        printf("Vecteur b: \n");
        for (i=0;i<n;i++){
            b[i] =(double) tmp[i]/normeTmp;
            printf("%lf \n", b[i]);
        }
        my_free(&res);

    }
    my_free(&m1_temp);
    my_free(&m2_temp);

    my_free(&m1_local);
    my_free(&m2_local);
    my_free(&res_local);
    MPI_Comm_free(&commutateur_colonne);
    MPI_Comm_free(&commutateur_ligne);

    MPI_Finalize();

    return 0;
}
