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

int main(int argc, char *argv[]){
    int rank; 		/* Mon rang*/
    int p;          /*Le nombre de processus*/
    int n = 12; 		/*La taille de la matrice carrée de connectivité du net*/
    int n_local;    /*La taille de la matrice locale pour chaque processeur*/
    int tag1= 1;
    int tag2 = 2;
    int tag3 = 3;

    MPI_Status status;  /*valeur de retour pour le récepteur*/
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
    /**La mission du processus 0 est de créer les deux matrices m1 et m2,
    * les diviser en blocs puis
    * envoyer les blocs aux autres processus en fonction de leur rang.
    * A la fin, il récupère les matrice res_local
    * et reconstruit la matrice du résultat global */

    //Création de 2 matrices m1 et m2
        double **m1;
        allocation(&m1,n,n);
        double **m2;
        allocation(&m2,n,n);
        int i,j;
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
            m1[i][j] = 0;
            m2[i][j] = 0;
            }
        }

        //Initialisation
        m1[0][0] = 0;m1[0][1] = 0;m1[0][2] = 1;m1[0][3] = 0;m1[0][4] = 0;m1[0][5] = 1;m1[0][6] = 0;m1[0][7] = 0;m1[0][8] = 0;m1[0][9] = 1;m1[0][10] = 0;m1[0][11] = 0;
        m1[1][0] = 1;m1[1][1] = 0;m1[1][2] = 0;m1[1][3] = 0;m1[1][4] = 0;m1[1][5] = 0;m1[1][6] = 0;m1[1][7] = 0;m1[1][8] = 1;m1[1][9] = 0;m1[1][10] = 0;m1[1][11] = 0;
        m1[2][0] = 0;m1[2][1] = 1;m1[2][2] = 1;m1[2][3] = 1;m1[2][4] = 1;m1[2][5] = 0;m1[2][6] = 0;m1[2][7] = 0;m1[2][8] = 0;m1[2][9] = 0;m1[2][10] = 0;m1[2][11] = 0;
        m1[3][0] = 0;m1[3][1] = 0;m1[3][2] = 1;m1[3][3] = 0;m1[3][4] = 0;m1[3][5] = 0;m1[3][6] = 0;m1[3][7] = 0;m1[3][8] = 0;m1[3][9] = 0;m1[3][10] = 1;m1[3][11] = 1;
        m1[4][0] = 1;m1[4][1] = 0;m1[4][2] = 1;m1[4][3] = 0;m1[4][4] = 0;m1[4][5] = 0;m1[4][6] = 0;m1[4][7] = 1;m1[4][8] = 1;m1[4][9] = 0;m1[4][10] = 0;m1[4][11] = 0;
        m1[5][0] = 0;m1[5][1] = 0;m1[5][2] = 0;m1[5][3] = 1;m1[5][4] = 1;m1[5][5] = 0;m1[5][6] = 0;m1[5][7] = 0;m1[5][8] = 0;m1[5][9] = 0;m1[5][10] = 1;m1[5][11] = 0;
        m1[6][0] = 0;m1[6][1] = 1;m1[6][2] = 0;m1[6][3] = 0;m1[6][4] = 0;m1[6][5] = 0;m1[6][6] = 1;m1[6][7] = 0;m1[6][8] = 1;m1[6][9] = 1;m1[6][10] = 0;m1[6][11] = 0;
        m1[7][0] = 1;m1[7][1] = 0;m1[7][2] = 0;m1[7][3] = 0;m1[7][4] = 1;m1[7][5] = 0;m1[7][6] = 0;m1[7][7] = 0;m1[7][8] = 0;m1[7][9] = 0;m1[7][10] = 0;m1[7][11] = 0;
        m1[8][0] = 0;m1[8][1] = 1;m1[8][2] = 0;m1[8][3] = 0;m1[8][4] = 0;m1[8][5] = 1;m1[8][6] = 0;m1[8][7] = 0;m1[8][8] = 0;m1[8][9] = 0;m1[8][10] = 0;m1[8][11] = 0;
        m1[9][0] = 0;m1[9][1] = 1;m1[9][2] = 0;m1[9][3] = 1;m1[9][4] = 0;m1[9][5] = 0;m1[9][6] = 0;m1[9][7] = 0;m1[9][8] = 0;m1[9][9] = 0;m1[9][10] = 0;m1[9][11] = 1;
        

        m2[0][0] = 0;m2[0][1] = 0;m2[0][2] = 0;m2[0][3] = 0;m2[0][4] = 1;m2[0][5] = 0;m2[0][6] = 0;m2[0][7] = 0;m2[0][8] = 0;m2[0][9] = 0;m2[0][10] = 0;m2[0][11] = 0;
        m2[1][0] = 0;m2[1][1] = 0;m2[1][2] = 0;m2[1][3] = 0;m2[1][4] = 0;m2[1][5] = 1;m2[1][6] = 1;m2[1][7] = 1;m2[1][8] = 0;m2[1][9] = 0;m2[1][10] = 1;m2[1][11] = 0;
        m2[2][0] = 0;m2[2][1] = 0;m2[2][2] = 0;m2[2][3] = 1;m2[2][4] = 0;m2[2][5] = 0;m2[2][6] = 0;m2[2][7] = 0;m2[2][8] = 0;m2[2][9] = 0;m2[2][10] = 1;m2[2][11] = 1;
        m2[3][0] = 0;m2[3][1] = 0;m2[3][2] = 0;m2[3][3] = 0;m2[3][4] = 0;m2[3][5] = 0;m2[3][6] = 0;m2[3][7] = 0;m2[3][8] = 0;m2[3][9] = 0;m2[3][10] = 0;m2[3][11] = 0;
        m2[4][0] = 1;m2[4][1] = 0;m2[4][2] = 0;m2[4][3] = 0;m2[4][4] = 0;m2[4][5] = 1;m2[4][6] = 0;m2[4][7] = 0;m2[4][8] = 0;m2[4][9] = 1;m2[4][10] = 0;m2[4][11] = 0;
        m2[5][0] = 0;m2[5][1] = 1;m2[5][2] = 0;m2[5][3] = 0;m2[5][4] = 0;m2[5][5] = 0;m2[5][6] = 0;m2[5][7] = 0;m2[5][8] = 0;m2[5][9] = 0;m2[5][10] = 0;m2[5][11] = 0;
        m2[6][0] = 0;m2[6][1] = 0;m2[6][2] = 1;m2[6][3] = 0;m2[6][4] = 1;m2[6][5] = 0;m2[6][6] = 0;m2[6][7] = 0;m2[6][8] = 0;m2[6][9] = 1;m2[6][10] = 0;m2[6][11] = 1;
        m2[7][0] = 0;m2[7][1] = 0;m2[7][2] = 0;m2[7][3] = 0;m2[7][4] = 0;m2[7][5] = 0;m2[7][6] = 0;m2[7][7] = 0;m2[7][8] = 0;m2[7][9] = 0;m2[7][10] = 1;m2[7][11] = 0;
        m2[8][0] = 0;m2[8][1] = 0;m2[8][2] = 0;m2[8][3] = 1;m2[8][4] = 0;m2[8][5] = 0;m2[8][6] = 0;m2[8][7] = 0;m2[8][8] = 1;m2[8][9] = 0;m2[8][10] = 0;m2[8][11] = 0;
        m2[9][0] = 1;m2[9][1] = 0;m2[9][2] = 0;m2[9][3] = 0;m2[9][4] = 0;m2[9][5] = 0;m2[9][6] = 0;m2[9][7] = 0;m2[9][8] = 1;m2[9][9] = 0;m2[9][10] = 0;m2[9][11] = 0;
        // Afficher les matrices m1 et m2
        printf("Matrice m1 \n");
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
            printf("%lf ",m1[i][j]);
            }
            printf("\n");
        }
        printf("Matrice m2 \n");
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
            printf("%lf ",m2[i][j]);
            }
            printf("\n");
        }

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
                m1_local[i][j] = m1[i][j];
                m2_local[i][j] = m2[i][j];
            }
        }
		//libérer la mémoire
		my_free(&m1);
        my_free(&m2);
    }
    else {
        //printf("coucou from processus %d \n", rank);

        // recevoir les matrices extraites stocké dans m1_local et m2_local
        MPI_Recv(&(m1_local[0][0]), n_local*n_local, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &status);
        MPI_Recv(&(m2_local[0][0]), n_local*n_local, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD, &status);
        //printf("Matrice extraite reçue pour le procesus %d \n", rank);


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


        /*********VERIFICATION: UNCOMMENT POUR AFFICHER m1_temp ET m2_temp REÇUS à chaque itération K******

        printf("This is processus %d. Matrice m1_temp \n", rank);
        for (i=0;i<n_local;i++){
            for (j=0;j<n_local;j++){
            printf("%lf ",m1_temp[i][j]);
            }
            printf("\n");
        }
        printf("This is processus %d. Matrice m2_temp\n", rank);
        for (i=0;i<n_local;i++){
            for (j=0;j<n_local;j++){
            printf("%lf ",m2_temp[i][j]);
            }
            printf("\n");
        }
        /************************************/
        produitPartielMM(m1_temp,m2_temp,res_local,n_local);
    }
    MPI_Barrier(MPI_COMM_WORLD);


    if (rank!=0){
        //envoyer les résultats locaux vers le processus 0
        MPI_Send(&(res_local[0][0]),n_local*n_local, MPI_DOUBLE,0,tag3,MPI_COMM_WORLD);
    }
    else{
        double **res;
        allocation(&res,n,n);
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
        //Afficher le résultat final
        printf("Matrice res = m1*m2 :\n");
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
            printf("%lf ",res[i][j]);
            }
            printf("\n");
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
