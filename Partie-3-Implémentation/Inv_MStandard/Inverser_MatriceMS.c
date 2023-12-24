#include<stdio.h>
#include<stdlib.h>

double** initialiserMatrice(int n){
    int i,j;
    double **M = (double**)malloc(n*sizeof(double*));
    for (i = 0; i < n; ++i) {
        M[i] = (double*)malloc(n * sizeof(double));
    }
    for (i = 0; i < n; i++) {
        for ( j = 0; j < n; j++) {
            M[i][j] = 0;
        }
    }
    return M;
}

void libererMatrice(int n, double** M){
    
    for (int i = 0; i < n; ++i) {
        free(M[i]);
    }
    free(M);
}

void afficherMatrice(int n, double **M){
    int i,j;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            printf("%lf ",M[i][j]);
        }
        printf("\n");
    }
    printf("\n ---------------------------------- \n");
}

double** multiplicationMatrices(int n, double** A, double** B){
    int i, j, k;
    double** M=initialiserMatrice(n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                M[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return M;
}

double** transposerMatrices(int n ,double** M){
    int i,j;
    double** Mt=initialiserMatrice(n);

    for (i = 0; i < n; i++) {
        for ( j = 0; j < n; j++) {
            Mt[i][j] = M[j][i];
        }
    }
    return Mt;
}


void decomposerMatrices(int n, double** M, double** B, double** C, double** Ct, double** D){
    
    int i,j;

    for(i=0; i<n/2 ; i++){
        for(j=0; j<n/2 ; j++){
            B[i][j]=M[i][j];        
        }    
    }
    
    for(i=0;i<n/2;i++){
        for(j=0;j<n/2;j++){
            C[i][j]=M[i+n/2][j];        
        }    
    }
    // Prendre Ct directement de la matrice. ( - d'appel de fonction Transpose)
    for(i=0;i<n/2;i++){
        for(j=0;j<n/2;j++){
            Ct[i][j]=M[i][j+n/2];        
        }    
    }

    for(i=0;i<n/2;i++){
        for(j=0;j<n/2;j++){
            D[i][j]=M[i+n/2][j+n/2];        
        }    
    }
}


double** soustractionMatrices(int n, double** A, double** B) {
    
    double** M=initialiserMatrice(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            M[i][j] = A[i][j] - B[i][j];
        }
    }
    return M;
}


double** additionMatrices(int n, double** A, double** B) {
    
    double** M=initialiserMatrice(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            M[i][j] = A[i][j] + B[i][j];
        }
    }
    return M;
}


double** changer_signe_Matrices(int n,double** B) {
    
    double** M=initialiserMatrice(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            M[i][j] = - B[i][j];
        }
    }
    return M;
}


double** concatenerMatrices(int n, double** A, double** B, double** C, double** D){
    
    int i,j;
    double** M=initialiserMatrice(n*2);
    
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            M[i][j]=A[i][j];        
        }    
    }
    
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            M[i][j+n]=B[i][j];   
        }    
    }

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            M[i+n][j]=C[i][j];        
        }    
    }

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            M[i+n][j+n]=D[i][j];        
        }    
    }
    return M;
}


double** inverserMatrices(int n, double** M){

    double** B=initialiserMatrice(n/2);
    double** D=initialiserMatrice(n/2);
    double** C=initialiserMatrice(n/2);
    double** Ct=initialiserMatrice(n/2);
    
    int taille=n/2;

    decomposerMatrices(n,M,B,C,Ct,D);

    if(taille==1){

        B[0][0]=1/B[0][0];

        double** CBinv=multiplicationMatrices(taille,C,B);
        double** BinvCt=transposerMatrices(taille,CBinv);        
        double** S=multiplicationMatrices(taille,C,BinvCt);
        S=soustractionMatrices(taille,D,S);
        
        double** Sinv=initialiserMatrice(1);
        Sinv[0][0]=1/S[0][0];
        double** SinvCBinv=multiplicationMatrices(taille,Sinv,CBinv);
        double** CBinv_tSinv=transposerMatrices(taille,SinvCBinv);
        double** BinvCtSinvCBinv=multiplicationMatrices(taille,BinvCt,SinvCBinv);
        double** I=additionMatrices(taille,B,BinvCtSinvCBinv);

        // --------------- Assemblage Final 2x2 -------------------//
        CBinv_tSinv=changer_signe_Matrices(taille,CBinv_tSinv);
        SinvCBinv=changer_signe_Matrices(taille,SinvCBinv);
        double** Minv=concatenerMatrices(taille,I,CBinv_tSinv,SinvCBinv,Sinv);
        libererMatrice(taille,B);
        libererMatrice(taille,C);
        libererMatrice(taille,Ct);
        libererMatrice(taille,D);
        libererMatrice(taille,CBinv);
        libererMatrice(taille,BinvCt);
        libererMatrice(taille,S);
        libererMatrice(taille,Sinv);
        libererMatrice(taille,SinvCBinv);
        libererMatrice(taille,CBinv_tSinv);
        libererMatrice(taille,BinvCtSinvCBinv);

        return Minv; 
    }
    
    double** Binv=inverserMatrices(taille,B); //inverse de B 
    
    double** CBinv=multiplicationMatrices(taille,C,Binv); // C*inverse de B
    double** BinvCt=transposerMatrices(taille,CBinv); //Transposé C*inverse de B

    double** S= multiplicationMatrices(taille,CBinv,Ct); 
    S=soustractionMatrices(taille,D,S); 

    double** Sinv=inverserMatrices(taille,S); //Inverse de S
    double** SinvCBinv=multiplicationMatrices(taille,Sinv,CBinv); // Inverse S * C * inverse de B 
    double** CBinv_tSinv=transposerMatrices(taille,SinvCBinv); // Transpose(Inverse S * C * inverse de B)
    double** BinvCtSinvCBinv=multiplicationMatrices(taille,BinvCt,SinvCBinv); //Transposer C*inverse de B * Inverse S * C * inverse de B 

    double** I=additionMatrices(taille,Binv,BinvCtSinvCBinv); // Partie Top-Left Matrice
    
    // --------------- Assemblage Final -------------------//
    CBinv_tSinv=changer_signe_Matrices(taille,CBinv_tSinv); 
    SinvCBinv=changer_signe_Matrices(taille,SinvCBinv);
    double** Minv=concatenerMatrices(taille,I,CBinv_tSinv,SinvCBinv,Sinv);
    
    // ------------ Liberation des matrices -------------//
    libererMatrice(taille,B);
    libererMatrice(taille,C);
    libererMatrice(taille,Ct);
    libererMatrice(taille,D);
    libererMatrice(taille,Binv);
    libererMatrice(taille,CBinv);
    libererMatrice(taille,BinvCt);
    libererMatrice(taille,S);
    libererMatrice(taille,Sinv);
    libererMatrice(taille,SinvCBinv);
    libererMatrice(taille,CBinv_tSinv);
    libererMatrice(taille,BinvCtSinvCBinv);
    libererMatrice(taille,I);

    return Minv;
}


int main(){

    FILE *fichier_A;
    int n, i, j;

    //Chargement des fichiers

    fichier_A = fopen("matrice.txt", "r");
    if (fichier_A == NULL){
        fprintf(stderr, "Erreur : fichier introuvable dans le repertoire courant\n");
        return EXIT_FAILURE;
    }

    // Lecture de la taille de la matrice carré
    fscanf(fichier_A, "%d", &n); 
    
    // Initialisation de la matrice
    double **A = initialiserMatrice(n);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            fscanf(fichier_A, "%lf", &A[i][j]);
        }
    }
    
    // -------------------------------------------------- Esapce de Test -------------------------------------------------
    
    double** AtAinv=inverserMatrices(n,multiplicationMatrices(n,transposerMatrices(n,A),A));
    double** Ainv=multiplicationMatrices(n,AtAinv,transposerMatrices(n,A));
    
    printf("Ainv : \n");
    afficherMatrice(n,Ainv);

    // Libération de mémoire et fermeture des fichiers 
    libererMatrice(n,A);
    libererMatrice(n,AtAinv);
    libererMatrice(n,Ainv);

    fclose(fichier_A);

    return EXIT_SUCCESS;
}