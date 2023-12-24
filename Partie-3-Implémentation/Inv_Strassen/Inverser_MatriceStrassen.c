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
    printf("\n---------------------------------- \n");
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

// ------------------- Debut Strassen ----------------- //

double** decomposerStrassen(double **M, int n, int lig, int col) {
	int i,j,r=lig,c=col;
	double** M1 = initialiserMatrice(n/2);	
	
	for(i = 0;i < n/2; i++) {
		c=col;
   	 	for(j = 0; j < n/2; j++) {
        	M1[i][j] = M[r][c];
			c++;
    		}
		r++;
	}
	return M1;
}

void concatenerStrassen(double** A, double** M, int row, int col, int n){
	int i,j,r=row,c=col;
	for(i=0;i<n;i++){
		c=col;
		for(j=0;j<n;j++){
			M[r][c]=A[i][j];	
			c++;	
		}
		r++;
	}
}

double** multiplicationStrassen(double** A, double** B, int n){
	double ** result = initialiserMatrice(n);
	if(n>1) {

		double** a11=decomposerStrassen(A,n,0,0);
		double** a12=decomposerStrassen(A,n,0,(n/2));
		double** a21=decomposerStrassen(A,n,(n/2),0);
		double** a22=decomposerStrassen(A,n,(n/2),(n/2));	
		double** b11=decomposerStrassen(B,n,0,0);
		double** b12=decomposerStrassen(B,n,0,n/2);
		double** b21=decomposerStrassen(B,n,n/2,0);
		double** b22=decomposerStrassen(B,n,n/2,n/2);
		
        // ------------------------ Appel recursive ---------------------

		double** m1= multiplicationStrassen(additionMatrices(n/2,a11,a22),additionMatrices(n/2,b11,b22),n/2);
		double** m2= multiplicationStrassen(additionMatrices(n/2,a21,a22),b11,n/2);
		double** m3= multiplicationStrassen(a11,soustractionMatrices(n/2,b12,b22),n/2);
		double** m4= multiplicationStrassen(a22,soustractionMatrices(n/2,b21,b11),n/2);
		double** m5= multiplicationStrassen(additionMatrices(n/2,a11,a12),b22,n/2);
		double** m6= multiplicationStrassen(soustractionMatrices(n/2,a21,a11),additionMatrices(n/2,b11,b12),n/2);
		double** m7= multiplicationStrassen(soustractionMatrices(n/2,a12,a22),additionMatrices(n/2,b21,b22),n/2);

        // ------------------------ Assemblage -------------------------

		double** c11 = additionMatrices(n/2,soustractionMatrices(n/2,additionMatrices(n/2,m1,m4),m5),m7);
		double** c12 = additionMatrices(n/2,m3,m5);
		double** c21 = additionMatrices(n/2,m2,m4);
		double** c22 = additionMatrices(n/2,soustractionMatrices(n/2,additionMatrices(n/2,m1,m3),m2),m6);
		
        concatenerStrassen(c11,result,0,0,n/2);
		concatenerStrassen(c12,result,0,n/2,n/2);
		concatenerStrassen(c21,result,n/2,0,n/2);
		concatenerStrassen(c22,result,n/2,n/2,n/2);
	} 
	else {
		result[0][0]=A[0][0]*B[0][0]; //Taille 1
	}
	return result;
}
// ------------------- FIN Strassen ----------------- //

// ------------------ Partie Inverse --------------- //

double** inverserMatrices(int n, double** M){

    double** B=initialiserMatrice(n/2);
    double** D=initialiserMatrice(n/2);
    double** C=initialiserMatrice(n/2);
    double** Ct=initialiserMatrice(n/2);
    
    int taille=n/2;

    decomposerMatrices(n,M,B,C,Ct,D);

    if(taille==1){

        B[0][0]=1/B[0][0];

        double** CBinv=multiplicationStrassen(C,B,taille);
        double** BinvCt=transposerMatrices(taille,CBinv);        
        double** S=multiplicationStrassen(C,BinvCt,taille);
        S=soustractionMatrices(taille,D,S);
        
        double** Sinv=initialiserMatrice(1); 
        Sinv[0][0]=1/S[0][0];

        double** SinvCBinv=multiplicationStrassen(Sinv,CBinv,taille);
        double** CBinv_tSinv=transposerMatrices(taille,SinvCBinv);
        double** BinvCtSinvCBinv=multiplicationStrassen(BinvCt,SinvCBinv,taille);
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
    
    double** CBinv=multiplicationStrassen(C,Binv,taille); // C*inverse de B
    double** BinvCt=transposerMatrices(taille,CBinv); //Transposé C*inverse de B

    double** S= multiplicationStrassen(CBinv,Ct,taille); 
    S=soustractionMatrices(taille,D,S);

    double** Sinv=inverserMatrices(taille,S); //Inverse de S
    double** SinvCBinv=multiplicationStrassen(Sinv,CBinv,taille); // Inverse S * C * inverse de B 
    double** CBinv_tSinv=transposerMatrices(taille,SinvCBinv); // Transpose(Inverse S * C * inverse de B)
    double** BinvCtSinvCBinv=multiplicationStrassen(BinvCt,SinvCBinv,taille); //Transposer C*inverse de B * Inverse S * C * inverse de B 

    double** I=additionMatrices(taille,Binv,BinvCtSinvCBinv); // Partie Top-Left Matrice inversé
    
    // --------------- Assemblage Final -------------------//
    CBinv_tSinv=changer_signe_Matrices(taille,CBinv_tSinv); 
    SinvCBinv=changer_signe_Matrices(taille,SinvCBinv);
    double** Minv=concatenerMatrices(taille,I,CBinv_tSinv,SinvCBinv,Sinv);
    
    // ---------- Liberation des matrices ----------------//
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

    
    // Initialisation de la matrice
    fscanf(fichier_A, "%d", &n);
    double **A = initialiserMatrice(n);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            fscanf(fichier_A, "%lf", &A[i][j]);
        }
    }
    
    // -------------------------------------------------- Esapce de Test -------------------------------------------------
    
    double** AtAinv=inverserMatrices(n,multiplicationStrassen(transposerMatrices(n,A),A,n));
    double** Ainv=multiplicationStrassen(AtAinv,transposerMatrices(n,A),n);
    
    printf("Ainv : \n");
    afficherMatrice(n,Ainv);
    // Libération de mémoire et fermeture des fichiers 
    
    libererMatrice(n,A);
    libererMatrice(n,AtAinv);
    libererMatrice(n,Ainv);

    fclose(fichier_A);


    return EXIT_SUCCESS;
}