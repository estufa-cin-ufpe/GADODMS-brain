#ifndef EKF_P
#define EKF_p



#include <math.h>
#define Npontos 10
#define Nsta 1
#define Mobs 1

typedef struct ekf_s
{
    double x[Nsta];
    //double xp[Nsta];    /* state vector */
    //double z[Mobs];   /* observation vector */

    double P[Nsta*Nsta];  /* prediction error covariance */
    double Q[Nsta*Nsta];  /* process noise covariance */
    double R[Mobs*Mobs];  /* measurement error covariance */

    double G[Nsta*Mobs];  /* Kalman gain; a.k.a. K */

    double F[Nsta*Nsta];  /* Jacobian of process model */
    double H[Mobs*Nsta];  /* Jacobian of measurement model */

    double Ht[Nsta*Mobs]; /* transpose of measurement Jacobian */
    double Ft[Nsta*Nsta]; /* transpose of process Jacobian */
    double Pp[Nsta*Nsta]; /* P, post-prediction, pre-update */

    double hx[Mobs];
    double fx[Nsta];

    double tmp0[Nsta*Nsta];
    double tmp1[Nsta*Mobs];
    double tmp2[Mobs*Nsta];
    double tmp3[Mobs*Mobs];
    double tmp4[Mobs*Mobs];
    double tmp5[Mobs];

} EKF;

EKF pontos[Npontos];

int choldc1(double * a, double * p, int n) {
    int i,j,k;
    double sum;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = a[i*n+j];
            for (k = i - 1; k >= 0; k--) {
                sum -= a[i*n+k] * a[j*n+k];
            }
            if (i == j) {
                if (sum <= 0) {
                    return 1; /* error */
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j*n+i] = sum / p[i];
            }
        }
    }

    return 0; /* success */
}

int choldcsl(double * A, double * a, double * p, int n)
{
    int i,j,k; double sum;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            a[i*n+j] = A[i*n+j];
    if (choldc1(a, p, n)) return 1;
    for (i = 0; i < n; i++) {
        a[i*n+i] = 1 / p[i];
        for (j = i + 1; j < n; j++) {
            sum = 0;
            for (k = i; k < j; k++) {
                sum -= a[j*n+k] * a[k*n+i];
            }
            a[j*n+i] = sum / p[j];
        }
    }

    return 0; /* success */
}


int cholsl(double * A, double * a, double * p, int n)
{
    int i,j,k;
    if (choldcsl(A,a,p,n)) return 1;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            a[i*n+j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
        a[i*n+i] *= a[i*n+i];
        for (k = i + 1; k < n; k++) {
            a[i*n+i] += a[k*n+i] * a[k*n+i];
        }
        for (j = i + 1; j < n; j++) {
            for (k = j; k < n; k++) {
                a[i*n+j] += a[k*n+i] * a[k*n+j];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            a[i*n+j] = a[j*n+i];
        }
    }

    return 0; /* success */
}

void zeros(double * a, int m, int n)
{
    int j;
    for (j=0; j<m*n; ++j)
        a[j] = 0;
}


/* C <- A * B */
void mulmat(double * a, double * b, double * c, int arows, int acols, int bcols)
{
    int i, j,l;

    for(i=0; i<arows; ++i)
        for(j=0; j<bcols; ++j) {
            c[i*bcols+j] = 0;
            for(l=0; l<acols; ++l)
                c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
        }
}

void mulvec(double * a, double * x, double * y, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i) {
        y[i] = 0;
        for(j=0; j<n; ++j)
            y[i] += x[j] * a[i*n+j];
    }
}

void transpose(double * a, double * at, int m, int n)
{
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j) {
            at[j*m+i] = a[i*n+j];
        }
}

/* A <- A + B */
void accum(double * a, double * b, int m, int n)
{
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] += b[i*n+j];
}

/* C <- A + B */
void add(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] + b[j];
}


/* C <- A - B */
void sub(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] - b[j];
}

void negate(double * a, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] = -a[i*n+j];
}

void mat_addeye(double * a, int n)
{
    int i;
    for (i=0; i<n; ++i)
        a[i*n+i] += 1;
}



void inicializar(int index){

    for(int i=0;i<Nsta;i++){
        pontos[index].P[i*Nsta+1] = 1;
        for(int j=0; j<Nsta; j++){
            pontos[index].F[i*Nsta+j] = 0;
            
            pontos[index].Pp[i*Nsta+j] = 0;
            pontos[index].Q[i*Nsta+j] = 0;
            pontos[index].Ft[i*Nsta+j] = 0;
            pontos[index].tmp0[i*Nsta+j] = 0;

        }
    }

    for(int i=0;i<Mobs;i++){
        for(int j=0; j<Mobs; j++){
            pontos[index].R[i*Mobs+j] = 0;
            pontos[index].tmp3[i*Mobs+j] = 0;
            pontos[index].tmp4[i*Mobs+j] = 0;
        }
    }

    for(int i=0;i<Nsta;i++){
        for(int j=0; j<Mobs; j++){
            pontos[index].G[i*Mobs+j] = 0;
            pontos[index].Ht[i*Mobs+j] = 0;
            pontos[index].tmp1[i*Mobs+j] = 0;
        }
    }

    for(int i=0;i<Mobs;i++){
        for(int j=0; j<Nsta; j++){
            pontos[index].H[i*Nsta+j] = 0;
            pontos[index].tmp2[i*Nsta+j] = 0;
        }
        pontos[index].tmp5[i]=0;
    }

    for(int i=0;i<Nsta;i++){
        pontos[index].x[i]=0;
    }
    

    for(int i=0;i<Nsta;i++){
        pontos[index].Q[i*Nsta+i] = 0.59;
        pontos[index].P[i*Nsta+i] = 15;
    }

    for(int i=0; i<Mobs;i++){
        pontos[index].R[i*Mobs+i] = 0.96;
    }

    //x = {0.0};

}

void modelar(int index){
    //std::cout << "delta = " << delta << std::endl;
    //Serial.print("                                      clk: ");
    //Serial.println(delta);
    pontos[index].F[0] = 1;;
    /* pontos[index].F[1] = 0;
    pontos[index].F[2] = 0;

    pontos[index].F[3] = 0;
    pontos[index].F[4] = 1;
    pontos[index].F[5] = 0;

    pontos[index].F[6] = 0;
    pontos[index].F[7] = 0;
    pontos[index].F[8] = 1;*/

    pontos[index].H[0] = 1;
    /*pontos[index].H[1] = 0;
    pontos[index].H[2] = 0;

    pontos[index].H[3] = 0;
    pontos[index].H[4] = 1;
    pontos[index].H[5] = 0;

    pontos[index].H[6] = 0;
    pontos[index].H[7] = 0;
    pontos[index].H[8] = 1;*/
   


    transpose(pontos[index].F, pontos[index].Ft, Nsta, Nsta);
    transpose(pontos[index].H, pontos[index].Ht, Mobs, Nsta);
}



int passo(double* z, int index){
    //std::cout << "entrou em passso" << std::endl;

    /* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
       mulmat(pontos[index].F, pontos[index].P, pontos[index].tmp0, Nsta, Nsta, Nsta);
       transpose(pontos[index].F, pontos[index].Ft, Nsta, Nsta);
       mulmat(pontos[index].tmp0, pontos[index].Ft, pontos[index].Pp, Nsta, Nsta, Nsta);
       accum(pontos[index].Pp, pontos[index].Q, Nsta, Nsta);

       /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
       transpose(pontos[index].H, pontos[index].Ht, Mobs, Nsta);
       mulmat(pontos[index].Pp, pontos[index].Ht, pontos[index].tmp1, Nsta, Nsta, Mobs);
       mulmat(pontos[index].H, pontos[index].Pp, pontos[index].tmp2, Mobs, Nsta, Nsta);
       mulmat(pontos[index].tmp2, pontos[index].Ht, pontos[index].tmp3, Mobs, Nsta, Mobs);
       accum(pontos[index].tmp3, pontos[index].R, Mobs, Mobs);
       if (cholsl(pontos[index].tmp3, pontos[index].tmp4, pontos[index].tmp5, Mobs)) return 1;
       mulmat(pontos[index].tmp1, pontos[index].tmp4, pontos[index].G, Nsta, Mobs, Mobs);

       //Serial.print("G: ");
       //Serial.print(G[0]);
       //Serial.print("     ");
       //Serial.println(G[1]);
       

       /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
       mulvec(pontos[index].H, pontos[index].x, pontos[index].hx, Mobs, Nsta);//////////////////
       sub(z, pontos[index].hx, pontos[index].tmp5, Mobs);
       mulvec(pontos[index].G, pontos[index].tmp5, pontos[index].tmp2, Nsta, Mobs);
       mulvec(pontos[index].F, pontos[index].x, pontos[index].fx, Nsta, Nsta);/////////////////
       add(pontos[index].fx, pontos[index].tmp2, pontos[index].x, Nsta);

       /* P_k = (I - G_k H_k) P_k */
       mulmat(pontos[index].G, pontos[index].H, pontos[index].tmp0, Nsta, Mobs, Nsta);
       negate(pontos[index].tmp0, Nsta, Nsta);
       mat_addeye(pontos[index].tmp0, Nsta);
       mulmat(pontos[index].tmp0, pontos[index].Pp, pontos[index].P, Nsta, Nsta, Nsta);
}

double getX(int index, int index_ponto){
    return pontos[index_ponto].x[index];
}

#endif
