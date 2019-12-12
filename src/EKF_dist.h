#ifndef EKF_P
#define EKF_p



#include <math.h>
#define Nvacas 10
#define Nsta_dist 2
#define Mobs_dist 2

typedef struct ekf_s_dist
{
    double x[Nsta_dist];
    //double xp[Nsta_dist];    /* state vector */
    //double z[Mobs_dist];   /* observation vector */

    double P[Nsta_dist*Nsta_dist];  /* prediction error covariance */
    double Q[Nsta_dist*Nsta_dist];  /* process noise covariance */
    double R[Mobs_dist*Mobs_dist];  /* measurement error covariance */

    double G[Nsta_dist*Mobs_dist];  /* Kalman gain; a.k.a. K */

    double F[Nsta_dist*Nsta_dist];  /* Jacobian of process model */
    double H[Mobs_dist*Nsta_dist];  /* Jacobian of measurement model */

    double Ht[Nsta_dist*Mobs_dist]; /* transpose_dist of measurement Jacobian */
    double Ft[Nsta_dist*Nsta_dist]; /* transpose_dist of process Jacobian */
    double Pp[Nsta_dist*Nsta_dist]; /* P, post-prediction, pre-update */

    double hx[Mobs_dist];
    double fx[Nsta_dist];

    double tmp0[Nsta_dist*Nsta_dist];
    double tmp1[Nsta_dist*Mobs_dist];
    double tmp2[Mobs_dist*Nsta_dist];
    double tmp3[Mobs_dist*Mobs_dist];
    double tmp4[Mobs_dist*Mobs_dist];
    double tmp5[Mobs_dist];

} EKF_dist;

EKF_dist vacas[Nvacas];

int choldc1_dist(double * a, double * p, int n) {
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

int choldcsl_dist(double * A, double * a, double * p, int n)
{
    int i,j,k; double sum;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            a[i*n+j] = A[i*n+j];
    if (choldc1_dist(a, p, n)) return 1;
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


int cholsl_dist(double * A, double * a, double * p, int n)
{
    int i,j,k;
    if (choldcsl_dist(A,a,p,n)) return 1;
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

void zeros_dist(double * a, int m, int n)
{
    int j;
    for (j=0; j<m*n; ++j)
        a[j] = 0;
}


/* C <- A * B */
void mulmat_dist(double * a, double * b, double * c, int arows, int acols, int bcols)
{
    int i, j,l;

    for(i=0; i<arows; ++i)
        for(j=0; j<bcols; ++j) {
            c[i*bcols+j] = 0;
            for(l=0; l<acols; ++l)
                c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
        }
}

void mulvec_dist(double * a, double * x, double * y, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i) {
        y[i] = 0;
        for(j=0; j<n; ++j)
            y[i] += x[j] * a[i*n+j];
    }
}

void transpose_dist(double * a, double * at, int m, int n)
{
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j) {
            at[j*m+i] = a[i*n+j];
        }
}

/* A <- A + B */
void accum_dist(double * a, double * b, int m, int n)
{
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] += b[i*n+j];
}

/* C <- A + B */
void add_dist(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] + b[j];
}


/* C <- A - B */
void sub_dist(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] - b[j];
}

void negate_dist(double * a, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] = -a[i*n+j];
}

void mat_add_disteye(double * a, int n)
{
    int i;
    for (i=0; i<n; ++i)
        a[i*n+i] += 1;
}



void inicializar_dist(int index){

    for(int i=0;i<Nsta_dist;i++){
        vacas[index].P[i*Nsta_dist+1] = 1;
        for(int j=0; j<Nsta_dist; j++){
            vacas[index].F[i*Nsta_dist+j] = 0;

            vacas[index].Pp[i*Nsta_dist+j] = 0;
            vacas[index].Q[i*Nsta_dist+j] = 0;
            vacas[index].Ft[i*Nsta_dist+j] = 0;
            vacas[index].tmp0[i*Nsta_dist+j] = 0;

        }
    }

    for(int i=0;i<Mobs_dist;i++){
        for(int j=0; j<Mobs_dist; j++){
            vacas[index].R[i*Mobs_dist+j] = 0;
            vacas[index].tmp3[i*Mobs_dist+j] = 0;
            vacas[index].tmp4[i*Mobs_dist+j] = 0;
        }
    }

    for(int i=0;i<Nsta_dist;i++){
        for(int j=0; j<Mobs_dist; j++){
            vacas[index].G[i*Mobs_dist+j] = 0;
            vacas[index].Ht[i*Mobs_dist+j] = 0;
            vacas[index].tmp1[i*Mobs_dist+j] = 0;
        }
    }

    for(int i=0;i<Mobs_dist;i++){
        for(int j=0; j<Nsta_dist; j++){
            vacas[index].H[i*Nsta_dist+j] = 0;
            vacas[index].tmp2[i*Nsta_dist+j] = 0;
        }
        vacas[index].tmp5[i]=0;
    }

    for(int i=0;i<Nsta_dist;i++){
        vacas[index].x[i]=0;
    }


    for(int i=0;i<Nsta_dist;i++){
        vacas[index].Q[i*Nsta_dist+i] = 0.15;
        vacas[index].P[i*Nsta_dist+i] = 1500.0;
    }

    for(int i=0; i<Mobs_dist;i++){
        vacas[index].R[i*Mobs_dist+i] = 70.0;
    }

    //x = {0.0};

}

void modelar_dist(int index){
    //std::cout << "delta = " << delta << std::endl;
    //Serial.print("                                      clk: ");
    //Serial.println(delta);
    vacas[index].F[0] = 1;
    vacas[index].F[1] = 0;
    vacas[index].F[2] = 0;
    vacas[index].F[3] = 1;
    /* vacas[index].F[1] = 0;
    vacas[index].F[2] = 0;

    vacas[index].F[3] = 0;
    vacas[index].F[4] = 1;
    vacas[index].F[5] = 0;

    vacas[index].F[6] = 0;
    vacas[index].F[7] = 0;
    vacas[index].F[8] = 1;*/

    vacas[index].H[0] = 1;
    vacas[index].H[1] = 0;
    vacas[index].H[2] = 0;
    vacas[index].H[3] = 1;
    /*vacas[index].H[1] = 0;
    vacas[index].H[2] = 0;

    vacas[index].H[3] = 0;
    vacas[index].H[4] = 1;
    vacas[index].H[5] = 0;

    vacas[index].H[6] = 0;
    vacas[index].H[7] = 0;
    vacas[index].H[8] = 1;*/



    transpose_dist(vacas[index].F, vacas[index].Ft, Nsta_dist, Nsta_dist);
    transpose_dist(vacas[index].H, vacas[index].Ht, Mobs_dist, Nsta_dist);
}



int passo_dist(double* z, int index){
    //std::cout << "entrou em passso" << std::endl;

    /* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
       mulmat_dist(vacas[index].F, vacas[index].P, vacas[index].tmp0, Nsta_dist, Nsta_dist, Nsta_dist);
       transpose_dist(vacas[index].F, vacas[index].Ft, Nsta_dist, Nsta_dist);
       mulmat_dist(vacas[index].tmp0, vacas[index].Ft, vacas[index].Pp, Nsta_dist, Nsta_dist, Nsta_dist);
       accum_dist(vacas[index].Pp, vacas[index].Q, Nsta_dist, Nsta_dist);

       /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
       transpose_dist(vacas[index].H, vacas[index].Ht, Mobs_dist, Nsta_dist);
       mulmat_dist(vacas[index].Pp, vacas[index].Ht, vacas[index].tmp1, Nsta_dist, Nsta_dist, Mobs_dist);
       mulmat_dist(vacas[index].H, vacas[index].Pp, vacas[index].tmp2, Mobs_dist, Nsta_dist, Nsta_dist);
       mulmat_dist(vacas[index].tmp2, vacas[index].Ht, vacas[index].tmp3, Mobs_dist, Nsta_dist, Mobs_dist);
       accum_dist(vacas[index].tmp3, vacas[index].R, Mobs_dist, Mobs_dist);
       if (cholsl_dist(vacas[index].tmp3, vacas[index].tmp4, vacas[index].tmp5, Mobs_dist)) return 1;
       mulmat_dist(vacas[index].tmp1, vacas[index].tmp4, vacas[index].G, Nsta_dist, Mobs_dist, Mobs_dist);

       //Serial.print("G: ");
       //Serial.print(G[0]);
       //Serial.print("     ");
       //Serial.println(G[1]);


       /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
       mulvec_dist(vacas[index].H, vacas[index].x, vacas[index].hx, Mobs_dist, Nsta_dist);//////////////////
       sub_dist(z, vacas[index].hx, vacas[index].tmp5, Mobs_dist);
       mulvec_dist(vacas[index].G, vacas[index].tmp5, vacas[index].tmp2, Nsta_dist, Mobs_dist);
       mulvec_dist(vacas[index].F, vacas[index].x, vacas[index].fx, Nsta_dist, Nsta_dist);/////////////////
       add_dist(vacas[index].fx, vacas[index].tmp2, vacas[index].x, Nsta_dist);

       /* P_k = (I - G_k H_k) P_k */
       mulmat_dist(vacas[index].G, vacas[index].H, vacas[index].tmp0, Nsta_dist, Mobs_dist, Nsta_dist);
       negate_dist(vacas[index].tmp0, Nsta_dist, Nsta_dist);
       mat_add_disteye(vacas[index].tmp0, Nsta_dist);
       mulmat_dist(vacas[index].tmp0, vacas[index].Pp, vacas[index].P, Nsta_dist, Nsta_dist, Nsta_dist);
}

double getX_dist(int index, int index_ponto){
    return vacas[index_ponto].x[index];
}

#endif
