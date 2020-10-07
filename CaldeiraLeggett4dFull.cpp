// ==============================================================================
//
//  CaldeiraLeggett4dFull.cpp
//  QTR
//
//  Created by Albert Lu on 1/30/20.
//  alu@tacc.utexas.edu
//
//  Last modified on 1/30/20
//
//  Note:
//
// ==============================================================================

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <omp.h>
#include <vector>
#include <parallel/algorithm>
#include <new>

#include "Constants.h"
#include "Containers.h"
#include "Error.h"
#include "Log.h"
#include "Parameters.h"
#include "CaldeiraLeggett4d.h"
 
using namespace QTR_NS;
using std::vector;
using std::cout;
using std::endl;
using std::nothrow;

/* ------------------------------------------------------------------------------- */

// DEFINE POTENTIAL TYPE

#if defined CLPOT_DW1

#define WAVEFUNCTION(x1,x2,x3,x4) Wavefunction_DW1(x1,x2,x3,x4)
#define POTENTIAL(x1,x2) Potential_DW1(x1,x2)
#define POTENTIAL_X1(x1,x2) Vx1_DW1(x1,x2)
#define POTENTIAL_X2(x1,x2) Vx2_DW1(x1,x2)
#define POTENTIAL_XXX1(x1, x2) Vxxx1_DW1(x1,x2)
#define POTENTIAL_XXX2(x1, x2) Vxxx2_DW1(x1,x2)
#define POTNAME "DoubleWell-1"

#elif defined CLPOT_DW2

#define WAVEFUNCTION(x1,x2,x3,x4) Wavefunction_DW2(x1,x2,x3,x4)
#define POTENTIAL(x1,x2) Potential_DW2(x1,x2)
#define POTENTIAL_X1(x1,x2) Vx1_DW2(x1,x2)
#define POTENTIAL_X2(x1,x2) Vx2_DW2(x1,x2)
#define POTENTIAL_XXX1(x1, x2) Vxxx1_DW2(x1,x2)
#define POTENTIAL_XXX2(x1, x2) Vxxx2_DW2(x1,x2)
#define POTNAME "DoubleWell-2"

#elif defined CLPOT_DW3

#define WAVEFUNCTION(x1,x2,x3,x4) Wavefunction_DW3(x1,x2,x3,x4)
#define POTENTIAL(x1,x2) Potential_DW3(x1,x2)
#define POTENTIAL_X1(x1,x2) Vx1_DW3(x1,x2)
#define POTENTIAL_X2(x1,x2) Vx2_DW3(x1,x2)
#define POTENTIAL_XXX1(x1, x2) Vxxx1_DW3(x1,x2)
#define POTENTIAL_XXX2(x1, x2) Vxxx2_DW3(x1,x2)
#define POTNAME "DoubleWell-3"

#endif

/* ------------------------------------------------------------------------------- */

CaldeiraLeggett4d::CaldeiraLeggett4d(class QTR *q)
{
    qtr = q;
    err = qtr->error;
    log = qtr->log;
    parameters = qtr->parameters;
    init();
} 
/* ------------------------------------------------------------------------------- */

CaldeiraLeggett4d::~CaldeiraLeggett4d()
{     
    return;
}
/* ------------------------------------------------------------------------------- */

void CaldeiraLeggett4d::init()
{
    log->log("\n\n[CaldeiraLeggett4d] INIT starts ...\n");
    log->log("\n\n[CaldeiraLeggett4d] Using Full Grid\n");
    log->log("\n\n[CaldeiraLeggett4d] Potential type: %s\n", POTNAME);

    // General parameters
    I = {0,1}; // sqrt(-1)
    PI_INV = 1.0 / PI; // 1/PI
    xZERO = {0,0}; // complex zero
    DIMENSIONS = parameters->scxd_dimensions;
    PERIOD = parameters->scxd_period;
    SORT_PERIOD = parameters->scxd_sortperiod;
    PRINT_PERIOD = parameters->scxd_printperiod;
    PRINT_WAVEFUNC_PERIOD = parameters->scxd_printwavefuncperiod;
    TIME = parameters->scxd_Tf;
    QUIET = parameters->quiet;
    TIMING = parameters->timing;
    isTrans = parameters->scxd_isTrans;
    isCorr = parameters->scxd_isAcf;
    isPrintEdge = parameters->scxd_isPrintEdge;
    isPrintDensity = parameters->scxd_isPrintDensity;
    isPrintWavefunc = parameters->scxd_isPrintWavefunc;

    // Fine grid
    H.resize(DIMENSIONS);
    Hi.resize(DIMENSIONS);
    Hisq.resize(DIMENSIONS);
    S.resize(DIMENSIONS);  
    kk = parameters->scxd_k;
    H[0] = parameters->scxd_h1;
    H[1] = parameters->scxd_h2;
    H[2] = parameters->scxd_h3;
    H[3] = parameters->scxd_h4;

    for (int i = 0; i < DIMENSIONS; i ++)  {
        Hi[i] = 1 / H[i];
        Hisq[i] = 1 / pow(H[i],2);
        S[i] = kk * Hisq[i];
    }

    // Domain size and # grids
    Box.resize(DIMENSIONS * 2);
    Box[0] = parameters->scxd_xi1;
    Box[1] = parameters->scxd_xf1;
    Box[2] = parameters->scxd_xi2;
    Box[3] = parameters->scxd_xf2;
    Box[4] = parameters->scxd_xi3;
    Box[5] = parameters->scxd_xf3; 
    Box[6] = parameters->scxd_xi4;
    Box[7] = parameters->scxd_xf4; 
    BoxShapeFine.resize(DIMENSIONS);

    // Parameters
    hb = parameters->scxd_hb;
    m  = parameters->scxd_m;
    kb = parameters->scxd_kb;
    temp = parameters->scxd_temp;
    gamma = parameters->scxd_gamma;
    k0 = parameters->scxd_k0;
    sigma = parameters->scxd_sigma;
    lambda = parameters->scxd_lambda;
    trans_x0 = parameters->scxd_trans_x0;
    HBSQ_INV = 1.0 / (hb * hb);
    PIHBSQ_INV = 1.0 / (PI * PI * hb * hb);

    // Wavefunction parameters
    Wave0.resize(DIMENSIONS);
    Wave0[0] = parameters->scxd_x01;
    Wave0[1] = parameters->scxd_x02;
    Wave0[2] = parameters->scxd_x03;
    Wave0[3] = parameters->scxd_x04;

    A.resize(DIMENSIONS);
    A[0] = parameters->scxd_a1;
    A[1] = parameters->scxd_a2;
    A[2] = parameters->scxd_a3;
    A[3] = parameters->scxd_a4;

    log->log("[CaldeiraLeggett4d] INIT done.\n\n");
}
/* ------------------------------------------------------------------------------- */

void CaldeiraLeggett4d::Evolve()
{
    log->log("[CaldeiraLeggett4d] Evolve starts ...\n");

    // Files
    FILE *pfile;

    // Variables 
    int count;
    int n1, n2, n3, n4;
    unsigned long int index;
    double sum;
    double density;
    double coeff1 = 0.0;
    double coeff2 = 0.0;
    double norm;   // normalization factor
    bool b1, b2, b3, b4, b5, b6, b7;

    double TolHdX1 = 2.0 * TolHd * H[0];
    double TolHdX2 = 2.0 * TolHd * H[1];
    double TolHdX3 = 2.0 * TolHd * H[2];
    double TolHdX4 = 2.0 * TolHd * H[3];
    double TolLdX1 = 2.0 * TolLd * H[0];
    double TolLdX2 = 2.0 * TolLd * H[1];
    double TolLdX3 = 2.0 * TolLd * H[2];
    double TolLdX4 = 2.0 * TolLd * H[3];

    // Timing variables
    double t_0_begin, t_0_end;
    double t_1_begin, t_1_end;
    double t_0_elapsed = 0.0;
    double t_1_elapsed = 0.0;

    // Constants
    double val;
    double kh0m = kk / (H[0] * m);
    double kh1m = kk / (H[1] * m);
    double k2h2 = kk / H[2];
    double k2h3 = kk / H[3];
    double i2h2 = 1.0 / (2.0 * H[2]);
    double i2h3 = 1.0 / (2.0 * H[3]);
    double kgamma = kk * gamma;
    double khbsq2h2 = kk * hb * hb / 24.0 / (H[2] * H[2] * H[2]);
    double khbsq2h3 = kk * hb * hb / 24.0 / (H[3] * H[3] * H[3]);
    double mkT2h2sq = m * kb * temp / (H[2] * H[2]);
    double mkT2h3sq = m * kb * temp / (H[3] * H[3]);

    double corr;
    double corr_0;

    //  2d Grid vector and indices
    VectorXi grid;
    int g1, g2, g3, g4;
    double xx1, xx2, xx3, xx4;
    double f0, kk0;
    double f1p1, f1m1, f1p2, f1m2;
    double f2p1, f2m1, f2p2, f2m2;
    double f3p1, f3m1, f3p2, f3m2;
    double f4p1, f4m1, f4p2, f4m2;
    double kk1p1, kk1m1, kk1p2, kk1m2;
    double kk2p1, kk2m1, kk2p2, kk2m2;
    double kk3p1, kk3m1, kk3p2, kk3m2;
    double kk4p1, kk4m1, kk4p2, kk4m2;

    // Vector iterater
    vector<unsigned long int>::iterator it;

    // PF_trans
    double pftrans;
    vector<double> PF_trans;
    PF_trans.push_back(0.0);

    // Compute coarse grid shape
    t_1_begin = omp_get_wtime();
    log->log("[CaldeiraLeggett4d] Fine grid size = (");

    GRIDS_TOT = 1;
    for (int i = 0; i < DIMENSIONS; i ++)  {
       BoxShapeFine[i] = (int)( (Box[2 * i + 1] - Box[2 * i]) / H[i] ) + 1;
       GRIDS_TOT *= BoxShapeFine[i];

        if ( i < DIMENSIONS - 1 )
            log->log("%d, ", BoxShapeFine[i]);
        else
            log->log("%d)=", BoxShapeFine[i]);
    }
    log->log("%lu\n",GRIDS_TOT);

    M1 = BoxShapeFine[1] * BoxShapeFine[2] * BoxShapeFine[3];
    M2 = BoxShapeFine[2] * BoxShapeFine[3];
    M3 = BoxShapeFine[3];
    W1 = BoxShapeFine[1];
    W2 = BoxShapeFine[3];
    O1 = BoxShapeFine[0] * BoxShapeFine[1];
    O2 = BoxShapeFine[2] * BoxShapeFine[3];

    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    if ( !QUIET && TIMING ) log->log("[CaldeiraLeggett4d] Elapsed time (make bounding box) = %lf sec\n\n", t_1_elapsed); 

    // Initialize containers
    t_1_begin = omp_get_wtime();

    log->log("[CaldeiraLeggett4d] Initializing containers ...\n");

    double **KK1 = new double*[O1];

    for (unsigned long int  i = 0; i < O1; i++)  {
        KK1[i] = new(nothrow) double[O2];

        if (!KK1[i])
            log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "KK1");
    }

    double **KK2 = new double*[O1];

    for (unsigned long int  i = 0; i < O1; i++)  {
        KK2[i] = new(nothrow) double[O2];

        if (!KK2[i])
            log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "KK2");
    }

    double **KK3 = new double*[O1];

    for (unsigned long int  i = 0; i < O1; i++)  {
        KK3[i] = new(nothrow) double[O2];

        if (!KK3[i])
            log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "KK3");
    }

    double **KK4 = new double*[O1];

    for (unsigned long int  i = 0; i < O1; i++)  {
        KK4[i] = new(nothrow) double[O2];

        if (!KK4[i])
            log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "KK4");
    }

    double **F = new double*[O1];

    for (unsigned long int  i = 0; i < O1; i++)  {
        F[i] = new(nothrow) double[O2];

        if (!F[i])
            log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "F");
    }

    double **FF = new double*[O1];

    for (unsigned long int  i = 0; i < O1; i++)  {
        FF[i] = new(nothrow) double[O2];

        if (!FF[i])
            log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "FF");
    }

    double **PF = new double*[O1];

    for (unsigned long int  i = 0; i < O1; i++)  {
        PF[i] = new(nothrow) double[O2];

        if (!PF[i])
            log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "PF");
    }

    #pragma omp parallel for
    for (int i1 = 0; i1 < BoxShapeFine[0]; i1 ++)  {
        for (int i2 = 0; i2 < BoxShapeFine[1]; i2 ++)  {
            for (int i3 = 0; i3 < BoxShapeFine[2]; i3 ++)  {
                for (int i4 = 0; i4 < BoxShapeFine[3]; i4 ++)  {
                    F[i1*W1+i2][i3*W2+i4] = 0.0;
                    PF[i1*W1+i2][i3*W2+i4] = 0.0;
                    FF[i1*W1+i2][i3*W2+i4] = 0.0;
                    KK1[i1*W1+i2][i3*W2+i4] = 0.0;
                    KK2[i1*W1+i2][i3*W2+i4] = 0.0;
                    KK3[i1*W1+i2][i3*W2+i4] = 0.0;
                    KK4[i1*W1+i2][i3*W2+i4] = 0.0;
                }
            }
        }
    }
    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    if ( !QUIET && TIMING ) log->log("[CaldeiraLeggett4d] Elapsed time (initializing containers) = %lf sec\n\n", t_1_elapsed); 

    // .........................................................................................

    // Initialize wavefunction
    t_1_begin = omp_get_wtime();
    log->log("[CaldeiraLeggett4d] Initializing wavefunction ...\n");  

    #pragma omp parallel for
    for (int i1 = 1; i1 < BoxShapeFine[0]-1; i1 ++)  {
        for (int i2 = 1; i2 < BoxShapeFine[1]-1; i2 ++)  {
            for (int i3 = 1; i3 < BoxShapeFine[2]-1; i3 ++)  {
                for (int i4 = 1; i4 < BoxShapeFine[3]-1; i4 ++)  {
                    F[i1*W1+i2][i3*W2+i4] = WAVEFUNCTION(Box[0]+i1*H[0], Box[2]+i2*H[1], Box[4]+i3*H[2], Box[6]+i4*H[3]);
                }
            }
        }
    }
    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    if (!QUIET && TIMING) log->log("[CaldeiraLeggett4d] Elapsed time (initializing wavefunction) = %lf sec\n\n", t_1_elapsed); 

    // Normalization
    t_1_begin = omp_get_wtime();
    norm = 0.0;

    #pragma omp parallel for reduction (+:norm)
    for (int i1 = 0; i1 <  BoxShapeFine[0]; i1 ++)  {
        for (int i2 = 0; i2 < BoxShapeFine[1]; i2 ++)  {
            for (int i3 = 0; i3 < BoxShapeFine[2]; i3 ++)  {
                for (int i4 = 0; i4 < BoxShapeFine[3]; i4 ++)  {
                    norm += F[i1*W1+i2][i3*W2+i4];
                }
            }
        }
    }
    norm *= H[0] * H[1] * H[2] * H[3];
    log->log("[CaldeiraLeggett4d] Normalization factor = %.16e\n",norm);
    norm = 1.0 / norm;

    #pragma omp parallel for
    for (int i1 = 0; i1 < BoxShapeFine[0]; i1 ++)  {
        for (int i2 = 0; i2 < BoxShapeFine[1]; i2 ++)  {
            for (int i3 = 0; i3 < BoxShapeFine[2]; i3 ++)  {
                for (int i4 = 0; i4 < BoxShapeFine[3]; i4 ++)  {
                    F[i1*W1+i2][i3*W2+i4] = norm * F[i1*W1+i2][i3*W2+i4];
                    PF[i1*W1+i2][i3*W2+i4] = F[i1*W1+i2][i3*W2+i4];
                }
            }
        }
    }
    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    if (!QUIET && TIMING) log->log("[CaldeiraLeggett4d] Elapsed time (Wavefunction normalization) = %lf sec\n\n", t_1_elapsed); 

    // Initial density

    VectorXd F0;
    VectorXd Ft;

    if ( isCorr )  {

        t_1_begin = omp_get_wtime();
        
        F0.resize(O1);
        Ft.resize(O1);

        #pragma omp parallel for
        for (unsigned long int i1 = 0; i1 < O1; i1 ++)  {
            F0[i1] = 0.0;
            Ft[i1] = 0.0;
        }

        #pragma omp parallel for private(density)
        for (int i1 = 2; i1 < BoxShapeFine[0] - 2; i1 ++)  {
            for (int i2 = 2; i2 < BoxShapeFine[1] - 2; i2 ++)  {
                density = 0.0;
                for (int i3 = 2; i3 < BoxShapeFine[2]-2; i3 ++)  {
                    for (int i4 = 2; i4 < BoxShapeFine[3]-2; i4 ++)  {      
                        density += PF[i1*W1+i2][i3*W2+i4];
                    }
                }
                F0[i1*W1+i2] = density * H[2] * H[3];
            }
        }      
        corr_0 = 0.0;
        for (int i1 = 2; i1 < BoxShapeFine[0] - 2; i1 ++)  {
            for (int i2 = 2; i2 < BoxShapeFine[1] - 2; i2 ++)  {
                corr_0 += F0[i1*W1+i2] * F0[i1*W1+i2];
            }
        }
        corr_0 = corr_0 * H[0] * H[1];

        log->log("[CaldeiraLeggett4d] corr_0 = %.16e\n",corr_0);
        log->log("[CaldeiraLeggett4d] Time %lf, Corr = %.16e\n",0.0,1.0);
    }

    // .........................................................................................

    // Time iteration 

    log->log("=======================================================\n\n"); 
    log->log("[CaldeiraLeggett4d] Time interation starts ...\n"); 
    log->log("[CaldeiraLeggett4d] Number of steps = %d\n\n", (int)(TIME / kk)); 
    log->log("=======================================================\n\n"); 

    for (unsigned long int tt = 0; tt < (unsigned long int)(TIME / kk); tt ++)
    {
        t_0_begin = omp_get_wtime(); 
            
        if ( tt % PRINT_PERIOD == 0 && isPrintDensity )  {
            pfile = fopen ("density.dat","a");
            fprintf(pfile, "%lu %lf %lu\n", tt, tt * kk, (unsigned long int)(BoxShapeFine[0]*BoxShapeFine[1]));

            for (int i1 = 0; i1 < BoxShapeFine[0]; i1 ++)  {
                for (int i2 = 0; i2 < BoxShapeFine[1]; i2 ++)  {

                    density = 0.0;

                    for (int i3 = 0; i3 < BoxShapeFine[2]; i3 ++)  {
                        for (int i4 = 0; i4 < BoxShapeFine[3]; i4 ++)  {
                                density += PF[i1*W1+i2][i3*W2+i4];
                        }
                    }
                    xx1 = Box[0] + i1 * H[0];
                    xx2 = Box[2] + i2 * H[1];
                    fprintf(pfile, "%.4f %.4f %.16e\n", xx1, xx2, density);
                }
            }
            fclose(pfile);
        }

        if ( true )
        {
            // RK4-1
            #pragma omp parallel
            {
                #pragma omp single nowait
                {
                    t_1_begin = omp_get_wtime();
                }
                #pragma omp for private(xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2) schedule(runtime)
                for (int i1 = 2; i1 < BoxShapeFine[0]-2; i1 ++)  {
                    for (int i2 = 2; i2 < BoxShapeFine[1]-2; i2 ++)  {
                        for (int i3 = 2; i3 < BoxShapeFine[2]-2; i3 ++)  {
                            for (int i4 = 2; i4 < BoxShapeFine[3]-2; i4 ++)  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                xx3 = Box[4] + i3 * H[2];
                                xx4 = Box[6] + i4 * H[3];
                                f0 = F[i1*W1+i2][i3*W2+i4];
                                f1p1 = F[(i1+1)*W1+i2][i3*W2+i4];
                                f1m1 = F[(i1-1)*W1+i2][i3*W2+i4];
                                f2p1 = F[i1*W1+(i2+1)][i3*W2+i4];
                                f2m1 = F[i1*W1+(i2-1)][i3*W2+i4];
                                f3p1 = F[i1*W1+i2][(i3+1)*W2+i4];
                                f3m1 = F[i1*W1+i2][(i3-1)*W2+i4];
                                f4p1 = F[i1*W1+i2][i3*W2+(i4+1)];
                                f4m1 = F[i1*W1+i2][i3*W2+(i4-1)];
                                f1p2 = F[(i1+2)*W1+i2][i3*W2+i4];
                                f1m2 = F[(i1-2)*W1+i2][i3*W2+i4];
                                f2p2 = F[i1*W1+(i2+2)][i3*W2+i4];
                                f2m2 = F[i1*W1+(i2-2)][i3*W2+i4];
                                f3p2 = F[i1*W1+i2][(i3+2)*W2+i4];
                                f3m2 = F[i1*W1+i2][(i3-2)*W2+i4];
                                f4p2 = F[i1*W1+i2][i3*W2+(i4+2)];
                                f4m2 = F[i1*W1+i2][i3*W2+(i4-2)];

                                KK1[i1*W1+i2][i3*W2+i4] = (-kh0m) * xx3 * (-f1p2/12.0 + 2/3.0*f1p1 - 2/3.0*f1m1 + f1m2/12.0) +
                                                        (-kh1m) * xx4 * (-f2p2/12.0 + 2/3.0*f2p1 - 2/3.0*f2m1 + f2m2/12.0) + 
                                                        k2h2 * POTENTIAL_X1(xx1, xx2) * (-f3p2/12.0 + 2/3.0*f3p1 - 2/3.0*f3m1 + f3m2/12.0) +
                                                        k2h3 * POTENTIAL_X2(xx1, xx2) * (-f4p2/12.0 + 2/3.0*f4p1 - 2/3.0*f4m1 + f4m2/12.0) +
                                                        (-khbsq2h2) * POTENTIAL_XXX1(xx1, xx2) * (0.5*f3p2 - f3p1 + f3m1 - 0.5*f3m2) +
                                                        (-khbsq2h3) * POTENTIAL_XXX2(xx1, xx2) * (0.5*f4p2 - f4p1 + f4m1 - 0.5*f4m2) +
                                                        kgamma * (f0 + i2h2 * xx3 * (f3p1 - f3m1) + mkT2h2sq * (f3p1 + f3m1 - 2*f0)) +
                                                        kgamma * (f0 + i2h3 * xx4 * (f4p1 - f4m1) + mkT2h3sq * (f4p1 + f4m1 - 2*f0));

                                FF[i1*W1+i2][i3*W2+i4] = F[i1*W1+i2][i3*W2+i4] + KK1[i1*W1+i2][i3*W2+i4] / 6.0;
                            }
                        }
                    }
                }
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-21: CASE 2 KK1) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-2
                #pragma omp for private(xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk3p1,kk3m1,kk4p1,kk4m1,kk1p2,kk1m2,kk2p2,kk2m2,kk3p2,kk3m2,kk4p2,kk4m2) schedule(runtime)
                for (int i1 = 2; i1 < BoxShapeFine[0]-2; i1 ++)  {
                    for (int i2 = 2; i2 < BoxShapeFine[1]-2; i2 ++)  {
                        for (int i3 = 2; i3 < BoxShapeFine[2]-2; i3 ++)  {
                            for (int i4 = 2; i4 < BoxShapeFine[3]-2; i4 ++)  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                xx3 = Box[4] + i3 * H[2];
                                xx4 = Box[6] + i4 * H[3];
                                f0 = F[i1*W1+i2][i3*W2+i4];
                                f1p1 = F[(i1+1)*W1+i2][i3*W2+i4];
                                f1m1 = F[(i1-1)*W1+i2][i3*W2+i4];
                                f2p1 = F[i1*W1+(i2+1)][i3*W2+i4];
                                f2m1 = F[i1*W1+(i2-1)][i3*W2+i4];
                                f3p1 = F[i1*W1+i2][(i3+1)*W2+i4];
                                f3m1 = F[i1*W1+i2][(i3-1)*W2+i4]; 
                                f4p1 = F[i1*W1+i2][i3*W2+(i4+1)];
                                f4m1 = F[i1*W1+i2][i3*W2+(i4-1)];          
                                f1p2 = F[(i1+2)*W1+i2][i3*W2+i4];
                                f1m2 = F[(i1-2)*W1+i2][i3*W2+i4];
                                f2p2 = F[i1*W1+(i2+2)][i3*W2+i4];
                                f2m2 = F[i1*W1+(i2-2)][i3*W2+i4];
                                f3p2 = F[i1*W1+i2][(i3+2)*W2+i4];
                                f3m2 = F[i1*W1+i2][(i3-2)*W2+i4]; 
                                f4p2 = F[i1*W1+i2][i3*W2+(i4+2)];
                                f4m2 = F[i1*W1+i2][i3*W2+(i4-2)];
                                kk0 = KK1[i1*W1+i2][i3*W2+i4];
                                kk1p1 = KK1[(i1+1)*W1+i2][i3*W2+i4];
                                kk1m1 = KK1[(i1-1)*W1+i2][i3*W2+i4];
                                kk2p1 = KK1[i1*W1+(i2+1)][i3*W2+i4];
                                kk2m1 = KK1[i1*W1+(i2-1)][i3*W2+i4];
                                kk3p1 = KK1[i1*W1+i2][(i3+1)*W2+i4];
                                kk3m1 = KK1[i1*W1+i2][(i3-1)*W2+i4];
                                kk4p1 = KK1[i1*W1+i2][i3*W2+(i4+1)];
                                kk4m1 = KK1[i1*W1+i2][i3*W2+(i4-1)];
                                kk1p2 = KK1[(i1+2)*W1+i2][i3*W2+i4];
                                kk1m2 = KK1[(i1-2)*W1+i2][i3*W2+i4];
                                kk2p2 = KK1[i1*W1+(i2+2)][i3*W2+i4];
                                kk2m2 = KK1[i1*W1+(i2-2)][i3*W2+i4];
                                kk3p2 = KK1[i1*W1+i2][(i3+2)*W2+i4];
                                kk3m2 = KK1[i1*W1+i2][(i3-2)*W2+i4];
                                kk4p2 = KK1[i1*W1+i2][i3*W2+(i4+2)];
                                kk4m2 = KK1[i1*W1+i2][i3*W2+(i4-2)];

                                KK2[i1*W1+i2][i3*W2+i4] = (-kh0m) * xx3 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                                        (-kh1m) * xx4 * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) + 
                                                        k2h2 * POTENTIAL_X1(xx1, xx2) * (-1/12.0*(f3p2+0.5*kk3p2) + 2/3.0*(f3p1+0.5*kk3p1) - 2/3.0*(f3m1+0.5*kk3m1) + 1/12.0*(f3m2+0.5*kk3m2)) +
                                                        k2h3 * POTENTIAL_X2(xx1, xx2) * (-1/12.0*(f4p2+0.5*kk4p2) + 2/3.0*(f4p1+0.5*kk4p1) - 2/3.0*(f4m1+0.5*kk4m1) + 1/12.0*(f4m2+0.5*kk4m2)) +
                                                        (-khbsq2h2) * POTENTIAL_XXX1(xx1, xx2) * (0.5*(f3p2+0.5*kk3p2) - (f3p1+0.5*kk3p1) + (f3m1+0.5*kk3m1) - 0.5*(f3m2+0.5*kk3m2)) +
                                                        (-khbsq2h3) * POTENTIAL_XXX2(xx1, xx2) * (0.5*(f4p2+0.5*kk4p2) - (f4p1+0.5*kk4p1) + (f4m1+0.5*kk4m1) - 0.5*(f4m2+0.5*kk4m2)) +
                                                        kgamma * (f0 + 0.5*kk0 + i2h2 * xx3 * (f3p1 + 0.5*kk3p1 - f3m1 - 0.5*kk3m1) + mkT2h2sq * (f3p1 + 0.5*kk3p1 + f3m1 + 0.5*kk3m1 - 2*f0 - kk0)) +
                                                        kgamma * (f0 + 0.5*kk0 + i2h3 * xx4 * (f4p1 + 0.5*kk4p1 - f4m1 - 0.5*kk4m1) + mkT2h3sq * (f4p1 + 0.5*kk4p1 + f4m1 + 0.5*kk4m1 - 2*f0 - kk0));

                                FF[i1*W1+i2][i3*W2+i4] += KK2[i1*W1+i2][i3*W2+i4] / 3.0;
                            }
                        }
                    }
                }
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-22: CASE 2 KK2) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-3
                #pragma omp for private(xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk3p1,kk3m1,kk4p1,kk4m1,kk1p2,kk1m2,kk2p2,kk2m2,kk3p2,kk3m2,kk4p2,kk4m2) schedule(runtime)
                for (int i1 = 2; i1 < BoxShapeFine[0]-2; i1 ++)  {
                    for (int i2 = 2; i2 < BoxShapeFine[1]-2; i2 ++)  {
                        for (int i3 = 2; i3 < BoxShapeFine[2]-2; i3 ++)  {
                            for (int i4 = 2; i4 < BoxShapeFine[3]-2; i4 ++)  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                xx3 = Box[4] + i3 * H[2];
                                xx4 = Box[6] + i4 * H[3];
                                f0 = F[i1*W1+i2][i3*W2+i4];
                                f1p1 = F[(i1+1)*W1+i2][i3*W2+i4];
                                f1m1 = F[(i1-1)*W1+i2][i3*W2+i4];
                                f2p1 = F[i1*W1+(i2+1)][i3*W2+i4];
                                f2m1 = F[i1*W1+(i2-1)][i3*W2+i4];
                                f3p1 = F[i1*W1+i2][(i3+1)*W2+i4];
                                f3m1 = F[i1*W1+i2][(i3-1)*W2+i4]; 
                                f4p1 = F[i1*W1+i2][i3*W2+(i4+1)];
                                f4m1 = F[i1*W1+i2][i3*W2+(i4-1)];          
                                f1p2 = F[(i1+2)*W1+i2][i3*W2+i4];
                                f1m2 = F[(i1-2)*W1+i2][i3*W2+i4];
                                f2p2 = F[i1*W1+(i2+2)][i3*W2+i4];
                                f2m2 = F[i1*W1+(i2-2)][i3*W2+i4];
                                f3p2 = F[i1*W1+i2][(i3+2)*W2+i4];
                                f3m2 = F[i1*W1+i2][(i3-2)*W2+i4]; 
                                f4p2 = F[i1*W1+i2][i3*W2+(i4+2)];
                                f4m2 = F[i1*W1+i2][i3*W2+(i4-2)];
                                kk0 = KK2[i1*W1+i2][i3*W2+i4];
                                kk1p1 = KK2[(i1+1)*W1+i2][i3*W2+i4];
                                kk1m1 = KK2[(i1-1)*W1+i2][i3*W2+i4];
                                kk2p1 = KK2[i1*W1+(i2+1)][i3*W2+i4];
                                kk2m1 = KK2[i1*W1+(i2-1)][i3*W2+i4];
                                kk3p1 = KK2[i1*W1+i2][(i3+1)*W2+i4];
                                kk3m1 = KK2[i1*W1+i2][(i3-1)*W2+i4];
                                kk4p1 = KK2[i1*W1+i2][i3*W2+(i4+1)];
                                kk4m1 = KK2[i1*W1+i2][i3*W2+(i4-1)];
                                kk1p2 = KK2[(i1+2)*W1+i2][i3*W2+i4];
                                kk1m2 = KK2[(i1-2)*W1+i2][i3*W2+i4];
                                kk2p2 = KK2[i1*W1+(i2+2)][i3*W2+i4];
                                kk2m2 = KK2[i1*W1+(i2-2)][i3*W2+i4];
                                kk3p2 = KK2[i1*W1+i2][(i3+2)*W2+i4];
                                kk3m2 = KK2[i1*W1+i2][(i3-2)*W2+i4];
                                kk4p2 = KK2[i1*W1+i2][i3*W2+(i4+2)];
                                kk4m2 = KK2[i1*W1+i2][i3*W2+(i4-2)];

                                KK3[i1*W1+i2][i3*W2+i4] = (-kh0m) * xx3 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                                        (-kh1m) * xx4 * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) + 
                                                        k2h2 * POTENTIAL_X1(xx1, xx2) * (-1/12.0*(f3p2+0.5*kk3p2) + 2/3.0*(f3p1+0.5*kk3p1) - 2/3.0*(f3m1+0.5*kk3m1) + 1/12.0*(f3m2+0.5*kk3m2)) +
                                                        k2h3 * POTENTIAL_X2(xx1, xx2) * (-1/12.0*(f4p2+0.5*kk4p2) + 2/3.0*(f4p1+0.5*kk4p1) - 2/3.0*(f4m1+0.5*kk4m1) + 1/12.0*(f4m2+0.5*kk4m2)) +
                                                        (-khbsq2h2) * POTENTIAL_XXX1(xx1, xx2) * (0.5*(f3p2+0.5*kk3p2) - (f3p1+0.5*kk3p1) + (f3m1+0.5*kk3m1) - 0.5*(f3m2+0.5*kk3m2)) +
                                                        (-khbsq2h3) * POTENTIAL_XXX2(xx1, xx2) * (0.5*(f4p2+0.5*kk4p2) - (f4p1+0.5*kk4p1) + (f4m1+0.5*kk4m1) - 0.5*(f4m2+0.5*kk4m2)) +
                                                        kgamma * (f0 + 0.5*kk0 + i2h2 * xx3 * (f3p1 + 0.5*kk3p1 - f3m1 - 0.5*kk3m1) + mkT2h2sq * (f3p1 + 0.5*kk3p1 + f3m1 + 0.5*kk3m1 - 2*f0 - kk0)) +
                                                        kgamma * (f0 + 0.5*kk0 + i2h3 * xx4 * (f4p1 + 0.5*kk4p1 - f4m1 - 0.5*kk4m1) + mkT2h3sq * (f4p1 + 0.5*kk4p1 + f4m1 + 0.5*kk4m1 - 2*f0 - kk0));

                                FF[i1*W1+i2][i3*W2+i4] += KK3[i1*W1+i2][i3*W2+i4] / 3.0;
                            }
                        }
                    }
                }
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-23: CASE 2 KK3) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-4
                #pragma omp for private(xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk3p1,kk3m1,kk4p1,kk4m1,kk1p2,kk1m2,kk2p2,kk2m2,kk3p2,kk3m2,kk4p2,kk4m2) schedule(runtime)
                for (int i1 = 2; i1 < BoxShapeFine[0]-2; i1 ++)  {
                    for (int i2 = 2; i2 < BoxShapeFine[1]-2; i2 ++)  {
                        for (int i3 = 2; i3 < BoxShapeFine[2]-2; i3 ++)  {
                            for (int i4 = 2; i4 < BoxShapeFine[3]-2; i4 ++)  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                xx3 = Box[4] + i3 * H[2];
                                xx4 = Box[6] + i4 * H[3];
                                f0 = F[i1*W1+i2][i3*W2+i4];
                                f1p1 = F[(i1+1)*W1+i2][i3*W2+i4];
                                f1m1 = F[(i1-1)*W1+i2][i3*W2+i4];
                                f2p1 = F[i1*W1+(i2+1)][i3*W2+i4];
                                f2m1 = F[i1*W1+(i2-1)][i3*W2+i4];
                                f3p1 = F[i1*W1+i2][(i3+1)*W2+i4];
                                f3m1 = F[i1*W1+i2][(i3-1)*W2+i4]; 
                                f4p1 = F[i1*W1+i2][i3*W2+(i4+1)];
                                f4m1 = F[i1*W1+i2][i3*W2+(i4-1)];          
                                f1p2 = F[(i1+2)*W1+i2][i3*W2+i4];
                                f1m2 = F[(i1-2)*W1+i2][i3*W2+i4];
                                f2p2 = F[i1*W1+(i2+2)][i3*W2+i4];
                                f2m2 = F[i1*W1+(i2-2)][i3*W2+i4];
                                f3p2 = F[i1*W1+i2][(i3+2)*W2+i4];
                                f3m2 = F[i1*W1+i2][(i3-2)*W2+i4]; 
                                f4p2 = F[i1*W1+i2][i3*W2+(i4+2)];
                                f4m2 = F[i1*W1+i2][i3*W2+(i4-2)];
                                kk0 = KK3[i1*W1+i2][i3*W2+i4];
                                kk1p1 = KK3[(i1+1)*W1+i2][i3*W2+i4];
                                kk1m1 = KK3[(i1-1)*W1+i2][i3*W2+i4];
                                kk2p1 = KK3[i1*W1+(i2+1)][i3*W2+i4];
                                kk2m1 = KK3[i1*W1+(i2-1)][i3*W2+i4];
                                kk3p1 = KK3[i1*W1+i2][(i3+1)*W2+i4];
                                kk3m1 = KK3[i1*W1+i2][(i3-1)*W2+i4];
                                kk4p1 = KK3[i1*W1+i2][i3*W2+(i4+1)];
                                kk4m1 = KK3[i1*W1+i2][i3*W2+(i4-1)];
                                kk1p2 = KK3[(i1+2)*W1+i2][i3*W2+i4];
                                kk1m2 = KK3[(i1-2)*W1+i2][i3*W2+i4];
                                kk2p2 = KK3[i1*W1+(i2+2)][i3*W2+i4];
                                kk2m2 = KK3[i1*W1+(i2-2)][i3*W2+i4];
                                kk3p2 = KK3[i1*W1+i2][(i3+2)*W2+i4];
                                kk3m2 = KK3[i1*W1+i2][(i3-2)*W2+i4];
                                kk4p2 = KK3[i1*W1+i2][i3*W2+(i4+2)];
                                kk4m2 = KK3[i1*W1+i2][i3*W2+(i4-2)];

                                KK4[i1*W1+i2][i3*W2+i4] = (-kh0m) * xx3 * (-1/12.0*(f1p2+kk1p2) + 2/3.0*(f1p1+kk1p1) - 2/3.0*(f1m1+kk1m1) + 1/12.0*(f1m2+kk1m2)) +
                                                        (-kh1m) * xx4 * (-1/12.0*(f2p2+kk2p2) + 2/3.0*(f2p1+kk2p1) - 2/3.0*(f2m1+kk2m1) + 1/12.0*(f2m2+kk2m2)) + 
                                                        k2h2 * POTENTIAL_X1(xx1, xx2) * (-1/12.0*(f3p2+kk3p2) + 2/3.0*(f3p1+kk3p1) - 2/3.0*(f3m1+kk3m1) + 1/12.0*(f3m2+kk3m2)) +
                                                        k2h3 * POTENTIAL_X2(xx1, xx2) * (-1/12.0*(f4p2+kk4p2) + 2/3.0*(f4p1+kk4p1) - 2/3.0*(f4m1+kk4m1) + 1/12.0*(f4m2+kk4m2)) +
                                                        (-khbsq2h2) * POTENTIAL_XXX1(xx1, xx2) * (0.5*(f3p2+kk3p2) - (f3p1+kk3p1) + (f3m1+kk3m1) - 0.5*(f3m2+kk3m2)) + 
                                                        (-khbsq2h3) * POTENTIAL_XXX2(xx1, xx2) * (0.5*(f4p2+kk4p2) - (f4p1+kk4p1) + (f4m1+kk4m1) - 0.5*(f4m2+kk4m2)) + 
                                                        kgamma * (f0 + kk0 + i2h2 * xx3 * (f3p1 + kk3p1 - f3m1 - kk3m1) + mkT2h2sq * (f3p1 + kk3p1 + f3m1 + kk3m1 - 2*f0 - 2*kk0)) +
                                                        kgamma * (f0 + kk0 + i2h3 * xx4 * (f4p1 + kk4p1 - f4m1 - kk4m1) + mkT2h3sq * (f4p1 + kk4p1 + f4m1 + kk4m1 - 2*f0 - 2*kk0));

                                FF[i1*W1+i2][i3*W2+i4] += KK4[i1*W1+i2][i3*W2+i4] / 6.0;
                            }
                        }
                    }                
                }
                
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-24: CASE 2 KK4) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }
            } // OMP PARALLEL
        } 
        // .........................................................................................

        // FF(t+1) Normailzed & go on

        t_1_begin = omp_get_wtime();
        norm = 0.0;

        #pragma omp parallel for reduction (+:norm)
        for (int i1 = 0; i1 < BoxShapeFine[0]; i1 ++)  {
            for (int i2 = 0; i2 < BoxShapeFine[1]; i2 ++)  {
                for (int i3 = 0; i3 < BoxShapeFine[2]; i3 ++)  {
                    for (int i4 = 0; i4 < BoxShapeFine[3]; i4 ++)  {
                        norm += FF[i1*W1+i2][i3*W2+i4];
                    }
                }
            }
        }
        norm *= H[0] * H[1] * H[2] * H[3];

        if ( (tt + 1) % PERIOD == 0 )
            log->log("[CaldeiraLeggett4d] Normalization factor = %.16e\n",norm);

        norm = 1.0 / norm; 

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-1-1 Norm) = %lf sec\n", t_1_elapsed);
        t_1_begin = omp_get_wtime();

        #pragma omp parallel for private(val)
        for (int i1 = 0; i1 < BoxShapeFine[0]; i1 ++)  {
            for (int i2 = 0; i2 < BoxShapeFine[1]; i2 ++)  {
                for (int i3 = 0; i3 < BoxShapeFine[2]; i3 ++)  {
                    for (int i4 = 0; i4 < BoxShapeFine[3]; i4 ++)  {
                        val = norm * FF[i1*W1+i2][i3*W2+i4];
                        FF[i1*W1+i2][i3*W2+i4] = val;
                        F[i1*W1+i2][i3*W2+i4] = val;
                        PF[i1*W1+i2][i3*W2+i4] = val;
                    }
                }
            }
        }

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-1-2 FF) = %lf sec\n", t_1_elapsed); 

        if ( (tt + 1) % PERIOD == 0 )
        {
            // REPORT MEASUREMENTS
            // ----------------------------------------------------------------------------
            // Compute Transmittance
            // ----------------------------------------------------------------------------

            if (isTrans)  {
            
                t_1_begin = omp_get_wtime();

                // position for x >= trans_x0
                idx_x0 = (int)std::round(std::round((trans_x0-Box[0])/H[0]));
                pftrans = 0.0;

                #pragma omp parallel for reduction (+:pftrans)
                for (int i1 = idx_x0; i1 < BoxShapeFine[0]; i1 ++)  {
                //for (int i1 = idx_x0 + 1; i1 < BoxShapeFine[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShapeFine[1]; i2 ++)  {
                        for (int i3 = 0; i3 < BoxShapeFine[2]; i3 ++)  {
                            for (int i4 = 0; i4 < BoxShapeFine[3]; i4 ++)  {
                                pftrans+=PF[i1*W1+i2][i3*W2+i4];
                            }
                        }
                    }
                }
                pftrans *= H[0] * H[1] * H[2] * H[3];
                PF_trans.push_back(pftrans);
                log->log("[CaldeiraLeggett4d] Time %lf, Trans = %.16e\n", ( tt + 1 ) * kk, pftrans);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (omp-x-2 trans) = %lf sec\n", t_1_elapsed); 
            }

            if (isCorr)  {

                // Compute Density at time t
                #pragma omp parallel for private(density)
                for (int i1 = 2; i1 < BoxShapeFine[0] - 2; i1 ++)  {
                    for (int i2 = 2; i2 < BoxShapeFine[1] - 2; i2 ++)  {
                        density = 0.0;
                        for (int i3 = 2; i3 < BoxShapeFine[2]-2; i3 ++)  {
                            for (int i4 = 2; i4 < BoxShapeFine[3]-2; i4 ++)  {      
                                density += PF[i1*W1+i2][i3*W2+i4];
                            }
                        }
                        Ft[i1*W1+i2] = density * H[2] * H[3];
                    }
                }      

                // Compute correlation
                corr = 0.0;
                for (int i1 = 2; i1 < BoxShapeFine[0] - 2; i1 ++)  {
                    for (int i2 = 2; i2 < BoxShapeFine[1] - 2; i2 ++)  {
                        corr += F0[i1*W1+i2] * Ft[i1*W1+i2];
                    }
                }
                corr *= H[0] * H[1];
                log->log("[CaldeiraLeggett4d] Time %lf, Corr = %.16e\n", ( tt + 1 ) * kk, corr/corr_0);
            }
        }

        // Reset

        t_1_begin = omp_get_wtime();

        #pragma omp parallel for
        for (int i1 = 0; i1 < BoxShapeFine[0]; i1++)  {
            for (int i2 = 0; i2 < BoxShapeFine[1]; i2++)  {
                for (int i3 = 0; i3 < BoxShapeFine[2]; i3++)  {
                    for (int i4 = 0; i4 < BoxShapeFine[3]; i4++)  {
                        FF[i1*W1+i2][i3*W2+i4] = 0.0;
                        KK1[i1*W1+i2][i3*W2+i4] = 0.0;
                        KK2[i1*W1+i2][i3*W2+i4] = 0.0;
                        KK3[i1*W1+i2][i3*W2+i4] = 0.0;
                        KK4[i1*W1+i2][i3*W2+i4] = 0.0;
                    }
                }
            }
        }

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if ( !QUIET && TIMING ) log->log("Elapsed time (omp-e-8: reset) = %lf sec\n", t_1_elapsed);  

        if ( (tt + 1) % PERIOD == 0 )
        {   
            t_0_end = omp_get_wtime();
            t_0_elapsed = t_0_end - t_0_begin;

            if ( !QUIET )  {

                log->log("[CaldeiraLeggett4d] Step: %lu, Elapsed time: %lf sec\n", tt + 1, t_0_elapsed);
                log->log("[CaldeiraLeggett4d] BoxPosition: [%.2f,%.2f][%.2f,%.2f][%.2f,%.2f][%.2f,%.2f]\n",Box[0],Box[1],Box[2],Box[3],Box[4],Box[5],Box[6],Box[7]);
                log->log("[CaldeiraLeggett4d] BoxShape: [%d,%d,%d,%d]\n",BoxShapeFine[0],BoxShapeFine[1],BoxShapeFine[2],BoxShapeFine[3]);
                log->log("\n........................................................\n\n");
            }
        }         
    } // Time iteration 

    for (unsigned long int i = 0; i < O1; i++)  {
        delete F[i];
        delete FF[i];
        delete PF[i];
        delete KK1[i];
        delete KK2[i];
        delete KK3[i];
        delete KK4[i];
    }
    delete F;
    delete FF;
    delete PF;
    delete KK1;
    delete KK2;
    delete KK3;
    delete KK4;

    log->log("[CaldeiraLeggett4d] Evolve done.\n");
}
/* =============================================================================== */

/* CL4DPOT_DW1 */

inline double CaldeiraLeggett4d::Wavefunction_DW1(double x1, double x2, double x3, double x4)
{
    return PIHBSQ_INV * exp(-2.0 * (A[0]*pow(x1-Wave0[0],2)+A[1]*pow(x2-Wave0[1],2))) * exp( -0.5 * HBSQ_INV * (pow(x3-Wave0[2],2)/A[2]+pow(x4-Wave0[3],2)/A[3]) );
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Potential_DW1(double x1, double x2)
{
    return 0.025 * (x1 * x1 - 1.0) * (x1 * x1 - 1.0) + 0.5 * k0 * (1.0 - sigma * exp(- lambda * x1 * x1)) * x2 * x2;    
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vx1_DW1(double x1, double x2)
{
    return 0.1 * x1 * (x1 * x1 - 1.0) + k0 * sigma * lambda * x1 * x2 * x2 * exp(-lambda * x1 * x1);
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vx2_DW1(double x1, double x2)
{
    return k0 * x2 * (1.0 - sigma * exp(-lambda * x1 * x1));
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vxxx1_DW1(double x1, double x2)
{
    return 0.6 * x1 + 2.0 * k0 * sigma * lambda * lambda * x1 * x2 * x2 * (2.0 * lambda * x1 * x1 - 3.0) * exp(-lambda * x1 * x1);
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vxxx2_DW1(double x1, double x2)
{
    return 0.0;
}
/* =============================================================================== */

/* CL4DPOT_DW2 */

inline double CaldeiraLeggett4d::Wavefunction_DW2(double x1, double x2, double x3, double x4)
{
    return PIHBSQ_INV * exp(-2.0 * (A[0]*pow(x1-Wave0[0],2)+A[1]*pow(x2-Wave0[1],2))) * exp( -0.5 * HBSQ_INV * (pow(x3-Wave0[2],2)/A[2]+pow(x4-Wave0[3],2)/A[3]) );
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Potential_DW2(double x1, double x2)
{
    return 0.007 * x1 * x1 * x1 * x1 - 0.01 * x1 * x1 + 0.5 * k0 * (1.0 - sigma * exp(- lambda * x1 * x1)) * x2 * x2;    
    //return 0.007 * (x1 * x1 * x1 * x1) - 0.01 * x1 * x1 + 0.5 * k0 * x2 * x2;
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vx1_DW2(double x1, double x2)
{
    return 0.028 * x1 * x1 * x1 - 0.02 * x1 + k0 * sigma * lambda * x1 * x2 * x2 * exp(-lambda * x1 * x1);
    //return 0.028 * x1 * x1 * x1 - 0.02 * x1;
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vx2_DW2(double x1, double x2)
{
    return k0 * x2 * (1.0 - sigma * exp(-lambda * x1 * x1));
    //return k0 * x2;
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vxxx1_DW2(double x1, double x2)
{
    return 0.168 * x1 + 2.0 * k0 * sigma * lambda * lambda * x1 * x2 * x2 * (2.0 * lambda * x1 * x1 - 3.0) * exp(-lambda * x1 * x1);
    //return 0.168 * x1;
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vxxx2_DW2(double x1, double x2)
{
    return 0.0;
}
/* =============================================================================== */

/* CL4DPOT_DW3 */

inline double CaldeiraLeggett4d::Wavefunction_DW3(double x1, double x2, double x3, double x4)
{
    return PIHBSQ_INV * exp(-2.0 * (A[0]*pow(x1-Wave0[0],2)+A[1]*pow(x2-Wave0[1],2))) * exp( -0.5 * HBSQ_INV * (pow(x3-Wave0[2],2)/A[2]+pow(x4-Wave0[3],2)/A[3]) );
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Potential_DW3(double x1, double x2)
{
    return 3.0 * x1 * x1 * x1 * x1 - 8.0 * x1 * x1 + x2 * x2 - 0.1 * x1 * x2;    
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vx1_DW3(double x1, double x2)
{
    return 12.0 * x1 * x1 * x1 - 16.0 * x1 - 0.1 * x2;
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vx2_DW3(double x1, double x2)
{
    return 2.0 * x2 - 0.1 * x1;
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vxxx1_DW3(double x1, double x2)
{
    return 72.0 * x1;
}
/* ------------------------------------------------------------------------------- */

inline double CaldeiraLeggett4d::Vxxx2_DW3(double x1, double x2)
{
    return 0.0;
}
/* ------------------------------------------------------------------------------- */

VectorXi CaldeiraLeggett4d::IdxToGrid(unsigned long int idx)
{
    int x4 = (int)( idx % M3 );
    int x3 = (int)(( idx % M2 ) / M3);
    int x2 = (int)(( idx % M1 ) / M2);
    int x1 = (int)( idx / M1 );

    VectorXi grid;
    grid.resize(DIMENSIONS);
    grid << x1, x2, x3, x4;

    return grid;
}
/* ------------------------------------------------------------------------------- */

inline unsigned long int CaldeiraLeggett4d::GridToIdx(int x1, int x2, int x3, int x4)
{
    return (unsigned long int)(x1 * M1 + x2 * M2 + x3 * M3 + x4);
}
/* ------------------------------------------------------------------------------- */
