// ==============================================================================
//
//  CaldeiraLeggett4dTr.cpp
//  QTR
//
//  Created by Albert Lu on 1/30/20.
//  alu@tacc.utexas.edu
//
//  Last modified on 3/4/20
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

#if defined CL4DPOT_DW1

#define WAVEFUNCTION(x1,x2,x3,x4) Wavefunction_DW1(x1,x2,x3,x4)
#define POTENTIAL(x1,x2) Potential_DW1(x1,x2)
#define POTENTIAL_X1(x1,x2) Vx1_DW1(x1,x2)
#define POTENTIAL_X2(x1,x2) Vx2_DW1(x1,x2)
#define POTENTIAL_XXX1(x1, x2) Vxxx1_DW1(x1,x2)
#define POTENTIAL_XXX2(x1, x2) Vxxx2_DW1(x1,x2)
#define POTNAME "DoubleWell-1"

#elif defined CL4DPOT_DW2

#define WAVEFUNCTION(x1,x2,x3,x4) Wavefunction_DW2(x1,x2,x3,x4)
#define POTENTIAL(x1,x2) Potential_DW2(x1,x2)
#define POTENTIAL_X1(x1,x2) Vx1_DW2(x1,x2)
#define POTENTIAL_X2(x1,x2) Vx2_DW2(x1,x2)
#define POTENTIAL_XXX1(x1, x2) Vxxx1_DW2(x1,x2)
#define POTENTIAL_XXX2(x1, x2) Vxxx2_DW2(x1,x2)
#define POTNAME "DoubleWell-2"

#elif defined CL4DPOT_DW3

#define WAVEFUNCTION(x1,x2,x3,x4) Wavefunction_DW3(x1,x2,x3,x4)
#define POTENTIAL(x1,x2) Potential_DW3(x1,x2)
#define POTENTIAL_X1(x1,x2) Vx1_DW3(x1,x2)
#define POTENTIAL_X2(x1,x2) Vx2_DW3(x1,x2)
#define POTENTIAL_XXX1(x1, x2) Vxxx1_DW3(x1,x2)
#define POTENTIAL_XXX2(x1, x2) Vxxx2_DW3(x1,x2)
#define POTNAME "DoubleWell-3"

#else // Default CL4DPOT_DW3

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
    log->log("\n\n[CaldeiraLeggett4d] Using Truncated Grid\n");
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

    // Grid size
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
    // Coarse grid
    CFactor = parameters->scxd_cfactor;
    C.resize(DIMENSIONS);
    C[0] = H[0] * CFactor; 
    C[1] = H[1] * CFactor; 
    C[2] = H[2] * CFactor; 
    C[3] = H[3] * CFactor; 

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

    BoxBoundary.resize(DIMENSIONS * 2);
    BoxBoundary[0] = parameters->scxd_bi1;
    BoxBoundary[1] = parameters->scxd_bf1;
    BoxBoundary[2] = parameters->scxd_bi2;
    BoxBoundary[3] = parameters->scxd_bf2;
    BoxBoundary[4] = parameters->scxd_bi3;
    BoxBoundary[5] = parameters->scxd_bf3;
    BoxBoundary[6] = parameters->scxd_bi4;
    BoxBoundary[7] = parameters->scxd_bf4;

    BoxDynamic.resize(DIMENSIONS * 2);
    BoxDynamicNew.resize(DIMENSIONS * 2);
    BoxShapeFine.resize(DIMENSIONS);
    BoxShapeFineNew.resize(DIMENSIONS);
    BoxShapeCoarse.resize(DIMENSIONS);
    DistToEdge.resize(DIMENSIONS * 2);

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

    // Truncate parameters
    TolH = parameters->scxd_TolH;    // Tolerance of probability density for Zero point Cutoff
    TolL = parameters->scxd_TolL;    // Tolerance of probability density for Edge point
    TolHd = parameters->scxd_TolHd;  // Tolerance of probability first diff for Zero point Cutoff
    TolLd = parameters->scxd_TolLd;  // Tolerance of probability density for Edge point
    ExReduce = parameters->scxd_ExReduce; //Extrapolation reduce factor
    ExLimit = parameters->scxd_ExLimit;   //Extrapolation counts limit

    // Domain buffer
    skin = std::max(parameters->scxd_skin, std::max(3,DIMENSIONS+ExLimit));
    log->log("[CaldeiraLeggett4d] Set skin thickness to %d.\n\n",skin);

    log->log("[CaldeiraLeggett4d] INIT done.\n\n");
}
/* ------------------------------------------------------------------------------- */

void CaldeiraLeggett4d::Evolve()
{
    #pragma omp declare reduction (merge : MeshIndexLU, MeshIndexFull, MeshValueD : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

    log->log("[CaldeiraLeggett4d] Evolve starts ...\n");

    // Files
    FILE *pfile;

    // Variables 
    int count;
    int n1, n2, n3, n4;
    unsigned long int index;
    unsigned long int ta_size, tb_size;
    double coeff1 = 0.0;
    double coeff2 = 0.0;
    double norm;
    double sum;
    double density;
    bool b1, b2, b3, b4, b5, b6, b7, b8;

    // Adaptive Box
    int x1_max, x2_max, x3_max, x4_max;
    int x1_min, x2_min, x3_min, x4_min;
    bool isEnterSkin;
    bool isTouchedBoundary = false;
    bool isAdaptiveBoxChanged;

    // Correlation
    double corr;
    double corr_0;
    double x1_min_0, x2_min_0;
    double x1_max_0, x2_max_0;
    double x1_min_t, x2_min_t;
    double x1_max_t, x2_max_t;
    double ibox_x1_min, ibox_x1_max;
    double ibox_x2_min, ibox_x2_max;
    int ibox_shape_1, ibox_shape_2;
    int shift_x1_0, shift_x1_t;
    int shift_x2_0, shift_x2_t;

    // Timing variables
    double t_0_begin, t_0_end;
    double t_1_begin, t_1_end;
    double t_0_elapsed = 0.0;
    double t_1_elapsed = 0.0;

    // Constants
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
    double TolHdX1 = 2.0 * TolHd * H[0];
    double TolHdX2 = 2.0 * TolHd * H[1];
    double TolHdX3 = 2.0 * TolHd * H[2];
    double TolHdX4 = 2.0 * TolHd * H[3];
    double TolLdX1 = 2.0 * TolLd * H[0];
    double TolLdX2 = 2.0 * TolLd * H[1];
    double TolLdX3 = 2.0 * TolLd * H[2];
    double TolLdX4 = 2.0 * TolLd * H[3];

    // temporary index container
    MeshIndexLU tmpVec;
    MeshIndexFull tmpVecFull;
    MeshValueD tmpVecD;

    // Boundary layer container for extrapolation loop
    MeshIndexLU ExBD;    
 
    // 2d Grid vector and indices
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

    // Extrapolation 
    int min_dir;
    bool isFirstExtrp;
    double val, val_min;
    double val_min_abs;
    vector<double> ExTBL;
    int Excount;

    // PF_trans
    double pftrans;
    vector<double> PF_trans;
    PF_trans.push_back(0.0);

    // Neighborlist
    int nneigh = 0;
    vector<vector<int>> neighlist;
    vector<int> neighs(DIMENSIONS);

    for (int d = 1; d <= 4; d ++)  {
        for (n1 = -d; n1 <= d; n1 ++)  {
            for (n2 = -(d-abs(n1)); n2 <= d-abs(n1); n2 ++)  {
                for (n3 = -(d-abs(n1)-abs(n2)); n3 <= d-abs(n1)-abs(n2); n3 ++)  {

                    n4 = d - abs(n1) - abs(n2) - abs(n3);

                    if (n4 != 0)  {
                        neighs = {n1,n2,n3,n4};
                        neighlist.push_back(neighs);
                        neighs = {n1,n2,n3,-n4};
                        neighlist.push_back(neighs);
                        nneigh += 2;
                    }
                    else  {
                        neighs = {n1,n2,n3,0};
                        neighlist.push_back(neighs);
                        nneigh += 1;
                    }
                }
            }
        }
    }
    log->log("[CaldeiraLeggett4d] nneigh = %d\n", nneigh); 

    // Compute coarse grid shape
    t_1_begin = omp_get_wtime();
    log->log("[CaldeiraLeggett4d] Coarse grid size = (");

    GRIDS_TOT = 1;
    for (int i = 0; i < DIMENSIONS; i ++)  {
        BoxShapeCoarse[i] = (int)( (Box[2 * i + 1] - Box[2 * i]) / C[i] ) + 1;
        GRIDS_TOT *= BoxShapeCoarse[i];

        if ( i < DIMENSIONS - 1 )
            log->log("%d, ", BoxShapeCoarse[i]);
        else
            log->log("%d)=", BoxShapeCoarse[i]);
    }
    log->log("%lu\n",GRIDS_TOT);

    // Use coarse grid to compute truncated boundary
    double vTolH = C[0] * C[1] * C[2] * C[3] * TolH;

    x1_min = BoxShapeCoarse[0];
    x2_min = BoxShapeCoarse[1];
    x3_min = BoxShapeCoarse[2];
    x4_min = BoxShapeCoarse[3];
    x1_max = 0;
    x2_max = 0;
    x3_max = 0;
    x4_max = 0;

    #pragma omp parallel for reduction(min: x1_min,x2_min,x3_min,x4_min) reduction(max: x1_max,x2_max,x3_max,x4_max)
    for (int i1 = 0; i1 < BoxShapeCoarse[0]; i1 ++)  {
        for (int i2 = 0; i2 < BoxShapeCoarse[1]; i2 ++)  {
            for (int i3 = 0; i3 < BoxShapeCoarse[2]; i3 ++)  {
                for (int i4 = 0; i4 < BoxShapeCoarse[3]; i4 ++)  {
                    if ( std::fabs(WAVEFUNCTION(Box[0]+i1*C[0], Box[2]+i2*C[1], Box[4]+i3*C[2], Box[6]+i4*C[3])) >= vTolH )  {
                        if ( i1 < x1_min ) x1_min = i1;
                        if ( i1 > x1_max ) x1_max = i1;
                        if ( i2 < x2_min ) x2_min = i2;
                        if ( i2 > x2_max ) x2_max = i2;
                        if ( i3 < x3_min ) x3_min = i3;
                        if ( i3 > x3_max ) x3_max = i3;
                        if ( i4 < x4_min ) x4_min = i4;
                        if ( i4 > x4_max ) x4_max = i4;
                    }
                }
            }
        }
    }
    log->log("[CaldeiraLeggett4d] Coarse TA Range [%d, %d][%d, %d][%d, %d][%d, %d]\n", x1_min, x1_max, x2_min, x2_max, x3_min, x3_max, x4_min, x4_max);

    if ( x1_min==0 || x2_min==0 || x3_min==0 || x4_min==0 || x1_max==BoxShapeCoarse[0]-1 || x2_max==BoxShapeCoarse[1]-1 || x3_max==BoxShapeCoarse[2]-1 || x4_max==BoxShapeCoarse[3]-1 )  {
        log->log("[CaldeiraLeggett4d] Coarse TA range touches boundary");
        exit(1);
    }

    BoxDynamic[0] = Box[0] + (x1_min-1) * C[0] - 2 * skin * H[0];
    BoxDynamic[1] = Box[0] + (x1_max+1) * C[0] + 2 * skin * H[0];
    BoxDynamic[2] = Box[2] + (x2_min-1) * C[1] - 2 * skin * H[1];
    BoxDynamic[3] = Box[2] + (x2_max+1) * C[1] + 2 * skin * H[1];
    BoxDynamic[4] = Box[4] + (x3_min-1) * C[2] - 2 * skin * H[2];
    BoxDynamic[5] = Box[4] + (x3_max+1) * C[2] + 2 * skin * H[2];
    BoxDynamic[6] = Box[6] + (x4_min-1) * C[3] - 2 * skin * H[3];
    BoxDynamic[7] = Box[6] + (x4_max+1) * C[3] + 2 * skin * H[3];

    // Compute fine grid shape
    log->log("[CaldeiraLeggett4d] Initital fine grid size = (");

    GRIDS_TOT = 1;
    for (int i = 0; i < DIMENSIONS; i ++)  {
        BoxShapeFine[i] = (int)std::round((BoxDynamic[2*i+1]-BoxDynamic[2*i])/H[i]) + 1;
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
    W1_0 = W1;

    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    if ( !QUIET && TIMING ) log->log("[CaldeiraLeggett4d] Elapsed time (make bounding box) = %lf sec\n\n", t_1_elapsed); 

    // Initialize containers
    t_0_begin = omp_get_wtime();

    log->log("[CaldeiraLeggett4d] Initializing containers ...\n");

    bool **TAMask = new bool*[O1];

    for (unsigned long int  i = 0; i < O1; i++)  {
        TAMask[i] = new(nothrow) bool[O2];

        if (!TAMask[i])
            log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "TAMask");
    }

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
                    TAMask[i1*W1+i2][i3*W2+i4] = 0;
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
                    F[i1*W1+i2][i3*W2+i4] = WAVEFUNCTION(BoxDynamic[0]+i1*H[0], BoxDynamic[2]+i2*H[1], BoxDynamic[4]+i3*H[2], BoxDynamic[6]+i4*H[3]);
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
    if (!QUIET && TIMING) log->log("[CaldeiraLeggett4d] Elapsed time (initializing wavefunction) = %lf sec\n\n", t_1_elapsed); 

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

        x1_min_0 = BoxDynamic[0];
        x1_max_0 = BoxDynamic[1];
        x2_min_0 = BoxDynamic[2];
        x2_max_0 = BoxDynamic[3];

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

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if (!QUIET && TIMING) log->log("[CaldeiraLeggett4d] Elapsed time (initial density) = %lf sec\n\n", t_1_elapsed); 
    }

    // .........................................................................................

    // Initial truncation & edge point check

    t_1_begin = omp_get_wtime();

    log->log("[CaldeiraLeggett4d] Initial truncation ...\n");

    ta_size = 0;

    // Truncation and TA
    #pragma omp parallel for reduction(+: ta_size) private(b1, b2, b3, b4, b5)
    for (int i1 = 1; i1 < BoxShapeFine[0] - 1; i1 ++)  {
        for (int i2 = 1; i2 < BoxShapeFine[1] - 1; i2 ++)  {
            for (int i3 = 1; i3 < BoxShapeFine[2] - 1; i3 ++)  {
                for (int i4 = 1; i4 < BoxShapeFine[3] - 1; i4 ++)  {

                    b1 = std::abs(PF[i1*W1+i2][i3*W2+i4]) < TolH;
                    b2 = std::abs(PF[(i1+1)*W1+i2][i3*W2+i4] - PF[(i1-1)*W1+i2][i3*W2+i4]) < TolHdX1;
                    b3 = std::abs(PF[i1*W1+(i2+1)][i3*W2+i4] - PF[i1*W1+(i2-1)][i3*W2+i4]) < TolHdX2;
                    b4 = std::abs(PF[i1*W1+i2][(i3+1)*W2+i4] - PF[i1*W1+i2][(i3-1)*W2+i4]) < TolHdX3;
                    b5 = std::abs(PF[i1*W1+i2][i3*W2+(i4+1)] - PF[i1*W1+i2][i3*W2+(i4-1)]) < TolHdX4;
            
                    if ( b1 && b2 && b3 && b4 && b5 )
                        F[i1*W1+i2][i3*W2+i4] = 0.0;
                    else {
                        TAMask[i1*W1+i2][i3*W2+i4] = 1;
                        ta_size += 1;
                    }
                }
            }
        }
    }   
    log->log("[CaldeiraLeggett4d] TA size = %lu\n", ta_size);
    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    if (!QUIET && TIMING) log->log("[CaldeiraLeggett4d] Elapsed time (truncation) = %lf sec\n\n", t_1_elapsed); 

    // TA box
    t_1_begin = omp_get_wtime();

    #pragma omp parallel for reduction(min: x1_min, x2_min, x3_min, x4_min) reduction(max: x1_max, x2_max, x3_max, x4_max)
    for (int i1 = 1; i1 < BoxShapeFine[0] - 1; i1 ++)  {
        for (int i2 = 1; i2 < BoxShapeFine[1] - 1; i2 ++)  {
            for (int i3 = 1; i3 < BoxShapeFine[2] - 1; i3 ++)  {
                for (int i4 = 1; i4 < BoxShapeFine[3] - 1; i4 ++)  {
                    if ( TAMask[i1*W1+i2][i3*W2+i4] )  {
                        if (i1 < x1_min) x1_min = i1;
                        if (i1 > x1_max) x1_max = i1;
                        if (i2 < x2_min) x2_min = i2;
                        if (i2 > x2_max) x2_max = i2;
                        if (i3 < x3_min) x3_min = i3;
                        if (i3 > x3_max) x3_max = i3;
                        if (i4 < x4_min) x4_min = i4;
                        if (i4 > x4_max) x4_max = i4;
                    }
                }
            }
        }
    }
    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    if (!QUIET && TIMING) log->log("[CaldeiraLeggett4d] Elapsed time (TA range) = %lf sec\n\n", t_1_elapsed); 

    // `````````````````````````````````````````````````````````````````
    // TB
    t_1_begin = omp_get_wtime();
    tmpVec.clear();

    #pragma omp parallel for reduction(merge: tmpVec)
    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
            for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                    if (TAMask[i1*W1+i2][i3*W2+i4] == 1 && (!TAMask[(i1-1)*W1+i2][i3*W2+i4] || !TAMask[(i1+1)*W1+i2][i3*W2+i4] || !TAMask[i1*W1+(i2-1)][i3*W2+i4] || !TAMask[i1*W1+(i2+1)][i3*W2+i4] || !TAMask[i1*W1+i2][(i3-1)*W2+i4] || !TAMask[i1*W1+i2][(i3+1)*W2+i4] || !TAMask[i1*W1+i2][i3*W2+(i4-1)] || !TAMask[i1*W1+i2][i3*W2+(i4+1)]))
                        tmpVec.push_back( i1*M1+i2*M2+i3*M3+i4 );       
                }
            }
        }
    }
    tmpVec.swap(TB);
    tmpVec.clear();
    tb_size = TB.size();
    log->log("[CaldeiraLeggett4d] TB size = %d\n", tb_size);
    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    if (!QUIET && TIMING) log->log("[CaldeiraLeggett4d] Elapsed time (TB) = %lf sec\n\n", t_1_elapsed); 

    // `````````````````````````````````````````````````````````````````

    // TA expansion
    t_1_begin = omp_get_wtime();

    #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, g3, g4)
    for (unsigned long int i = 0; i < TB.size(); i++)
    {
        g1 = (int)(TB[i] / M1);
        g2 = (int)((TB[i] % M1) / M2);
        g3 = (int)((TB[i] % M2) / M3);
        g4 = (int)(TB[i] % M3);

        if ( g1+1 != BoxShapeFine[0]-3 && !TAMask[(g1+1)*W1+g2][g3*W2+g4])
            tmpVec.push_back((g1+1)*M1+g2*M2+g3*M3+g4);
        if ( g1-1 != 2 && !TAMask[(g1-1)*W1+g2][g3*W2+g4])
            tmpVec.push_back((g1-1)*M1+g2*M2+g3*M3+g4);
        if ( g2+1 != BoxShapeFine[1]-3 && !TAMask[g1*W1+(g2+1)][g3*W2+g4])
            tmpVec.push_back(g1*M1+(g2+1)*M2+g3*M3+g4);
        if ( g2-1 != 2 && !TAMask[g1*W1+(g2-1)][g3*W2+g4])
            tmpVec.push_back(g1*M1+(g2-1)*M2+g3*M3+g4);
        if ( g3+1 != BoxShapeFine[2]-3 && !TAMask[g1*W1+g2][(g3+1)*W2+g4])
            tmpVec.push_back(g1*M1+g2*M2+(g3+1)*M3+g4);
        if ( g3-1 != 2 && !TAMask[g1*W1+g2][(g3-1)*W2+g4])
            tmpVec.push_back(g1*M1+g2*M2+(g3-1)*M3+g4);
        if ( g4+1 != BoxShapeFine[3]-3 && !TAMask[g1*W1+g2][g3*W2+(g4+1)])
            tmpVec.push_back(g1*M1+g2*M2+g3*M3+g4+1);
        if ( g4-1 != 2 && !TAMask[g1*W1+g2][g3*W2+(g4-1)])
            tmpVec.push_back(g1*M1+g2*M2+g3*M3+g4-1);
    }
    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    if (!QUIET && TIMING) log->log("[CaldeiraLeggett4d] Elapsed time (TA expansion) = %lf sec\n\n", t_1_elapsed); 

    t_1_begin = omp_get_wtime();

    #pragma omp parallel for reduction(min: x1_min, x2_min, x3_min, x4_min) reduction(max: x1_max, x2_max, x3_max, x4_max) private(g1, g2, g3, g4)
    for (unsigned long int i = 0; i < tmpVec.size(); i ++)  {

        g1 = (int)(tmpVec[i]/M1);
        g2 = (int)((tmpVec[i]%M1)/M2);
        g3 = (int)((tmpVec[i]%M2)/M3);
        g4 = (int)(tmpVec[i]%M3);
        TAMask[g1*W1+g2][g3*W2+g4] = 1;

        // Update TA box
        x1_min = (g1 < x1_min) ? g1 : x1_min;
        x2_min = (g2 < x2_min) ? g2 : x2_min;
        x3_min = (g3 < x3_min) ? g3 : x3_min;
        x4_min = (g4 < x4_min) ? g4 : x4_min;
        x1_max = (g1 > x1_max) ? g1 : x1_max;
        x2_max = (g2 > x2_max) ? g2 : x2_max;
        x3_max = (g3 > x3_max) ? g3 : x3_max;
        x4_max = (g4 > x4_max) ? g4 : x4_max;
    }
    tmpVec.clear();

    ta_size = 0;
    #pragma omp parallel for reduction(+: ta_size)
    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
            for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                    if (TAMask[i1*W1+i2][i3*W2+i4]) {
                        ta_size += 1; 
                    }
                }
            }
        }
    }
    log->log("[CaldeiraLeggett4d] TA size = %lu, TB size = %lu\n", ta_size, tb_size);

    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    if (!QUIET && TIMING)  log->log("[CaldeiraLeggett4d] Elapsed time (TA,TB range) = %lf sec\n\n", t_1_elapsed);
    // .........................................................................................

    // Time iteration 

    log->log("=======================================================\n\n"); 
    log->log("[CaldeiraLeggett4d] Time interation starts ...\n"); 
    log->log("[CaldeiraLeggett4d] Number of steps = %d\n\n", (int)(TIME / kk)); 
    log->log("=======================================================\n\n"); 

    for (unsigned long int tt = 0; tt < (unsigned long int)(TIME / kk); tt ++)
    {
        t_0_begin = omp_get_wtime(); 
        Excount = 0;

        if ( isPrintWavefunc && tt % PRINT_WAVEFUNC_PERIOD == 0 )
        {
            pfile = fopen("wave.dat","a");
            fprintf(pfile, "%lu %lu\n", tt, ta_size);
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                        for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                            if (TAMask[i1*W1+i2][i3*W2+i4])
                                fprintf(pfile, "%.4f %.4f %.4f %.4f %.8e\n", BoxDynamic[0]+i1*H[0], BoxDynamic[2]+i2*H[1], BoxDynamic[4]+i3*H[2], BoxDynamic[6]+i4*H[3], F[i1*W1+i2][i3*W2+i4]);
                        }
                    }  
                }
            }
            fclose(pfile);
        }

        if ( tt % PRINT_PERIOD == 0 && isPrintEdge )  {
            pfile = fopen ("edge.dat","a");
            fprintf(pfile, "%lu %lf %lu\n", tt, tt * kk, TB.size());

            for (unsigned long int i = 0; i < TB.size(); i++)
            {
                g1 = (int)(TB[i] / M1);
                g2 = (int)((TB[i] % M1) / M2);
                g3 = (int)((TB[i] % M2) / M3);
                g4 = (int)(TB[i] % M3);
                xx1 = BoxDynamic[0] + g1 * H[0];
                xx2 = BoxDynamic[2] + g2 * H[1];
                xx3 = BoxDynamic[4] + g3 * H[2];
                xx4 = BoxDynamic[6] + g4 * H[3];
                fprintf(pfile, "%.4f %.4f %.4e\n", xx1, xx2,PF[g1*W1+g2][g3*W2+g4]);          
            }
            fclose(pfile);
        }
            
        if ( tt % PRINT_PERIOD == 0 && isPrintDensity )  {
            pfile = fopen ("density.dat","a");
            fprintf(pfile, "%lu %lf %lu\n", tt, tt * kk, (unsigned long int)((x1_max-x1_min+1)*(x2_max-x2_min+1)));

            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    density = 0.0;
                    
                    for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                        for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                            if (TAMask[i1*W1+i2][i3*W2+i4])  {
                                density += PF[i1*W1+i2][i3*W2+i4];
                            }
                        }
                    }
                    xx1 = BoxDynamic[0] + i1 * H[0];
                    xx2 = BoxDynamic[2] + i2 * H[1];
                    fprintf(pfile, "%.4f %.4f %.16e\n", xx1, xx2, density);
                }
            }
            fclose(pfile);
        }

        // Check if TB of f is higher than TolL
        t_1_begin = omp_get_wtime();

        // TBL = Index of Extrapolating Edge points.
        // TBL_P = Index history of TBL of this iteration. To Prevent from extrapolating the same points multiple times.

        TBL.clear();
        tmpVec.clear();

        #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, g3, g4, b1, b2, b3, b4, b5, b6, b7)
        for (unsigned long int i = 0; i < TB.size(); i++)
        {
            g1 = (int)(TB[i] / M1);
            g2 = (int)((TB[i] % M1) / M2);
            g3 = (int)((TB[i] % M2) / M3);
            g4 = (int)(TB[i] % M3);
            b1 = PF[g1*W1+g2][g3*W2+g4] >= TolL;
            b2 = std::abs(PF[(g1+1)*W1+g2][g3*W2+g4]*TAMask[(g1+1)*W1+g2][g3*W2+g4] - PF[(g1-1)*W1+g2][g3*W2+g4]*TAMask[(g1-1)*W1+g2][g3*W2+g4]) >= TolLdX1;
            b3 = std::abs(PF[g1*W1+(g2+1)][g3*W2+g4]*TAMask[g1*W1+(g2+1)][g3*W2+g4] - PF[g1*W1+(g2-1)][g3*W2+g4]*TAMask[g1*W1+(g2-1)][g3*W2+g4]) >= TolLdX2;
            b4 = std::abs(PF[g1*W1+g2][(g3+1)*W2+g4]*TAMask[g1*W1+g2][(g3+1)*W2+g4] - PF[g1*W1+g2][(g3-1)*W2+g4]*TAMask[g1*W1+g2][(g3-1)*W2+g4]) >= TolLdX3;
            b5 = std::abs(PF[g1*W1+g2][g3*W2+(g4+1)]*TAMask[g1*W1+g2][g3*W2+(g4+1)] - PF[g1*W1+g2][g3*W2+(g4-1)]*TAMask[g1*W1+g2][g3*W2+(g4-1)]) >= TolLdX4;

            // Not in DBi2
            b6 = g1 > 2 && g2 > 2 && g3 > 2 && g4 > 2;
            b7 = g1 < BoxShapeFine[0] - 3 && g2 < BoxShapeFine[1] - 3 && g3 < BoxShapeFine[2] - 3 && g4 < BoxShapeFine[3] - 3;

            if ( (b1 || b2 || b3 || b4 || b5) && b6 && b7 )
                tmpVec.push_back(TB[i]);
        }
        tmpVec.swap(TBL);
        tmpVec.clear();
        TBL_P = TBL;

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-a-1: TBL) = %lf sec\n", t_1_elapsed);   

        isExtrapolate = false;
        isFirstExtrp = true;
        // .........................................................................................

        // CASE 1: Truncating with extrapolation

        while ( TBL.size() != 0 && Excount < ExLimit )
        {
            isExtrapolate = true;

            // Extrapolation
            // .............................................................................................

            t_1_begin = omp_get_wtime();

            // Avoid unexpected arrangement of TBL
            __gnu_parallel::sort(TBL.begin(),TBL.end());
            it = std::unique (TBL.begin(), TBL.end()); 
            TBL.resize(std::distance(TBL.begin(),it)); 

            // Find extrapolation target
            // ExFF: Index of Extrapolated points
            ExFF.clear();

            for ( it = TBL.begin(); it != TBL.end(); it ++ )  {

                index = std::distance( TBL.begin(), it );

                g1 = (int)(TBL[index] / M1);
                g2 = (int)((TBL[index] % M1) / M2);
                g3 = (int)((TBL[index] % M2) / M3);
                g4 = (int)(TBL[index] % M3);

                if ( g1-1 != 2 && F[(g1-1)*W1+g2][g3*W2+g4] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1-1,g2,g3,g4));
                }
                if ( g1+1 != BoxShapeFine[0]-3 && F[(g1+1)*W1+g2][g3*W2+g4] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1+1,g2,g3,g4));
                }
                if ( g2-1 != 2 && F[g1*W1+(g2-1)][g3*W2+g4] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1,g2-1,g3,g4));
                }
                if ( g2+1 != BoxShapeFine[1]-3 && F[g1*W1+(g2+1)][g3*W2+g4] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1,g2+1,g3,g4));
                }
                if ( g3-1 != 2 && F[g1*W1+g2][(g3-1)*W2+g4] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1,g2,g3-1,g4));
                }
                if ( g3+1 != BoxShapeFine[2]-3 && F[g1*W1+g2][(g3+1)*W2+g4] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1,g2,g3+1,g4));
                }
                if ( g4-1 != 2 && F[g1*W1+g2][g3*W2+(g4-1)] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1,g2,g3,g4-1));
                }
                if ( g4+1 != BoxShapeFine[3]-3 && F[g1*W1+g2][g3*W2+(g4+1)] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1,g2,g3,g4+1));
                }
            }

            // ExFF & TBL set difference
            tmpVec.resize(ExFF.size() + TBL.size());
            __gnu_parallel::sort(TBL.begin(), TBL.end());
            __gnu_parallel::sort(ExFF.begin(), ExFF.end());
            it = std::set_difference( ExFF.begin(), ExFF.end(), TBL.begin(), TBL.end(), tmpVec.begin() );
            tmpVec.resize(it - tmpVec.begin()); 
            tmpVec.swap(ExFF);
            tmpVec.clear();

            // Find unique elements
            __gnu_parallel::sort(ExFF.begin(),ExFF.end());
            it = std::unique (ExFF.begin(), ExFF.end()); 
            ExFF.resize(std::distance(ExFF.begin(),it));

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-b-1: ExFF) = %lf sec\n", t_1_elapsed);   

            // .....................................................................

            // Find the direction of Outer to Edge points

            t_1_begin = omp_get_wtime();
            ExTBL.clear();
            it = ExFF.begin();
            
            while ( it != ExFF.end() )  {

                index = std::distance( ExFF.begin(), it );
                g1 = (int)(ExFF[index] / M1);
                g2 = (int)((ExFF[index] % M1) / M2);
                g3 = (int)((ExFF[index] % M2) / M3);
                g4 = (int)(ExFF[index] % M3);                
                sum = 0.0;
                count = 0;
                isEmpty = true;
                val_min_abs = 100000000;
                val_min = 100000000;
                min_dir = -1;

                if ( F[(g1-1)*W1+g2][g3*W2+g4] != 0.0 )  {

                    if ( std::abs(F[(g1-1)*W1+g2][g3*W2+g4]) < val_min_abs &&  F[(g1-2)*W1+g2][g3*W2+g4] != 0.0 )  {
                        val_min_abs = std::abs(F[(g1-1)*W1+g2][g3*W2+g4]);
                        val_min = F[(g1-1)*W1+g2][g3*W2+g4];
                        min_dir = 0;
                    }
                    if ( F[(g1-2)*W1+g2][g3*W2+g4] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[(g1-1)*W1+g2][g3*W2+g4]) - std::log(F[(g1-2)*W1+g2][g3*W2+g4]) );

                        if ( !(isnan(val) || isnan(-val)) && !(isinf(val) || isinf(-val)) )  
                        {
                            sum += val;
                            count += 1;
                            isEmpty = false;
                        }
                    }
                }

                if ( F[(g1+1)*W1+g2][g3*W2+g4] != 0.0 )  {

                    if ( std::abs(F[(g1+1)*W1+g2][g3*W2+g4]) < val_min_abs && F[(g1+2)*W1+g2][g3*W2+g4] != 0.0  )  {
                        val_min_abs = std::abs(F[(g1+1)*W1+g2][g3*W2+g4]);
                        val_min = F[(g1+1)*W1+g2][g3*W2+g4];
                        min_dir = 0;
                    }
                    if ( F[(g1+2)*W1+g2][g3*W2+g4] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[(g1+1)*W1+g2][g3*W2+g4]) - std::log(F[(g1+2)*W1+g2][g3*W2+g4]) );

                        if ( !(isnan(val) || isnan(-val)) && !(isinf(val) || isinf(-val)) ) 
                        {
                            sum += val;
                            count += 1;
                            isEmpty = false;
                        }
                    }
                }

                if ( F[g1*W1+(g2-1)][g3*W2+g4] != 0.0 )  {

                    if ( std::abs(F[g1*W1+(g2-1)][g3*W2+g4]) < val_min_abs && F[g1*W1+(g2-2)][g3*W2+g4] != 0.0  )  {

                        val_min_abs = std::abs(F[g1*W1+(g2-1)][g3*W2+g4]);
                        val_min = F[g1*W1+(g2-1)][g3*W2+g4];
                        min_dir = 1;
                    }
                    if ( F[g1*W1+(g2-2)][g3*W2+g4] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[g1*W1+(g2-1)][g3*W2+g4]) - std::log(F[g1*W1+(g2-2)][g3*W2+g4]) );

                        if ( !(isnan(val) || isnan(-val)) && !(isinf(val) || isinf(-val)) ) 
                        {
                            sum += val;
                            count += 1;
                            isEmpty = false;
                        }
                    }
                }

                if ( F[g1*W1+(g2+1)][g3*W2+g4] != 0.0 )  {

                    if ( std::abs(F[g1*W1+(g2+1)][g3*W2+g4]) < val_min_abs && F[g1*W1+(g2+2)][g3*W2+g4] != 0.0 )  {

                        val_min_abs = std::abs(F[g1*W1+(g2+1)][g3*W2+g4]);
                        val_min = F[g1*W1+(g2+1)][g3*W2+g4];
                        min_dir = 1;
                    }
                    if ( F[g1*W1+(g2+2)][g3*W2+g4] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[g1*W1+(g2+1)][g3*W2+g4]) - std::log(F[g1*W1+(g2+2)][g3*W2+g4]) );

                        if ( !(isnan(val) || isnan(-val)) && !(isinf(val) || isinf(-val)) ) 
                        {
                            sum += val;
                            count += 1;
                            isEmpty = false;
                        }
                    }
                }

                if ( F[g1*W1+g2][(g3-1)*W2+g4] != 0.0 )  {

                    if ( std::abs(F[g1*W1+g2][(g3-1)*W2+g4]) < val_min_abs && F[g1*W1+g2][(g3-2)*W2+g4] != 0.0  )  {

                        val_min_abs = std::abs(F[g1*W1+g2][(g3-1)*W2+g4]);
                        val_min = F[g1*W1+g2][(g3-1)*W2+g4];
                        min_dir = 2;
                    }
                    if ( F[g1*W1+g2][(g3-2)*W2+g4] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[g1*W1+g2][(g3-1)*W2+g4]) - std::log(F[g1*W1+g2][(g3-2)*W2+g4]) );

                        if ( !(isnan(val) || isnan(-val)) && !(isinf(val) || isinf(-val)) ) 
                        {
                            sum += val;
                            count += 1;
                            isEmpty = false;
                        }
                    }
                }

                if ( F[g1*W1+g2][(g3+1)*W2+g4] != 0.0 )  {

                    if ( std::abs(F[g1*W1+g2][(g3+1)*W2+g4]) < val_min_abs && F[g1*W1+g2][(g3+2)*W2+g4] != 0.0 )  {

                        val_min_abs = std::abs(F[g1*W1+g2][(g3+1)*W2+g4]);
                        val_min = F[g1*W1+g2][(g3+1)*W2+g4];
                        min_dir = 2;
                    }
                    if ( F[g1*W1+g2][(g3+2)*W2+g4] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[g1*W1+g2][(g3+1)*W2+g4]) - std::log(F[g1*W1+g2][(g3+2)*W2+g4]) );

                        if ( !(isnan(val) || isnan(-val)) && !(isinf(val) || isinf(-val)) ) 
                        {
                            sum += val;
                            count += 1;
                            isEmpty = false;
                        }
                    }
                }

                if ( F[g1*W1+g2][g3*W2+(g4-1)] != 0.0 )  {

                    if ( std::abs(F[g1*W1+g2][g3*W2+(g4-1)]) < val_min_abs && F[g1*W1+g2][g3*W2+(g4-2)] != 0.0  )  {

                        val_min_abs = std::abs(F[g1*W1+g2][g3*W2+(g4-1)]);
                        val_min = F[g1*W1+g2][g3*W2+(g4-1)];
                        min_dir = 3;
                    }
                    if ( F[g1*W1+g2][g3*W2+(g4-2)] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[g1*W1+g2][g3*W2+(g4-1)]) - std::log(F[g1*W1+g2][g3*W2+(g4-2)]) );

                        if ( !(isnan(val) || isnan(-val)) && !(isinf(val) || isinf(-val)) ) 
                        {
                            sum += val;
                            count += 1;
                            isEmpty = false;
                        }
                    }
                }

                if ( F[g1*W1+g2][g3*W2+(g4+1)] != 0.0 )  {

                    if ( std::abs(F[g1*W1+g2][g3*W2+(g4+1)]) < val_min_abs && F[g1*W1+g2][g3*W2+(g4+2)] != 0.0 )  {

                        val_min_abs = std::abs(F[g1*W1+g2][g3*W2+(g4+1)]);
                        val_min = F[g1*W1+g2][g3*W2+(g4+1)];
                        min_dir = 3;
                    }
                    if ( F[g1*W1+g2][g3*W2+(g4+2)] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[g1*W1+g2][g3*W2+(g4+1)]) - std::log(F[g1*W1+g2][g3*W2+(g4+2)]) );

                        if ( !(isnan(val) || isnan(-val)) && !(isinf(val) || isinf(-val)) ) 
                        {
                            sum += val;
                            count += 1;
                            isEmpty = false;
                        }
                    }
                }

                if ( isEmpty )
                {
                    it = ExFF.erase(it);   

                }  else  {

                    // Assume the probability of outer points always smaller than edge point.
                    // if larger, then set the edge point with smallest P to outer
                    // point instead of using the extrapolation result.
                    if ( std::abs(sum / count) > val_min_abs )  {
                        ExTBL.push_back(val_min * exp(-ExReduce * H[min_dir]));
                    }
                    else  {
                        ExTBL.push_back(sum / count);
                    }
                    ++it;
                }
            }

            for ( unsigned long int i = 0; i < ExFF.size(); i++ )  {
                g1 = (int)(ExFF[i] / M1);
                g2 = (int)((ExFF[i] % M1) / M2);
                g3 = (int)((ExFF[i] % M2) / M3);
                g4 = (int)(ExFF[i] % M3);  
                F[g1*W1+g2][g3*W2+g4] = ExTBL[i];
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-b-2: ExFF) = %lf sec\n", t_1_elapsed);  

            // ............................................................................................. Extrapolation

            if ( isFirstExtrp )  {

                // Check Extending nonzero Area

                t_1_begin = omp_get_wtime();

                #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4)
                for (unsigned long int i = 0; i < ExFF.size(); i++)  {
                    g1 = (int)(ExFF[i] / M1);
                    g2 = (int)((ExFF[i] % M1) / M2);
                    g3 = (int)((ExFF[i] % M2) / M3);
                    g4 = (int)(ExFF[i] % M3); 

                    if ( !TAMask[g1*W1+g2][g3*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4));
                    if ( g1+1 != BoxShapeFine[0]-3 && !TAMask[(g1+1)*W1+g2][g3*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                    if ( g1-1 != 2 && !TAMask[(g1-1)*W1+g2][g3*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));
                    if ( g2+1 != BoxShapeFine[1]-3 && !TAMask[g1*W1+(g2+1)][g3*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                    if ( g2-1 != 2 && !TAMask[g1*W1+(g2-1)][g3*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));
                    if ( g3+1 != BoxShapeFine[2]-3 && !TAMask[g1*W1+g2][(g3+1)*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                    if ( g3-1 != 2 && !TAMask[g1*W1+g2][(g3-1)*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));
                    if ( g4+1 != BoxShapeFine[3]-3 && !TAMask[g1*W1+g2][g3*W2+(g4+1)] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                    if ( g4-1 != 2 && !TAMask[g1*W1+g2][g3*W2+(g4-1)] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
                }

                #pragma omp parallel for reduction(min: x1_min, x2_min, x3_min, x4_min) reduction(max: x1_max, x2_max, x3_max, x4_max) private(g1, g2, g3, g4)
                for (unsigned long int i = 0; i < tmpVec.size(); i ++)  {

                    g1 = (int)(tmpVec[i] / M1);
                    g2 = (int)((tmpVec[i] % M1) / M2);
                    g3 = (int)((tmpVec[i] % M2) / M3);
                    g4 = (int)(tmpVec[i] % M3); 
                    TAMask[g1*W1+g2][g3*W2+g4] = 1;
                    
                    // Update TA box
                    x1_min = (g1 < x1_min) ? g1 : x1_min;
                    x2_min = (g2 < x2_min) ? g2 : x2_min;
                    x3_min = (g3 < x3_min) ? g3 : x3_min;
                    x4_min = (g4 < x4_min) ? g4 : x4_min;
                    x1_max = (g1 > x1_max) ? g1 : x1_max;
                    x2_max = (g2 > x2_max) ? g2 : x2_max;
                    x3_max = (g3 > x3_max) ? g3 : x3_max;
                    x4_max = (g4 > x4_max) ? g4 : x4_max;
                }
                tmpVec.clear();
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-c-1: CASE 1 TA) = %lf sec\n", t_1_elapsed); 
        
                // Runge–Kutta 4
                #pragma omp parallel 
                {
                    // RK4-1
                    #pragma omp single nowait
                    {
                        t_1_begin = omp_get_wtime();
                    }
                    #pragma omp for private(xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2) schedule(runtime)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                                for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                    if ( TAMask[i1*W1+i2][i3*W2+i4] )  {
                                        xx1 = BoxDynamic[0] + i1 * H[0];
                                        xx2 = BoxDynamic[2] + i2 * H[1];
                                        xx3 = BoxDynamic[4] + i3 * H[2];
                                        xx4 = BoxDynamic[6] + i4 * H[3];
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
                    }
                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-11: CASE 1 KK1) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-2
                    #pragma omp for private(xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk3p1,kk3m1,kk4p1,kk4m1,kk1p2,kk1m2,kk2p2,kk2m2,kk3p2,kk3m2,kk4p2,kk4m2) schedule(runtime)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                                for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                    if ( TAMask[i1*W1+i2][i3*W2+i4] )  {
                                        xx1 = BoxDynamic[0] + i1 * H[0];
                                        xx2 = BoxDynamic[2] + i2 * H[1];
                                        xx3 = BoxDynamic[4] + i3 * H[2];
                                        xx4 = BoxDynamic[6] + i4 * H[3];
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
                    }
                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-12: CASE 1 KK2) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-3
                    #pragma omp for private(xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk3p1,kk3m1,kk4p1,kk4m1,kk1p2,kk1m2,kk2p2,kk2m2,kk3p2,kk3m2,kk4p2,kk4m2) schedule(runtime)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                                for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                    if ( TAMask[i1*W1+i2][i3*W2+i4] )  {
                                        xx1 = BoxDynamic[0] + i1 * H[0];
                                        xx2 = BoxDynamic[2] + i2 * H[1];
                                        xx3 = BoxDynamic[4] + i3 * H[2];
                                        xx4 = BoxDynamic[6] + i4 * H[3];
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
                    }
                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-13: CASE 1 KK3) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-4
                    #pragma omp for private(xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk3p1,kk3m1,kk4p1,kk4m1,kk1p2,kk1m2,kk2p2,kk2m2,kk3p2,kk3m2,kk4p2,kk4m2) schedule(runtime)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                                for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                    if ( TAMask[i1*W1+i2][i3*W2+i4] )  {
                                        xx1 = BoxDynamic[0] + i1 * H[0];
                                        xx2 = BoxDynamic[2] + i2 * H[1];
                                        xx3 = BoxDynamic[4] + i3 * H[2];
                                        xx4 = BoxDynamic[6] + i4 * H[3];
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
                    }
                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-14: CASE 1 KK4) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }
                } // OMP PARALLEL
            }
            else  
            {
                // Extrapolation loop when multiple expanding occured

                t_1_begin = omp_get_wtime();

                #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4)
                for (unsigned long int i = 0; i < ExFF.size(); i++)  {
                    g1 = (int)(ExFF[i] / M1);
                    g2 = (int)((ExFF[i] % M1) / M2);
                    g3 = (int)((ExFF[i] % M2) / M3);
                    g4 = (int)(ExFF[i] % M3); 

                    if ( !TAMask[g1*W1+g2][g3*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4));
                    if ( g1+1 != BoxShapeFine[0]-3 && !TAMask[(g1+1)*W1+g2][g3*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                    if ( g1-1 != 2 && !TAMask[(g1-1)*W1+g2][g3*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));
                    if ( g2+1 != BoxShapeFine[1]-3 && !TAMask[g1*W1+(g2+1)][g3*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                    if ( g2-1 != 2 && !TAMask[g1*W1+(g2-1)][g3*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));
                    if ( g3+1 != BoxShapeFine[2]-3 && !TAMask[g1*W1+g2][(g3+1)*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                    if ( g3-1 != 2 && !TAMask[g1*W1+g2][(g3-1)*W2+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));
                    if ( g4+1 != BoxShapeFine[3]-3 && !TAMask[g1*W1+g2][g3*W2+(g4+1)] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                    if ( g4-1 != 2 && !TAMask[g1*W1+g2][g3*W2+(g4-1)] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
                }

                #pragma omp parallel for reduction(min: x1_min, x2_min, x3_min, x4_min) reduction(max: x1_max, x2_max, x3_max, x4_max) private(g1, g2, g3, g4)
                for (unsigned long int i = 0; i < tmpVec.size(); i ++)  {

                    g1 = (int)(tmpVec[i] / M1);
                    g2 = (int)((tmpVec[i] % M1) / M2);
                    g3 = (int)((tmpVec[i] % M2) / M3);
                    g4 = (int)(tmpVec[i] % M3); 
                    TAMask[g1*W1+g2][g3*W2+g4] = 1;

                    // Update TA box
                    x1_min = (g1 < x1_min) ? g1 : x1_min;
                    x2_min = (g2 < x2_min) ? g2 : x2_min;
                    x3_min = (g3 < x3_min) ? g3 : x3_min;
                    x4_min = (g4 < x4_min) ? g4 : x4_min;
                    x1_max = (g1 > x1_max) ? g1 : x1_max;
                    x2_max = (g2 > x2_max) ? g2 : x2_max;
                    x3_max = (g3 > x3_max) ? g3 : x3_max;
                    x4_max = (g4 > x4_max) ? g4 : x4_max;
                }
                tmpVec.clear();
                ExBD.clear();

                #pragma omp parallel for reduction(merge: ExBD) private(g1,g2,g3,g4,n1,n2,n3,n4)
                for (unsigned long int i = 0; i < ExFF.size(); i++)
                {
                    g1 = (int)(ExFF[i] / M1);
                    g2 = (int)((ExFF[i] % M1) / M2);
                    g3 = (int)((ExFF[i] % M2) / M3);
                    g4 = (int)(ExFF[i] % M3); 
                    ExBD.push_back(ExFF[i]);

                    for (int j = 0; j < nneigh; j ++)  {
                        n1 = neighlist[j][0];
                        n2 = neighlist[j][1];
                        n3 = neighlist[j][2];
                        n4 = neighlist[j][3];

                        if (TAMask[(g1+n1)*W1+(g2+n2)][(g3+n3)*W2+(g4+n4)])  {
                            ExBD.push_back(GridToIdx(g1+n1,g2+n2,g3+n3,g4+n4));
                        }
                    }
                }
                // Find unique elements (ExBD)
                __gnu_parallel::sort(ExBD.begin(),ExBD.end());
                it = std::unique (ExBD.begin(), ExBD.end()); 
                ExBD.resize(std::distance(ExBD.begin(),it));

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-cx-1: CASE 1 ExBD) = %lf sec\n", t_1_elapsed); 

                // Runge–Kutta 4

                t_1_begin = omp_get_wtime();

                // RK4-1
                #pragma omp parallel for private(g1,g2,g3,g4,xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2) schedule(runtime)
                for (unsigned long int i = 0; i < ExBD.size(); i++)  {

                    g1 = (int)(ExBD[i] / M1);
                    g2 = (int)((ExBD[i] % M1) / M2);
                    g3 = (int)((ExBD[i] % M2) / M3);
                    g4 = (int)(ExBD[i] % M3);
                    xx1 = BoxDynamic[0] + g1 * H[0];
                    xx2 = BoxDynamic[2] + g2 * H[1];
                    xx3 = BoxDynamic[4] + g3 * H[2];
                    xx4 = BoxDynamic[6] + g4 * H[3];
                    f0 = F[g1*W1+g2][g3*W2+g4];
                    f1p1 = F[(g1+1)*W1+g2][g3*W2+g4];
                    f1m1 = F[(g1-1)*W1+g2][g3*W2+g4];
                    f2p1 = F[g1*W1+(g2+1)][g3*W2+g4];
                    f2m1 = F[g1*W1+(g2-1)][g3*W2+g4];
                    f3p1 = F[g1*W1+g2][(g3+1)*W2+g4];
                    f3m1 = F[g1*W1+g2][(g3-1)*W2+g4];
                    f4p1 = F[g1*W1+g2][g3*W2+(g4+1)];
                    f4m1 = F[g1*W1+g2][g3*W2+(g4-1)];
                    f1p2 = F[(g1+2)*W1+g2][g3*W2+g4];
                    f1m2 = F[(g1-2)*W1+g2][g3*W2+g4];
                    f2p2 = F[g1*W1+(g2+2)][g3*W2+g4];
                    f2m2 = F[g1*W1+(g2-2)][g3*W2+g4];
                    f3p2 = F[g1*W1+g2][(g3+2)*W2+g4];
                    f3m2 = F[g1*W1+g2][(g3-2)*W2+g4];
                    f4p2 = F[g1*W1+g2][g3*W2+(g4+2)];
                    f4m2 = F[g1*W1+g2][g3*W2+(g4-2)];

                    KK1[g1*W1+g2][g3*W2+g4] = (-kh0m) * xx3 * (-f1p2/12.0 + 2/3.0*f1p1 - 2/3.0*f1m1 + f1m2/12.0) +
                                            (-kh1m) * xx4 * (-f2p2/12.0 + 2/3.0*f2p1 - 2/3.0*f2m1 + f2m2/12.0) + 
                                            k2h2 * POTENTIAL_X1(xx1, xx2) * (-f3p2/12.0 + 2/3.0*f3p1 - 2/3.0*f3m1 + f3m2/12.0) +
                                            k2h3 * POTENTIAL_X2(xx1, xx2) * (-f4p2/12.0 + 2/3.0*f4p1 - 2/3.0*f4m1 + f4m2/12.0) +
                                            (-khbsq2h2) * POTENTIAL_XXX1(xx1, xx2) * (0.5*f3p2 - f3p1 + f3m1 - 0.5*f3m2) +
                                            (-khbsq2h3) * POTENTIAL_XXX2(xx1, xx2) * (0.5*f4p2 - f4p1 + f4m1 - 0.5*f4m2) +
                                            kgamma * (f0 + i2h2 * xx3 * (f3p1 - f3m1) + mkT2h2sq * (f3p1 + f3m1 - 2*f0)) +
                                            kgamma * (f0 + i2h3 * xx4 * (f4p1 - f4m1) + mkT2h3sq * (f4p1 + f4m1 - 2*f0));

                    FF[g1*W1+g2][g3*W2+g4] = F[g1*W1+g2][g3*W2+g4] + KK1[g1*W1+g2][g3*W2+g4] / 6.0;
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-11: CASE 1 KK1) = %lf sec\n", t_1_elapsed);
                t_1_begin = omp_get_wtime();

                // RK4-2
                #pragma omp parallel for private(g1,g2,g3,g4,xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk3p1,kk3m1,kk4p1,kk4m1,kk1p2,kk1m2,kk2p2,kk2m2,kk3p2,kk3m2,kk4p2,kk4m2) schedule(runtime)
                for (int i = 0; i < ExBD.size(); i++)  {

                    g1 = (int)(ExBD[i] / M1);
                    g2 = (int)((ExBD[i] % M1) / M2);
                    g3 = (int)((ExBD[i] % M2) / M3);
                    g4 = (int)(ExBD[i] % M3);
                    xx1 = BoxDynamic[0] + g1 * H[0];
                    xx2 = BoxDynamic[2] + g2 * H[1];
                    xx3 = BoxDynamic[4] + g3 * H[2];
                    xx4 = BoxDynamic[6] + g4 * H[3];
                    f0 = F[g1*W1+g2][g3*W2+g4];
                    f1p1 = F[(g1+1)*W1+g2][g3*W2+g4];
                    f1m1 = F[(g1-1)*W1+g2][g3*W2+g4];
                    f2p1 = F[g1*W1+(g2+1)][g3*W2+g4];
                    f2m1 = F[g1*W1+(g2-1)][g3*W2+g4];
                    f3p1 = F[g1*W1+g2][(g3+1)*W2+g4];
                    f3m1 = F[g1*W1+g2][(g3-1)*W2+g4]; 
                    f4p1 = F[g1*W1+g2][g3*W2+(g4+1)];
                    f4m1 = F[g1*W1+g2][g3*W2+(g4-1)];          
                    f1p2 = F[(g1+2)*W1+g2][g3*W2+g4];
                    f1m2 = F[(g1-2)*W1+g2][g3*W2+g4];
                    f2p2 = F[g1*W1+(g2+2)][g3*W2+g4];
                    f2m2 = F[g1*W1+(g2-2)][g3*W2+g4];
                    f3p2 = F[g1*W1+g2][(g3+2)*W2+g4];
                    f3m2 = F[g1*W1+g2][(g3-2)*W2+g4]; 
                    f4p2 = F[g1*W1+g2][g3*W2+(g4+2)];
                    f4m2 = F[g1*W1+g2][g3*W2+(g4-2)];
                    kk0 = KK1[g1*W1+g2][g3*W2+g4];
                    kk1p1 = KK1[(g1+1)*W1+g2][g3*W2+g4];
                    kk1m1 = KK1[(g1-1)*W1+g2][g3*W2+g4];
                    kk2p1 = KK1[g1*W1+(g2+1)][g3*W2+g4];
                    kk2m1 = KK1[g1*W1+(g2-1)][g3*W2+g4];
                    kk3p1 = KK1[g1*W1+g2][(g3+1)*W2+g4];
                    kk3m1 = KK1[g1*W1+g2][(g3-1)*W2+g4];
                    kk4p1 = KK1[g1*W1+g2][g3*W2+(g4+1)];
                    kk4m1 = KK1[g1*W1+g2][g3*W2+(g4-1)];
                    kk1p2 = KK1[(g1+2)*W1+g2][g3*W2+g4];
                    kk1m2 = KK1[(g1-2)*W1+g2][g3*W2+g4];
                    kk2p2 = KK1[g1*W1+(g2+2)][g3*W2+g4];
                    kk2m2 = KK1[g1*W1+(g2-2)][g3*W2+g4];
                    kk3p2 = KK1[g1*W1+g2][(g3+2)*W2+g4];
                    kk3m2 = KK1[g1*W1+g2][(g3-2)*W2+g4];
                    kk4p2 = KK1[g1*W1+g2][g3*W2+(g4+2)];
                    kk4m2 = KK1[g1*W1+g2][g3*W2+(g4-2)];

                    KK2[g1*W1+g2][g3*W2+g4] = (-kh0m) * xx3 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                              (-kh1m) * xx4 * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) + 
                                              k2h2 * POTENTIAL_X1(xx1, xx2) * (-1/12.0*(f3p2+0.5*kk3p2) + 2/3.0*(f3p1+0.5*kk3p1) - 2/3.0*(f3m1+0.5*kk3m1) + 1/12.0*(f3m2+0.5*kk3m2)) +
                                              k2h3 * POTENTIAL_X2(xx1, xx2) * (-1/12.0*(f4p2+0.5*kk4p2) + 2/3.0*(f4p1+0.5*kk4p1) - 2/3.0*(f4m1+0.5*kk4m1) + 1/12.0*(f4m2+0.5*kk4m2)) +
                                              (-khbsq2h2) * POTENTIAL_XXX1(xx1, xx2) * (0.5*(f3p2+0.5*kk3p2) - (f3p1+0.5*kk3p1) + (f3m1+0.5*kk3m1) - 0.5*(f3m2+0.5*kk3m2)) +
                                              (-khbsq2h3) * POTENTIAL_XXX2(xx1, xx2) * (0.5*(f4p2+0.5*kk4p2) - (f4p1+0.5*kk4p1) + (f4m1+0.5*kk4m1) - 0.5*(f4m2+0.5*kk4m2)) +
                                              kgamma * (f0 + 0.5*kk0 + i2h2 * xx3 * (f3p1 + 0.5*kk3p1 - f3m1 - 0.5*kk3m1) + mkT2h2sq * (f3p1 + 0.5*kk3p1 + f3m1 + 0.5*kk3m1 - 2*f0 - kk0)) +
                                              kgamma * (f0 + 0.5*kk0 + i2h3 * xx4 * (f4p1 + 0.5*kk4p1 - f4m1 - 0.5*kk4m1) + mkT2h3sq * (f4p1 + 0.5*kk4p1 + f4m1 + 0.5*kk4m1 - 2*f0 - kk0));

                    FF[g1*W1+g2][g3*W2+g4] += KK2[g1*W1+g2][g3*W2+g4] / 3.0;
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-12: CASE 1 KK2) = %lf sec\n", t_1_elapsed);
                t_1_begin = omp_get_wtime();

                // RK4-3
                #pragma omp parallel for private(g1,g2,g3,g4,xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk3p1,kk3m1,kk4p1,kk4m1,kk1p2,kk1m2,kk2p2,kk2m2,kk3p2,kk3m2,kk4p2,kk4m2) schedule(runtime)
                for (int i = 0; i < ExBD.size(); i++)  {

                    g1 = (int)(ExBD[i] / M1);
                    g2 = (int)((ExBD[i] % M1) / M2);
                    g3 = (int)((ExBD[i] % M2) / M3);
                    g4 = (int)(ExBD[i] % M3);
                    xx1 = BoxDynamic[0] + g1 * H[0];
                    xx2 = BoxDynamic[2] + g2 * H[1];
                    xx3 = BoxDynamic[4] + g3 * H[2];
                    xx4 = BoxDynamic[6] + g4 * H[3];
                    f0 = F[g1*W1+g2][g3*W2+g4];
                    f1p1 = F[(g1+1)*W1+g2][g3*W2+g4];
                    f1m1 = F[(g1-1)*W1+g2][g3*W2+g4];
                    f2p1 = F[g1*W1+(g2+1)][g3*W2+g4];
                    f2m1 = F[g1*W1+(g2-1)][g3*W2+g4];
                    f3p1 = F[g1*W1+g2][(g3+1)*W2+g4];
                    f3m1 = F[g1*W1+g2][(g3-1)*W2+g4]; 
                    f4p1 = F[g1*W1+g2][g3*W2+(g4+1)];
                    f4m1 = F[g1*W1+g2][g3*W2+(g4-1)];          
                    f1p2 = F[(g1+2)*W1+g2][g3*W2+g4];
                    f1m2 = F[(g1-2)*W1+g2][g3*W2+g4];
                    f2p2 = F[g1*W1+(g2+2)][g3*W2+g4];
                    f2m2 = F[g1*W1+(g2-2)][g3*W2+g4];
                    f3p2 = F[g1*W1+g2][(g3+2)*W2+g4];
                    f3m2 = F[g1*W1+g2][(g3-2)*W2+g4]; 
                    f4p2 = F[g1*W1+g2][g3*W2+(g4+2)];
                    f4m2 = F[g1*W1+g2][g3*W2+(g4-2)];
                    kk0 = KK2[g1*W1+g2][g3*W2+g4];
                    kk1p1 = KK2[(g1+1)*W1+g2][g3*W2+g4];
                    kk1m1 = KK2[(g1-1)*W1+g2][g3*W2+g4];
                    kk2p1 = KK2[g1*W1+(g2+1)][g3*W2+g4];
                    kk2m1 = KK2[g1*W1+(g2-1)][g3*W2+g4];
                    kk3p1 = KK2[g1*W1+g2][(g3+1)*W2+g4];
                    kk3m1 = KK2[g1*W1+g2][(g3-1)*W2+g4];
                    kk4p1 = KK2[g1*W1+g2][g3*W2+(g4+1)];
                    kk4m1 = KK2[g1*W1+g2][g3*W2+(g4-1)];
                    kk1p2 = KK2[(g1+2)*W1+g2][g3*W2+g4];
                    kk1m2 = KK2[(g1-2)*W1+g2][g3*W2+g4];
                    kk2p2 = KK2[g1*W1+(g2+2)][g3*W2+g4];
                    kk2m2 = KK2[g1*W1+(g2-2)][g3*W2+g4];
                    kk3p2 = KK2[g1*W1+g2][(g3+2)*W2+g4];
                    kk3m2 = KK2[g1*W1+g2][(g3-2)*W2+g4];
                    kk4p2 = KK2[g1*W1+g2][g3*W2+(g4+2)];
                    kk4m2 = KK2[g1*W1+g2][g3*W2+(g4-2)];

                    KK3[g1*W1+g2][g3*W2+g4] = (-kh0m) * xx3 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                              (-kh1m) * xx4 * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) + 
                                              k2h2 * POTENTIAL_X1(xx1, xx2) * (-1/12.0*(f3p2+0.5*kk3p2) + 2/3.0*(f3p1+0.5*kk3p1) - 2/3.0*(f3m1+0.5*kk3m1) + 1/12.0*(f3m2+0.5*kk3m2)) +
                                              k2h3 * POTENTIAL_X2(xx1, xx2) * (-1/12.0*(f4p2+0.5*kk4p2) + 2/3.0*(f4p1+0.5*kk4p1) - 2/3.0*(f4m1+0.5*kk4m1) + 1/12.0*(f4m2+0.5*kk4m2)) +
                                              (-khbsq2h2) * POTENTIAL_XXX1(xx1, xx2) * (0.5*(f3p2+0.5*kk3p2) - (f3p1+0.5*kk3p1) + (f3m1+0.5*kk3m1) - 0.5*(f3m2+0.5*kk3m2)) +
                                              (-khbsq2h3) * POTENTIAL_XXX2(xx1, xx2) * (0.5*(f4p2+0.5*kk4p2) - (f4p1+0.5*kk4p1) + (f4m1+0.5*kk4m1) - 0.5*(f4m2+0.5*kk4m2)) +
                                              kgamma * (f0 + 0.5*kk0 + i2h2 * xx3 * (f3p1 + 0.5*kk3p1 - f3m1 - 0.5*kk3m1) + mkT2h2sq * (f3p1 + 0.5*kk3p1 + f3m1 + 0.5*kk3m1 - 2*f0 - kk0)) +
                                              kgamma * (f0 + 0.5*kk0 + i2h3 * xx4 * (f4p1 + 0.5*kk4p1 - f4m1 - 0.5*kk4m1) + mkT2h3sq * (f4p1 + 0.5*kk4p1 + f4m1 + 0.5*kk4m1 - 2*f0 - kk0));

                    FF[g1*W1+g2][g3*W2+g4] += KK3[g1*W1+g2][g3*W2+g4] / 3.0;
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-13: CASE 1 KK3) = %lf sec\n", t_1_elapsed);
                t_1_begin = omp_get_wtime();

                // RK4-4
                #pragma omp for private(g1,g2,g3,g4,xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk3p1,kk3m1,kk4p1,kk4m1,kk1p2,kk1m2,kk2p2,kk2m2,kk3p2,kk3m2,kk4p2,kk4m2) schedule(runtime)
                for (int i = 0; i < ExBD.size(); i++)  {

                    g1 = (int)(ExBD[i] / M1);
                    g2 = (int)((ExBD[i] % M1) / M2);
                    g3 = (int)((ExBD[i] % M2) / M3);
                    g4 = (int)(ExBD[i] % M3);
                    xx1 = BoxDynamic[0] + g1 * H[0];
                    xx2 = BoxDynamic[2] + g2 * H[1];
                    xx3 = BoxDynamic[4] + g3 * H[2];
                    xx4 = BoxDynamic[6] + g4 * H[3];
                    f0 = F[g1*W1+g2][g3*W2+g4];
                    f1p1 = F[(g1+1)*W1+g2][g3*W2+g4];
                    f1m1 = F[(g1-1)*W1+g2][g3*W2+g4];
                    f2p1 = F[g1*W1+(g2+1)][g3*W2+g4];
                    f2m1 = F[g1*W1+(g2-1)][g3*W2+g4];
                    f3p1 = F[g1*W1+g2][(g3+1)*W2+g4];
                    f3m1 = F[g1*W1+g2][(g3-1)*W2+g4]; 
                    f4p1 = F[g1*W1+g2][g3*W2+(g4+1)];
                    f4m1 = F[g1*W1+g2][g3*W2+(g4-1)];          
                    f1p2 = F[(g1+2)*W1+g2][g3*W2+g4];
                    f1m2 = F[(g1-2)*W1+g2][g3*W2+g4];
                    f2p2 = F[g1*W1+(g2+2)][g3*W2+g4];
                    f2m2 = F[g1*W1+(g2-2)][g3*W2+g4];
                    f3p2 = F[g1*W1+g2][(g3+2)*W2+g4];
                    f3m2 = F[g1*W1+g2][(g3-2)*W2+g4]; 
                    f4p2 = F[g1*W1+g2][g3*W2+(g4+2)];
                    f4m2 = F[g1*W1+g2][g3*W2+(g4-2)];
                    kk0 = KK3[g1*W1+g2][g3*W2+g4];
                    kk1p1 = KK3[(g1+1)*W1+g2][g3*W2+g4];
                    kk1m1 = KK3[(g1-1)*W1+g2][g3*W2+g4];
                    kk2p1 = KK3[g1*W1+(g2+1)][g3*W2+g4];
                    kk2m1 = KK3[g1*W1+(g2-1)][g3*W2+g4];
                    kk3p1 = KK3[g1*W1+g2][(g3+1)*W2+g4];
                    kk3m1 = KK3[g1*W1+g2][(g3-1)*W2+g4];
                    kk4p1 = KK3[g1*W1+g2][g3*W2+(g4+1)];
                    kk4m1 = KK3[g1*W1+g2][g3*W2+(g4-1)];
                    kk1p2 = KK3[(g1+2)*W1+g2][g3*W2+g4];
                    kk1m2 = KK3[(g1-2)*W1+g2][g3*W2+g4];
                    kk2p2 = KK3[g1*W1+(g2+2)][g3*W2+g4];
                    kk2m2 = KK3[g1*W1+(g2-2)][g3*W2+g4];
                    kk3p2 = KK3[g1*W1+g2][(g3+2)*W2+g4];
                    kk3m2 = KK3[g1*W1+g2][(g3-2)*W2+g4];
                    kk4p2 = KK3[g1*W1+g2][g3*W2+(g4+2)];
                    kk4m2 = KK3[g1*W1+g2][g3*W2+(g4-2)];

                    KK4[g1*W1+g2][g3*W2+g4] = (-kh0m) * xx3 * (-1/12.0*(f1p2+kk1p2) + 2/3.0*(f1p1+kk1p1) - 2/3.0*(f1m1+kk1m1) + 1/12.0*(f1m2+kk1m2)) +
                                              (-kh1m) * xx4 * (-1/12.0*(f2p2+kk2p2) + 2/3.0*(f2p1+kk2p1) - 2/3.0*(f2m1+kk2m1) + 1/12.0*(f2m2+kk2m2)) + 
                                              k2h2 * POTENTIAL_X1(xx1, xx2) * (-1/12.0*(f3p2+kk3p2) + 2/3.0*(f3p1+kk3p1) - 2/3.0*(f3m1+kk3m1) + 1/12.0*(f3m2+kk3m2)) +
                                              k2h3 * POTENTIAL_X2(xx1, xx2) * (-1/12.0*(f4p2+kk4p2) + 2/3.0*(f4p1+kk4p1) - 2/3.0*(f4m1+kk4m1) + 1/12.0*(f4m2+kk4m2)) +
                                              (-khbsq2h2) * POTENTIAL_XXX1(xx1, xx2) * (0.5*(f3p2+kk3p2) - (f3p1+kk3p1) + (f3m1+kk3m1) - 0.5*(f3m2+kk3m2)) + 
                                              (-khbsq2h3) * POTENTIAL_XXX2(xx1, xx2) * (0.5*(f4p2+kk4p2) - (f4p1+kk4p1) + (f4m1+kk4m1) - 0.5*(f4m2+kk4m2)) + 
                                              kgamma * (f0 + kk0 + i2h2 * xx3 * (f3p1 + kk3p1 - f3m1 - kk3m1) + mkT2h2sq * (f3p1 + kk3p1 + f3m1 + kk3m1 - 2*f0 - 2*kk0)) +
                                              kgamma * (f0 + kk0 + i2h3 * xx4 * (f4p1 + kk4p1 - f4m1 - kk4m1) + mkT2h3sq * (f4p1 + kk4p1 + f4m1 + kk4m1 - 2*f0 - 2*kk0));

                    FF[g1*W1+g2][g3*W2+g4] += KK4[g1*W1+g2][g3*W2+g4] / 6.0;

                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-14: CASE 1 KK4) = %lf sec\n", t_1_elapsed);
                t_1_begin = omp_get_wtime();
            }

            // Check Multiple Expanding 
            // TBL = index of FF that FF(TBL) is higher than TolL

            t_1_begin = omp_get_wtime();
            TBL.clear();

            #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4,b1,b2,b3,b4,b5,b6,b7)
            for (unsigned long int i = 0; i < ExFF.size(); i++)
            {
                g1 = (int)(ExFF[i] / M1);
                g2 = (int)((ExFF[i] % M1) / M2);
                g3 = (int)((ExFF[i] % M2) / M3);
                g4 = (int)(ExFF[i] % M3); 
                b1 = std::abs(FF[g1*W1+g2][g3*W2+g4]) >= TolH;
                b2 = std::abs(FF[(g1+1)*W1+g2][g3*W2+g4] - FF[(g1-1)*W1+g2][g3*W2+g4]) >= TolHdX1;
                b3 = std::abs(FF[g1*W1+(g2+1)][g3*W2+g4] - FF[g1*W1+(g2-1)][g3*W2+g4]) >= TolHdX2;
                b4 = std::abs(FF[g1*W1+g2][(g3+1)*W2+g4] - FF[g1*W1+g2][(g3-1)*W2+g4]) >= TolHdX3;
                b5 = std::abs(FF[g1*W1+g2][g3*W2+(g4+1)] - FF[g1*W1+g2][g3*W2+(g4-1)]) >= TolHdX4;
                b6 = g1 > 2 && g2 > 2 && g3 > 2 && g4 > 2;
                b7 = g1 < BoxShapeFine[0] - 3 && g2 < BoxShapeFine[1] - 3 && g3 < BoxShapeFine[2] - 3 && g4 < BoxShapeFine[3] - 3;

                if (  ( b1 || b2 || b3 || b4 || b5 ) && b6 && b7 )  {
                    tmpVec.push_back(ExFF[i]);
                }
            }        
            tmpVec.swap(TBL);
            tmpVec.clear(); 

            // TBL & TBL_P set difference
            tmpVec.resize(TBL_P.size() + TBL.size());
            __gnu_parallel::sort(TBL.begin(), TBL.end());
            __gnu_parallel::sort(TBL_P.begin(), TBL_P.end());
            it=std::set_difference( TBL.begin(), TBL.end(), TBL_P.begin(), TBL_P.end(), tmpVec.begin() );
            tmpVec.resize(it - tmpVec.begin()); 
            tmpVec.swap(TBL);
            tmpVec.clear();

            // Combine TBL and TBL_P
            TBL_P.reserve(TBL_P.size() + TBL.size());
            TBL_P.insert(TBL_P.end(), TBL.begin(), TBL.end());

            // Find unique elements
            __gnu_parallel::sort(TBL_P.begin(),TBL_P.end());
            it = std::unique (TBL_P.begin(), TBL_P.end()); 
            TBL_P.resize(std::distance(TBL_P.begin(),it));

            // Update isFirstExtrp
            isFirstExtrp = (TBL.size() == 0) ? 1 : 0;
            Excount += 1;
            if (Excount == ExLimit) TBL.clear();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;

            if (!QUIET && TIMING) log->log("Elapsed time (omp-c-3 CASE 1 TBL) = %lf sec\n", t_1_elapsed); 
        }
        // .........................................................................................

        // CASE 2: Truncating without extrapolation

        if ( !isExtrapolate )
        {
            // RK4-1
            #pragma omp parallel
            {
                #pragma omp single nowait
                {
                    t_1_begin = omp_get_wtime();
                }
                #pragma omp for private(xx1,xx2,xx3,xx4,f0,f1p1,f1m1,f2p1,f2m1,f3p1,f3m1,f4p1,f4m1,f1p2,f1m2,f2p2,f2m2,f3p2,f3m2,f4p2,f4m2) schedule(runtime)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                if ( TAMask[i1*W1+i2][i3*W2+i4] )  {

                                    xx1 = BoxDynamic[0] + i1 * H[0];
                                    xx2 = BoxDynamic[2] + i2 * H[1];
                                    xx3 = BoxDynamic[4] + i3 * H[2];
                                    xx4 = BoxDynamic[6] + i4 * H[3];
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
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                if ( TAMask[i1*W1+i2][i3*W2+i4] )  {
                                    
                                    xx1 = BoxDynamic[0] + i1 * H[0];
                                    xx2 = BoxDynamic[2] + i2 * H[1];
                                    xx3 = BoxDynamic[4] + i3 * H[2];
                                    xx4 = BoxDynamic[6] + i4 * H[3];
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
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                if ( TAMask[i1*W1+i2][i3*W2+i4] )  {

                                    xx1 = BoxDynamic[0] + i1 * H[0];
                                    xx2 = BoxDynamic[2] + i2 * H[1];
                                    xx3 = BoxDynamic[4] + i3 * H[2];
                                    xx4 = BoxDynamic[6] + i4 * H[3];
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
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                if ( TAMask[i1*W1+i2][i3*W2+i4] )  {

                                    xx1 = BoxDynamic[0] + i1 * H[0];
                                    xx2 = BoxDynamic[2] + i2 * H[1];
                                    xx3 = BoxDynamic[4] + i3 * H[2];
                                    xx4 = BoxDynamic[6] + i4 * H[3];
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

                                    KK4[i1][i2] = (-kh0m) * xx3 * (-1/12.0*(f1p2+kk1p2) + 2/3.0*(f1p1+kk1p1) - 2/3.0*(f1m1+kk1m1) + 1/12.0*(f1m2+kk1m2)) +
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
        for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
            for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                    for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                        if (TAMask[i1*W1+i2][i3*W2+i4])
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
        for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
            for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                    for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                        val = norm * FF[i1*W1+i2][i3*W2+i4] * TAMask[i1*W1+i2][i3*W2+i4];
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
                idx_x0 = (int)std::round(std::round((trans_x0-BoxDynamic[0])/H[0]));
                pftrans = 0.0;

                #pragma omp parallel for reduction (+:pftrans)
                for (int i1 = idx_x0; i1 < BoxShapeFine[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShapeFine[1]; i2 ++)  {
                        for (int i3 = 0; i3 < BoxShapeFine[2]; i3 ++)  {
                            for (int i4 = 0; i4 < BoxShapeFine[3]; i4 ++)  {
                                pftrans += PF[i1*W1+i2][i3*W2+i4] * TAMask[i1*W1+i2][i3*W2+i4];
                            }
                        }
                    }
                }
                pftrans *= H[0] * H[1] * H[2] * H[3];
                PF_trans.push_back(pftrans);
                log->log("[CaldeiraLeggett4d] trans_x0 = %lf\n", idx_x0 * H[0] + BoxDynamic[0]);
                log->log("[CaldeiraLeggett4d] Time %lf, Trans = %.16e\n", ( tt + 1 ) * kk, pftrans);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (omp-x-2 trans) = %lf sec\n", t_1_elapsed); 
            }

            if (isCorr)  {

                Ft.resize(O1);

                #pragma omp parallel for
                for (unsigned long int i1 = 0; i1 < O1; i1 ++)  {
                    Ft[i1] = 0.0;
                }

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

                // Find out the overlapped region
                x1_min_t = BoxDynamic[0];
                x1_max_t = BoxDynamic[1];
                x2_min_t = BoxDynamic[2];
                x2_max_t = BoxDynamic[3];

                ibox_x1_min = std::max(x1_min_t,x1_min_0);
                ibox_x1_max = std::min(x1_max_t,x1_max_0);
                ibox_x2_min = std::max(x2_min_t,x2_min_0);
                ibox_x2_max = std::min(x2_max_t,x2_max_0);
                  
                ibox_shape_1 = (int)std::round((ibox_x1_max - ibox_x1_min)/H[0]) + 1;
                ibox_shape_2 = (int)std::round((ibox_x2_max - ibox_x2_min)/H[1]) + 1;
   
                shift_x1_0 = (int)std::round((ibox_x1_min - x1_min_0) / H[0]);
                shift_x1_t = (int)std::round((ibox_x1_min - x1_min_t) / H[0]);
                shift_x2_0 = (int)std::round((ibox_x2_min - x2_min_0) / H[1]);
                shift_x2_t = (int)std::round((ibox_x2_min - x2_min_t) / H[1]);

                // Compute correlation
                corr = 0.0;

                for (int i1 = 0; i1 < ibox_shape_1; i1++)  {
                    for (int i2 = 0; i2 < ibox_shape_2; i2++)  {
                        corr += F0[(i1+shift_x1_0)*W1_0+(i2+shift_x2_0)] * Ft[(i1+shift_x1_t)*W1+(i2+shift_x2_t)];
                    }
                }
                corr *= H[0] * H[1];
                log->log("[CaldeiraLeggett4d] Time %lf, Corr = %.16e\n", ( tt + 1 ) * kk, corr/corr_0);
            }
        }

        // Truncation and TA

        tmpVec.clear();

        t_1_begin = omp_get_wtime();
        #pragma omp parallel 
        {
            #pragma omp for reduction(merge: tmpVec) private(b1, b2, b3, b4, b5)
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                        for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                            if ( TAMask[i1*W1+i2][i3*W2+i4] )  {
                                b1 = std::abs(PF[i1*W1+i2][i3*W2+i4]) < TolH;
                                b2 = std::abs(PF[(i1+1)*W1+i2][i3*W2+i4]*TAMask[(i1+1)*W1+i2][i3*W2+i4] - PF[(i1-1)*W1+i2][i3*W2+i4]*TAMask[(i1-1)*W1+i2][i3*W2+i4]) < TolHdX1;
                                b3 = std::abs(PF[i1*W1+(i2+1)][i3*W2+i4]*TAMask[i1*W1+(i2+1)][i3*W2+i4] - PF[i1*W1+(i2-1)][i3*W2+i4]*TAMask[i1*W1+(i2-1)][i3*W2+i4]) < TolHdX2;
                                b4 = std::abs(PF[i1*W1+i2][(i3+1)*W2+i4]*TAMask[i1*W1+i2][(i3+1)*W2+i4] - PF[i1*W1+i2][(i3-1)*W2+i4]*TAMask[i1*W1+i2][(i3-1)*W2+i4]) < TolHdX3;
                                b5 = std::abs(PF[i1*W1+i2][i3*W2+(i4+1)]*TAMask[i1*W1+i2][i3*W2+(i4+1)] - PF[i1*W1+i2][i3*W2+(i4-1)]*TAMask[i1*W1+i2][i3*W2+(i4-1)]) < TolHdX4;

                                if (b1 && b2 && b3 && b4 && b5) {
                                    F[i1*W1+i2][i3*W2+i4] = 0.0;
                                    tmpVec.push_back(i1*M1+i2*M2+i3*M3+i4);
                                }
                            }
                        }
                    }
                }
            }

            #pragma omp for private(g1,g2,g3,g4)
            for (int i = 0; i < tmpVec.size(); i++)  {
                g1 = (int)(tmpVec[i] / M1);
                g2 = (int)((tmpVec[i] % M1) / M2);
                g3 = (int)((tmpVec[i] % M2) / M3);
                g4 = (int)(tmpVec[i] % M3); 
                TAMask[g1*W1+g2][g3*W2+g4] = 0;
                PF[g1*W1+g2][g3*W2+g4] = 0.0;
            }
        }
        tmpVec.clear();

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3 TA) = %lf sec\n", t_1_elapsed);
        t_1_begin = omp_get_wtime();

        // Rebuild TA box
        #pragma omp parallel for reduction(min: x1_min, x2_min, x3_min, x4_min) reduction(max: x1_max, x2_max, x3_max, x4_max)
        for (int i1 = 1; i1 < BoxShapeFine[0] - 1; i1 ++)  {
            for (int i2 = 1; i2 < BoxShapeFine[1] - 1; i2 ++)  {
                for (int i3 = 1; i3 < BoxShapeFine[2] - 1; i3 ++)  {
                    for (int i4 = 1; i4 < BoxShapeFine[3] - 1; i4 ++)  {
                        if (TAMask[i1*W1+i2][i3*W2+i4])  {
                            if (i1 < x1_min)  x1_min = i1;
                            if (i1 > x1_max)  x1_max = i1;
                            if (i2 < x2_min)  x2_min = i2;
                            if (i2 > x2_max)  x2_max = i2;
                            if (i3 < x3_min)  x3_min = i3;
                            if (i3 > x3_max)  x3_max = i3;
                            if (i4 < x4_min)  x4_min = i4;
                            if (i4 > x4_max)  x4_max = i4;
                        }
                    }
                }
            }
        }
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin; 
        if (!QUIET && TIMING) log->log("Elapsed time (omp-y-1 ta) = %lf sec\n", t_1_elapsed); 

        // Check if TA enters the skin region
        isEnterSkin = x1_min < skin || x2_min < skin || x3_min < skin || x4_min < skin || BoxShapeFine[0]-1-x1_max < skin || BoxShapeFine[1]-1-x2_max < skin || BoxShapeFine[2]-1-x3_max < skin || BoxShapeFine[3]-1-x4_max < skin;

        // Rebuild adaptive box
        if ( isEnterSkin )  {

            t_1_begin = omp_get_wtime();

            // Check if the adaptive box touches the boundary
            b1 = BoxDynamic[0] + (x1_min - 2 * skin) * H[0] < BoxBoundary[0];
            b2 = BoxDynamic[2] + (x2_min - 2 * skin) * H[1] < BoxBoundary[2];
            b3 = BoxDynamic[4] + (x3_min - 2 * skin) * H[2] < BoxBoundary[4];
            b4 = BoxDynamic[6] + (x4_min - 2 * skin) * H[3] < BoxBoundary[6];
            b5 = BoxDynamic[0] + (x1_max + 2 * skin) * H[0] > BoxBoundary[1];
            b6 = BoxDynamic[2] + (x2_max + 2 * skin) * H[1] > BoxBoundary[3];
            b7 = BoxDynamic[4] + (x3_max + 2 * skin) * H[2] > BoxBoundary[5];
            b8 = BoxDynamic[6] + (x4_max + 2 * skin) * H[3] > BoxBoundary[7];

            isTouchedBoundary = b1 || b2 || b3 || b4 || b5 || b6 || b7 || b8;

            DistToEdge[0] = b1 ? (int)(std::round( (BoxDynamic[0] + x1_min * H[0] - BoxBoundary[0]) / H[0])) : 2 * skin;
            DistToEdge[2] = b2 ? (int)(std::round( (BoxDynamic[2] + x2_min * H[1] - BoxBoundary[2]) / H[1])) : 2 * skin;
            DistToEdge[4] = b3 ? (int)(std::round( (BoxDynamic[4] + x3_min * H[2] - BoxBoundary[4]) / H[2])) : 2 * skin;
            DistToEdge[6] = b4 ? (int)(std::round( (BoxDynamic[6] + x4_min * H[3] - BoxBoundary[6]) / H[3])) : 2 * skin;
            DistToEdge[1] = b5 ? (int)(std::round(-(BoxDynamic[0] + x1_max * H[0] - BoxBoundary[1]) / H[0])) : 2 * skin;
            DistToEdge[3] = b6 ? (int)(std::round(-(BoxDynamic[2] + x2_max * H[1] - BoxBoundary[3]) / H[1])) : 2 * skin;
            DistToEdge[5] = b7 ? (int)(std::round(-(BoxDynamic[4] + x3_max * H[2] - BoxBoundary[5]) / H[2])) : 2 * skin;
            DistToEdge[7] = b8 ? (int)(std::round(-(BoxDynamic[6] + x4_max * H[3] - BoxBoundary[7]) / H[3])) : 2 * skin;

            BoxDynamicNew[0] = b1 ? BoxBoundary[0] : BoxDynamic[0] + (x1_min - 2 * skin) * H[0];
            BoxDynamicNew[2] = b2 ? BoxBoundary[2] : BoxDynamic[2] + (x2_min - 2 * skin) * H[1];
            BoxDynamicNew[4] = b3 ? BoxBoundary[4] : BoxDynamic[4] + (x3_min - 2 * skin) * H[2];
            BoxDynamicNew[6] = b4 ? BoxBoundary[6] : BoxDynamic[6] + (x4_min - 2 * skin) * H[3];
            BoxDynamicNew[1] = b5 ? BoxBoundary[1] : BoxDynamic[0] + (x1_max + 2 * skin) * H[0];
            BoxDynamicNew[3] = b6 ? BoxBoundary[3] : BoxDynamic[2] + (x2_max + 2 * skin) * H[1];
            BoxDynamicNew[5] = b7 ? BoxBoundary[5] : BoxDynamic[4] + (x3_max + 2 * skin) * H[2];
            BoxDynamicNew[7] = b8 ? BoxBoundary[7] : BoxDynamic[6] + (x4_max + 2 * skin) * H[3];

            b1 = BoxDynamicNew[0] != BoxDynamic[0]; 
            b2 = BoxDynamicNew[1] != BoxDynamic[1];
            b3 = BoxDynamicNew[2] != BoxDynamic[2]; 
            b4 = BoxDynamicNew[3] != BoxDynamic[3];
            b5 = BoxDynamicNew[4] != BoxDynamic[4]; 
            b6 = BoxDynamicNew[5] != BoxDynamic[5];
            b7 = BoxDynamicNew[6] != BoxDynamic[6]; 
            b8 = BoxDynamicNew[7] != BoxDynamic[7];                      

            isAdaptiveBoxChanged = b1 || b2 || b3 || b4 || b5 || b6 || b7 || b8;

            if ( isAdaptiveBoxChanged )  {

                if ( isTouchedBoundary )
                    log->log("[CaldeiraLeggett4d] Adaptive box touches boundary\n");

                log->log("[CaldeiraLeggett4d] BoxDynamic: [%.6f,%.6f][%.6f,%.6f][%.6f,%.6f][%.6f,%.6f]\n",BoxDynamic[0],BoxDynamic[1],BoxDynamic[2],BoxDynamic[3],BoxDynamic[4],BoxDynamic[5],BoxDynamic[6],BoxDynamic[7]);
                log->log("[CaldeiraLeggett4d] BoxDynamicNew: [%.6f,%.6f][%.6f,%.6f][%.6f,%.6f][%.6f,%.6f]\n",BoxDynamicNew[0],BoxDynamicNew[1],BoxDynamicNew[2],BoxDynamicNew[3],BoxDynamicNew[4],BoxDynamicNew[5],BoxDynamicNew[6],BoxDynamicNew[7]);
                log->log("[CaldeiraLeggett4d] New fine grid size = (");

                GRIDS_TOT = 1;
                for (int i = 0; i < DIMENSIONS; i ++)  {
                    BoxShapeFineNew[i] = (int)std::round((BoxDynamicNew[2*i+1]-BoxDynamicNew[2*i])/H[i]) + 1;
                    GRIDS_TOT *= BoxShapeFineNew[i];

                    if ( i < DIMENSIONS - 1 )
                        log->log("%d, ", BoxShapeFineNew[i]);
                    else
                        log->log("%d)=", BoxShapeFineNew[i]);
                }
                log->log("%lu\n", GRIDS_TOT);

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (omp-z-1 rebuild f) = %lf sec\n", t_1_elapsed); 

                // Rebuilt TAMask
                t_1_begin = omp_get_wtime();
                #pragma omp parallel for reduction(merge: tmpVecFull) reduction(merge: tmpVecD)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {    
                                if ( TAMask[i1*W1+i2][i3*W2+i4] )  {
                                    tmpVecFull.push_back({i1,i2,i3,i4});
                                    tmpVecD.push_back(FF[i1*W1+i2][i3*W2+i4]);
                                }
                            }
                        }
                    }
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (omp-z-3 rebuild tamask) = %lf sec\n", t_1_elapsed); 

                // Rebuild FF, PF, KK1, and KK2
                t_1_begin = omp_get_wtime();

                for (unsigned long int i = 0; i < O1; i++)  {
                    delete F[i];
                    delete FF[i];
                    delete PF[i];
                    delete KK1[i];
                    delete KK2[i];
                    delete KK3[i];
                    delete KK4[i];
                    delete TAMask[i];                
                }
                delete F;
                delete FF;
                delete PF;
                delete KK1;
                delete KK2;
                delete KK3;
                delete KK4;
                delete TAMask;

                M1 = (unsigned long int)(BoxShapeFineNew[1] * BoxShapeFineNew[2] * BoxShapeFineNew[3]);
                M2 = (unsigned long int)(BoxShapeFineNew[2] * BoxShapeFineNew[3]);
                M3 = (unsigned long int)(BoxShapeFineNew[3]);
                W1 = (unsigned long int)(BoxShapeFineNew[1]);
                W2 = (unsigned long int)(BoxShapeFineNew[3]);
                O1 = (unsigned long int)(BoxShapeFineNew[0] * BoxShapeFineNew[1]);
                O2 = (unsigned long int)(BoxShapeFineNew[2] * BoxShapeFineNew[3]);

                TAMask = new bool*[O1];

                for (unsigned long int i = 0; i < O1; i++)  {
                    TAMask[i] = new(nothrow) bool[O2];

                    if (!TAMask[i])
                        log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "TAMask");
                }

                KK1 = new double*[O1];

                for (unsigned long int i = 0; i < O1; i++)  {
                    KK1[i] = new(nothrow) double[O2];

                    if (!KK1[i])
                        log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "KK1");
                }

                KK2 = new double*[O1];

                for (unsigned long int i = 0; i < O1; i++)  {
                    KK2[i] = new(nothrow) double[O2];

                    if (!KK2[i])
                        log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "KK2");
                }

                KK3 = new double*[O1];

                for (unsigned long int i = 0; i < O1; i++)  {
                    KK3[i] = new(nothrow) double[O2];

                    if (!KK3[i])
                        log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "KK3");
                }

                KK4 = new double*[O1];

                for (unsigned long int i = 0; i < O1; i++)  {
                    KK4[i] = new(nothrow) double[O2];

                    if (!KK4[i])
                        log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "KK4");
                }

                F = new double*[O1];

                for (unsigned long int i = 0; i < O1; i++)  {
                    F[i] = new(nothrow) double[O2];

                    if (!F[i])
                        log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "F");
                }

                FF = new double*[O1];

                for (unsigned long int i = 0; i < O1; i++)  {
                    FF[i] = new(nothrow) double[O2];

                    if (!FF[i])
                        log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "FF");
                }

                PF = new double*[O1];

                for (unsigned long int i = 0; i < O1; i++)  {
                    PF[i] = new(nothrow) double[O2];

                    if (!PF[i])
                        log->log("[CaldeiraLeggett4d] Memory allocation %s failed\n", "PF");
                }

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (omp-z-3 re-allocation memory) = %lf sec\n", t_1_elapsed); 

                // Rebuild arrays
                t_1_begin = omp_get_wtime();

                #pragma omp parallel for
                for (int i1 = 0; i1 < BoxShapeFineNew[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShapeFineNew[1]; i2 ++)  {
                        for (int i3 = 0; i3 < BoxShapeFineNew[2]; i3 ++)  {
                            for (int i4 = 0; i4 < BoxShapeFineNew[3]; i4 ++)  {   
                                F[i1*W1+i2][i3*W2+i4] = 0.0;         
                                PF[i1*W1+i2][i3*W2+i4] = 0.0;
                                FF[i1*W1+i2][i3*W2+i4] = 0.0;
                                KK1[i1*W1+i2][i3*W2+i4] = 0.0;
                                KK2[i1*W1+i2][i3*W2+i4] = 0.0;
                                KK3[i1*W1+i2][i3*W2+i4] = 0.0;
                                KK4[i1*W1+i2][i3*W2+i4] = 0.0;
                                TAMask[i1*W1+i2][i3*W2+i4] = 0;
                            }
                        }
                    }
                } 

                #pragma omp parallel for private(g1, g2, g3, g4)
                for (unsigned long int i = 0; i < tmpVecFull.size(); i ++)  {

                    g1 = (int)(tmpVecFull[i][0] - x1_min + DistToEdge[0]);
                    g2 = (int)(tmpVecFull[i][1] - x2_min + DistToEdge[2]);
                    g3 = (int)(tmpVecFull[i][2] - x3_min + DistToEdge[4]);
                    g4 = (int)(tmpVecFull[i][3] - x4_min + DistToEdge[6]);

                    TAMask[g1*W1+g2][g3*W2+g4] = 1;
                    PF[g1*W1+g2][g3*W2+g4]=tmpVecD[i];
                    F[g1*W1+g2][g3*W2+g4]=tmpVecD[i];
                }
                // Free tmpVecFull
                tmpVecFull.clear();
                tmpVecD.clear();

                // Update TA x_max and x_min
                x1_max = x1_max - x1_min + DistToEdge[0];
                x2_max = x2_max - x2_min + DistToEdge[2];
                x3_max = x3_max - x3_min + DistToEdge[4];
                x4_max = x4_max - x4_min + DistToEdge[6];
                x1_min = DistToEdge[0];
                x2_min = DistToEdge[2];
                x3_min = DistToEdge[4];
                x4_min = DistToEdge[6];

                // Update Box information
                BoxShapeFine = BoxShapeFineNew;
                BoxDynamic = BoxDynamicNew;                  

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (omp-z-4 rebuild arrays) = %lf sec\n", t_1_elapsed); 

            } // isAdaptiveBoxChange
        } // isEnterSkin

        // Applying absorbing bounday conditions

        if ( isTouchedBoundary )  {

            // X1_MIN
            if ( (BoxDynamic[0] == BoxBoundary[0]) && x1_min < 3 )  {

                t_1_begin = omp_get_wtime();
                log->log("[CaldeiraLeggett4d] x1_min to boundary is less than 3\n");
                log->log("[CaldeiraLeggett4d] BoxBoundary[0] = %lf, dist = %d\n", BoxBoundary[0], x1_min);

                #pragma omp parallel for
                for (int i1 = 0; i1 < 3; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                if (TAMask[i1*W1+i2][i3*W2+i4])  {
                                    TAMask[i1*W1+i2][i3*W2+i4] = 0;
                                    F[i1*W1+i2][i3*W2+i4] = 0.0;
                                    PF[i1*W1+i2][i3*W2+i4] = 0.0;
                                }
                            }
                        }
                    }
                }
                x1_min = 3;
                log->log("[CaldeiraLeggett4d] New x1_min = %d\n", x1_min);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (abs-1) = %lf sec\n", t_1_elapsed); 
            }

            // X1_MAX
            if ( (BoxDynamic[1] == BoxBoundary[1]) && (int)(std::round(BoxBoundary[1] - (BoxDynamic[0] + x1_max * H[0]))/H[0]) < 3)  {

                t_1_begin = omp_get_wtime();
                log->log("[CaldeiraLeggett4d] x1_max to boundary is less than 3\n");
                log->log("[CaldeiraLeggett4d] BoxBoundary[1] = %lf, dist = %d\n",BoxBoundary[1], (int)(std::round(BoxBoundary[1] - (BoxDynamic[0] + x1_max * H[0]))/H[0]));

                #pragma omp parallel for
                for (int i1 = BoxShapeFine[0] - 3; i1 < BoxShapeFine[0]; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                if (TAMask[i1*W1+i2][i3*W2+i4])  {
                                    TAMask[i1*W1+i2][i3*W2+i4] = 0;
                                    F[i1*W1+i2][i3*W2+i4] = 0.0;
                                    PF[i1*W1+i2][i3*W2+i4] = 0.0;
                                }
                            }
                        }
                    }
                }
                x1_max = BoxShapeFine[0] - 4;
                log->log("[CaldeiraLeggett4d] New x1_max = %d\n", x1_max);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (abs-2) = %lf sec\n", t_1_elapsed); 
            }

            // X2_MIN
            if ( (BoxDynamic[2] == BoxBoundary[2]) && x2_min < 3 )  {

                t_1_begin = omp_get_wtime();
                log->log("[CaldeiraLeggett4d] x2_min to boundary is less than 3\n");
                log->log("[CaldeiraLeggett4d] BoxBoundary[2] = %lf, dist = %d\n", BoxBoundary[2], x2_min);

                #pragma omp parallel for
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = 0; i2 < 3; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                if (TAMask[i1*W1+i2][i3*W2+i4])  {
                                    TAMask[i1*W1+i2][i3*W2+i4] = 0;
                                    F[i1*W1+i2][i3*W2+i4] = 0.0;
                                    PF[i1*W1+i2][i3*W2+i4] = 0.0;
                                }
                            }
                        }
                    }
                }
                x2_min = 3;
                log->log("[CaldeiraLeggett4d] New x2_min = %d\n", x2_min);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (abs-3) = %lf sec\n", t_1_elapsed); 
            }

            // X2_MAX
            if ( (BoxDynamic[3] == BoxBoundary[3]) && (int)(std::round(BoxBoundary[3] - (BoxDynamic[2] + x2_max * H[1]))/H[1]) < 3)  {

                t_1_begin = omp_get_wtime();
                log->log("[CaldeiraLeggett4d] x2_max to boundary is less than 3\n");
                log->log("[CaldeiraLeggett4d] BoxBoundary[3] = %lf, dist = %d\n",BoxBoundary[3], (int)(std::round(BoxBoundary[3] - (BoxDynamic[2] + x2_max * H[1]))/H[1]));

                #pragma omp parallel for
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = BoxShapeFine[1] - 3; i2 < BoxShapeFine[1]; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                if (TAMask[i1*W1+i2][i3*W2+i4])  {
                                    TAMask[i1*W1+i2][i3*W2+i4] = 0;
                                    F[i1*W1+i2][i3*W2+i4] = 0.0;
                                    PF[i1*W1+i2][i3*W2+i4] = 0.0;
                                }
                            }
                        }
                    }
                }
                x2_max = BoxShapeFine[1] - 4;
                log->log("[CaldeiraLeggett4d] New x2_max = %d\n", x2_max);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (abs-4) = %lf sec\n", t_1_elapsed); 
            }

            // X3_MIN
            if ( (BoxDynamic[4] == BoxBoundary[4]) && x3_min < 3 )  {

                t_1_begin = omp_get_wtime();
                log->log("[CaldeiraLeggett4d] x3_min to boundary is less than 3\n");
                log->log("[CaldeiraLeggett4d] BoxBoundary[4] = %lf, dist = %d\n", BoxBoundary[4], x3_min);

                #pragma omp parallel for
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = 0; i3 < 3; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                if (TAMask[i1*W1+i2][i3*W2+i4])  {
                                    TAMask[i1*W1+i2][i3*W2+i4] = 0;
                                    F[i1*W1+i2][i3*W2+i4] = 0.0;
                                    PF[i1*W1+i2][i3*W2+i4] = 0.0;
                                }
                            }
                        }
                    }
                }
                x3_min = 3;
                log->log("[CaldeiraLeggett4d] New x3_min = %d\n", x3_min);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (abs-5) = %lf sec\n", t_1_elapsed); 
            }

            // X3_MAX
            if ( (BoxDynamic[5] == BoxBoundary[5]) && (int)(std::round(BoxBoundary[5] - (BoxDynamic[4] + x3_max * H[2]))/H[2]) < 3)  {

                t_1_begin = omp_get_wtime();
                log->log("[CaldeiraLeggett4d] x3_max to boundary is less than 3\n");
                log->log("[CaldeiraLeggett4d] BoxBoundary[5] = %lf, dist = %d\n", BoxBoundary[5], (int)(std::round(BoxBoundary[5] - (BoxDynamic[4] + x3_max * H[2]))/H[2]));

                #pragma omp parallel for
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = BoxShapeFine[2] - 3; i3 < BoxShapeFine[2]; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                if (TAMask[i1*W1+i2][i3*W2+i4])  {
                                    TAMask[i1*W1+i2][i3*W2+i4] = 0;
                                    F[i1*W1+i2][i3*W2+i4] = 0.0;
                                    PF[i1*W1+i2][i3*W2+i4] = 0.0;
                                }
                            }
                        }
                    }
                }
                x3_max = BoxShapeFine[2] - 4;
                log->log("[CaldeiraLeggett4d] New x3_max = %d\n", x3_max);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (abs-6) = %lf sec\n", t_1_elapsed); 
            }

            // X4_MIN
            if ( (BoxDynamic[6] == BoxBoundary[6]) && x4_min < 3 )  {

                t_1_begin = omp_get_wtime();
                log->log("[CaldeiraLeggett4d] x4_min to boundary is less than 3\n");
                log->log("[CaldeiraLeggett4d] BoxBoundary[6] = %lf, dist = %d\n", BoxBoundary[6], x4_min);

                #pragma omp parallel for
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = 0; i4 < 3; i4 ++)  {
                                if (TAMask[i1*W1+i2][i3*W2+i4])  {
                                    TAMask[i1*W1+i2][i3*W2+i4] = 0;
                                    F[i1*W1+i2][i3*W2+i4] = 0.0;
                                    PF[i1*W1+i2][i3*W2+i4] = 0.0;
                                }
                            }
                        }
                    }
                }
                x4_min = 3;
                log->log("[CaldeiraLeggett4d] New x4_min = %d\n", x4_min);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (abs-7) = %lf sec\n", t_1_elapsed); 
            }

            // X4_MAX
            if ( (BoxDynamic[7] == BoxBoundary[7]) && (int)(std::round(BoxBoundary[7] - (BoxDynamic[6] + x4_max * H[3]))/H[3]) < 3)  {

                t_1_begin = omp_get_wtime();
                log->log("[CaldeiraLeggett4d] x4_max to boundary is less than 3\n");
                log->log("[CaldeiraLeggett4d] BoxBoundary[7] = %lf, dist = %d\n", BoxBoundary[7], (int)(std::round(BoxBoundary[7] - (BoxDynamic[6] + x4_max * H[3]))/H[3]));

                #pragma omp parallel for
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = BoxShapeFine[3] - 3; i4 < BoxShapeFine[3]; i4 ++)  {
                                if (TAMask[i1*W1+i2][i3*W2+i4])  {
                                    TAMask[i1*W1+i2][i3*W2+i4] = 0;
                                    F[i1*W1+i2][i3*W2+i4] = 0.0;
                                    PF[i1*W1+i2][i3*W2+i4] = 0.0;
                                }
                            }
                        }
                    }
                }
                x4_max = BoxShapeFine[3] - 4;
                log->log("[CaldeiraLeggett4d] New x4_max = %d\n", x4_max);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (abs-8) = %lf sec\n", t_1_elapsed); 
            }
        } // isTouchedBoundary

        // TB
        t_1_begin = omp_get_wtime();
        tmpVec.clear();

        #pragma omp parallel for reduction(merge: tmpVec)
        for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
            for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                    for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                        if (TAMask[i1*W1+i2][i3*W2+i4])  {
                            if (!TAMask[(i1-1)*W1+i2][i3*W2+i4] || !TAMask[(i1+1)*W1+i2][i3*W2+i4] || !TAMask[i1*W1+(i2-1)][i3*W2+i4] || !TAMask[i1*W1+(i2+1)][i3*W2+i4] || !TAMask[i1*W1+i2][(i3-1)*W2+i4] || !TAMask[i1*W1+i2][(i3+1)*W2+i4] || !TAMask[i1*W1+i2][i3*W2+(i4-1)] || !TAMask[i1*W1+i2][i3*W2+(i4+1)])
                                tmpVec.push_back(i1*M1+i2*M2+i3*M3+i4);
                        }
                    }
                }
            }
        }
        tmpVec.swap(TB);
        tmpVec.clear();

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-5 TB) = %lf sec\n", t_1_elapsed);
        t_1_begin = omp_get_wtime();

        // `````````````````````````````````````````````````````````````````

        // TA expansion
        #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, g3, g4)
        for (unsigned long int i = 0; i < TB.size(); i++)
        {
            g1 = (int)(TB[i] / M1);
            g2 = (int)((TB[i] % M1) / M2);
            g3 = (int)((TB[i] % M2) / M3);
            g4 = (int)(TB[i] % M3); 

            if ( g1+1 != BoxShapeFine[0]-3 && !TAMask[(g1+1)*W1+g2][g3*W2+g4])
                tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
            if ( g1-1 != 2 && !TAMask[(g1-1)*W1+g2][g3*W2+g4])
                tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));
            if ( g2+1 != BoxShapeFine[1]-3 && !TAMask[g1*W1+(g2+1)][g3*W2+g4])
                tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
            if ( g2-1 != 2 && !TAMask[g1*W1+(g2-1)][g3*W2+g4])
                tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));
            if ( g3+1 != BoxShapeFine[2]-3 && !TAMask[g1*W1+g2][(g3+1)*W2+g4])
                tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
            if ( g3-1 != 2 && !TAMask[g1*W1+g2][(g3-1)*W2+g4])
                tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));
            if ( g4+1 != BoxShapeFine[3]-3 && !TAMask[g1*W1+g2][g3*W2+(g4+1)])
                tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
            if ( g4-1 != 2 && !TAMask[g1*W1+g2][g3*W2+(g4-1)])
                tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
        }

        #pragma omp parallel for reduction(min: x1_min, x2_min, x3_min, x4_min) reduction(max: x1_max, x2_max, x3_max, x4_max) private(g1, g2, g3, g4)
        for (unsigned long int i = 0; i < tmpVec.size(); i ++)  {

            g1 = (int)(tmpVec[i] / M1);
            g2 = (int)((tmpVec[i] % M1) / M2);
            g3 = (int)((tmpVec[i] % M2) / M3);
            g4 = (int)(tmpVec[i] % M3); 
            TAMask[g1*W1+g2][g3*W2+g4] = 1;

            // Update TA box
            x1_min = (g1 < x1_min) ? g1 : x1_min;
            x2_min = (g2 < x2_min) ? g2 : x2_min;
            x3_min = (g3 < x3_min) ? g3 : x3_min;
            x4_min = (g4 < x4_min) ? g4 : x4_min;
            x1_max = (g1 > x1_max) ? g1 : x1_max;
            x2_max = (g2 > x2_max) ? g2 : x2_max;
            x3_max = (g3 > x3_max) ? g3 : x3_max;
            x4_max = (g4 > x4_max) ? g4 : x4_max;
        }
        tmpVec.clear();

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-6 TAEX) = %lf sec\n", t_1_elapsed);

        // Reset
        t_1_begin = omp_get_wtime();

        #pragma omp parallel for
        for (int i1 = std::max(x1_min-ExLimit-DIMENSIONS, 1); i1 <= std::min(x1_max+ExLimit+DIMENSIONS, BoxShapeFine[0]-1); i1 ++)  {
            for (int i2 = std::max(x2_min-ExLimit-DIMENSIONS, 1); i2 <= std::min(x2_max+ExLimit+DIMENSIONS, BoxShapeFine[1]-1); i2 ++)  {
                for (int i3 = std::max(x3_min-ExLimit-DIMENSIONS, 1); i3 <= std::min(x3_max+ExLimit+DIMENSIONS, BoxShapeFine[2]-1); i3 ++)  {
                    for (int i4 = std::max(x4_min-ExLimit-DIMENSIONS, 1); i4 <= std::min(x4_max+ExLimit+DIMENSIONS, BoxShapeFine[3]-1); i4 ++)  {
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

                log->log("[CaldeiraLeggett4d] Step: %d, Elapsed time: %lf sec\n", tt + 1, t_0_elapsed);
                ta_size = 0;
                #pragma omp parallel for reduction(+: ta_size) schedule(runtime)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        for (int i3 = x3_min; i3 <= x3_max; i3 ++)  {
                            for (int i4 = x4_min; i4 <= x4_max; i4 ++)  {
                                if ( TAMask[i1*W1+i2][i3*W2+i4] )
                                    ta_size += 1;    
                            }
                        }
                    }
                }
                tb_size = TB.size();
                log->log("[CaldeiraLeggett4d] TA size = %lu, TB size = %lu\n", ta_size, tb_size);
                log->log("[CaldeiraLeggett4d] BoxPosition: [%.2f,%.2f][%.2f,%.2f][%.2f,%.2f][%.2f,%.2f]\n",BoxDynamic[0],BoxDynamic[1],BoxDynamic[2],BoxDynamic[3],BoxDynamic[4],BoxDynamic[5],BoxDynamic[6],BoxDynamic[7]);
                log->log("[CaldeiraLeggett4d] BoxShape: [%d,%d,%d,%d]\n",BoxShapeFine[0],BoxShapeFine[1],BoxShapeFine[2],BoxShapeFine[3]);
                log->log("[CaldeiraLeggett4d] TA Range [%d, %d][%d, %d][%d, %d][%d, %d]\n", x1_min, x1_max, x2_min, x2_max, x3_min, x3_max, x4_min, x4_max);
                log->log("[CaldeiraLeggett4d] TA / total grids = %lf\n", ( ta_size * 1.0 ) / GRIDS_TOT);
                log->log("[CaldeiraLeggett4d] ExCount = %d ExLimit = %d\n", Excount, ExLimit);
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
        delete TAMask[i];                
    }
    delete F;
    delete FF;
    delete PF;
    delete KK1;
    delete KK2;
    delete KK3;
    delete KK4;
    delete TAMask;

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
