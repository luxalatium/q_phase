// ==============================================================================
//
//  WignerMoyal2d.cpp
//  QTR
//
//  Created by Albert Lu on 9/27/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 4/4/20
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
#include "WignerMoyal2d.h"

using namespace QTR_NS;
using std::vector;
using std::cout;
using std::endl;
using std::nothrow;

/* ------------------------------------------------------------------------------- */

// DEFINE POTENTIAL TYPE

#if defined WM2DPOT_DW1

#define WAVEFUNCTION(x1,x2) Wavefunction_DW1(x1,x2)
#define POTENTIAL(x1,x2) Potential_DW1(x1,x2)
#define POTENTIAL_X(x1,x2) Vx_DW1(x1,x2)
#define POTENTIAL_XXX(x1,x2) Vxxx_DW1(x1,x2)
#define POTNAME "DoubleWell-1"

#elif defined WM2DPOT_DW2

#define WAVEFUNCTION(x1,x2) Wavefunction_DW2(x1,x2)
#define POTENTIAL(x1,x2) Potential_DW2(x1,x2)
#define POTENTIAL_X(x1,x2) Vx_DW2(x1,x2)
#define POTENTIAL_XXX(x1,x2) Vxxx_DW2(x1,x2)
#define POTNAME "MSW-1"

#elif defined WM2DPOT_DW3

#define WAVEFUNCTION(x1,x2) Wavefunction_DW3(x1,x2)
#define POTENTIAL(x1,x2) Potential_DW3(x1,x2)
#define POTENTIAL_X(x1,x2) Vx_DW3(x1,x2)
#define POTENTIAL_XXX(x1,x2) Vxxx_DW3(x1,x2)
#define POTNAME "MSW-2"

#elif defined WM2DPOT_DW4

#define WAVEFUNCTION(x1,x2) Wavefunction_DW4(x1,x2)
#define POTENTIAL(x1,x2) Potential_DW4(x1,x2)
#define POTENTIAL_X(x1,x2) Vx_DW4(x1,x2)
#define POTENTIAL_XXX(x1,x2) Vxxx_DW4(x1,x2)
#define POTNAME "MSW-3"

#elif defined WM2DPOT_DW5

#define WAVEFUNCTION(x1,x2) Wavefunction_DW5(x1,x2)
#define POTENTIAL(x1,x2) Potential_DW5(x1,x2)
#define POTENTIAL_X(x1,x2) Vx_DW5(x1,x2)
#define POTENTIAL_XXX(x1,x2) Vxxx_DW5(x1,x2)
#define POTNAME "EB"

#elif defined WM2DPOT_DW6

#define WAVEFUNCTION(x1,x2) Wavefunction_DW6(x1,x2)
#define POTENTIAL(x1,x2) Potential_DW6(x1,x2)
#define POTENTIAL_X(x1,x2) Vx_DW6(x1,x2)
#define POTENTIAL_XXX(x1,x2) Vxxx_DW6(x1,x2)
#define POTNAME "GB"

#elif defined WM2DPOT_DW7

#define WAVEFUNCTION(x1,x2) Wavefunction_DW7(x1,x2)
#define POTENTIAL(x1,x2) Potential_DW7(x1,x2)
#define POTENTIAL_X(x1,x2) Vx_DW7(x1,x2)
#define POTENTIAL_XXX(x1,x2) Vxxx_DW7(x1,x2)
#define POTNAME "Flat"

#elif defined WM2DPOT_DW8

#define WAVEFUNCTION(x1,x2) Wavefunction_DW8(x1,x2)
#define POTENTIAL(x1,x2) Potential_DW8(x1,x2)
#define POTENTIAL_X(x1,x2) Vx_DW8(x1,x2)
#define POTENTIAL_XXX(x1,x2) Vxxx_DW8(x1,x2)
#define POTNAME "MSW-4"

#else // Default  WM2DPOT_DW4

#define WAVEFUNCTION(x1,x2) Wavefunction_DW4(x1,x2)
#define POTENTIAL(x1,x2) Potential_DW4(x1,x2)
#define POTENTIAL_X(x1,x2) Vx_DW4(x1,x2)
#define POTENTIAL_XXX(x1,x2) Vxxx_DW4(x1,x2)
#define POTNAME "MSW-3"

#endif

// Define finite-difference order

#if defined FD_O2

#define EDGE 2
#define FDORDER "2"
#define FDO2

#elif defined FD_O4

#define EDGE 3
#define FDORDER "4"
#define FDO4

#elif defined FD_O6

#define EDGE 4
#define FDORDER "6"
#define FDO6

#else

#define EDGE 2
#define FDORDER "2"
#define FDO2

#endif

/* ------------------------------------------------------------------------------- */

WignerMoyal2d::WignerMoyal2d(class QTR *q)
{
    qtr = q;
    err = qtr->error;
    log = qtr->log;
    parameters = qtr->parameters;
    init();
} 
/* ------------------------------------------------------------------------------- */

WignerMoyal2d::~WignerMoyal2d()
{     
    return;
}
/* ------------------------------------------------------------------------------- */

void WignerMoyal2d::init()
{
    log->log("\n\n[WignerMoyal2d] INIT starts ...\n");
    log->log("\n\n[WignerMoyal2d] Potential type: %s\n", POTNAME);
    log->log("\n\n[WignerMoyal2d] Finite-difference order: %s\n", FDORDER);

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
    BoxShape.resize(DIMENSIONS);

    GRIDS_TOT = 1;
    log->log("[WignerMoyal2d] Number of grids = (");

    for (int i = 0; i < DIMENSIONS; i ++)  {

        BoxShape[i] = (int)( (Box[2 * i + 1] - Box[2 * i]) / H[i] ) + 1;
        GRIDS_TOT *= BoxShape[i];

        if ( i < DIMENSIONS - 1 )
            log->log("%d, ", BoxShape[i]);
        else
            log->log("%d)\n", BoxShape[i]);
    }
    M1 = BoxShape[1];
    W1 = BoxShape[1];
    O1 = BoxShape[0] * BoxShape[1];

    // Parameters
    hb = parameters->scxd_hb;
    m  = parameters->scxd_m;
    kb = parameters->scxd_kb;
    V0 = parameters->scxd_V0;
    alpha = parameters->scxd_alpha;
    lambda = parameters->scxd_lambda;
    skin = parameters->scxd_skin;

    // Wavefunction parameters
    Wave0.resize(DIMENSIONS);
    Wave0[0] = parameters->scxd_x01;
    Wave0[1] = parameters->scxd_x02;

    A.resize(DIMENSIONS);
    A[0] = parameters->scxd_a1;
    A[1] = parameters->scxd_a2;

    // Truncate parameters
    isFullGrid = parameters->scxd_isFullGrid;
    TolH = parameters->scxd_TolH;    // Tolerance of probability density for Zero point Cutoff
    TolL = parameters->scxd_TolL;    // Tolerance of probability density for Edge point
    TolHd = parameters->scxd_TolHd;  // Tolerance of probability first diff for Zero point Cutoff
    TolLd = parameters->scxd_TolLd;  // Tolerance of probability density for Edge point
    ExReduce = parameters->scxd_ExReduce; //Extrapolation reduce factor
    ExLimit = parameters->scxd_ExLimit;   //Extrapolation counts limit

    // Transition position
    trans_x0 = parameters->scxd_trans_x0;
    idx_x0 = (int) std::round( ( trans_x0 - Box[0] ) / H[0] );

    log->log("[WignerMoyal2d] INIT done.\n\n");
}
/* ------------------------------------------------------------------------------- */

void WignerMoyal2d::Evolve()
{
    #pragma omp declare reduction (merge : MeshIndex : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

    log->log("[WignerMoyal2d] Evolve starts ...\n");

    // Files
    FILE *pfile;

    // Variables 
    int index;
    int count;
    int n1, n2;
    int ta_size, tb_size;
    int x1_max, x1_min;  // TA box range
    int x2_max, x2_min;  // TA box range
    double sum;
    double coeff1 = 0.0;
    double coeff2 = 0.0;
    double norm;   // normalization factor
    double density;
    double corr;
    double corr_0;
    bool b1, b2, b3, b4, b5;

    // Timing variables
    double t_0_begin, t_0_end;
    double t_1_begin, t_1_end;
    double t_0_elapsed = 0.0;
    double t_1_elapsed = 0.0;

    // Core computation time (RK4, normalization, initialization, etc)
    double t_full = 0.0;
    double t_truncate = 0.0;

    // Overhead time (truncate)
    double t_overhead = 0.0;

    // Constants
    double kh0m = kk / (H[0] * m);
    double k2h1 = kk / H[1];
    double khbsq2h1 = kk * hb * hb / 24.0 / (H[1] * H[1] * H[1]);
    double TolHdX1 = 2.0 * TolHd * H[0];
    double TolHdX2 = 2.0 * TolHd * H[1];
    double TolLdX1 = 2.0 * TolLd * H[0];
    double TolLdX2 = 2.0 * TolLd * H[1];

    // temporary index container
    MeshIndex tmpVec; 

    // Boundary layer container for extrapolation loop
    MeshIndex ExBD;     
 
    //  2d Grid vector and indices
    VectorXi grid;
    int g1, g2;
    double xx1, xx2;
    double f0, kk0;
    double f1p1, f1m1, f1p2, f1m2, f1p3, f1m3, f1p4, f1m4;
    double f2p1, f2m1, f2p2, f2m2, f2p3, f2m3, f2p4, f2m4;
    double kk1p1, kk1m1, kk1p2, kk1m2, kk1p3, kk1m3, kk1p4, kk1m4;
    double kk2p1, kk2m1, kk2p2, kk2m2, kk2p3, kk2m3, kk2p4, kk2m4;

    // Vector iterater
    vector<int>::iterator it;

    // Extrapolation 
    int min_dir;
    bool isFirstExtrp;
    double val, val_min;
    double val_min_abs;
    vector<double> ExTBL;
    int Excount;

    // Neighborlist
    int nneigh = 0;
    vector<vector<int>> neighlist;
    vector<int> neighs(DIMENSIONS);

    // PF_trans
    double pftrans;
    vector<double> PF_trans;
    PF_trans.push_back(0.0);

    log->log("[WignerMoyal2d] Initializing containers ...\n");

    // Initialize containers

    t_0_begin = omp_get_wtime();

    bool *TAMask;

    if ( !isFullGrid ) 
        TAMask = new bool[O1];
    
    double *F = new double[O1];
    double *FF = new double[O1];
    double *PF = new double[O1];
    double *KK1 = new double[O1];
    double *KK2 = new double[O1];
    double *KK3 = new double[O1];
    double *KK4 = new double[O1];

    double *F0;
    double *Ft;

    if ( isCorr )  {
        F0 = new double[BoxShape[0]];
        Ft = new double[BoxShape[0]];

        #pragma omp parallel for
        for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
            F0[i1] = 0.0;
            Ft[i1] = 0.0;
        }
    }

    #pragma omp parallel for
    for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
        for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
            F[i1*W1+i2] = 0.0;
            PF[i1*W1+i2] = 0.0;
            FF[i1*W1+i2] = 0.0;
            KK1[i1*W1+i2] = 0.0;
            KK2[i1*W1+i2] = 0.0;
            KK3[i1*W1+i2] = 0.0;
            KK4[i1*W1+i2] = 0.0;
        }
    }

    if ( !isFullGrid )  {

        t_1_begin = omp_get_wtime();
        
        #pragma omp parallel for
        for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
            for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                TAMask[i1*W1+i2] = 0;
            }
        }

        nneigh = 0;

        for (int d = 1; d <= 4; d ++)  {
            for (n1 = -d; n1 <= d; n1 ++)  {

                n2 = d - abs(n1);

                if (n2 != 0)  {
                    neighs = {n1,n2};
                    neighlist.push_back(neighs);
                    neighs = {n1,-n2};
                    neighlist.push_back(neighs);
                    nneigh += 2;
                }
                else  {
                    neighs = {n1,0};
                    neighlist.push_back(neighs);
                    nneigh += 1;
                }
            }
        }
        log->log("[WignerMoyal2d] nneigh = %d\n", nneigh); 
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_overhead += t_1_elapsed;
    }
    t_0_end = omp_get_wtime();
    t_0_elapsed = t_0_end - t_0_begin;
    t_full += t_0_elapsed;

    if ( !isFullGrid )
        t_truncate += t_0_elapsed - t_1_elapsed; // subtract overhead

    if (!QUIET && TIMING) log->log("[WignerMoyal2d] Elapsed time (initializing containers) = %lf sec\n\n", t_0_elapsed); 

    // .........................................................................................

    // Initialize wavefunction
    t_1_begin = omp_get_wtime();
    log->log("[WignerMoyal2d] Initializing wavefunction ...\n");  

    #pragma omp parallel for
    for (int i1 = EDGE; i1 < BoxShape[0] - EDGE ; i1 ++)  {
        for (int i2 = EDGE; i2 < BoxShape[1] - EDGE ; i2 ++)  {
            F[i1*W1+i2] = WAVEFUNCTION(Box[0]+i1*H[0], Box[2]+i2*H[1]);
        }
    }

    // Normalization
    norm = 0.0;

    #pragma omp parallel for reduction (+:norm)
    for (int i1 = EDGE; i1 <  BoxShape[0] - EDGE; i1 ++)  {
        for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  { 
            norm += F[i1*W1+i2];
        }
    }
    norm *= H[0] * H[1];
    log->log("[WignerMoyal2d] Normalization factor = %.16e\n",norm);
    norm = 1.0 / norm;

    #pragma omp parallel for
    for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
        for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
            F[i1*W1+i2] = norm * F[i1*W1+i2];
            PF[i1*W1+i2] = F[i1*W1+i2];
        }
    }

    // Initial density
    if ( isCorr )   {
        #pragma omp parallel for private(density)
        for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
            density = 0.0;
            for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                density += PF[i1*W1+i2]; 
            }
            F0[i1] = density * H[1];
        }

        corr_0 = 0.0;
        for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
            corr_0 += F0[i1] * F0[i1];
        }

        corr_0 = corr_0 * H[0];

        log->log("[WignerMoyal2d] corr_0 = %.16e\n",corr_0);
        log->log("[WignerMoyal2d] Time %lf, Corr = %.16e\n",0.0,1.0);
    }

    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    t_full += t_1_elapsed;
    t_truncate += t_1_elapsed;
    if (!QUIET && TIMING) log->log("[WignerMoyal2d] Elapsed time (initializing wavefunction) = %lf sec\n\n", t_1_elapsed); 

    // .........................................................................................

    // Initial truncation & edge point check

    if ( !isFullGrid )
    {
        t_1_begin = omp_get_wtime();

        log->log("[WignerMoyal2d] Initial truncation ...\n");

        ta_size = 0;

        // Truncation and TA
        #pragma omp parallel for reduction(+: ta_size) private(b1, b2, b3)
        for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
            for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {

                b1 = std::abs(PF[i1*W1+i2]) < TolH;
                b2 = std::abs(PF[(i1+1)*W1+i2] - PF[(i1-1)*W1+i2]) < TolHdX1;
                b3 = std::abs(PF[i1*W1+(i2+1)] - PF[i1*W1+(i2-1)]) < TolHdX2;
        
                if (b1 && b2 && b3)
                    F[i1*W1+i2] = 0.0;
                else {
                    TAMask[i1*W1+i2] = 1;
                    ta_size += 1;
                }
            }
        }   
        log->log("[WignerMoyal2d] TA size = %d\n", ta_size);

        // TA box
        x1_min = BoxShape[0];
        x2_min = BoxShape[1];
        x1_max = 0;
        x2_max = 0;

        #pragma omp parallel for reduction(min: x1_min, x2_min) reduction(max: x1_max, x2_max)
        for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
            for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                if (TAMask[i1*W1+i2])  {
                    if (i1 < x1_min)  x1_min = i1;
                    if (i1 > x1_max)  x1_max = i1;
                    if (i2 < x2_min)  x2_min = i2;
                    if (i2 > x2_max)  x2_max = i2;
                }
            }
        }
        // `````````````````````````````````````````````````````````````````
        // TB
        tmpVec.clear();

        #pragma omp parallel for reduction(merge: tmpVec) private(b1, b2, b3, b4) schedule(runtime)
        for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
            for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                if (TAMask[i1*W1+i2] == 1 && (!TAMask[(i1-1)*W1+i2] || !TAMask[(i1+1)*W1+i2] || !TAMask[i1*W1+(i2-1)] || !TAMask[i1*W1+(i2+1)]))
                    tmpVec.push_back(i1*M1+i2);
            }
        }
        tmpVec.swap(TB);
        tmpVec.clear();
        tb_size = TB.size();
        log->log("[WignerMoyal2d] TB size = %d\n", tb_size);

        // `````````````````````````````````````````````````````````````````

        // TA expansion

        #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2) schedule(runtime)
        for (int i = 0; i < TB.size(); i++)
        {
            g1 = (int)(TB[i] / M1);
            g2 = (int)(TB[i] % M1);

            if ( g1+1 != BoxShape[0]-EDGE-1 && !TAMask[(g1+1)*W1+g2])
                tmpVec.push_back(GridToIdx(g1+1,g2));
            if ( g1-1 != EDGE && !TAMask[(g1-1)*W1+g2])
                tmpVec.push_back(GridToIdx(g1-1,g2));
            if ( g2+1 != BoxShape[1]-EDGE-1 && !TAMask[g1*W1+(g2+1)])
                tmpVec.push_back(GridToIdx(g1,g2+1));
            if ( g2-1 != EDGE && !TAMask[g1*W1+(g2-1)])
                tmpVec.push_back(GridToIdx(g1,g2-1));
        }

        for (int i = 0; i < tmpVec.size(); i ++)  {

            g2 = (int)(tmpVec[i] % M1);
            g1 = (int)(tmpVec[i] / M1);
            TAMask[g1*W1+g2] = 1;

            // Update TA box
            x1_min = (g1 < x1_min) ? g1 : x1_min;
            x2_min = (g2 < x2_min) ? g2 : x2_min;
            x1_max = (g1 > x1_max) ? g1 : x1_max;
            x2_max = (g2 > x2_max) ? g2 : x2_max;
        }
        tmpVec.clear();

        ta_size = 0;
        #pragma omp parallel for reduction(+: ta_size) schedule(runtime)
        for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
            for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                if (TAMask[i1*W1+i2])
                    ta_size += 1;    
            }
        }
        log->log("[WignerMoyal2d] TA size = %d, TB size = %d\n", ta_size, tb_size);
        log->log("[WignerMoyal2d] TA Range [%d, %d][%d, %d]\n", x1_min, x1_max, x2_min, x2_max);

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_overhead += t_1_elapsed;

        if (!QUIET && TIMING)  {
            log->log("[WignerMoyal2d] Elapsed time (initial truncation) = %lf sec\n\n", t_1_elapsed);
            log->log("[WignerMoyal2d] Initialization core computation time: %lf sec\n", t_truncate);            
            log->log("[WignerMoyal2d] Initialization overhead: %lf sec\n", t_overhead); 
        }
    }
    else  // Full grid approach
    {
        log->log("[WignerMoyal2d] Initialization core computation time: %lf sec\n", t_full);   
    }
    // .........................................................................................

    // Time iteration 

    log->log("=======================================================\n\n"); 
    log->log("[WignerMoyal2d] Time interation starts ...\n"); 
    log->log("[WignerMoyal2d] Number of steps = %d\n\n", (int)(TIME / kk)); 
    log->log("=======================================================\n\n"); 

    for (int tt = 0; tt < (int)(TIME / kk); tt ++)
    {
        t_0_begin = omp_get_wtime(); 
        Excount = 0;

        if ( isPrintWavefunc && tt % PRINT_WAVEFUNC_PERIOD == 0 )
        {
            pfile = fopen("wave.dat","a");

            if ( !isFullGrid )  {
                fprintf(pfile, "%d %d\n", tt, ta_size);
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])
                            fprintf(pfile, "%d %d %.8e\n", i1, i2, F[i1*W1+i2]);  
                    }
                }
            }
            else  {
                fprintf(pfile, "%d %d\n", tt, GRIDS_TOT);
                for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                        fprintf(pfile, "%d %d %.8e\n", i1, i2, F[i1*W1+i2]);
                    }
                }
            }
            fclose(pfile);
        }

        if ( tt % PRINT_PERIOD == 0 )
        {
            if ( isPrintEdge  && !isFullGrid )  {

                pfile = fopen ("edge.dat","a");
                fprintf(pfile, "%d %lf %lu\n", tt, tt * kk, TB.size());

                for (int i = 0; i < TB.size(); i++)
                {
                    g2 = (int)(TB[i] % M1);
                    g1 = (int)(TB[i] / M1);
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];
                    fprintf(pfile, "%d %d %lf %lf\n", g1, g2, xx1, xx2);          
                }
                fclose(pfile);
            }
            if ( isPrintDensity && !isFullGrid )  {

                pfile = fopen ("density.dat","a");
                fprintf(pfile, "%d %lf %d\n", tt, tt * kk, ta_size);

                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            xx1 = Box[0] + i1 * H[0];
                            xx2 = Box[2] + i2 * H[1];
                            fprintf(pfile, "%d %d %lf %lf %.16e\n", i1, i2, xx1, xx2, PF[i1*W1+i2]);
                        }
                    }
                }
                fclose(pfile);
            }
            if ( isPrintDensity && isFullGrid )  {

                pfile = fopen ("density.dat","a");
                fprintf(pfile, "%d %lf %d\n", tt, tt * kk, BoxShape[0] * BoxShape[1] );
                //fprintf(pfile, "%d %lf %d\n", tt, tt * kk, BoxShape[0] );

                for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {

                    //density = 0.0;
                    for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                        g1 = i1;
                        g2 = i2;
                        xx1 = Box[0] + g1 * H[0];
                        xx2 = Box[2] + g2 * H[1];
                        fprintf(pfile, "%d %d %lf %lf %.16e\n", g1, g2, xx1, xx2, PF[g1*W1+g2]);
                        //density += PF[i1*W1+i2];
                    }
                    //xx1 = Box[0] + i1 * H[0];
                    //fprintf(pfile, "%d %lf %.16e\n", i1, xx1, density);
                }
                fclose(pfile);  
            }
        }

        // Check if TB of f is higher than TolL
        
        if ( !isFullGrid )
        {
            t_1_begin = omp_get_wtime();
            t_truncate = 0.0;
            t_overhead = 0.0;

            // TBL = Index of Extrapolating Edge points.
            // TBL_P = Index history of TBL of this iteration. To Prevent from extrapolating the same points multiple times.

            TBL.clear();
            tmpVec.clear();

            #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, b1, b2, b3, b4, b5)
            for (int i = 0; i < TB.size(); i++)
            {
                g1 = (int)(TB[i] / M1);
                g2 = (int)(TB[i] % M1);
                b1 = std::abs(PF[g1*W1+g2]) >= TolL;
                b2 = std::abs(PF[(g1+1)*W1+g2]*TAMask[(g1+1)*W1+g2] - PF[(g1-1)*W1+g2]*TAMask[(g1-1)*W1+g2]) >= TolLdX1;
                b3 = std::abs(PF[g1*W1+(g2+1)]*TAMask[g1*W1+(g2+1)] - PF[g1*W1+(g2-1)]*TAMask[g1*W1+(g2-1)]) >= TolLdX2;

                // Not in DBi2
                b4 = g1 > EDGE && g2 > EDGE;
                b5 = g1 < BoxShape[0]-EDGE-1 && g2 < BoxShape[1]-EDGE-1;

                if ( (b1 || b2 || b3 ) && b4 && b5 )
                    tmpVec.push_back(TB[i]);
            }
            tmpVec.swap(TBL);
            tmpVec.clear();
            TBL_P = TBL;

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-a-1: TBL) = %lf sec\n", t_1_elapsed);   
            //if (!QUIET) log->log("TBL size = %d TBL_P size = %d\n", TBL.size(), TBL_P.size());
        }
        else  
        {
            t_full = 0.0;
        }
        isExtrapolate = false;
        isFirstExtrp = true;
        // .........................................................................................

        // CASE 1: Truncating with extrapolation

        while ( TBL.size() != 0 && !isFullGrid && Excount < ExLimit)
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
                g2 = (int)(TBL[index] % M1);

                if ( g1-1 != EDGE && F[(g1-1)*W1+g2] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1-1,g2));
                }
                if ( g1+1 != BoxShape[0]-EDGE-1 && F[(g1+1)*W1+g2] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1+1,g2));
                }
                if ( g2-1 != EDGE && F[g1*W1+(g2-1)] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1,g2-1));
                }
                if ( g2+1 != BoxShape[1]-EDGE-1 && F[g1*W1+(g2+1)] == 0.0 )  {
                    ExFF.push_back(GridToIdx(g1,g2+1));
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
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-b-1: ExFF) = %lf sec\n", t_1_elapsed);   

            // .....................................................................

            // Find the direction of Outer to Edge points

            t_1_begin = omp_get_wtime();
            ExTBL.clear();
            it = ExFF.begin();
            
            while ( it != ExFF.end() )  {

                index = std::distance( ExFF.begin(), it );
                g1 = (int)(ExFF[index] / M1);
                g2 = (int)(ExFF[index] % M1);
                sum = 0.0;
                count = 0;
                isEmpty = true;
                val_min_abs = 100000000;
                val_min = 100000000;
                min_dir = -1;

                if ( F[(g1-1)*W1+g2] != 0.0 )  {

                    if ( std::abs(F[(g1-1)*W1+g2]) < val_min_abs && F[(g1-2)*W1+g2] != 0.0 )  {
                        val_min_abs = std::abs(F[(g1-1)*W1+g2]);
                        val_min = F[(g1-1)*W1+g2];
                        min_dir = 0;
                    }
                    if ( F[(g1-2)*W1+g2] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[(g1-1)*W1+g2]) - std::log(F[(g1-2)*W1+g2]) );

                        if ( !(isnan(val) || isnan(-val)) && !(isinf(val) || isinf(-val)) )  
                        {
                            sum += val;
                            count += 1;
                            isEmpty = false;
                        }
                    }
                }

                if ( F[(g1+1)*W1+g2] != 0.0 )  {

                    if ( std::abs(F[(g1+1)*W1+g2]) < val_min_abs && F[(g1+2)*W1+g2] != 0.0  )  {
                        val_min_abs = std::abs(F[(g1+1)*W1+g2]);
                        val_min = F[(g1+1)*W1+g2];
                        min_dir = 0;
                    }
                    if ( F[(g1+2)*W1+g2] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[(g1+1)*W1+g2]) - std::log(F[(g1+2)*W1+g2]) );

                        if ( !(isnan(val) || isnan(-val)) && !(isinf(val) || isinf(-val)) ) 
                        {
                            sum += val;
                            count += 1;
                            isEmpty = false;
                        }
                    }
                }

                if ( F[g1*W1+(g2-1)] != 0.0 )  {

                    if ( std::abs(F[g1*W1+(g2-1)]) < val_min_abs && F[g1*W1+(g2-2)] != 0.0  )  {

                        val_min_abs = std::abs(F[g1*W1+(g2-1)]);
                        val_min = F[g1*W1+(g2-1)];
                        min_dir = 1;
                    }
                    if ( F[g1*W1+(g2-2)] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[g1*W1+(g2-1)]) - std::log(F[g1*W1+(g2-2)]) );

                        if ( !(isnan(val) || isnan(-val)) && !(isinf(val) || isinf(-val)) ) 
                        {
                            sum += val;
                            count += 1;
                            isEmpty = false;
                        }
                    }
                }

                if ( F[g1*W1+(g2+1)] != 0.0 )  {

                    if ( std::abs(F[g1*W1+(g2+1)]) < val_min_abs && F[g1*W1+(g2+2)] != 0.0 )  {

                        val_min_abs = std::abs(F[g1*W1+(g2+1)]);
                        val_min = F[g1*W1+(g2+1)];
                        min_dir = 1;
                    }
                    if ( F[g1*W1+(g2+2)] != 0.0 )  {

                        val = exp( 2.0 * std::log(F[g1*W1+(g2+1)]) - std::log(F[g1*W1+(g2+2)]) );

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

            for ( int i = 0; i < ExFF.size(); i++ )  {
                g1 = (int)(ExFF[i] / M1);
                g2 = (int)(ExFF[i] % M1);
                F[g1*W1+g2] = ExTBL[i];
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-b-2: ExFF) = %lf sec\n", t_1_elapsed);  

            // ............................................................................................. Extrapolation

            if ( isFirstExtrp )  {

                // Check Extending nonzero Area

                t_1_begin = omp_get_wtime();

                #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2)
                for (int i = 0; i < ExFF.size(); i++)  {
                    g1 = (int)(ExFF[i] / M1);
                    g2 = (int)(ExFF[i] % M1);

                    if ( !TAMask[g1*W1+g2] )
                        tmpVec.push_back(GridToIdx(g1,g2));
                    if ( g1+1 != BoxShape[0]-EDGE-1 && !TAMask[(g1+1)*W1+g2] )
                        tmpVec.push_back(GridToIdx(g1+1,g2));
                    if ( g1-1 != EDGE && !TAMask[(g1-1)*W1+g2] )
                        tmpVec.push_back(GridToIdx(g1-1,g2));
                    if ( g2+1 != BoxShape[1]-EDGE-1 && !TAMask[g1*W1+(g2+1)] )
                        tmpVec.push_back(GridToIdx(g1,g2+1));
                    if ( g2-1 != EDGE && !TAMask[g1*W1+(g2-1)] )
                        tmpVec.push_back(GridToIdx(g1,g2-1));
                }

                #pragma omp parallel for reduction(min: x1_min,x2_min) reduction(max: x1_max,x2_max) private(g1, g2)
                for (int i = 0; i < tmpVec.size(); i ++)  {

                    g1 = (int)(tmpVec[i] / M1);
                    g2 = (int)(tmpVec[i] % M1);
                    TAMask[g1*W1+g2] = 1;
                    
                    // Update TA box
                    x1_min = (g1 < x1_min) ? g1 : x1_min;
                    x2_min = (g2 < x2_min) ? g2 : x2_min;
                    x1_max = (g1 > x1_max) ? g1 : x1_max;
                    x2_max = (g2 > x2_max) ? g2 : x2_max;
                }
                tmpVec.clear();

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-c-1: CASE 1 TA) = %lf sec\n", t_1_elapsed); 
        
                // Runge–Kutta 4
                #pragma omp parallel 
                {
                    // RK4-1
                    #pragma omp single nowait
                    {
                        t_1_begin = omp_get_wtime();
                    }
                    #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4) schedule(runtime)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            if (TAMask[i1*W1+i2])  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                f0 = F[i1*W1+i2];
                                f1p1 = F[(i1+1)*W1+i2];
                                f1m1 = F[(i1-1)*W1+i2];
                                f2p1 = F[i1*W1+(i2+1)];
                                f2m1 = F[i1*W1+(i2-1)];
                                f1p2 = F[(i1+2)*W1+i2];
                                f1m2 = F[(i1-2)*W1+i2];
                                f2p2 = F[i1*W1+(i2+2)];
                                f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)                             
                                f1p3 = F[(i1+3)*W1+i2];
                                f1m3 = F[(i1-3)*W1+i2];
                                f2p3 = F[i1*W1+(i2+3)];
                                f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                                f1p4 = F[(i1+4)*W1+i2];
                                f1m4 = F[(i1-4)*W1+i2];
                                f2p4 = F[i1*W1+(i2+4)];
                                f2m4 = F[i1*W1+(i2-4)];
#endif
                                KK1[i1*W1+i2] = -kh0m * xx2 * (-f1p2/12.0 + 2/3.0*f1p1 - 2/3.0*f1m1 + f1m2/12.0) + 
                                                k2h1 * POTENTIAL_X(xx1, xx2) * (-f2p2/12.0 + 2/3.0*f2p1 - 2/3.0*f2m1 + f2m2/12.0) -khbsq2h1 * POTENTIAL_XXX(xx1, xx2) * 
#if defined(FDO4)
                                                (-f2p3/8.0+f2p2-13.0*f2p1/8.0+13*f2m1/8.0-f2m2+f2m3/8.0);
#elif defined(FDO6)                                               
                                                (7.0/240.0*f2p4-3.0/10.0*f2p3+169.0/120.0*f2p2-61.0/30.0*f2p1+61.0/30.0*f2m1-169.0/120.0*f2m2+3.0/10.0*f2m3-7.0/240.0*f2m4);
#else
                                                (0.5*f2p2-f2p1+f2m1-0.5*f2m2);
#endif
                                FF[i1*W1+i2] = F[i1*W1+i2] + KK1[i1*W1+i2] / 6.0;
                            }
                        }
                    }
                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-11: CASE 1 KK1) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-2
                    #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            if (TAMask[i1*W1+i2])  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                f0 = F[i1*W1+i2];
                                f1p1 = F[(i1+1)*W1+i2];
                                f1m1 = F[(i1-1)*W1+i2];
                                f2p1 = F[i1*W1+(i2+1)];
                                f2m1 = F[i1*W1+(i2-1)];
                                f1p2 = F[(i1+2)*W1+i2];
                                f1m2 = F[(i1-2)*W1+i2];
                                f2p2 = F[i1*W1+(i2+2)];
                                f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)                             
                                f1p3 = F[(i1+3)*W1+i2];
                                f1m3 = F[(i1-3)*W1+i2];
                                f2p3 = F[i1*W1+(i2+3)];
                                f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                                f1p4 = F[(i1+4)*W1+i2];
                                f1m4 = F[(i1-4)*W1+i2];
                                f2p4 = F[i1*W1+(i2+4)];
                                f2m4 = F[i1*W1+(i2-4)];
#endif
                                kk0 = KK1[i1*W1+i2];
                                kk1p1 = KK1[(i1+1)*W1+i2];
                                kk1m1 = KK1[(i1-1)*W1+i2];
                                kk2p1 = KK1[i1*W1+(i2+1)];
                                kk2m1 = KK1[i1*W1+(i2-1)];
                                kk1p2 = KK1[(i1+2)*W1+i2];
                                kk1m2 = KK1[(i1-2)*W1+i2];
                                kk2p2 = KK1[i1*W1+(i2+2)];
                                kk2m2 = KK1[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)  
                                kk1p3 = KK1[(i1+3)*W1+i2];
                                kk1m3 = KK1[(i1-3)*W1+i2];
                                kk2p3 = KK1[i1*W1+(i2+3)];
                                kk2m3 = KK1[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                                kk1p4 = KK1[(i1+4)*W1+i2];
                                kk1m4 = KK1[(i1-4)*W1+i2];
                                kk2p4 = KK1[i1*W1+(i2+4)];
                                kk2m4 = KK1[i1*W1+(i2-4)];
#endif
                                KK2[i1*W1+i2] = -kh0m * xx2 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                                k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) -khbsq2h1 * POTENTIAL_XXX(xx1, xx2) *
#if defined(FDO4)
                                                (-(f2p3+0.5*kk2p3)/8.0+(f2p2+0.5*kk2p2)-13.0*(f2p1+0.5*kk2p1)/8.0+13.0*(f2m1+0.5*kk2m1)/8.0-(f2m2+0.5*kk2m2)+(f2m3+0.5*kk2m3)/8.0);
#elif defined(FDO6)
                                                (7.0/240.0*(f2p4+0.5*kk2p4)-3.0/10.0*(f2p3+0.5*kk2p3)+169.0/120.0*(f2p2+0.5*kk2p2)-61.0/30.0*(f2p1+0.5*kk2p1)+61.0/30.0*(f2m1+0.5*kk2m1)-169.0/120.0*(f2m2+0.5*kk2m2)+3.0/10.0*(f2m3+0.5*kk2m3)-7.0/240.0*(f2m4+0.5*kk2m4));
#else
                                                (0.5*(f2p2+0.5*kk2p2)-(f2p1+0.5*kk2p1)+(f2m1+0.5*kk2m1)-0.5*(f2m2+0.5*kk2m2));
#endif
                                FF[i1*W1+i2] += KK2[i1*W1+i2] / 3.0;
                            }
                        }
                    }
                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-12: CASE 1 KK2) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-3
                    #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            if (TAMask[i1*W1+i2])  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                f0 = F[i1*W1+i2];
                                f1p1 = F[(i1+1)*W1+i2];
                                f1m1 = F[(i1-1)*W1+i2];
                                f2p1 = F[i1*W1+(i2+1)];
                                f2m1 = F[i1*W1+(i2-1)];
                                f1p2 = F[(i1+2)*W1+i2];
                                f1m2 = F[(i1-2)*W1+i2];
                                f2p2 = F[i1*W1+(i2+2)];
                                f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                                f1p3 = F[(i1+3)*W1+i2];
                                f1m3 = F[(i1-3)*W1+i2];
                                f2p3 = F[i1*W1+(i2+3)];
                                f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6)                                 
                                f1p4 = F[(i1+4)*W1+i2];
                                f1m4 = F[(i1-4)*W1+i2];
                                f2p4 = F[i1*W1+(i2+4)];
                                f2m4 = F[i1*W1+(i2-4)];
#endif
                                kk0 = KK2[i1*W1+i2];
                                kk1p1 = KK2[(i1+1)*W1+i2];
                                kk1m1 = KK2[(i1-1)*W1+i2];
                                kk2p1 = KK2[i1*W1+(i2+1)];
                                kk2m1 = KK2[i1*W1+(i2-1)];
                                kk1p2 = KK2[(i1+2)*W1+i2];
                                kk1m2 = KK2[(i1-2)*W1+i2];
                                kk2p2 = KK2[i1*W1+(i2+2)];
                                kk2m2 = KK2[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                                kk1p3 = KK2[(i1+3)*W1+i2];
                                kk1m3 = KK2[(i1-3)*W1+i2];
                                kk2p3 = KK2[i1*W1+(i2+3)];
                                kk2m3 = KK2[i1*W1+(i2-3)];
#endif
#if defined(FDO6)  
                                kk1p4 = KK2[(i1+4)*W1+i2];
                                kk1m4 = KK2[(i1-4)*W1+i2];
                                kk2p4 = KK2[i1*W1+(i2+4)];
                                kk2m4 = KK2[i1*W1+(i2-4)];
#endif
                                KK3[i1*W1+i2] = -kh0m * xx2 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                                k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) -khbsq2h1 * POTENTIAL_XXX(xx1, xx2) *
#if defined(FDO4)
                                                (-(f2p3+0.5*kk2p3)/8.0+(f2p2+0.5*kk2p2)-13.0*(f2p1+0.5*kk2p1)/8.0+13.0*(f2m1+0.5*kk2m1)/8.0-(f2m2+0.5*kk2m2)+(f2m3+0.5*kk2m3)/8.0);
#elif defined(FDO6)
                                                (7.0/240.0*(f2p4+0.5*kk2p4)-3.0/10.0*(f2p3+0.5*kk2p3)+169.0/120.0*(f2p2+0.5*kk2p2)-61.0/30.0*(f2p1+0.5*kk2p1)+61.0/30.0*(f2m1+0.5*kk2m1)-169.0/120.0*(f2m2+0.5*kk2m2)+3.0/10.0*(f2m3+0.5*kk2m3)-7.0/240.0*(f2m4+0.5*kk2m4));
#else
                                                (0.5*(f2p2+0.5*kk2p2) - (f2p1+0.5*kk2p1) + (f2m1+0.5*kk2m1) - 0.5*(f2m2+0.5*kk2m2));
#endif
                                FF[i1*W1+i2] += KK3[i1*W1+i2] / 3.0;
                            }
                        }            
                    }
                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-13: CASE 1 KK3) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-4
                    #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)
                    for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                            if (TAMask[i1*W1+i2])  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                f0 = F[i1*W1+i2];
                                f1p1 = F[(i1+1)*W1+i2];
                                f1m1 = F[(i1-1)*W1+i2];
                                f2p1 = F[i1*W1+(i2+1)];
                                f2m1 = F[i1*W1+(i2-1)];
                                f1p2 = F[(i1+2)*W1+i2];
                                f1m2 = F[(i1-2)*W1+i2];
                                f2p2 = F[i1*W1+(i2+2)];
                                f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                                f1p3 = F[(i1+3)*W1+i2];
                                f1m3 = F[(i1-3)*W1+i2];
                                f2p3 = F[i1*W1+(i2+3)];
                                f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6)
                                f1p4 = F[(i1+4)*W1+i2];
                                f1m4 = F[(i1-4)*W1+i2];
                                f2p4 = F[i1*W1+(i2+4)];
                                f2m4 = F[i1*W1+(i2-4)];
#endif
                                kk0 = KK3[i1*W1+i2];
                                kk1p1 = KK3[(i1+1)*W1+i2];
                                kk1m1 = KK3[(i1-1)*W1+i2];
                                kk2p1 = KK3[i1*W1+(i2+1)];
                                kk2m1 = KK3[i1*W1+(i2-1)];
                                kk1p2 = KK3[(i1+2)*W1+i2];
                                kk1m2 = KK3[(i1-2)*W1+i2];
                                kk2p2 = KK3[i1*W1+(i2+2)];
                                kk2m2 = KK3[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                                kk1p3 = KK3[(i1+3)*W1+i2];
                                kk1m3 = KK3[(i1-3)*W1+i2];
                                kk2p3 = KK3[i1*W1+(i2+3)];
                                kk2m3 = KK3[i1*W1+(i2-3)];
#endif
#if defined(FDO6)
                                kk1p4 = KK3[(i1+4)*W1+i2];
                                kk1m4 = KK3[(i1-4)*W1+i2];
                                kk2p4 = KK3[i1*W1+(i2+4)];
                                kk2m4 = KK3[i1*W1+(i2-4)];
#endif
                                KK4[i1*W1+i2] = -kh0m * xx2 * (-1/12.0*(f1p2+kk1p2) + 2/3.0*(f1p1+kk1p1) - 2/3.0*(f1m1+kk1m1) + 1/12.0*(f1m2+kk1m2)) + 
                                                k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+kk2p2) + 2/3.0*(f2p1+kk2p1) - 2/3.0*(f2m1+kk2m1) + 1/12.0*(f2m2+kk2m2)) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) * 
#if defined(FDO4)
                                                (-(f2p3+kk2p3)/8.0+(f2p2+kk2p2)-13.0*(f2p1+kk2p1)/8.0+13.0*(f2m1+kk2m1)/8.0-(f2m2+kk2m2)+(f2m3+kk2m3)/8.0);
#elif defined(FDO6)
                                                (7.0/240.0*(f2p4+kk2p4)-3.0/10.0*(f2p3+kk2p3)+169.0/120.0*(f2p2+kk2p2)-61.0/30.0*(f2p1+kk2p1)+61.0/30.0*(f2m1+kk2m1)-169.0/120.0*(f2m2+kk2m2)+3.0/10.0*(f2m3+kk2m3)-7.0/240.0*(f2m4+kk2m4));
#else
                                                (0.5*(f2p2+kk2p2)-(f2p1+kk2p1)+(f2m1+kk2m1)-0.5*(f2m2+kk2m2));
#endif
                                FF[i1*W1+i2] += KK4[i1*W1+i2] / 6.0;
                           }
                        }
                    }
                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-14: CASE 1 KK4) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }
                } // OMP PARALLEL
            }
            else  
            {
                // Extrapolation loop when multiple expanding occured

                t_1_begin = omp_get_wtime();

                #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2)
                for (int i = 0; i < ExFF.size(); i++)  {

                    g1 = (int)(ExFF[i] / M1);
                    g2 = (int)(ExFF[i] % M1);

                    if ( !TAMask[g1*W1+g2] )
                        tmpVec.push_back(GridToIdx(g1,g2));
                    if ( g1+1 != BoxShape[0]-EDGE-1 && !TAMask[(g1+1)*W1+g2] )
                        tmpVec.push_back(GridToIdx(g1+1,g2));
                    if ( g1-1 != EDGE && !TAMask[(g1-1)*W1+g2] )
                        tmpVec.push_back(GridToIdx(g1-1,g2));
                    if ( g2+1 != BoxShape[1]-EDGE-1 && !TAMask[g1*W1+(g2+1)] )
                        tmpVec.push_back(GridToIdx(g1,g2+1));
                    if ( g2-1 != EDGE && !TAMask[g1*W1+(g2-1)] )
                        tmpVec.push_back(GridToIdx(g1,g2-1));
                }

                for (int i = 0; i < tmpVec.size(); i ++)  {

                    g1 = (int)(tmpVec[i] / M1); 
                    g2 = (int)(tmpVec[i] % M1);
                    TAMask[g1*W1+g2] = 1;

                    // Update TA box
                    x1_min = (g1 < x1_min) ? g1 : x1_min;
                    x2_min = (g2 < x2_min) ? g2 : x2_min;
                    x1_max = (g1 > x1_max) ? g1 : x1_max;
                    x2_max = (g2 > x2_max) ? g2 : x2_max;
                }
                tmpVec.clear();
                ExBD.clear();

                #pragma omp parallel for reduction(merge: ExBD) private(g1, g2, n1, n2)
                for (int i = 0; i < ExFF.size(); i++)
                {
                    g1 = (int)(ExFF[i] / M1);
                    g2 = (int)(ExFF[i] % M1);
                    ExBD.push_back(ExFF[i]);

                    for (int j = 0; j < nneigh; j ++)  {

                        n1 = neighlist[j][0];
                        n2 = neighlist[j][1];

                        if (TAMask[(g1+n1)*W1+(g2+n2)])
                            ExBD.push_back(GridToIdx(g1+n1,g2+n2));
                    }
                }

                // Find unique elements (ExBD)
                __gnu_parallel::sort(ExBD.begin(),ExBD.end());
                it = std::unique (ExBD.begin(), ExBD.end()); 
                ExBD.resize(std::distance(ExBD.begin(),it));

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-cx-1: CASE 1 ExBD) = %lf sec\n", t_1_elapsed); 

                // Runge–Kutta 4

                t_1_begin = omp_get_wtime();

                // RK4-1
                #pragma omp for private(g1,g2,xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4) schedule(runtime)
                for (int i = 0; i < ExBD.size(); i++)  {
                    g1 = (int)(ExBD[i] / M1);
                    g2 = (int)(ExBD[i] % M1);
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];
                    f0 = F[g1*W1+g2];
                    f1p1 = F[(g1+1)*W1+g2];
                    f1m1 = F[(g1-1)*W1+g2];
                    f2p1 = F[g1*W1+(g2+1)];
                    f2m1 = F[g1*W1+(g2-1)];
                    f1p2 = F[(g1+2)*W1+g2];
                    f1m2 = F[(g1-2)*W1+g2];
                    f2p2 = F[g1*W1+(g2+2)];
                    f2m2 = F[g1*W1+(g2-2)];
#if defined(FDO4) || defined(FDO6)     
                    f1p3 = F[(g1+3)*W1+g2];
                    f1m3 = F[(g1-3)*W1+g2];
                    f2p3 = F[g1*W1+(g2+3)];
                    f2m3 = F[g1*W1+(g2-3)];
#endif
#if defined(FDO6) 
                    f1p4 = F[(g1+4)*W1+g2];
                    f1m4 = F[(g1-4)*W1+g2];
                    f2p4 = F[g1*W1+(g2+4)];
                    f2m4 = F[g1*W1+(g2-4)];
#endif
                    KK1[g1*W1+g2] = -kh0m * xx2 * (-f1p2/12.0 + 2/3.0*f1p1 - 2/3.0*f1m1 + f1m2/12.0) + 
                                    k2h1 * POTENTIAL_X(xx1, xx2) * (-f2p2/12.0 + 2/3.0*f2p1 - 2/3.0*f2m1 + f2m2/12.0) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) * 
#if defined(FDO4)
                                    (-f2p3/8.0+f2p2-13.0*f2p1/8.0+13*f2m1/8.0-f2m2+f2m3/8.0);
#elif defined(FDO6)                                               
                                    (7.0/240.0*f2p4-3.0/10.0*f2p3+169.0/120.0*f2p2-61.0/30.0*f2p1+61.0/30.0*f2m1-169.0/120.0*f2m2+3.0/10.0*f2m3-7.0/240.0*f2m4);
#else
                                    (0.5*f2p2-f2p1+f2m1-0.5*f2m2);
#endif

                    FF[g1*W1+g2] = F[g1*W1+g2] + KK1[g1*W1+g2] / 6.0;
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-11: CASE 1 KK1) = %lf sec\n", t_1_elapsed);
                t_1_begin = omp_get_wtime();

                // RK4-2
                #pragma omp for private(g1,g2,xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)

                for (int i = 0; i < ExBD.size(); i++)  {
                    g1 = (int)(ExBD[i] / M1);
                    g2 = (int)(ExBD[i] % M1);
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];
                    f0 = F[g1*W1+g2];
                    f1p1 = F[(g1+1)*W1+g2];
                    f1m1 = F[(g1-1)*W1+g2];
                    f2p1 = F[g1*W1+(g2+1)];
                    f2m1 = F[g1*W1+(g2-1)];
                    f1p2 = F[(g1+2)*W1+g2];
                    f1m2 = F[(g1-2)*W1+g2];
                    f2p2 = F[g1*W1+(g2+2)];
                    f2m2 = F[g1*W1+(g2-2)];
#if defined(FDO4) || defined(FDO6)   
                    f1p3 = F[(g1+3)*W1+g2];
                    f1m3 = F[(g1-3)*W1+g2];
                    f2p3 = F[g1*W1+(g2+3)];
                    f2m3 = F[g1*W1+(g2-3)];
#endif
#if defined(FDO6) 
                    f1p4 = F[(g1+4)*W1+g2];
                    f1m4 = F[(g1-4)*W1+g2];
                    f2p4 = F[g1*W1+(g2+4)];
                    f2m4 = F[g1*W1+(g2-4)];
#endif
                    kk0 = KK1[g1*W1+g2];
                    kk1p1 = KK1[(g1+1)*W1+g2];
                    kk1m1 = KK1[(g1-1)*W1+g2];
                    kk2p1 = KK1[g1*W1+(g2+1)];
                    kk2m1 = KK1[g1*W1+(g2-1)];
                    kk1p2 = KK1[(g1+2)*W1+g2];
                    kk1m2 = KK1[(g1-2)*W1+g2];
                    kk2p2 = KK1[g1*W1+(g2+2)];
                    kk2m2 = KK1[g1*W1+(g2-2)];
#if defined(FDO4) || defined(FDO6)  
                    kk1p3 = KK1[(g1+3)*W1+g2];
                    kk1m3 = KK1[(g1-3)*W1+g2];
                    kk2p3 = KK1[g1*W1+(g2+3)];
                    kk2m3 = KK1[g1*W1+(g2-3)];
#endif
#if defined(FDO6) 
                    kk1p4 = KK1[(g1+4)*W1+g2];
                    kk1m4 = KK1[(g1-4)*W1+g2];
                    kk2p4 = KK1[g1*W1+(g2+4)];
                    kk2m4 = KK1[g1*W1+(g2-4)];
#endif
                    KK2[g1*W1+g2] = -kh0m * xx2 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                    k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) * 
#if defined(FDO4)
                                    (-(f2p3+0.5*kk2p3)/8.0+(f2p2+0.5*kk2p2)-13.0*(f2p1+0.5*kk2p1)/8.0+13.0*(f2m1+0.5*kk2m1)/8.0-(f2m2+0.5*kk2m2)+(f2m3+0.5*kk2m3)/8.0);
#elif defined(FDO6)
                                    (7.0/240.0*(f2p4+0.5*kk2p4)-3.0/10.0*(f2p3+0.5*kk2p3)+169.0/120.0*(f2p2+0.5*kk2p2)-61.0/30.0*(f2p1+0.5*kk2p1)+61.0/30.0*(f2m1+0.5*kk2m1)-169.0/120.0*(f2m2+0.5*kk2m2)+3.0/10.0*(f2m3+0.5*kk2m3)-7.0/240.0*(f2m4+0.5*kk2m4));
#else
                                    (0.5*(f2p2+0.5*kk2p2)-(f2p1+0.5*kk2p1)+(f2m1+0.5*kk2m1)-0.5*(f2m2+0.5*kk2m2));
#endif

                    FF[g1*W1+g2] += KK2[g1*W1+g2] / 3.0;
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-12: CASE 1 KK2) = %lf sec\n", t_1_elapsed);
                t_1_begin = omp_get_wtime();

                // RK4-3
                #pragma omp for private(g1,g2,xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)
                for (int i = 0; i < ExBD.size(); i++)  {
                    g1 = (int)(ExBD[i] / M1);
                    g2 = (int)(ExBD[i] % M1);
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];
                    f0 = F[g1*W1+g2];
                    f1p1 = F[(g1+1)*W1+g2];
                    f1m1 = F[(g1-1)*W1+g2];
                    f2p1 = F[g1*W1+(g2+1)];
                    f2m1 = F[g1*W1+(g2-1)];
                    f1p2 = F[(g1+2)*W1+g2];
                    f1m2 = F[(g1-2)*W1+g2];
                    f2p2 = F[g1*W1+(g2+2)];
                    f2m2 = F[g1*W1+(g2-2)];
#if defined(FDO4) || defined(FDO6)
                    f1p3 = F[(g1+3)*W1+g2];
                    f1m3 = F[(g1-3)*W1+g2];
                    f2p3 = F[g1*W1+(g2+3)];
                    f2m3 = F[g1*W1+(g2-3)];
#endif
#if defined(FDO6)     
                    f1p4 = F[(g1+4)*W1+g2];
                    f1m4 = F[(g1-4)*W1+g2];
                    f2p4 = F[g1*W1+(g2+4)];
                    f2m4 = F[g1*W1+(g2-4)];
#endif
                    kk0 = KK2[g1*W1+g2];
                    kk1p1 = KK2[(g1+1)*W1+g2];
                    kk1m1 = KK2[(g1-1)*W1+g2];
                    kk2p1 = KK2[g1*W1+(g2+1)];
                    kk2m1 = KK2[g1*W1+(g2-1)];
                    kk1p2 = KK2[(g1+2)*W1+g2];
                    kk1m2 = KK2[(g1-2)*W1+g2];
                    kk2p2 = KK2[g1*W1+(g2+2)];
                    kk2m2 = KK2[g1*W1+(g2-2)];
#if defined(FDO4) || defined(FDO6)
                    kk1p3 = KK2[(g1+3)*W1+g2];
                    kk1m3 = KK2[(g1-3)*W1+g2];
                    kk2p3 = KK2[g1*W1+(g2+3)];
                    kk2m3 = KK2[g1*W1+(g2-3)];
#endif
#if defined(FDO6) 
                    kk1p4 = KK2[(g1+4)*W1+g2];
                    kk1m4 = KK2[(g1-4)*W1+g2];
                    kk2p4 = KK2[g1*W1+(g2+4)];
                    kk2m4 = KK2[g1*W1+(g2-4)];
#endif
                    KK3[g1*W1+g2] = -kh0m * xx2 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                    k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) *        
#if defined(FDO4)
                                    (-(f2p3+0.5*kk2p3)/8.0+(f2p2+0.5*kk2p2)-13.0*(f2p1+0.5*kk2p1)/8.0+13.0*(f2m1+0.5*kk2m1)/8.0-(f2m2+0.5*kk2m2)+(f2m3+0.5*kk2m3)/8.0);
#elif defined(FDO6)
                                    (7.0/240.0*(f2p4+0.5*kk2p4)-3.0/10.0*(f2p3+0.5*kk2p3)+169.0/120.0*(f2p2+0.5*kk2p2)-61.0/30.0*(f2p1+0.5*kk2p1)+61.0/30.0*(f2m1+0.5*kk2m1)-169.0/120.0*(f2m2+0.5*kk2m2)+3.0/10.0*(f2m3+0.5*kk2m3)-7.0/240.0*(f2m4+0.5*kk2m4));
#else
                                    (0.5*(f2p2+0.5*kk2p2) - (f2p1+0.5*kk2p1) + (f2m1+0.5*kk2m1) - 0.5*(f2m2+0.5*kk2m2));
#endif
                    FF[g1*W1+g2] += KK3[g1*W1+g2] / 3.0;    
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-13: CASE 1 KK3) = %lf sec\n", t_1_elapsed);
                t_1_begin = omp_get_wtime();

                // RK4-4
                #pragma omp for private(g1,g2,xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)
                for (int i = 0; i < ExBD.size(); i++)  {
                    g1 = (int)(ExBD[i] / M1);
                    g2 = (int)(ExBD[i] % M1);
                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];
                    f0 = F[g1*W1+g2];
                    f1p1 = F[(g1+1)*W1+g2];
                    f1m1 = F[(g1-1)*W1+g2];
                    f2p1 = F[g1*W1+(g2+1)];
                    f2m1 = F[g1*W1+(g2-1)];
                    f1p2 = F[(g1+2)*W1+g2];
                    f1m2 = F[(g1-2)*W1+g2];
                    f2p2 = F[g1*W1+(g2+2)];
                    f2m2 = F[g1*W1+(g2-2)];
#if defined(FDO4) || defined(FDO6)
                    f1p3 = F[(g1+3)*W1+g2];
                    f1m3 = F[(g1-3)*W1+g2];
                    f2p3 = F[g1*W1+(g2+3)];
                    f2m3 = F[g1*W1+(g2-3)];
#endif
#if defined(FDO6)
                    f1p4 = F[(g1+4)*W1+g2];
                    f1m4 = F[(g1-4)*W1+g2];
                    f2p4 = F[g1*W1+(g2+4)];
                    f2m4 = F[g1*W1+(g2-4)];
#endif
                    kk0 = KK3[g1*W1+g2];
                    kk1p1 = KK3[(g1+1)*W1+g2];
                    kk1m1 = KK3[(g1-1)*W1+g2];
                    kk2p1 = KK3[g1*W1+(g2+1)];
                    kk2m1 = KK3[g1*W1+(g2-1)];
                    kk1p2 = KK3[(g1+2)*W1+g2];
                    kk1m2 = KK3[(g1-2)*W1+g2];
                    kk2p2 = KK3[g1*W1+(g2+2)];
                    kk2m2 = KK3[g1*W1+(g2-2)];
#if defined(FDO4) || defined(FDO6)
                    kk1p3 = KK3[(g1+3)*W1+g2];
                    kk1m3 = KK3[(g1-3)*W1+g2];
                    kk2p3 = KK3[g1*W1+(g2+3)];
                    kk2m3 = KK3[g1*W1+(g2-3)];
#endif
#if defined(FDO6)
                    kk1p4 = KK3[(g1+4)*W1+g2];
                    kk1m4 = KK3[(g1-4)*W1+g2];
                    kk2p4 = KK3[g1*W1+(g2+4)];
                    kk2m4 = KK3[g1*W1+(g2-4)];
#endif
                    KK4[g1*W1+g2] = -kh0m * xx2 * (-1/12.0*(f1p2+kk1p2) + 2/3.0*(f1p1+kk1p1) - 2/3.0*(f1m1+kk1m1) + 1/12.0*(f1m2+kk1m2)) + 
                                    k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+kk2p2) + 2/3.0*(f2p1+kk2p1) - 2/3.0*(f2m1+kk2m1) + 1/12.0*(f2m2+kk2m2)) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) * 
#if defined(FDO4)
                                    (-(f2p3+kk2p3)/8.0+(f2p2+kk2p2)-13.0*(f2p1+kk2p1)/8.0+13.0*(f2m1+kk2m1)/8.0-(f2m2+kk2m2)+(f2m3+kk2m3)/8.0);
#elif defined(FDO6)
                                    (7.0/240.0*(f2p4+kk2p4)-3.0/10.0*(f2p3+kk2p3)+169.0/120.0*(f2p2+kk2p2)-61.0/30.0*(f2p1+kk2p1)+61.0/30.0*(f2m1+kk2m1)-169.0/120.0*(f2m2+kk2m2)+3.0/10.0*(f2m3+kk2m3)-7.0/240.0*(f2m4+kk2m4));
#else
                                    (0.5*(f2p2+kk2p2)-(f2p1+kk2p1)+(f2m1+kk2m1)-0.5*(f2m2+kk2m2));
#endif
                    FF[g1*W1+g2] += KK4[g1*W1+g2] / 6.0;
                }
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-14: CASE 1 KK4) = %lf sec\n", t_1_elapsed);
                t_1_begin = omp_get_wtime();
                ExBD.clear();
            }

            // Check Multiple Expanding 
            // TBL = index of FF that FF(TBL) is higher than TolL

            t_1_begin = omp_get_wtime();
            TBL.clear();

            #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2, b1, b2, b3, b4, b5)
            for (int i = 0; i < ExFF.size(); i++)
            {
                g1 = (int)(ExFF[i] / M1);
                g2 = (int)(ExFF[i] % M1);
                b1 = std::abs(FF[g1*W1+g2]) >= TolH;
                b2 = std::abs(FF[(g1+1)*W1+g2] - FF[(g1-1)*W1+g2]) >= TolHdX1;
                b3 = std::abs(FF[g1*W1+(g2+1)] - FF[g1*W1+(g2-1)]) >= TolHdX2;
                b4 = g1 > EDGE && g2 > EDGE;
                b5 = g1 < BoxShape[0] - EDGE - 1 && g2 < BoxShape[1] - EDGE - 1;

                if (  ( b1 || b2 || b3 ) && b4 && b5 )  {
                    tmpVec.push_back(ExFF[i]);
                }
            }        
            tmpVec.swap(TBL);
            tmpVec.clear(); 

            // TBL & TBL_P set difference
            tmpVec.resize(TBL_P.size() + TBL.size());
            __gnu_parallel::sort (TBL.begin(), TBL.end());
            __gnu_parallel::sort (TBL_P.begin(), TBL_P.end());
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
            t_overhead += t_1_elapsed;
            //if (!QUIET) log->log("TBL size = %d TBL_P size = %d\n", TBL.size(), TBL_P.size()); 
            if (!QUIET && TIMING) log->log("Elapsed time (omp-c-3 CASE 1 TBL) = %lf sec\n", t_1_elapsed); 
        }
        // .........................................................................................

        // CASE 2: Truncating without extrapolation

        if ( !isExtrapolate && !isFullGrid )
        {
            // RK4-1
            #pragma omp parallel
            {
                #pragma omp single nowait
                {
                    t_1_begin = omp_get_wtime();
                }
                #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4) schedule(runtime)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            xx1 = Box[0] + i1 * H[0];
                            xx2 = Box[2] + i2 * H[1];
                            f0 = F[i1*W1+i2];
                            f1p1 = F[(i1+1)*W1+i2];
                            f1m1 = F[(i1-1)*W1+i2];
                            f2p1 = F[i1*W1+(i2+1)];
                            f2m1 = F[i1*W1+(i2-1)];
                            f1p2 = F[(i1+2)*W1+i2];
                            f1m2 = F[(i1-2)*W1+i2];
                            f2p2 = F[i1*W1+(i2+2)];
                            f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)     
                            f1p3 = F[(i1+3)*W1+i2];
                            f1m3 = F[(i1-3)*W1+i2];
                            f2p3 = F[i1*W1+(i2+3)];
                            f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                            f1p4 = F[(i1+4)*W1+i2];
                            f1m4 = F[(i1-4)*W1+i2];
                            f2p4 = F[i1*W1+(i2+4)];
                            f2m4 = F[i1*W1+(i2-4)];
#endif
                            KK1[i1*W1+i2] = -kh0m * xx2 * (-f1p2/12.0 + 2/3.0*f1p1 - 2/3.0*f1m1 + f1m2/12.0) + 
                                            k2h1 * POTENTIAL_X(xx1, xx2) * (-f2p2/12.0 + 2/3.0*f2p1 - 2/3.0*f2m1 + f2m2/12.0) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) *                                         
#if defined(FDO4)
                                            (-f2p3/8.0+f2p2-13.0*f2p1/8.0+13*f2m1/8.0-f2m2+f2m3/8.0);
#elif defined(FDO6)                                               
                                            (7.0/240.0*f2p4-3.0/10.0*f2p3+169.0/120.0*f2p2-61.0/30.0*f2p1+61.0/30.0*f2m1-169.0/120.0*f2m2+3.0/10.0*f2m3-7.0/240.0*f2m4);
#else
                                            (0.5*f2p2-f2p1+f2m1-0.5*f2m2);
#endif                            
                            FF[i1*W1+i2] = F[i1*W1+i2] + KK1[i1*W1+i2] / 6.0;
                        }
                    }
                }
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-21: CASE 2 KK1) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-2
                #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            xx1 = Box[0] + i1 * H[0];
                            xx2 = Box[2] + i2 * H[1];
                            f0 = F[i1*W1+i2];
                            f1p1 = F[(i1+1)*W1+i2];
                            f1m1 = F[(i1-1)*W1+i2];
                            f2p1 = F[i1*W1+(i2+1)];
                            f2m1 = F[i1*W1+(i2-1)];
                            f1p2 = F[(i1+2)*W1+i2];
                            f1m2 = F[(i1-2)*W1+i2];
                            f2p2 = F[i1*W1+(i2+2)];
                            f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)   
                            f1p3 = F[(i1+3)*W1+i2];
                            f1m3 = F[(i1-3)*W1+i2];
                            f2p3 = F[i1*W1+(i2+3)];
                            f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                            f1p4 = F[(i1+4)*W1+i2];
                            f1m4 = F[(i1-4)*W1+i2];
                            f2p4 = F[i1*W1+(i2+4)];
                            f2m4 = F[i1*W1+(i2-4)];
#endif
                            kk0 = KK1[i1*W1+i2];
                            kk1p1 = KK1[(i1+1)*W1+i2];
                            kk1m1 = KK1[(i1-1)*W1+i2];
                            kk2p1 = KK1[i1*W1+(i2+1)];
                            kk2m1 = KK1[i1*W1+(i2-1)];
                            kk1p2 = KK1[(i1+2)*W1+i2];
                            kk1m2 = KK1[(i1-2)*W1+i2];
                            kk2p2 = KK1[i1*W1+(i2+2)];
                            kk2m2 = KK1[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)  
                            kk1p3 = KK1[(i1+3)*W1+i2];
                            kk1m3 = KK1[(i1-3)*W1+i2];
                            kk2p3 = KK1[i1*W1+(i2+3)];
                            kk2m3 = KK1[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                            kk1p4 = KK1[(i1+4)*W1+i2];
                            kk1m4 = KK1[(i1-4)*W1+i2];
                            kk2p4 = KK1[i1*W1+(i2+4)];
                            kk2m4 = KK1[i1*W1+(i2-4)];
#endif
                            KK2[i1*W1+i2] = -kh0m * xx2 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                            k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) *                                   
#if defined(FDO4)
                                            (-(f2p3+0.5*kk2p3)/8.0+(f2p2+0.5*kk2p2)-13.0*(f2p1+0.5*kk2p1)/8.0+13.0*(f2m1+0.5*kk2m1)/8.0-(f2m2+0.5*kk2m2)+(f2m3+0.5*kk2m3)/8.0);
#elif defined(FDO6)
                                            (7.0/240.0*(f2p4+0.5*kk2p4)-3.0/10.0*(f2p3+0.5*kk2p3)+169.0/120.0*(f2p2+0.5*kk2p2)-61.0/30.0*(f2p1+0.5*kk2p1)+61.0/30.0*(f2m1+0.5*kk2m1)-169.0/120.0*(f2m2+0.5*kk2m2)+3.0/10.0*(f2m3+0.5*kk2m3)-7.0/240.0*(f2m4+0.5*kk2m4));
#else
                                            (0.5*(f2p2+0.5*kk2p2)-(f2p1+0.5*kk2p1)+(f2m1+0.5*kk2m1)-0.5*(f2m2+0.5*kk2m2));
#endif
                            FF[i1*W1+i2] += KK2[i1*W1+i2] / 3.0;
                        }
                    }
                }
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-22: CASE 2 KK2) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-3
                #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            xx1 = Box[0] + i1 * H[0];
                            xx2 = Box[2] + i2 * H[1];
                            f0 = F[i1*W1+i2];
                            f1p1 = F[(i1+1)*W1+i2];
                            f1m1 = F[(i1-1)*W1+i2];
                            f2p1 = F[i1*W1+(i2+1)];
                            f2m1 = F[i1*W1+(i2-1)];
                            f1p2 = F[(i1+2)*W1+i2];
                            f1m2 = F[(i1-2)*W1+i2];
                            f2p2 = F[i1*W1+(i2+2)];
                            f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                            f1p3 = F[(i1+3)*W1+i2];
                            f1m3 = F[(i1-3)*W1+i2];
                            f2p3 = F[i1*W1+(i2+3)];
                            f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6)  
                            f1p4 = F[(i1+4)*W1+i2];
                            f1m4 = F[(i1-4)*W1+i2];
                            f2p4 = F[i1*W1+(i2+4)];
                            f2m4 = F[i1*W1+(i2-4)];
#endif
                            kk0 = KK2[i1*W1+i2];
                            kk1p1 = KK2[(i1+1)*W1+i2];
                            kk1m1 = KK2[(i1-1)*W1+i2];
                            kk2p1 = KK2[i1*W1+(i2+1)];
                            kk2m1 = KK2[i1*W1+(i2-1)];
                            kk1p2 = KK2[(i1+2)*W1+i2];
                            kk1m2 = KK2[(i1-2)*W1+i2];
                            kk2p2 = KK2[i1*W1+(i2+2)];
                            kk2m2 = KK2[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                            kk1p3 = KK2[(i1+3)*W1+i2];
                            kk1m3 = KK2[(i1-3)*W1+i2];
                            kk2p3 = KK2[i1*W1+(i2+3)];
                            kk2m3 = KK2[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                            kk1p4 = KK2[(i1+4)*W1+i2];
                            kk1m4 = KK2[(i1-4)*W1+i2];
                            kk2p4 = KK2[i1*W1+(i2+4)];
                            kk2m4 = KK2[i1*W1+(i2-4)];
#endif
                            KK3[i1*W1+i2] = -kh0m * xx2 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                            k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) *                           
#if defined(FDO4)
                                            (-(f2p3+0.5*kk2p3)/8.0+(f2p2+0.5*kk2p2)-13.0*(f2p1+0.5*kk2p1)/8.0+13.0*(f2m1+0.5*kk2m1)/8.0-(f2m2+0.5*kk2m2)+(f2m3+0.5*kk2m3)/8.0);
#elif defined(FDO6)
                                            (7.0/240.0*(f2p4+0.5*kk2p4)-3.0/10.0*(f2p3+0.5*kk2p3)+169.0/120.0*(f2p2+0.5*kk2p2)-61.0/30.0*(f2p1+0.5*kk2p1)+61.0/30.0*(f2m1+0.5*kk2m1)-169.0/120.0*(f2m2+0.5*kk2m2)+3.0/10.0*(f2m3+0.5*kk2m3)-7.0/240.0*(f2m4+0.5*kk2m4));
#else
                                            (0.5*(f2p2+0.5*kk2p2) - (f2p1+0.5*kk2p1) + (f2m1+0.5*kk2m1) - 0.5*(f2m2+0.5*kk2m2));
#endif
                            FF[i1*W1+i2] += KK3[i1*W1+i2] / 3.0;
                        }
                    }
                }
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-23: CASE 2 KK3) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-4
                #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            xx1 = Box[0] + i1 * H[0];
                            xx2 = Box[2] + i2 * H[1];
                            f0 = F[i1*W1+i2];
                            f1p1 = F[(i1+1)*W1+i2];
                            f1m1 = F[(i1-1)*W1+i2];
                            f2p1 = F[i1*W1+(i2+1)];
                            f2m1 = F[i1*W1+(i2-1)];
                            f1p2 = F[(i1+2)*W1+i2];
                            f1m2 = F[(i1-2)*W1+i2];
                            f2p2 = F[i1*W1+(i2+2)];
                            f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                            f1p3 = F[(i1+3)*W1+i2];
                            f1m3 = F[(i1-3)*W1+i2];
                            f2p3 = F[i1*W1+(i2+3)];
                            f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6)
                            f1p4 = F[(i1+4)*W1+i2];
                            f1m4 = F[(i1-4)*W1+i2];
                            f2p4 = F[i1*W1+(i2+4)];
                            f2m4 = F[i1*W1+(i2-4)];
#endif
                            kk0 = KK3[i1*W1+i2];
                            kk1p1 = KK3[(i1+1)*W1+i2];
                            kk1m1 = KK3[(i1-1)*W1+i2];
                            kk2p1 = KK3[i1*W1+(i2+1)];
                            kk2m1 = KK3[i1*W1+(i2-1)];
                            kk1p2 = KK3[(i1+2)*W1+i2];
                            kk1m2 = KK3[(i1-2)*W1+i2];
                            kk2p2 = KK3[i1*W1+(i2+2)];
                            kk2m2 = KK3[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                            kk1p3 = KK3[(i1+3)*W1+i2];
                            kk1m3 = KK3[(i1-3)*W1+i2];
                            kk2p3 = KK3[i1*W1+(i2+3)];
                            kk2m3 = KK3[i1*W1+(i2-3)];
#endif
#if defined(FDO6)
                            kk1p4 = KK3[(i1+4)*W1+i2];
                            kk1m4 = KK3[(i1-4)*W1+i2];
                            kk2p4 = KK3[i1*W1+(i2+4)];
                            kk2m4 = KK3[i1*W1+(i2-4)];
#endif
                            KK4[i1*W1+i2] = -kh0m * xx2 * (-1/12.0*(f1p2+kk1p2) + 2/3.0*(f1p1+kk1p1) - 2/3.0*(f1m1+kk1m1) + 1/12.0*(f1m2+kk1m2)) + 
                                            k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+kk2p2) + 2/3.0*(f2p1+kk2p1) - 2/3.0*(f2m1+kk2m1) + 1/12.0*(f2m2+kk2m2)) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) *                                         
#if defined(FDO4)
                                            (-(f2p3+kk2p3)/8.0+(f2p2+kk2p2)-13.0*(f2p1+kk2p1)/8.0+13.0*(f2m1+kk2m1)/8.0-(f2m2+kk2m2)+(f2m3+kk2m3)/8.0);
#elif defined(FDO6)
                                            (7.0/240.0*(f2p4+kk2p4)-3.0/10.0*(f2p3+kk2p3)+169.0/120.0*(f2p2+kk2p2)-61.0/30.0*(f2p1+kk2p1)+61.0/30.0*(f2m1+kk2m1)-169.0/120.0*(f2m2+kk2m2)+3.0/10.0*(f2m3+kk2m3)-7.0/240.0*(f2m4+kk2m4));
#else
                                            (0.5*(f2p2+kk2p2)-(f2p1+kk2p1)+(f2m1+kk2m1)-0.5*(f2m2+kk2m2));
#endif
                            FF[i1*W1+i2] += KK4[i1*W1+i2] / 6.0;
                        }
                    }                
                }
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-24: CASE 2 KK4) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }
            } // OMP PARALLEL
        } 
        else if ( !isExtrapolate && isFullGrid )
        {
            // .........................................................................................

            // CASE 3: Full grid

            // RK4-1
            #pragma omp parallel
            {
                #pragma omp single nowait
                {
                    t_1_begin = omp_get_wtime();
                }
                #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4) schedule(runtime)

                for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                        xx1 = Box[0] + i1 * H[0];
                        xx2 = Box[2] + i2 * H[1];
                        f0 = F[i1*W1+i2];
                        f1p1 = F[(i1+1)*W1+i2];
                        f1m1 = F[(i1-1)*W1+i2];
                        f2p1 = F[i1*W1+(i2+1)];
                        f2m1 = F[i1*W1+(i2-1)];
                        f1p2 = F[(i1+2)*W1+i2];
                        f1m2 = F[(i1-2)*W1+i2];
                        f2p2 = F[i1*W1+(i2+2)];
                        f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)     
                        f1p3 = F[(i1+3)*W1+i2];
                        f1m3 = F[(i1-3)*W1+i2];
                        f2p3 = F[i1*W1+(i2+3)];
                        f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                        f1p4 = F[(i1+4)*W1+i2];
                        f1m4 = F[(i1-4)*W1+i2];
                        f2p4 = F[i1*W1+(i2+4)];
                        f2m4 = F[i1*W1+(i2-4)];
#endif
                        KK1[i1*W1+i2] = -kh0m * xx2 * (-f1p2/12.0 + 2/3.0*f1p1 - 2/3.0*f1m1 + f1m2/12.0) + 
                                        k2h1 * POTENTIAL_X(xx1, xx2) * (-f2p2/12.0 + 2/3.0*f2p1 - 2/3.0*f2m1 + f2m2/12.0) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) *       
#if defined(FDO4)
                                        (-f2p3/8.0+f2p2-13.0*f2p1/8.0+13*f2m1/8.0-f2m2+f2m3/8.0);
#elif defined(FDO6)                                               
                                        (7.0/240.0*f2p4-3.0/10.0*f2p3+169.0/120.0*f2p2-61.0/30.0*f2p1+61.0/30.0*f2m1-169.0/120.0*f2m2+3.0/10.0*f2m3-7.0/240.0*f2m4);
#else
                                        (0.5*f2p2-f2p1+f2m1-0.5*f2m2);
#endif
                        FF[i1*W1+i2] = F[i1*W1+i2] + KK1[i1*W1+i2] / 6.0;
                    }
                }
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-31: CASE 3 KK1) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }
                // RK4-2
                #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)
                for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                        xx1 = Box[0] + i1 * H[0];
                        xx2 = Box[2] + i2 * H[1];
                        f0 = F[i1*W1+i2];
                        f1p1 = F[(i1+1)*W1+i2];
                        f1m1 = F[(i1-1)*W1+i2];
                        f2p1 = F[i1*W1+(i2+1)];
                        f2m1 = F[i1*W1+(i2-1)];
                        f1p2 = F[(i1+2)*W1+i2];
                        f1m2 = F[(i1-2)*W1+i2];
                        f2p2 = F[i1*W1+(i2+2)];
                        f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6) 
                        f1p3 = F[(i1+3)*W1+i2];
                        f1m3 = F[(i1-3)*W1+i2];
                        f2p3 = F[i1*W1+(i2+3)];
                        f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                        f1p4 = F[(i1+4)*W1+i2];
                        f1m4 = F[(i1-4)*W1+i2];
                        f2p4 = F[i1*W1+(i2+4)];
                        f2m4 = F[i1*W1+(i2-4)];
#endif
                        kk0 = KK1[i1*W1+i2];
                        kk1p1 = KK1[(i1+1)*W1+i2];
                        kk1m1 = KK1[(i1-1)*W1+i2];
                        kk2p1 = KK1[i1*W1+(i2+1)];
                        kk2m1 = KK1[i1*W1+(i2-1)];
                        kk1p2 = KK1[(i1+2)*W1+i2];
                        kk1m2 = KK1[(i1-2)*W1+i2];
                        kk2p2 = KK1[i1*W1+(i2+2)];
                        kk2m2 = KK1[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6) 
                        kk1p3 = KK1[(i1+3)*W1+i2];
                        kk1m3 = KK1[(i1-3)*W1+i2];
                        kk2p3 = KK1[i1*W1+(i2+3)];
                        kk2m3 = KK1[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                        kk1p4 = KK1[(i1+4)*W1+i2];
                        kk1m4 = KK1[(i1-4)*W1+i2];
                        kk2p4 = KK1[i1*W1+(i2+4)];
                        kk2m4 = KK1[i1*W1+(i2-4)];
#endif
                        KK2[i1*W1+i2] = -kh0m * xx2 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                        k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) * 
#if defined(FDO4)
                                        (-(f2p3+0.5*kk2p3)/8.0+(f2p2+0.5*kk2p2)-13.0*(f2p1+0.5*kk2p1)/8.0+13.0*(f2m1+0.5*kk2m1)/8.0-(f2m2+0.5*kk2m2)+(f2m3+0.5*kk2m3)/8.0);
#elif defined(FDO6)
                                        (7.0/240.0*(f2p4+0.5*kk2p4)-3.0/10.0*(f2p3+0.5*kk2p3)+169.0/120.0*(f2p2+0.5*kk2p2)-61.0/30.0*(f2p1+0.5*kk2p1)+61.0/30.0*(f2m1+0.5*kk2m1)-169.0/120.0*(f2m2+0.5*kk2m2)+3.0/10.0*(f2m3+0.5*kk2m3)-7.0/240.0*(f2m4+0.5*kk2m4));
#else
                                        (0.5*(f2p2+0.5*kk2p2)-(f2p1+0.5*kk2p1)+(f2m1+0.5*kk2m1)-0.5*(f2m2+0.5*kk2m2));
#endif
                        FF[i1*W1+i2] += KK2[i1*W1+i2] / 3.0;
                    }
                }
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-32: CASE 3 KK2) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-3
                #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)
                for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                        xx1 = Box[0] + i1 * H[0];
                        xx2 = Box[2] + i2 * H[1];
                        f0 = F[i1*W1+i2];
                        f1p1 = F[(i1+1)*W1+i2];
                        f1m1 = F[(i1-1)*W1+i2];
                        f2p1 = F[i1*W1+(i2+1)];
                        f2m1 = F[i1*W1+(i2-1)];
                        f1p2 = F[(i1+2)*W1+i2];
                        f1m2 = F[(i1-2)*W1+i2];
                        f2p2 = F[i1*W1+(i2+2)];
                        f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                        f1p3 = F[(i1+3)*W1+i2];
                        f1m3 = F[(i1-3)*W1+i2];
                        f2p3 = F[i1*W1+(i2+3)];
                        f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                        f1p4 = F[(i1+4)*W1+i2];
                        f1m4 = F[(i1-4)*W1+i2];
                        f2p4 = F[i1*W1+(i2+4)];
                        f2m4 = F[i1*W1+(i2-4)];
#endif
                        kk0 = KK2[i1*W1+i2];
                        kk1p1 = KK2[(i1+1)*W1+i2];
                        kk1m1 = KK2[(i1-1)*W1+i2];
                        kk2p1 = KK2[i1*W1+(i2+1)];
                        kk2m1 = KK2[i1*W1+(i2-1)];
                        kk1p2 = KK2[(i1+2)*W1+i2];
                        kk1m2 = KK2[(i1-2)*W1+i2];
                        kk2p2 = KK2[i1*W1+(i2+2)];
                        kk2m2 = KK2[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                        kk1p3 = KK2[(i1+3)*W1+i2];
                        kk1m3 = KK2[(i1-3)*W1+i2];
                        kk2p3 = KK2[i1*W1+(i2+3)];
                        kk2m3 = KK2[i1*W1+(i2-3)];
#endif
#if defined(FDO6) 
                        kk1p4 = KK2[(i1+4)*W1+i2];
                        kk1m4 = KK2[(i1-4)*W1+i2];
                        kk2p4 = KK2[i1*W1+(i2+4)];
                        kk2m4 = KK2[i1*W1+(i2-4)];
#endif
                        KK3[i1*W1+i2] = -kh0m * xx2 * (-1/12.0*(f1p2+0.5*kk1p2) + 2/3.0*(f1p1+0.5*kk1p1) - 2/3.0*(f1m1+0.5*kk1m1) + 1/12.0*(f1m2+0.5*kk1m2)) + 
                                        k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+0.5*kk2p2) + 2/3.0*(f2p1+0.5*kk2p1) - 2/3.0*(f2m1+0.5*kk2m1) + 1/12.0*(f2m2+0.5*kk2m2)) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) *           
#if defined(FDO4)
                                        (-(f2p3+0.5*kk2p3)/8.0+(f2p2+0.5*kk2p2)-13.0*(f2p1+0.5*kk2p1)/8.0+13.0*(f2m1+0.5*kk2m1)/8.0-(f2m2+0.5*kk2m2)+(f2m3+0.5*kk2m3)/8.0);
#elif defined(FDO6)
                                        (7.0/240.0*(f2p4+0.5*kk2p4)-3.0/10.0*(f2p3+0.5*kk2p3)+169.0/120.0*(f2p2+0.5*kk2p2)-61.0/30.0*(f2p1+0.5*kk2p1)+61.0/30.0*(f2m1+0.5*kk2m1)-169.0/120.0*(f2m2+0.5*kk2m2)+3.0/10.0*(f2m3+0.5*kk2m3)-7.0/240.0*(f2m4+0.5*kk2m4));
#else
                                        (0.5*(f2p2+0.5*kk2p2) - (f2p1+0.5*kk2p1) + (f2m1+0.5*kk2m1) - 0.5*(f2m2+0.5*kk2m2));
#endif
                        FF[i1*W1+i2] += KK3[i1*W1+i2] / 3.0;  
                    }
                }
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-33: CASE 3 KK3) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-4
                #pragma omp for private(xx1,xx2,f0,f1p1,f1m1,f2p1,f2m1,f1p2,f1m2,f2p2,f2m2,f1p3,f1m3,f2p3,f2m3,f1p4,f1m4,f2p4,f2m4,kk0,kk1p1,kk1m1,kk2p1,kk2m1,kk1p2,kk1m2,kk2p2,kk2m2,kk1p3,kk1m3,kk2p3,kk2m3,kk1p4,kk1m4,kk2p4,kk2m4) schedule(runtime)
                for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                        xx1 = Box[0] + i1 * H[0];
                        xx2 = Box[2] + i2 * H[1];
                        f0 = F[i1*W1+i2];
                        f1p1 = F[(i1+1)*W1+i2];
                        f1m1 = F[(i1-1)*W1+i2];
                        f2p1 = F[i1*W1+(i2+1)];
                        f2m1 = F[i1*W1+(i2-1)];
                        f1p2 = F[(i1+2)*W1+i2];
                        f1m2 = F[(i1-2)*W1+i2];
                        f2p2 = F[i1*W1+(i2+2)];
                        f2m2 = F[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                        f1p3 = F[(i1+3)*W1+i2];
                        f1m3 = F[(i1-3)*W1+i2];
                        f2p3 = F[i1*W1+(i2+3)];
                        f2m3 = F[i1*W1+(i2-3)];
#endif
#if defined(FDO6)
                        f1p4 = F[(i1+4)*W1+i2];
                        f1m4 = F[(i1-4)*W1+i2];
                        f2p4 = F[i1*W1+(i2+4)];
                        f2m4 = F[i1*W1+(i2-4)];
#endif
                        kk0 = KK3[i1*W1+i2];
                        kk1p1 = KK3[(i1+1)*W1+i2];
                        kk1m1 = KK3[(i1-1)*W1+i2];
                        kk2p1 = KK3[i1*W1+(i2+1)];
                        kk2m1 = KK3[i1*W1+(i2-1)];
                        kk1p2 = KK3[(i1+2)*W1+i2];
                        kk1m2 = KK3[(i1-2)*W1+i2];
                        kk2p2 = KK3[i1*W1+(i2+2)];
                        kk2m2 = KK3[i1*W1+(i2-2)];
#if defined(FDO4) || defined(FDO6)
                        kk1p3 = KK3[(i1+3)*W1+i2];
                        kk1m3 = KK3[(i1-3)*W1+i2];
                        kk2p3 = KK3[i1*W1+(i2+3)];
                        kk2m3 = KK3[i1*W1+(i2-3)];
#endif
#if defined(FDO6)
                        kk1p4 = KK3[(i1+4)*W1+i2];
                        kk1m4 = KK3[(i1-4)*W1+i2];
                        kk2p4 = KK3[i1*W1+(i2+4)];
                        kk2m4 = KK3[i1*W1+(i2-4)];
#endif
                        KK4[i1*W1+i2] = -kh0m * xx2 * (-1/12.0*(f1p2+kk1p2) + 2/3.0*(f1p1+kk1p1) - 2/3.0*(f1m1+kk1m1) + 1/12.0*(f1m2+kk1m2)) + 
                                        k2h1 * POTENTIAL_X(xx1, xx2) * (-1/12.0*(f2p2+kk2p2) + 2/3.0*(f2p1+kk2p1) - 2/3.0*(f2m1+kk2m1) + 1/12.0*(f2m2+kk2m2)) - khbsq2h1 * POTENTIAL_XXX(xx1, xx2) *                               
#if defined(FDO4)
                                        (-(f2p3+kk2p3)/8.0+(f2p2+kk2p2)-13.0*(f2p1+kk2p1)/8.0+13.0*(f2m1+kk2m1)/8.0-(f2m2+kk2m2)+(f2m3+kk2m3)/8.0);
#elif defined(FDO6)
                                        (7.0/240.0*(f2p4+kk2p4)-3.0/10.0*(f2p3+kk2p3)+169.0/120.0*(f2p2+kk2p2)-61.0/30.0*(f2p1+kk2p1)+61.0/30.0*(f2m1+kk2m1)-169.0/120.0*(f2m2+kk2m2)+3.0/10.0*(f2m3+kk2m3)-7.0/240.0*(f2m4+kk2m4));
#else
                                        (0.5*(f2p2+kk2p2)-(f2p1+kk2p1)+(f2m1+kk2m1)-0.5*(f2m2+kk2m2));
#endif
                        FF[i1*W1+i2] += KK4[i1*W1+i2] / 6.0; 
                    }
                }
                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-34: CASE 3 KK4) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }
            }
        }
        // .........................................................................................

        // FF(t+1) Normailzed & go on

        t_1_begin = omp_get_wtime();
        norm = 0.0; 

        if (!isFullGrid)  {

            #pragma omp parallel for reduction (+:norm) schedule(runtime)
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    if (TAMask[i1*W1+i2])
                        norm += FF[i1*W1+i2];
                }
            }
        }  
        else  {

            #pragma omp parallel for reduction (+:norm) schedule(runtime)
            for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                    norm += FF[i1*W1+i2];
                }
            }
        }
        norm *= H[0] * H[1];

        if ( (tt + 1) % PERIOD == 0 )
            log->log("[WignerMoyal2d] Normalization factor = %.16e\n",norm);

        norm = 1.0 / norm; 

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_full += t_1_elapsed;
        t_truncate += t_1_elapsed;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-1-1 Norm) = %lf sec\n", t_1_elapsed);
        t_1_begin = omp_get_wtime();

        if (!isFullGrid)  {
            if (tt % SORT_PERIOD == 0)  {
                #pragma omp parallel for schedule(runtime)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if ( i1 < skin || i1 >= BoxShape[0] - skin || i2 < skin || i2 >= BoxShape[1] - skin )
                            FF[i1*W1+i2] = FF[i1*W1+i2] * (i1 >= skin ? 1.0 : exp(-lambda*pow(skin-i1,2))) * (i2 >= skin ? 1.0 : exp(-lambda*pow(skin-i2,2))) * (i1 <= BoxShape[0] - skin ? 1.0 : exp(-lambda*pow(i1-BoxShape[0]+skin,2))) * (i2 <= BoxShape[1] - skin ? 1.0 : exp(-lambda*pow(i2-BoxShape[1]+skin,2)));
                    }
                }
            }
            #pragma omp parallel for private(val) schedule(runtime)
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    val = norm * FF[i1*W1+i2];
                    FF[i1*W1+i2] = val;
                    F[i1*W1+i2] = val;
                    PF[i1*W1+i2] = val;
                }
            }
        }  
        else  {
            if (tt % SORT_PERIOD == 0)  {
                #pragma omp parallel for schedule(runtime)
                for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                        if ( i1 < skin || i1 >= BoxShape[0] - skin || i2 < skin || i2 >= BoxShape[1] - skin )
                            FF[i1*W1+i2] = FF[i1*W1+i2] * (i1 >= skin ? 1.0 : exp(-lambda*pow(skin-i1,2))) * (i2 >= skin ? 1.0 : exp(-lambda*pow(skin-i2,2))) * (i1 <= BoxShape[0] - skin ? 1.0 : exp(-lambda*pow(i1-BoxShape[0]+skin,2))) * (i2 <= BoxShape[1] - skin ? 1.0 : exp(-lambda*pow(i2-BoxShape[1]+skin,2)));
                    }
                }
            }
            #pragma omp parallel for private(val) schedule(runtime)
            for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                    val = norm * FF[i1*W1+i2];
                    FF[i1*W1+i2] = val;
                    F[i1*W1+i2] = val;
                    PF[i1*W1+i2] = val;
                }
            }
        }
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_full += t_1_elapsed;
        t_truncate += t_1_elapsed;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-1-2 FF) = %lf sec\n", t_1_elapsed); 

        if ( (tt + 1) % PERIOD == 0 )
        {
            // REPORT MEASUREMENTS
            // ----------------------------------------------------------------------------
            // Compute Transmittance
            // ----------------------------------------------------------------------------

            if (isTrans)  {
            
                t_1_begin = omp_get_wtime();
                pftrans = 0.0;

                if (!isFullGrid)  {
                    #pragma omp parallel for reduction (+:pftrans)
                    for (int i1 = idx_x0; i1 <= x1_max; i1 ++)  {
                        for (int i2 = x2_min; i2 <= x2_max; i2 ++)
                            pftrans+=PF[i1*W1+i2];
                    }
                }
                else  {
                    #pragma omp parallel for reduction (+:pftrans)
                    for (int i1 = idx_x0; i1 < BoxShape[0]-EDGE; i1 ++)  {
                        for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)
                            pftrans+=PF[i1*W1+i2];
                    }
                }
                pftrans *= H[0] * H[1];
                PF_trans.push_back(pftrans);
                log->log("[WignerMoyal2d] idx_x0 = %d\n", idx_x0);
                log->log("[WignerMoyal2d] Time %lf, Trans = %.16e\n", ( tt + 1 ) * kk, pftrans);
                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                if (!QUIET && TIMING) log->log("Elapsed time (omp-x-2 trans) = %lf sec\n", t_1_elapsed); 
            }

            if (isCorr)  {

                // Density
                #pragma omp parallel for private(density)
                for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                    density = 0.0;
                    for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                        density += PF[i1*W1+i2]; 
                    }
                    Ft[i1] = density * H[1];
                }

                corr = 0.0;
                for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                    corr += Ft[i1] * F0[i1];
                }
                corr *= H[0];
                log->log("[WignerMoyal2d] Time %lf, Corr = %.16e\n", ( tt + 1 ) * kk, corr/corr_0);
            }
        }

        // Truncation and TA

        tmpVec.clear();

        if ( !isFullGrid )
        {
            t_1_begin = omp_get_wtime();
            #pragma omp parallel 
            {
                #pragma omp for reduction(merge: tmpVec) private(b1, b2, b3)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            b1 = std::abs(PF[i1*W1+i2]) < TolH;
                            b2 = std::abs(PF[(i1+1)*W1+i2]*TAMask[(i1+1)*W1+i2] - PF[(i1-1)*W1+i2]*TAMask[(i1-1)*W1+i2]) < TolHdX1;
                            b3 = std::abs(PF[i1*W1+(i2+1)]*TAMask[i1*W1+(i2+1)] - PF[i1*W1+(i2-1)]*TAMask[i1*W1+(i2-1)]) < TolHdX2;

                            if (b1 && b2 && b3) {
                                F[i1*W1+i2] = 0.0;
                                tmpVec.push_back(i1*M1+i2);
                            }
                        }
                    }
                } 

                #pragma omp for private(g1,g2)
                for (int i = 0; i < tmpVec.size(); i++)  {
                    g1 = (int)(tmpVec[i] / M1);
                    g2 = (int)(tmpVec[i] % M1);
                    TAMask[g1*W1+g2] = 0;
                    PF[g1*W1+g2] = 0.0;
                }
            }
            tmpVec.clear();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3 TA) = %lf sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            // Rebuild TA box
            if ( (tt + 1) %  SORT_PERIOD == 0 )  {

                x1_min = BoxShape[0];
                x2_min = BoxShape[1];
                x1_max = 0;
                x2_max = 0;

                #pragma omp parallel for reduction(min: x1_min, x2_min) reduction(max: x1_max, x2_max) schedule(runtime)
                for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
                        if (TAMask[i1*W1+i2])  {
                            if (i1 < x1_min)  x1_min = i1;
                            if (i1 > x1_max)  x1_max = i1;
                            if (i2 < x2_min)  x2_min = i2;
                            if (i2 > x2_max)  x2_max = i2;
                        }
                    }
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-4 TARB) = %lf sec\n", t_1_elapsed);

            // TB
            t_1_begin = omp_get_wtime();
            tmpVec.clear();

            #pragma omp parallel for reduction(merge: tmpVec) schedule(runtime)
            for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                    if (TAMask[i1*W1+i2])  {
                        if (!TAMask[(i1-1)*W1+i2] || !TAMask[(i1+1)*W1+i2] || !TAMask[i1*W1+(i2-1)] || !TAMask[i1*W1+(i2+1)])
                            tmpVec.push_back(i1*M1+i2);
                    }
                }
            }
            tmpVec.swap(TB);
            tmpVec.clear();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-5 TB) = %lf sec\n", t_1_elapsed);

            // `````````````````````````````````````````````````````````````````

            // TA expansion
            t_1_begin = omp_get_wtime();
            
            #pragma omp parallel for reduction(merge: tmpVec) private(g1, g2) schedule(runtime)
            for (int i = 0; i < TB.size(); i++)
            {
                g1 = (int)(TB[i] / M1);
                g2 = (int)(TB[i] % M1);

                if ( g1+1 != BoxShape[0]-EDGE-1 && !TAMask[(g1+1)*W1+g2])
                    tmpVec.push_back(GridToIdx(g1+1,g2));
                if ( g1-1 != EDGE && !TAMask[(g1-1)*W1+g2])
                    tmpVec.push_back(GridToIdx(g1-1,g2));
                if ( g2+1 != BoxShape[1]-EDGE-1 && !TAMask[g1*W1+(g2+1)])
                    tmpVec.push_back(GridToIdx(g1,g2+1));
                if ( g2-1 != EDGE && !TAMask[g1*W1+(g2-1)])
                    tmpVec.push_back(GridToIdx(g1,g2-1));
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-6 TAEX) = %lf sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            #pragma omp parallel for reduction(min: x1_min, x2_min) reduction(max: x1_max, x2_max) private(g1, g2)
            for (int i = 0; i < tmpVec.size(); i ++)  {

                g1 = (int)(tmpVec[i] / M1);
                g2 = (int)(tmpVec[i] % M1);
                TAMask[g1*W1+g2] = 1;

                // Update TA box
                x1_min = (g1 < x1_min) ? g1 : x1_min;
                x2_min = (g2 < x2_min) ? g2 : x2_min;
                x1_max = (g1 > x1_max) ? g1 : x1_max;
                x2_max = (g2 > x2_max) ? g2 : x2_max;
            }
            tmpVec.clear();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-7 TARB) = %lf sec\n", t_1_elapsed);
        }

        // Reset
        t_1_begin = omp_get_wtime();

        if (!isFullGrid)  {
            #pragma omp parallel for
            for (int i1 = std::max(x1_min-ExLimit-DIMENSIONS, 1); i1 <= std::min(x1_max+ExLimit+DIMENSIONS, BoxShape[0]-1); i1 ++)  {
                for (int i2 = std::max(x2_min-ExLimit-DIMENSIONS, 1); i2 <= std::min(x2_max+ExLimit+DIMENSIONS, BoxShape[1]-1); i2 ++)  {
            
                    FF[i1*W1+i2] = 0.0;
                    KK1[i1*W1+i2] = 0.0;
                    KK2[i1*W1+i2] = 0.0;
                    KK3[i1*W1+i2] = 0.0;
                    KK4[i1*W1+i2] = 0.0;
                }
            }
        }  
        else  {
            #pragma omp parallel for
            for (int i1 = EDGE; i1 < BoxShape[0] - EDGE; i1 ++)  {
                for (int i2 = EDGE; i2 < BoxShape[1] - EDGE; i2 ++)  {
            
                    FF[i1*W1+i2] = 0.0;
                    KK1[i1*W1+i2] = 0.0;
                    KK2[i1*W1+i2] = 0.0;
                    KK3[i1*W1+i2] = 0.0;
                    KK4[i1*W1+i2] = 0.0;
                }
            }
        }
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_full += t_1_elapsed;
        t_truncate += t_1_elapsed;
        if ( !QUIET && TIMING ) log->log("Elapsed time (omp-e-8: reset) = %lf sec\n", t_1_elapsed);  

        if ( (tt + 1) % PERIOD == 0 )
        {   
            t_0_end = omp_get_wtime();
            t_0_elapsed = t_0_end - t_0_begin;

            if ( !QUIET ) log->log("[WignerMoyal2d] Step: %d, Elapsed time: %lf sec\n", tt + 1, t_0_elapsed);

            if ( !isFullGrid && !QUIET )  {

                ta_size = 0;
                #pragma omp parallel for reduction(+: ta_size) schedule(runtime)
                for (int i1 = x1_min; i1 <= x1_max; i1 ++)  {
                    for (int i2 = x2_min; i2 <= x2_max; i2 ++)  {
                        if (TAMask[i1*W1+i2])
                            ta_size += 1;    
                    }
                }
                tb_size = TB.size();
                log->log("[WignerMoyal2d] Core computation time = %lf\n", t_truncate);
                log->log("[WignerMoyal2d] Overhead time = %lf\n", t_overhead);
                log->log("[WignerMoyal2d] TA size = %d, TB size = %d\n", ta_size, tb_size);
                log->log("[WignerMoyal2d] TA Range [%d, %d][%d, %d]\n", x1_min, x1_max, x2_min, x2_max);
                log->log("[WignerMoyal2d] TA / total grids = %lf\n", ( ta_size * 1.0 ) / GRIDS_TOT);
                log->log("[WignerMoyal2d] ExCount = %d , ExLimit = %d\n", Excount, ExLimit);
            }
            else if ( isFullGrid && !QUIET )  {

                log->log("[WignerMoyal2d] Core computation time = %lf\n", t_full);
            }
            if ( !QUIET ) log->log("\n........................................................\n\n");
        }         
    } // Time iteration 

    delete F;
    delete FF;
    delete PF;
    delete KK1;
    delete KK2;
    delete KK3;
    delete KK4;

    if ( !isFullGrid )
        delete TAMask;

    log->log("[WignerMoyal2d] Evolve done.\n");
}
/* =============================================================================== */

/* WM2DPOT_DW1 */

inline double WignerMoyal2d::Wavefunction_DW1(double x1, double x2)
{
    return PI_INV / hb * exp(-2 * A[0] * pow(x1-Wave0[0],2)) * exp((-1/(2 * hb * hb * A[1])) * pow(x2 - Wave0[1],2));
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Potential_DW1(double x1, double x2)
{
    return  0.007 * (x1 * x1 * x1 * x1) - 0.01 * (x1 * x1);
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vx_DW1(double x1, double x2)
{
    return 0.028 * (x1 * x1 * x1) - 0.02 * x1;
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vxxx_DW1(double x1, double x2)
{
    return 0.168 * x1;
}
/* =============================================================================== */

/* WM2DPOT_DW2 */

inline double WignerMoyal2d::Wavefunction_DW2(double x1, double x2)
{
    return PI_INV / hb * exp(-2 * A[0] * pow(x1-Wave0[0],2)) * exp((-1/(2 * hb * hb * A[1])) * pow(x2 - Wave0[1],2));
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Potential_DW2(double x1, double x2)
{
    return ( x1 > 1.12556 ) ? -0.015 : x1 * x1 * (0.1 - 0.09936666666667 * x1);
    //return ( x1 > 1.12556 ) ? -0.015 : ( x1 <= 1.0 ) ? x1 * x1 * (0.1 - 0.09936666666667 * x1) : 9.57277081 * x1 * x1 * x1 - 30.1305982 * x1 * x1 + 31.44478397 * x1 - 10.88632325;
    //return ( x1 > 1.12556 ) ? -0.015 : ( x1 <= 1.0 ) ? x1 * x1 * (0.1 - 0.09936666666667 * x1) : 359360.7772604*pow(x1,7)-2669484.49810168*pow(x1,6)+8492609.81504811*pow(x1,5)-14999490.30127248*pow(x1,4)+15883925.71316182*pow(x1,3)-10085242.53137483*x1*x1+3555001.50280088*x1-536680.47688888;
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vx_DW2(double x1, double x2)
{
    return ( x1 > 1.12556 ) ? 0.0 : x1 * (0.2 - 0.2981 * x1);
    //return ( x1 > 1.12556 ) ? 0.0 : ( x1 <= 1.0 ) ? x1 * (0.2 - 0.2981 * x1) : 3.0 * 9.57277081 * x1 * x1 - 2.0 * 30.1305982 * x1 + 31.44478397;
    //return ( x1 > 1.12556 ) ? 0.0 : ( x1 <= 1.0 ) ? x1 * (0.2 - 0.2981 * x1) : 7.0*359360.7772604*pow(x1,6)-6.0*2669484.49810168*pow(x1,5)+5.0*8492609.81504811*pow(x1,4)-4.0*14999490.30127248*pow(x1,3)+3.0*15883925.71316182*pow(x1,2)-2.0*10085242.53137483*x1+3555001.50280088;
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vxxx_DW2(double x1, double x2)
{
    return ( x1 > 1.12556 ) ? 0.0 : -0.5962;
    //return ( x1 > 1.12556 ) ? -0.5962 : -0.5962;
    //return ( x1 > 1.12556 ) ? 0.0 : ( x1 <= 1.0 ) ? -0.5962 : 6.0 * 9.57277081;
    //return ( x1 > 1.12556 ) ? 0.0 : ( x1 <= 1.0 ) ? -0.5962 : 210.0*359360.7772604*pow(x1,4)-120.0*2669484.49810168*pow(x1,3)+60.0*8492609.81504811*pow(x1,2)-24.0*14999490.30127248*x1+6.0*15883925.71316182;
}
/* =============================================================================== */

/* WM2DPOT_DW3 */

inline double WignerMoyal2d::Wavefunction_DW3(double x1, double x2)
{
    return PI_INV / hb * exp(-2 * A[0] * pow(x1-Wave0[0],2)) * exp((-1/(2 * hb * hb * A[1])) * pow(x2 - Wave0[1],2));
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Potential_DW3(double x1, double x2)
{
    return  0.04 * (x1+0.845154255) * (x1+0.845154255) / (8.0 * exp(8.5 * (x1+0.00422577127)) + 1.0);
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vx_DW3(double x1, double x2)
{
    return (x1 > 1.0) ? 0.0 : (exp(-17*x1) * (0.0009832118398463022 + 0.0011633519376437773 * x1) + exp(-8.5*x1) * (-0.02113273521439367 - 0.059656362811940186 * x1 - 0.041000531092843885 * x1 * x1)) / pow((0.12058979733189377*exp(-8.5*x1)+1.0),2);
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vxxx_DW3(double x1, double x2)
{
    return (x1 > 1.0) ? 0.0 : (exp(-8.5*x1) * (-0.5946830136226476-2.91615415600599*x1-2.962288371457973*x1*x1)
    + exp(-25.5*x1) * (-0.0600458851556218-0.10322143014460589*x1-0.04307729895994074*x1*x1)
    + exp(-17.0*x1)* (0.9613026350451455+2.4152598845824413*x1+1.42888701741097*x1*x1)) / pow(0.12058979733189377*(-8.5*x1)+1.0,4);
}
/* =============================================================================== */

/* WM2DPOT_DW4 */

inline double WignerMoyal2d::Wavefunction_DW4(double x1, double x2)
{
    return PI_INV / hb * exp(-2 * A[0] * pow(x1-Wave0[0],2)) * exp((-1/(2 * hb * hb * A[1])) * pow(x2 - Wave0[1],2));
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Potential_DW4(double x1, double x2)
{
    return  (x1 < sqrt(1.0/1.4)) ? 0.007 * (x1 * x1 * x1 * x1) - 0.01 * (x1 * x1) : -(0.1 / 28.0);
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vx_DW4(double x1, double x2)
{
    return (x1 < sqrt(1.0/1.4)) ? 0.028 * (x1 * x1 * x1) - 0.02 * x1 : 0.0;
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vxxx_DW4(double x1, double x2)
{
    return (x1 < sqrt(1.0/1.4)) ? 0.168 * x1 : 0.0;
}

/* =============================================================================== */

/* WM2DPOT_DW5 */

inline double WignerMoyal2d::Wavefunction_DW5(double x1, double x2)
{
    return PI_INV / hb * exp(-2 * A[0] * pow(x1-Wave0[0],2)) * exp((-1/(2 * hb * hb * A[1])) * pow(x2 - Wave0[1],2));
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Potential_DW5(double x1, double x2)
{
    return V0 / pow(cosh(alpha * x1),2);
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vx_DW5(double x1, double x2)
{
    return -2.0 * V0 * alpha * sinh(alpha * x1) / pow(cosh(alpha * x1),3);
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vxxx_DW5(double x1, double x2)
{
    return -8.0 * V0 * alpha * alpha * alpha * tanh(alpha * x1) * (3.0 * pow(tanh(alpha * x1),2) - 2.0) / pow(cosh(alpha * x1),2);
    //return 24.0 * alpha * alpha * alpha * V0 * pow(sinh(alpha * x1),3) / pow(cosh(alpha * x1),5) - 16.0 * alpha * alpha * alpha * V0 * sinh(alpha * x1) / pow(cosh(alpha * x1),3);
}
/* =============================================================================== */

/* WM2DPOT_DW6 */

inline double WignerMoyal2d::Wavefunction_DW6(double x1, double x2)
{
    return PI_INV / hb * exp(-2 * A[0] * pow(x1-Wave0[0],2)) * exp((-1/(2 * hb * hb * A[1])) * pow(x2 - Wave0[1],2));
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Potential_DW6(double x1, double x2)
{
    return V0 * exp(-alpha * x1 * x1);
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vx_DW6(double x1, double x2)
{
    return -2.0 * alpha * x1 * V0 * exp(-alpha * x1 * x1);
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vxxx_DW6(double x1, double x2)
{ 
    return V0 * (12.0 * alpha * alpha * x1 * exp(-alpha*x1*x1) - pow(2.0*alpha*x1,3) * exp(-alpha*x1*x1));
}
/* =============================================================================== */

/* WM2DPOT_DW7 */

inline double WignerMoyal2d::Wavefunction_DW7(double x1, double x2)
{
    return PI_INV / hb * exp(-2 * A[0] * pow(x1-Wave0[0],2)) * exp((-1/(2 * hb * hb * A[1])) * pow(x2 - Wave0[1],2));
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Potential_DW7(double x1, double x2)
{
    return 0.0;
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vx_DW7(double x1, double x2)
{
    return 0.0;
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Vxxx_DW7(double x1, double x2)
{
    return 0.0;
}
/* =============================================================================== */

/* WM2DPOT_DW8 */

inline double WignerMoyal2d::Wavefunction_DW8(double x1, double x2)
{
    return PI_INV / hb * exp(-2 * A[0] * pow(x1-Wave0[0],2)) * exp((-1/(2 * hb * hb * A[1])) * pow(x2 - Wave0[1],2));
}
/* ------------------------------------------------------------------------------- */

inline double WignerMoyal2d::Potential_DW8(double x1, double x2)
{
    return ( x1 > 5.0 ) ? -9.92083333 : x1 * x1 * (0.1 - 0.09936666666667 * x1);
}
/* ------------------------------------------------------------------------------- */
inline double WignerMoyal2d::Vx_DW8(double x1, double x2)
{
    return ( x1 > 5.0 ) ? 0.0 : x1 * (0.2 - 0.2981 * x1);
}
/* ------------------------------------------------------------------------------- */
    
inline double WignerMoyal2d::Vxxx_DW8(double x1, double x2)
{
    return ( x1 > 5.0 ) ? 0.0 : -0.5962;
}
/* ------------------------------------------------------------------------------- */

VectorXi WignerMoyal2d::IdxToGrid(int idx)
{
    int x1 = (int)(idx / M1);
    int x2 = (int)(idx % M1);

    VectorXi grid;
    grid.resize(DIMENSIONS);
    grid << x1, x2;

    return grid;
}
/* ------------------------------------------------------------------------------- */

inline int WignerMoyal2d::GridToIdx(int x1, int x2)
{
    return (int)(x1 * M1 + x2);
}
/* ------------------------------------------------------------------------------- */
