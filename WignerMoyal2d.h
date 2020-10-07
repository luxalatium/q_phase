// ==============================================================================
//
//  WignerMoyal2d.h
//  QTR
//
//  Created by Albert Lu on 9/27/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 4/1/20
//
//  Note:
//
// ==============================================================================

#ifndef QTR_WIGNERMOYAL2D_H
#define QTR_WIGNERMOYAL2D_H

#include <complex>

#include "Containers.h"
#include "Eigen.h"
#include "Pointers.h"

namespace QTR_NS {
    
    class WignerMoyal2d {
        
    public:
        WignerMoyal2d(class QTR *q);
        ~WignerMoyal2d();
  
        void                          Evolve();
        VectorXi                      IdxToGrid(int idx);
        inline int                    GridToIdx(int x1, int x2);

        inline double                 Wavefunction_DW1(double x1, double x2);
        inline double                 Potential_DW1(double x1, double x2);
        inline double                 Vx_DW1(double x1, double x2);
        inline double                 Vxxx_DW1(double x1, double x2);

        inline double                 Wavefunction_DW2(double x1, double x2);
        inline double                 Potential_DW2(double x1, double x2);
        inline double                 Vx_DW2(double x1, double x2);
        inline double                 Vxxx_DW2(double x1, double x2);

        inline double                 Wavefunction_DW3(double x1, double x2);
        inline double                 Potential_DW3(double x1, double x2);
        inline double                 Vx_DW3(double x1, double x2);
        inline double                 Vxxx_DW3(double x1, double x2);

        inline double                 Wavefunction_DW4(double x1, double x2);
        inline double                 Potential_DW4(double x1, double x2);
        inline double                 Vx_DW4(double x1, double x2);
        inline double                 Vxxx_DW4(double x1, double x2);

        inline double                 Wavefunction_DW5(double x1, double x2);
        inline double                 Potential_DW5(double x1, double x2);
        inline double                 Vx_DW5(double x1, double x2);
        inline double                 Vxxx_DW5(double x1, double x2);

        inline double                 Wavefunction_DW6(double x1, double x2);
        inline double                 Potential_DW6(double x1, double x2);
        inline double                 Vx_DW6(double x1, double x2);
        inline double                 Vxxx_DW6(double x1, double x2);

        inline double                 Wavefunction_DW7(double x1, double x2);
        inline double                 Potential_DW7(double x1, double x2);
        inline double                 Vx_DW7(double x1, double x2);
        inline double                 Vxxx_DW7(double x1, double x2);

        inline double                 Wavefunction_DW8(double x1, double x2);
        inline double                 Potential_DW8(double x1, double x2);
        inline double                 Vx_DW8(double x1, double x2);
        inline double                 Vxxx_DW8(double x1, double x2);

    private:

        void            init();
        QTR             *qtr;
        Error           *err;
        Log             *log;
        Parameters      *parameters;

        // General parameters
        std::complex<double>  I;      // sqrt(-1)
        std::complex<double>  xZERO;  // complex zero
        int             DIMENSIONS;
        int             PERIOD;
        int             SORT_PERIOD;
        int             PRINT_PERIOD;
        int             PRINT_WAVEFUNC_PERIOD;
        int             GRIDS_TOT;
        bool            QUIET;
        bool            TIMING;
        double          TIME;   
        double          PI_INV;  // 1/pi

        // Grid size
        double          kk;    // time resolution
        VectorXd        H;     // grid size
        VectorXd        Hi;    // inverse grid size
        VectorXd        Hisq;  // inverse grid size square     
        VectorXd        S;  

        // Domain size
        VectorXd        Box;
        VectorXi        BoxShape;
        int             M1;
        int             W1;
        int             O1;

        // Potential parameters
        int             idx_x0;
        int             skin; 
        double          trans_x0;        
        double          hb;
        double          m;
        double          kb;
        double          temp;
        double          gamma;
        double          alpha;
        double          V0;
        double          lambda;

        // Wavefunction
        VectorXd        Wave0;
        VectorXd        A;
        VectorXd        P;

        // Truncate parameters
        bool            isEmpty;
        bool            isFullGrid; 
        bool            isExtrapolate;         
        double          TolH;
        double          TolL;
        double          TolHd;
        double          TolLd;
        double          ExReduce;
        int             ExLimit;

        // Domains
        MeshIndex       TA;
        MeshIndex       TB;    // Truncation boundary
        MeshIndex       TBL;
        MeshIndex       TBL_P;          
        MeshIndex       DBi;   // Grid boundary
        MeshIndex       DBi2;  // Extrapolation-restricted area
        MeshIndex       ExFF;
        MeshIndex       ExFF2;

        // Output
        bool            isTrans;
        bool            isCorr;
        bool            isPrintEdge;
        bool            isPrintDensity;
        bool            isPrintWavefunc;
    };
}

#endif /* QTR_WIGNERMOYAL2D_H */
