// ==============================================================================
//
//  CaldeiraLeggett4d.h
//  QTR  i      
//
//  Created by Albert Lu on 1/30/20.
//  alu@tacc.utexas.edu
//
//  Last modified on 1/30/20
//
//  Note:
//
// ==============================================================================

#ifndef QTR_CALDEIRALEGGETT4D_H
#define QTR_CALDEIRALEGGETT4D_H

#include <complex>

#include "Containers.h"
#include "Eigen.h"
#include "Pointers.h"

namespace QTR_NS {
    
    class CaldeiraLeggett4d {
        
    public:
        CaldeiraLeggett4d(class QTR *q);
        ~CaldeiraLeggett4d();
  
        void                          Evolve();
        VectorXi                      IdxToGrid(unsigned long int idx);
        inline unsigned long int      GridToIdx(int x1, int x2, int x3, int x4);
        
        inline double                 Wavefunction_DW1(double x1, double x2, double x3, double x4);
        inline double                 Potential_DW1(double x1, double x2);
        inline double                 Vx1_DW1(double x1, double x2);
        inline double                 Vx2_DW1(double x1, double x2);
        inline double                 Vxxx1_DW1(double x1, double x2);
        inline double                 Vxxx2_DW1(double x1, double x2);

        inline double                 Wavefunction_DW2(double x1, double x2, double x3, double x4);
        inline double                 Potential_DW2(double x1, double x2);
        inline double                 Vx1_DW2(double x1, double x2);
        inline double                 Vx2_DW2(double x1, double x2);
        inline double                 Vxxx1_DW2(double x1, double x2);
        inline double                 Vxxx2_DW2(double x1, double x2);

        inline double                 Wavefunction_DW3(double x1, double x2, double x3, double x4);
        inline double                 Potential_DW3(double x1, double x2);
        inline double                 Vx1_DW3(double x1, double x2);
        inline double                 Vx2_DW3(double x1, double x2);
        inline double                 Vxxx1_DW3(double x1, double x2);
        inline double                 Vxxx2_DW3(double x1, double x2);

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
        unsigned long int GRIDS_TOT;
        bool            QUIET;
        bool            TIMING;
        double          TIME;   
        double          PI_INV;  // 1/pi
        double          HBSQ_INV; // (1/hb)^2
        double          PIHBSQ_INV; // (1/hb/pi)^2

        // Grid size
        double          kk;    // time resolution
        VectorXd        H;     // grid size
        VectorXd        C;     // coarse grid size
        VectorXd        Hi;    // inverse grid size
        VectorXd        Hisq;  // inverse grid size square     
        VectorXd        S;  
        VectorXd        DistToEdge;
        int             CFactor; // Coarse dx / fine dx
        int             skin; 

        // Domain size
        VectorXd        Box;
        VectorXd        BoxBoundary;
        VectorXd        BoxDynamic;
        VectorXd        BoxDynamicNew;
        VectorXi        BoxShapeFine;     // Fine-grid box shape
        VectorXi        BoxShapeFineNew;  // Fine-grid box shape
        VectorXi        BoxShapeCoarse;   // Coarse-grid box shape
        unsigned long int M1, M2, M3;
        unsigned long int W1, W2, W1_0;
        unsigned long int O1, O2;

        // Potential parameters
        int             idx_x0;  
        double          trans_x0;     
        double          hb;
        double          m;
        double          kb;
        double          temp;
        double          gamma;
        double          k0;
        double          sigma;
        double          lambda;

        // Wavefunction
        VectorXd        Wave0;
        VectorXd        A;
        VectorXd        P;

        // Truncate parameters
        bool            isEmpty;
        bool            isExtrapolate;         
        double          TolH;
        double          TolL;
        double          TolHd;
        double          TolLd;
        double          ExReduce;
        int             ExLimit;

        // Domains
        MeshIndexLU       TB;    // Truncation boundary
        MeshIndexLU       TBL;
        MeshIndexLU       TBL_P;          
        MeshIndexLU       ExFF;
        MeshIndexLU       ExFF2;

        // Output
        bool            isTrans;
        bool            isCorr;
        bool            isPrintEdge;
        bool            isPrintDensity;
        bool            isPrintWavefunc;
    };
}

#endif /* QTR_CALDEIRALEGGETT4D_H */
