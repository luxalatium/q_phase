// ==============================================================================
//
//  Job_Diosi2d.h
//  QTR
//
//  Created by Albert Lu on 10/8/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 10/8/19
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_DIOSI2D_H
#define QTR_JOB_DIOSI2D_H

#include "Job.h"
#include "Diosi2d.h"

namespace QTR_NS  {

    class JobDiosi2d: public Job    {
        
    public:
        JobDiosi2d(class QTR *);
        ~JobDiosi2d(void);
        
        void run(class QTR *);
        
    private:
        
        Log        *log;
        Parameters *parameters;
        QTR        *qtr;
        Diosi2d    *ds;
    };
}
#endif /* QTR_JOB_DIOSI2D_H */
