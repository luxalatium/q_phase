// ==============================================================================
//
//  Job_Diosi4d.h
//  QTR
//
//  Created by Albert Lu on 12/16/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 2/24/20
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_DIOSI4D_H
#define QTR_JOB_DIOSI4D_H

#include "Job.h"
#include "Diosi4d.h"

namespace QTR_NS  {

    class JobDiosi4d: public Job    {
        
    public:
        JobDiosi4d(class QTR *);
        ~JobDiosi4d(void);
        
        void run(class QTR *);
        
    private:
        
        Log        *log;
        Parameters *parameters;
        QTR        *qtr;
        Diosi4d    *ds;
    };
}
#endif /* QTR_JOB_DIOSI4D_H */

