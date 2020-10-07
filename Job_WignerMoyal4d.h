// ==============================================================================
//
//  Job_WignerMoyal4d.h
//  QTR
//
//  Created by Albert Lu on 12/16/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 12/16/19
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_WIGNERMOYAL4D_H
#define QTR_JOB_WIGNERMOYAL4D_H

#include "Job.h"
#include "WignerMoyal4d.h"

namespace QTR_NS  {

    class JobWignerMoyal4d: public Job    {
        
    public:
        JobWignerMoyal4d(class QTR *);
        ~JobWignerMoyal4d(void);
        
        void run(class QTR *);
        
    private:
        
        Log        *log;
        Parameters *parameters;
        QTR        *qtr;
        WignerMoyal4d *wm;
    };
}
#endif /* QTR_JOB_WIGNERMOYAL4D_H */

