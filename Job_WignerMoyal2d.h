// ==============================================================================
//
//  Job_WignerMoyal2d.h
//  QTR
//
//  Created by Albert Lu on 9/27/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 9/27/19
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_WIGNERMOYAL2D_H
#define QTR_JOB_WIGNERMOYAL2D_H

#include "Job.h"
#include "WignerMoyal2d.h"

namespace QTR_NS  {

    class JobWignerMoyal2d: public Job    {
        
    public:
        JobWignerMoyal2d(class QTR *);
        ~JobWignerMoyal2d(void);
        
        void run(class QTR *);
        
    private:
        
        Log        *log;
        Parameters *parameters;
        QTR        *qtr;
        WignerMoyal2d  *wm;
    };
}
#endif /* QTR_JOB_WIGNERMOYAL2D_H */
