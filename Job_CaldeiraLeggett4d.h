// ==============================================================================
//
//  Job_CaldeiraLeggett4d.h
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

#ifndef QTR_JOB_CALDEIRALEGGETT4D_H
#define QTR_JOB_CALDEIRALEGGETT4D_H

#include "Job.h"
#include "CaldeiraLeggett4d.h"

namespace QTR_NS  {

    class JobCaldeiraLeggett4d: public Job    {
        
    public:
        JobCaldeiraLeggett4d(class QTR *);
        ~JobCaldeiraLeggett4d(void);
        
        void run(class QTR *);
        
    private:
        
        Log        *log;
        Parameters *parameters;
        QTR        *qtr;
        CaldeiraLeggett4d *cl;
    };
}
#endif /* QTR_JOB_CALDEIRALEGGETT4D_H */

