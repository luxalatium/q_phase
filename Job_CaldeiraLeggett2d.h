// ==============================================================================
//
//  Job_CaldeiraLeggett2d.h
//  QTR
//
//  Created by Albert Lu on 10/5/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 10/5/19
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_CALDEIRALEGGETT2D_H
#define QTR_JOB_CALDEIRALEGGETT2D_H

#include "Job.h"
#include "CaldeiraLeggett2d.h"

namespace QTR_NS  {

    class JobCaldeiraLeggett2d: public Job    {
        
    public:
        JobCaldeiraLeggett2d(class QTR *);
        ~JobCaldeiraLeggett2d(void);
        
        void run(class QTR *);
        
    private:
        
        Log        *log;
        Parameters *parameters;
        QTR        *qtr;
        CaldeiraLeggett2d  *cl;
    };
}
#endif /* QTR_JOB_CALDEIRALEGGETT2D_H */
