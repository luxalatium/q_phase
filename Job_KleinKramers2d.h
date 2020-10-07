// ==============================================================================
//
//  Job_KleinKramers2d.h
//  QTR
//
//  Created by Albert Lu on 9/2/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 9/27/19
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_KLEINKRAMERS2D_H
#define QTR_JOB_KLEINKRAMERS2D_H

#include "Job.h"
#include "KleinKramers2d.h"

namespace QTR_NS  {

    class JobKleinKramers2d: public Job    {
        
    public:
        JobKleinKramers2d(class QTR *);
        ~JobKleinKramers2d(void);
        
        void run(class QTR *);
        
    private:
        
        Log        *log;
        Parameters *parameters;
        QTR        *qtr;
        KleinKramers2d *kk;
    };
}
#endif /* QTR_JOB_KLEINKRAMERS2D_H */
