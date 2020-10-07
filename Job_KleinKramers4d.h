// ==============================================================================
//
//  Job_KleinKramers4d.h
//  QTR
//
//  Created by Albert Lu on 10/28/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 10/28/19
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_KLEINKRAMERS4D_H
#define QTR_JOB_KLEINKRAMERS4D_H

#include "Job.h"
#include "KleinKramers4d.h"

namespace QTR_NS  {

    class JobKleinKramers4d: public Job    {
        
    public:
        JobKleinKramers4d(class QTR *);
        ~JobKleinKramers4d(void);
        
        void run(class QTR *);
        
    private:
        
        Log        *log;
        Parameters *parameters;
        QTR        *qtr;
        KleinKramers4d *kk;
    };
}
#endif /* QTR_JOB_KLEINKRAMERS4D_H */

