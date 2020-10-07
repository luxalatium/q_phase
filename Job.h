// ==============================================================================
//
//  Job.h
//  QTR
//
//  Created by Albert Lu on 8/31/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 2/21/20
//
//  Note:
//
// ==============================================================================

#ifndef QTR_JOB_H
#define QTR_JOB_H

#include "Qtr.h"

namespace QTR_NS {
    
    class Job {
        
    public:
        
        virtual ~Job() {};
        
        virtual void run(class QTR *) = 0;
        
        static Job *getJob(class QTR *);

        static const char TEST[];
        static const char WIGNERMOYAL2D[];
        static const char WIGNERMOYAL4D[];
        static const char KLEINKRAMERS2D[];
        static const char KLEINKRAMERS4D[];
        static const char CALDEIRALEGGETT2D[];
        static const char CALDEIRALEGGETT4D[];
        static const char DIOSI2D[];
        static const char DIOSI4D[];
    };
}
#endif /* QTR_JOB_H */
