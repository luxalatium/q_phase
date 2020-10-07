// ==============================================================================
//
//  Job_KleinKramers2d.cpp
//  QTR
//
//  Created by Albert Lu on 9/2/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 9/2/19
//
//  Note:
//
// ==============================================================================

#include <string>

#include "Job_KleinKramers2d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "KleinKramers2d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobKleinKramers2d::JobKleinKramers2d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    kk = new KleinKramers2d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobKleinKramers2d::~JobKleinKramers2d()
{
    delete kk;
}

/* ------------------------------------------------------------------------------- */

void JobKleinKramers2d::run(class QTR *qtr)
{     
    kk->Evolve();
    log->log("[Job_KleinKramers2d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

