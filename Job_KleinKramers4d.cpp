// ==============================================================================
//
//  Job_KleinKramers4d.cpp
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

#include <string>

#include "Job_KleinKramers4d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "KleinKramers4d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobKleinKramers4d::JobKleinKramers4d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    kk = new KleinKramers4d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobKleinKramers4d::~JobKleinKramers4d()
{
    delete kk;
}

/* ------------------------------------------------------------------------------- */

void JobKleinKramers4d::run(class QTR *qtr)
{     
    kk->Evolve();
    log->log("[Job_KleinKramers4d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

