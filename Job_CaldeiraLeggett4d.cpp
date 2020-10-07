// ==============================================================================
//
//  Job_CaldeiraLeggett4d.cpp
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

#include <string>

#include "Job_CaldeiraLeggett4d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "CaldeiraLeggett4d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobCaldeiraLeggett4d::JobCaldeiraLeggett4d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    cl = new CaldeiraLeggett4d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobCaldeiraLeggett4d::~JobCaldeiraLeggett4d()
{
    delete cl;
}

/* ------------------------------------------------------------------------------- */

void JobCaldeiraLeggett4d::run(class QTR *qtr)
{     
    cl->Evolve();
    log->log("[Job_CaldeiraLeggett4d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

