// ==============================================================================
//
//  Job_CaldeiraLeggett2d.cpp
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

#include <string>

#include "Job_CaldeiraLeggett2d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "CaldeiraLeggett2d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobCaldeiraLeggett2d::JobCaldeiraLeggett2d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    cl = new CaldeiraLeggett2d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobCaldeiraLeggett2d::~JobCaldeiraLeggett2d()
{
    delete cl;
}

/* ------------------------------------------------------------------------------- */

void JobCaldeiraLeggett2d::run(class QTR *qtr)
{     
    cl->Evolve();
    log->log("[Job_CaldeiraLeggett2d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

