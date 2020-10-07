// ==============================================================================
//
//  Job_Diosi4d.cpp
//  QTR
//
//  Created by Albert Lu on 12/16/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 2/24/19
//
//  Note:
//
// ==============================================================================

#include <string>

#include "Job_Diosi4d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Diosi4d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobDiosi4d::JobDiosi4d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    ds = new Diosi4d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobDiosi4d::~JobDiosi4d()
{
    delete ds;
}

/* ------------------------------------------------------------------------------- */

void JobDiosi4d::run(class QTR *qtr)
{     
    ds->Evolve();
    log->log("[Job_Diosi4d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

