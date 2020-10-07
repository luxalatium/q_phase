// ==============================================================================
//
//  Job_Diosi2d.cpp
//  QTR
//
//  Created by Albert Lu on 10/8/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 10/8/19
//
//  Note:
//
// ==============================================================================

#include <string>

#include "Job_Diosi2d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Diosi2d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobDiosi2d::JobDiosi2d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    ds = new Diosi2d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobDiosi2d::~JobDiosi2d()
{
    delete ds;
}

/* ------------------------------------------------------------------------------- */

void JobDiosi2d::run(class QTR *qtr)
{     
    ds->Evolve();
    log->log("[Job_Diosi2d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

