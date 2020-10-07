// ==============================================================================
//
//  Job_WignerMoyal4d.cpp
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

#include "Job_WignerMoyal4d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "WignerMoyal4d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobWignerMoyal4d::JobWignerMoyal4d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    wm = new WignerMoyal4d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobWignerMoyal4d::~JobWignerMoyal4d()
{
    delete wm;
}

/* ------------------------------------------------------------------------------- */

void JobWignerMoyal4d::run(class QTR *qtr)
{     
    wm->Evolve();
    log->log("[Job_WignerMoyal4d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

