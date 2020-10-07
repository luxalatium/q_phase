// ==============================================================================
//
//  Job_WignerMoyal2d.cpp
//  QTR
//
//  Created by Albert Lu on 9/27/19.
//  alu@tacc.utexas.edu
//
//  Last modified on 9/27/19
//
//  Note:
//
// ==============================================================================

#include <string>

#include "Job_WignerMoyal2d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "WignerMoyal2d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobWignerMoyal2d::JobWignerMoyal2d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    wm = new WignerMoyal2d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobWignerMoyal2d::~JobWignerMoyal2d()
{
    delete wm;
}

/* ------------------------------------------------------------------------------- */

void JobWignerMoyal2d::run(class QTR *qtr)
{     
    wm->Evolve();
    log->log("[Job_WignerMoyal2d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

