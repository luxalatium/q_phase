// ==============================================================================
//
//  Job.cpp
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

# include "Job.h"
# include "Parameters.h"

//# include "Job_Test.h"
//# include "Job_Imt2d.h"
//# include "Job_Imt3d.h"
//# include "Job_Imt4d.h"
//# include "Job_Scatter2d.h"
//# include "Job_Scatter3d.h"
//# include "Job_Scatter4d.h"
# include "Job_KleinKramers2d.h"
# include "Job_KleinKramers4d.h"
# include "Job_WignerMoyal2d.h"
# include "Job_WignerMoyal4d.h"
# include "Job_CaldeiraLeggett2d.h"
# include "Job_CaldeiraLeggett4d.h"
# include "Job_Diosi2d.h"
# include "Job_Diosi4d.h"

using namespace QTR_NS;

const char Job::TEST[] = "test";
//const char Job::IMT2D[] = "imt2d";
//const char Job::IMT3D[] = "imt3d";
//const char Job::IMT4D[] = "imt4d";
//const char Job::SCATTER2D[] = "scatter2d";
//const char Job::SCATTER3D[] = "scatter3d";
//const char Job::SCATTER4D[] = "scatter4d";
const char Job::KLEINKRAMERS2D[] = "kleinkramers2d";
const char Job::KLEINKRAMERS4D[] = "kleinkramers4d";
const char Job::WIGNERMOYAL2D[] = "wignermoyal2d";
const char Job::WIGNERMOYAL4D[] = "wignermoyal4d";
const char Job::CALDEIRALEGGETT2D[] = "caldeiraleggett2d";
const char Job::CALDEIRALEGGETT4D[] = "caldeiraleggett4d";
const char Job::DIOSI2D[] = "diosi2d";
const char Job::DIOSI4D[] = "diosi4d";

/* ------------------------------------------------------------------------------- */

Job *Job::getJob(class QTR *qtr) {
    
    Job *job = NULL;
    
    if (qtr->parameters->job == KLEINKRAMERS4D)
    {
        job = new JobKleinKramers4d(qtr);
    }
    else if (qtr->parameters->job == KLEINKRAMERS2D)
    {
        job = new JobKleinKramers2d(qtr);
    }
    else if (qtr->parameters->job == WIGNERMOYAL2D)
    {
        job = new JobWignerMoyal2d(qtr);
    }
    else if (qtr->parameters->job == WIGNERMOYAL4D)
    {
        job = new JobWignerMoyal4d(qtr);
    }
    else if (qtr->parameters->job == CALDEIRALEGGETT2D)
    {
        job = new JobCaldeiraLeggett2d(qtr);
    }
    else if (qtr->parameters->job == CALDEIRALEGGETT4D)
    {
        job = new JobCaldeiraLeggett4d(qtr);
    }
    else if (qtr->parameters->job == DIOSI2D)
    {
        job = new JobDiosi2d(qtr);
    }
    else if (qtr->parameters->job == DIOSI4D)
    {
        job = new JobDiosi4d(qtr);
    }
    return job;
}
/* ----------------------------------------------------------------- */
