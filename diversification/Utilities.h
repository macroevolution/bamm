/*
 *  Utilities.h
 *  BAMMt
 *
 *  Created by Dan Rabosky on 6/20/12.
 *  Copyright 2012 DLR. All rights reserved.
 *
 */

// comparison operator

#ifndef UTILS_H
#define UTILS_H

#include <math.h>

#include "event.h"
#include "TraitEvent.h"




class comp_history {
    
public:
	bool		operator()(BranchEvent* m1, BranchEvent* m2) const { return (*m1 < *m2); }
	bool		operator()(TraitBranchEvent* m1, TraitBranchEvent* m2) const { return (*m1 < *m2); }
	
};

/*

double safeExponentiation(double x) {
	
	if (x > 0.0)
		return 1.0;
	else if (x < -100.0)
		return 0.0;
	else
		return exp(x);
	
}
*/


#endif