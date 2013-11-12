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

#include "BranchEvent.h"
#include "TraitBranchEvent.h"


class comp_history
{

public:

    bool operator()(BranchEvent* m1, BranchEvent* m2) const;
    bool operator()(TraitBranchEvent* m1, TraitBranchEvent* m2) const;
};


inline bool comp_history::operator()(BranchEvent* m1, BranchEvent* m2) const
{
   return (*m1 < *m2);
}


inline bool comp_history::operator()(TraitBranchEvent* m1,
                                     TraitBranchEvent* m2) const
{
    return (*m1 < *m2);
}


#endif
