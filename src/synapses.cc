// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#include "MyNeuron.h"
#include "synapse_op.h"
#include "synapses.h"

namespace bdm {
BDM_REGISTER_OP(Synapsification_op, "Synapsification_op", kCpu);
}  // namespace bdm

int main(int argc, const char** argv) { return bdm::Simulate(argc, argv); }