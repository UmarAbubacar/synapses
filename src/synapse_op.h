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
#ifndef SYNAPSE_OP_H
#define SYNAPSE_OP_H

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

// Synapse behaviour
struct SynapseFormation : public Behavior {
  BDM_BEHAVIOR_HEADER(SynapseFormation, Behavior, 1);

  // This function is used to detect dendritic agents in the simulation.
  // It takes an agent and a pointer to a direction vector as arguments.
  // The function first identifies the mother neuron of the given agent.
  // Then, it iterates over all neighboring agents within a certain distance.
  // If a neighboring agent belongs to a different neuron, it updates the
  // closest neighbor information. Finally, it returns the unique identifier
  // (UID) of the closest neighbor agent.
  AgentUid DendriticDetector(Agent* agent, Real3* neighbours_direction) {
    // Cast the agent to a neurite element
    auto* neurite = bdm_static_cast<NeuriteElement*>(agent);
    // Get the active simulation and its execution context
    auto* sim = Simulation::GetActive();
    auto* ctxt = sim->GetExecutionContext();

    // Find the mother neuron of the neurite
    auto* mother = dynamic_cast<NeuriteElement*>(neurite)->GetMother().Get();
    basic_neuron* mother_neuron = dynamic_cast<basic_neuron*>(mother);
    while (mother_neuron == nullptr && mother != nullptr) {
      mother = dynamic_cast<NeuriteElement*>(mother)->GetMother().Get();
      mother_neuron = dynamic_cast<basic_neuron*>(mother);
    }

    // Get the UID of the mother neuron
    int mother_cell_uid;
    if (mother_neuron != nullptr) {
      mother_cell_uid = mother_neuron->GetUid();
    }

    // Initialize closest neighbor information
    int closest_neighbour_uid = -1;
    real_t closest_distance = std::numeric_limits<real_t>::max();

    // Define a lambda function to process each neighbor
    auto print_id_distance = L2F([&](Agent* a, real_t squared_distance) {
      // If the neighbor is a neurite element
      if (auto* neighbor_den = dynamic_cast<NeuriteElement*>(a)) {
        // Get the mother neuron of the neighbor
        auto* neighbor_mother =
            dynamic_cast<NeuriteElement*>(neighbor_den)->GetMother().Get();

        basic_neuron* neighbor_mother_neuron =
            dynamic_cast<basic_neuron*>(neighbor_mother);

        // If the mother is not a neuron, keep going up the tree
        while (neighbor_mother_neuron == nullptr &&
               neighbor_mother != nullptr) {
          neighbor_mother =
              dynamic_cast<NeuriteElement*>(neighbor_mother)->GetMother().Get();
          neighbor_mother_neuron = dynamic_cast<basic_neuron*>(neighbor_mother);
        }

        // Get the UID of the neighbor's mother neuron
        int neighbor_mother_cell_uid = neighbor_mother_neuron->GetUid();
        // If the neighbor belongs to a different neuron
        if (neighbor_mother_cell_uid != mother_cell_uid) {
          // Update the closest neighbor information if the neighbor is closer
          real_t distance = std::sqrt(squared_distance);
          if (distance < closest_distance && distance < 1) {
            closest_neighbour_uid = neighbor_den->GetUid().GetIndex();
            closest_distance = distance;
          }
        }
      }
    });
    // Process all neighbors within a certain distance
    ctxt->ForEachNeighbor(print_id_distance, *neurite,
                          25);  // 25 is the squared distance threshold

    // Return the UID of the closest neighbor
    return int(closest_neighbour_uid) != -1 ? AgentUid(closest_neighbour_uid)
                                            : AgentUid(-1);
  };

  // This function is executed for each agent in the simulation.
  void Run(Agent* agent) override {
    // Get the active simulation's resource manager
    auto* rm = Simulation::GetActive()->GetResourceManager();
    // Get the current simulation time step
    int time_step =
        Simulation::GetActive()->GetScheduler()->GetSimulatedSteps();

    // Cast the agent to a neurite element
    auto* neurite = bdm_static_cast<NeuriteElement*>(agent);

    // If the neurite has not formed a synapse yet
    if (!synapsed_) {
      // Initialize a direction vector for the closest neighbor
      Real3 neighbours_direction;
      // Detect the closest dendritic neighbor
      AgentUid closest_neighbour_uid =
          DendriticDetector(agent, &neighbours_direction);

      // If a closest neighbor was found
      if (closest_neighbour_uid != AgentUid(-1)) {
        // Create a synapse between the neurite and the closest neighbor
        CreateSynapseBetweenNeurites(
            neurite,
            // Get the closest neighbor agent from the resource manager and cast
            // it to a neurite element
            dynamic_cast<NeuriteElement*>(rm->GetAgent(closest_neighbour_uid)),
            0.0,       // The synapse strength
            1,         // The synapse type
            time_step  // The current simulation time step
        );
      }
    }
  }

  // A flag indicating whether the neurite has checked for a synapse
  bool synapsed_ = false;
};
// This function adds a SynapseFormation behavior to the given axon element.
// It first checks if the axon element is a neurite element.
// If it is, it adds a new SynapseFormation behavior to it.
inline void AddSynapseBehavior(Agent* axon_element) {
  // Cast the agent to a neurite element
  auto axon = dynamic_cast<NeuriteElement*>(axon_element);

  // Check if axon is nullptr before accessing its methods
  if (axon == nullptr) {
    return;
  }

  // Add a new SynapseFormation behavior to the neurite element
  axon->AddBehavior(new SynapseFormation());
}

// This struct represents a standalone operation in the simulation.
// It overrides the operator() method to add a SynapseFormation behavior
// to each agent in the simulation at a certain time step.
struct synapse_op : public StandaloneOperationImpl {
  // Macro to set the operation header
  BDM_OP_HEADER(synapse_op);

  // Override the operator() method
  void operator()() override {
    // Get the active simulation's resource manager
    auto* rm = Simulation::GetActive()->GetResourceManager();

    // call the function when simulation time is -3 the maximum value of the
    // simulation time. This is to ensure the operation is registered and runs
    // for the last time step.
    if (Simulation::GetActive()->GetScheduler()->GetSimulatedSteps() >
        500 - 3) {
      // For each agent in the simulation
      rm->ForEachAgent([&](Agent* agent) {
        // If the agent is not nullptr
        if (agent != nullptr) {
          // Add a SynapseFormation behavior to the agent
          AddSynapseBehavior(agent);
        }
      });
    }
  }
};

}  // namespace bdm

#endif  // SYNAPSE_OP_H