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
#ifndef MY_NEURON_H_
#define MY_NEURON_H_

#include <vector>
#include "biodynamo.h"
#include "core/agent/cell_division_event.h"
#include "core/behavior/behavior.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

enum States { alive, dead };

class basic_neuron;
class Synapses;

// our object extends the Cell object create the header with our new data member
class basic_neuron : public neuroscience::NeuronSoma {
  BDM_AGENT_HEADER(basic_neuron, neuroscience::NeuronSoma, 1);

 public:
  basic_neuron() {}
  explicit basic_neuron(const Real3& position) : Base(position) {}
  virtual ~basic_neuron() {}

  void AddSynapse(basic_neuron* target, real_t distance, int strength = 1,
                  int time = 0);
  const std::vector<Synapses>& GetSynapses() const { return synapses_; }

  int GetState() const { return state_; }
  void SetState(int state) { state_ = state; }

  int state_ = States::alive;
  std::vector<Synapses> synapses_;
};

// This class represents a Synapse in the simulation
class Synapses {
 public:
  // Default constructor
  Synapses()
      : source_(nullptr), target_(nullptr), distance_(-1.0), strength_(0) {}

  // Parameterized constructor
  Synapses(basic_neuron* source, basic_neuron* target, double distance = 0.0,
           int strength = 1, int time = 0)
      : source_(source),
        target_(target),
        distance_(distance),
        strength_(strength),
        time_(time) {}

  basic_neuron* GetSource() const { return source_; }
  basic_neuron* GetTarget() const { return target_; }
  double GetDistance() const { return distance_; }
  int GetStrength() const { return strength_; }
  int GetTime() const { return time_; }

  void IncreaseStrength(int amount = 1) { strength_ += amount; }

 private:
  basic_neuron* source_;
  basic_neuron* target_;
  double distance_;
  int strength_;
  int time_;
};

// This function adds a new Synapse to the current neuron.
// It takes a target neuron, the distance to the target, the strength of the
// synapse, and the time of synapse formation as arguments. It creates a new
// Synapse with these parameters and adds it to the neuron's list of synapses.
inline void basic_neuron::AddSynapse(basic_neuron* target, real_t distance,
                                     int strength, int time) {
  Synapses synapse(this, target, distance, strength, time);
  synapses_.push_back(synapse);
}

// This function finds the parent neuron of a given neurite element.
// It takes a neurite element as an argument and iterates up the tree of neurite
// elements until it finds a neuron. It returns the found neuron, or nullptr if
// no neuron was found.
inline basic_neuron* FindParentNeuron(NeuriteElement* neurite) {
  basic_neuron* mother_neuron = nullptr;

  while (neurite && !mother_neuron) {
    auto* mother = neurite->GetMother().Get();
    mother_neuron = dynamic_cast<basic_neuron*>(mother);
    neurite = dynamic_cast<NeuriteElement*>(mother);
  }

  return mother_neuron;
}

// This function checks if there is a synapse between two given neurons.
// It iterates over the synapses of each neuron and checks if the target of any
// synapse is the other neuron. If such a synapse is found, it returns true.
// Otherwise, it returns false.
inline bool hasSynapse(const basic_neuron* neuronA,
                       const basic_neuron* neuronB) {
  for (const auto& synapse : neuronA->GetSynapses()) {
    if (synapse.GetTarget() == neuronB) {
      return true;
    }
  }
  for (const auto& synapse : neuronB->GetSynapses()) {
    if (synapse.GetTarget() == neuronA) {
      return true;
    }
  }
  return false;
}

// This function creates a synapse between two given neurite elements.
// It first finds the parent neurons of the neurite elements.
// If both parent neurons are valid and there is no existing synapse between
// them, it adds a new synapse from the first neuron to the second. The synapse
// is characterized by its distance, strength, and the time when it was formed.
inline void CreateSynapseBetweenNeurites(NeuriteElement* neurite1,
                                         NeuriteElement* neurite2,
                                         real_t distance = 0.0,
                                         int strength = 1, int time = 0) {
  // Trace back to parent neurons
  basic_neuron* neuronA = FindParentNeuron(neurite1);
  basic_neuron* neuronB = FindParentNeuron(neurite2);

  // Check if both neurons are valid
  if (neuronA && neuronB) {
    // avoid duplicate synapses
    if (!hasSynapse(neuronA, neuronB)) {
      // Create a Synapses object
      neuronA->AddSynapse(neuronB, distance, strength, time);
    }
  } else {
    std::cerr << "Failed to find parent neurons for neurites!" << std::endl;
  }
}

// This struct is a custom hash function for pairs of values.
// It is used to create a unique hash for each pair by XORing the hashes of the
// first and second values of the pair. The PairHash struct is used in the
// export_connection_list function to create a unique hash for
// each pair of neuron UIDs. This is necessary because the adjacency matrix is
// stored in an std::unordered_map where each key is a pair of neuron UIDs and
// each value is the count of synapses between the pair of neurons.
struct PairHash {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2>& pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

// This function exports the adjacency matrix of all neurons in the simulation
// to a CSV file. The adjacency matrix represents the connections between
// neurons. Each row in the CSV file represents a connection from one neuron to
// another. The columns of the CSV file are the source neuron's UID, the target
// neuron's UID, the cell type of the source neuron, and the count of synapses
// between the source and target neurons.
inline void export_connection_list() {
  // Get the active simulation and its resource manager
  auto* sim = Simulation::GetActive();
  auto* rm = sim->GetResourceManager();

  // This will store the adjacency matrix.
  std::unordered_map<std::pair<int, int>, int, PairHash> adjacency;

  // This will store the cell types of each neuron.
  std::unordered_map<int, int> cell_types;

  // Populate the set of unique neuron UIDs and their cell types
  rm->ForEachAgent([&](Agent* agent) {
    auto* neuron = dynamic_cast<basic_neuron*>(agent);
    if (neuron && neuron->GetState() != States::dead) {
      int uid = neuron->GetUid();
      cell_types[uid] = neuron->GetState();

      for (const auto& synapse : neuron->GetSynapses()) {
        int target_uid = synapse.GetTarget()->GetUid();
        adjacency[std::make_pair(uid, target_uid)]++;
      }
    }
  });

  // Export to CSV
  std::ofstream file("connection_list.csv");
  if (!file.is_open()) {
    std::cout << "Failed to open file" << std::endl;
    return;
  }

  // Write CSV header
  file << "Source_UID,Target_UID,Cell_Type,Synapse_Count" << std::endl;

  // Write neurons with connections
  for (const auto& entry : adjacency) {
    file << entry.first.first << "," << entry.first.second << ","
         << cell_types[entry.first.first] << "," << entry.second << std::endl;
  }

  // Write neurons without connections
  for (const auto& entry : cell_types) {
    if (adjacency.find(std::make_pair(entry.first, entry.first)) ==
        adjacency.end()) {
      file << entry.first << "," << entry.first << "," << entry.second << ",0"
           << std::endl;
    }
  }

  file.close();
}

}  // namespace bdm

#endif  // MY_NEURON_H_
