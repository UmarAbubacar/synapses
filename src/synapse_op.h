#ifndef SYNAPSE_OP_H
#define SYNAPSE_OP_H

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

// Synapse behaviour
struct SynapseFormation : public Behavior {
  BDM_BEHAVIOR_HEADER(SynapseFormation, Behavior, 1);

  AgentUid DendriticDetector(Agent* agent, Real3* neighbours_direction) {
    auto* dendrite = bdm_static_cast<NeuriteElement*>(agent);
    auto* sim = Simulation::GetActive();
    auto* ctxt = sim->GetExecutionContext();

    auto* mother = dynamic_cast<NeuriteElement*>(dendrite)->GetMother().Get();
    MyNeuron* mother_neuron = dynamic_cast<MyNeuron*>(mother);
    while (mother_neuron == nullptr && mother != nullptr) {
      mother = dynamic_cast<NeuriteElement*>(mother)->GetMother().Get();
      mother_neuron = dynamic_cast<MyNeuron*>(mother);
    }

    int mother_cell_uid;
    if (mother_neuron != nullptr) {
      mother_cell_uid = mother_neuron->GetUid();
    }

    int closest_neighbour_uid = -1;
    Real3 closest_neighbour_direction;

    real_t closest_distance = std::numeric_limits<real_t>::max();

    auto print_id_distance = L2F([&](Agent* a, real_t squared_distance) {
      // if it is a neurite element
      if (auto* neighbor_den = dynamic_cast<NeuriteElement*>(a)) {
        // get the mother agent
        auto* neighbor_mother =
            dynamic_cast<NeuriteElement*>(neighbor_den)->GetMother().Get();

        MyNeuron* neighbor_mother_neuron =
            dynamic_cast<MyNeuron*>(neighbor_mother);

        // if it is not a neuron, keep going up the tree
        while (neighbor_mother_neuron == nullptr &&
               neighbor_mother != nullptr) {
          neighbor_mother =
              dynamic_cast<NeuriteElement*>(neighbor_mother)->GetMother().Get();
          neighbor_mother_neuron = dynamic_cast<MyNeuron*>(neighbor_mother);
        }

        int neighbor_mother_cell_uid = neighbor_mother_neuron->GetUid();
        if (neighbor_mother_cell_uid != mother_cell_uid) {
          closest_neighbour_direction +=
              neighbor_den->GetPosition() - dendrite->GetPosition();

          // NOTE: squared distance is the distance between the two dendrites
          real_t distance = std::sqrt(squared_distance);
          if (distance < closest_distance && distance < 1) {
            closest_neighbour_uid = neighbor_den->GetUid().GetIndex();
            closest_distance = distance;
          }
        }
      }
    });
    ctxt->ForEachNeighbor(print_id_distance, *dendrite, 25);  // 25 is the
                                                              // squared
                                                              // distance
                                                              // threshold
    *neighbours_direction = closest_neighbour_direction;
    return int(closest_neighbour_uid) != -1 ? AgentUid(closest_neighbour_uid)
                                            : AgentUid(-1);
  };

  void Run(Agent* agent) override {
    // auto* random = Simulation::GetActive()->GetRandom();
    auto* rm = Simulation::GetActive()->GetResourceManager();
    // auto* sparam = Simulation::GetActive()->GetParam()->Get<SimParam>();
    int time_step =
        Simulation::GetActive()->GetScheduler()->GetSimulatedSteps();

    auto* dendrite = bdm_static_cast<NeuriteElement*>(agent);

    if (!synapsed_) {
      Real3 neighbours_direction;
      AgentUid closest_neighbour_uid =
          DendriticDetector(agent, &neighbours_direction);

      if (closest_neighbour_uid != AgentUid(-1)) {
        CreateSynapseBetweenNeurites(
            dendrite,
            dynamic_cast<NeuriteElement*>(rm->GetAgent(closest_neighbour_uid)),
            0.0, 1, time_step);
      }
    }
  }

  bool synapsed_ = false;
};

inline void Synapsification(Agent* axon_element) {
  auto axon = dynamic_cast<NeuriteElement*>(axon_element);

  // Check if axon is nullptr before accessing its methods
  if (axon == nullptr) {
    return;
  }

  axon->AddBehavior(new SynapseFormation());
}

struct Synapsification_op : public StandaloneOperationImpl {
  BDM_OP_HEADER(Synapsification_op);

  void operator()() override {
    auto* rm = Simulation::GetActive()->GetResourceManager();

    // call the function when simulation time is -3 the maximum value of the
    // simulation time
    if (Simulation::GetActive()->GetScheduler()->GetSimulatedSteps() >
        500 - 3) {
      rm->ForEachAgent([&](Agent* agent) {
        if (agent != nullptr) {
          Synapsification(agent);
        }
      });
    }
  }
};

}  // namespace bdm

#endif  // SYNAPSE_OP_H