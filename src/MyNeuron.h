#ifndef MY_NEURON_H_
#define MY_NEURON_H_

#include <vector>
#include "biodynamo.h"
#include "core/agent/cell_division_event.h"
#include "core/behavior/behavior.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

enum States { progenitor, dead };

class MyNeuron;
class Synapses;

class MyNeuron
    : public neuroscience::NeuronSoma {  // our object extends the Cell object
  // create the header with our new data member
  BDM_AGENT_HEADER(MyNeuron, neuroscience::NeuronSoma, 1);

 public:
  MyNeuron() {}
  explicit MyNeuron(const Real3& position) : Base(position) {}
  virtual ~MyNeuron() {}

  /// If MyNeuron divides, the daughter has to initialize its attributes
  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);

    if (auto* mother = dynamic_cast<MyNeuron*>(event.existing_agent)) {
      cell_colour_ = mother->cell_colour_;
      this->SetMass(mother->GetMass());
    }
  }

  void AddSynapse(MyNeuron* target, real_t distance, int strength = 1,
                  int time = 0);
  const std::vector<Synapses>& GetSynapses() const { return synapses_; }

  void SetCellColour(int cell_colour) { cell_colour_ = cell_colour; }
  int GetCellColour() const { return cell_colour_; }

  int GetState() const { return state_; }
  void SetState(int state) { state_ = state; }

  int GetCellType() const { return cell_type_; }
  void SetCellType(int cell_type) { cell_type_ = cell_type; }

  int cell_colour_;
  int state_ = States::progenitor;
  int cell_type_;
  std::vector<Synapses> synapses_;
};

class Synapses {
 public:
  // Default constructor
  Synapses()
      : source_(nullptr), target_(nullptr), distance_(-1.0), strength_(0) {}

  // Parameterized constructor
  Synapses(MyNeuron* source, MyNeuron* target, double distance = 0.0,
           int strength = 1, int time = 0)
      : source_(source),
        target_(target),
        distance_(distance),
        strength_(strength),
        time_(time) {}

  MyNeuron* GetSource() const { return source_; }
  MyNeuron* GetTarget() const { return target_; }
  double GetDistance() const { return distance_; }
  int GetStrength() const { return strength_; }
  int GetTime() const { return time_; }

  void IncreaseStrength(int amount = 1) { strength_ += amount; }

 private:
  MyNeuron* source_;
  MyNeuron* target_;
  double distance_;
  int strength_;
  int time_;
};

inline void MyNeuron::AddSynapse(MyNeuron* target, real_t distance,
                                 int strength, int time) {
  Synapses synapse(this, target, distance, strength, time);
  synapses_.push_back(synapse);
}

inline MyNeuron* FindParentNeuron(NeuriteElement* neurite) {
  MyNeuron* mother_neuron = nullptr;
  auto* dendrite = neurite;

  while (dendrite && !mother_neuron) {
    auto* mother = dendrite->GetMother().Get();
    mother_neuron = dynamic_cast<MyNeuron*>(mother);
    dendrite = dynamic_cast<NeuriteElement*>(mother);
  }

  return mother_neuron;
}

inline bool hasSynapse(const MyNeuron* neuronA, const MyNeuron* neuronB) {
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

inline void CreateSynapseBetweenNeurites(NeuriteElement* neurite1,
                                         NeuriteElement* neurite2,
                                         real_t distance = 0.0,
                                         int strength = 1, int time = 0) {
  // Trace back to parent neurons
  MyNeuron* neuronA = FindParentNeuron(neurite1);
  MyNeuron* neuronB = FindParentNeuron(neurite2);

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

struct PairHash {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2>& pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

inline void export_adjacency_matrix_with_all_neurons() {
  auto* sim = Simulation::GetActive();
  auto* rm = sim->GetResourceManager();

  // This will store the adjacency matrix.
  std::unordered_map<std::pair<int, int>, int, PairHash> adjacency;

  // This will store the cell types of each neuron.
  std::unordered_map<int, int> cell_types;

  // Populate the set of unique neuron UIDs and their cell types
  rm->ForEachAgent([&](Agent* agent) {
    auto* neuron = dynamic_cast<MyNeuron*>(agent);
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
  std::ofstream file("adjacency_matrix_all.csv");
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
