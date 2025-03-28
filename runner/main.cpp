#include <G4BosonConstructor.hh>
#include <G4LeptonConstructor.hh>
#include <G4MesonConstructor.hh>
#include <G4BaryonConstructor.hh>
#include <G4IonConstructor.hh>
#include <G4ProcessManager.hh>
#include <G4StateManager.hh>
#include <G4RunManager.hh>
#include "G4ExcitationHandler.hh"
#include "G4ExceptionHandler.hh"

#include "G4LorentzVector.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"

#include "G4FermiBreakUpAN.hh"
#include "G4FermiBreakUpVI.hh"
#include <globals.hh>

#include "G4FermiNucleiProperties.hh"
#include <CLHEP/Units/PhysicalConstants.h>

#include <fstream>
#include <iostream>


void CalculateFragments(
  G4VFermiBreakUp& model,
  G4FermiAtomicMass mass,
  G4FermiChargeNumber charge,
  const std::string& dumpName,
  G4double step = 0.2,
  size_t tests = 1e4)
{
  std::vector<G4double> energyNucleonValues;
  std::vector<float> avgParts;
  for (G4double energyNucleon = 0.; energyNucleon <= 10.; energyNucleon += step) {
    size_t partsCounter = 0;
    auto additionalEnergy = energyNucleon * G4double(mass);
    for (size_t i = 0; i < tests; ++i) {
      auto vec = G4LorentzVector(
        0, 0, 0, G4FermiNucleiProperties::GetNuclearMass(mass, charge) + additionalEnergy);
      G4FragmentVector res;
      G4Fragment frag(G4int(mass), G4int(charge), vec);
      model.BreakFragment(&res, &frag);
      partsCounter += res.size();
    }
    avgParts.push_back(float(partsCounter) / float(tests));
    energyNucleonValues.push_back(energyNucleon);
  }

  std::ofstream out(dumpName);

  out << "energy,avg_count\n";
  for (size_t i = 0; i < avgParts.size(); ++i) {
    out << energyNucleonValues[i] << ',' << avgParts[i] << '\n';
  }

  std::cout << dumpName << ": done\n";
}

void CalculateMomentum(G4VFermiBreakUp& model, G4FermiAtomicMass mass, G4FermiChargeNumber charge,
                       const std::string& dumpName, G4double energy,
                       const G4Vector3D& momentum, size_t tests = 1e4)
{
  std::ofstream out(dumpName);
  auto vec = G4LorentzVector(
    momentum.x(), momentum.y(), momentum.z(),
    std::sqrt(std::pow(G4FermiNucleiProperties::GetNuclearMass(mass, charge) + energy, 2)
              + momentum.mag2()));
  out << vec / mass << '\n';
  std::vector<G4double> xComponent, yComponent, zComponent, magnitude;
  for (size_t i = 0; i < tests; ++i) {
    G4FragmentVector particles;
    G4Fragment frag(G4int(mass), G4int(charge), vec);
    model.BreakFragment(&particles, &frag);
    auto sum = G4LorentzVector();
    for (const auto particle : particles) {
      sum += particle->GetMomentum();
      out << particle->GetMomentum() / particle->GetA_asInt() << ' ';
    }
    out << '\n';
  }

  std::cout << dumpName << ": done\n";
}

#include <chrono>

int main()
{
  G4StateManager::GetStateManager()->SetExceptionHandler(new G4ExceptionHandler());

  G4BosonConstructor pCBos;
  pCBos.ConstructParticle();

  G4LeptonConstructor pCLept;
  pCLept.ConstructParticle();

  G4MesonConstructor pCMes;
  pCMes.ConstructParticle();

  G4BaryonConstructor pCBar;
  pCBar.ConstructParticle();

  G4IonConstructor pCIon;
  pCIon.ConstructParticle();

  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager(new G4ProcessManager(gion));

  G4StateManager::GetStateManager()->SetNewState(G4State_Init); // To let create ions
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ionTable = partTable->GetIonTable();
  partTable->SetReadiness();
  ionTable->CreateAllIon();
  ionTable->CreateAllIsomer();

  {
    auto model = G4FermiBreakUpAN();
    model.Initialise();
    auto begin = std::chrono::steady_clock::now();
    CalculateMomentum(model, 12_m, 6_c, "../Results/stat.data", 12 * 10 * CLHEP::GeV, {0, 0, 0});
    CalculateMomentum(model, 12_m, 6_c, "../Results/mov_x.data", 12 * 5 * CLHEP::MeV,
                      {12 * 10 * CLHEP::GeV, 0, 0});
    CalculateMomentum(model, 12_m, 6_c, "../Results/mov_y.data", 12 * 5 * CLHEP::MeV,
                      {0, 12 * 10 * CLHEP::GeV, 0});
    CalculateMomentum(model, 12_m, 6_c, "../Results/mov_z.data", 12 * 5 * CLHEP::MeV,
                      {0, 0, 12 * 10 * CLHEP::GeV});

    CalculateFragments(model, 12_m, 6_c, "../Results/C12.csv");

    CalculateFragments(model, 13_m, 6_c, "../Results/C13.csv");

    CalculateFragments(model, 12_m, 7_c, "../Results/N12.csv");

    CalculateFragments(model, 13_m, 7_c, "../Results/N13.csv");

    auto end = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;
  }

  {
    auto model = G4FermiBreakUpVI();
    model.Initialise();
    auto begin = std::chrono::steady_clock::now();
    CalculateMomentum(model, 12_m, 6_c, "../Results/oldstat.data", 12 * 10 * CLHEP::GeV, {0, 0, 0});
    CalculateMomentum(model, 12_m, 6_c, "../Results/oldmov_x.data", 12 * 5 * CLHEP::MeV,
                      {12 * 10 * CLHEP::GeV, 0, 0});
    CalculateMomentum(model, 12_m, 6_c, "../Results/oldmov_y.data", 12 * 5 * CLHEP::MeV,
                      {0, 12 * 10 * CLHEP::GeV, 0});
    CalculateMomentum(model, 12_m, 6_c, "../Results/oldmov_z.data", 12 * 5 * CLHEP::MeV,
                      {0, 0, 12 * 10 * CLHEP::GeV});

    CalculateFragments(model, 12_m, 6_c, "../Results/oldC12.csv");

    CalculateFragments(model, 13_m, 6_c, "../Results/oldC13.csv");

    CalculateFragments(model, 12_m, 7_c, "../Results/oldN12.csv");

    CalculateFragments(model, 13_m, 7_c, "../Results/oldN13.csv");

    auto end = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
  }
}
