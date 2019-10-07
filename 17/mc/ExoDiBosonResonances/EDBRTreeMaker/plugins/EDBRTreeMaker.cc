// system include files
#include <iostream>
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"  
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <DataFormats/JetReco/interface/Jet.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
// DeepAK8
#include "NNKit/FatJetNN/interface/FatJetNN.h"
#include "NNKit/FatJetNN/interface/FatJetNNHelper.h"
#include "EDBRChannels.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include <TFormula.h>
#define Pi 3.141593
using namespace std;
using namespace deepntuples; // DeepAK8
//
// class declaration
//
/*
struct sortPt
{
   bool operator()(TLorentzVector* s1, TLorentzVector* s2) const
   {
      return s1->Pt() >= s2->Pt();
   }
} mysortPt;
*/
class EDBRTreeMaker : public edm::EDAnalyzer {
public:
  explicit EDBRTreeMaker(const edm::ParameterSet&);
  ~EDBRTreeMaker();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
  virtual void endRun(const edm::Run&, const edm::EventSetup&) override;

  virtual bool looseJetID( const pat::Jet& j);
  virtual const reco::Candidate* findLastParticle(const reco::Candidate *particle,int IDpdg);
  virtual const reco::Candidate* findFirstParticle(const reco::Candidate *particle,int IDpdg);

  virtual bool tightJetID( const pat::Jet& j);
  virtual float dEtaInSeed( const pat::Electron* ele) ;
  virtual void initJetCorrFactors( void );
  virtual void addTypeICorr( edm::Event const & event );
  virtual double getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
  virtual double getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
  math::XYZTLorentzVector getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType);
  std::vector<std::string>                    jecAK8PayloadNames_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8_            ;
  std::vector<std::string>                    jecAK8PayloadNamesGroomed_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8Groomed_            ;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8GroomedSD_            ;
  std::vector<std::string>                    jecAK8puppiPayloadNames_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8puppi_            ;
  std::vector<std::string>                    jecAK8puppiPayloadNamesGroomed_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK8puppiGroomed_            ;
  std::vector<std::string>                    jecAK4PayloadNames_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK4_            ;
  std::vector<std::string> offsetCorrLabel_;
  boost::shared_ptr<FactorizedJetCorrector> jecOffset_;
  edm::Handle< double >  rho_;
  edm::InputTag  METsRawLabel_;
  edm::Handle<pat::METCollection>  METs_;
  edm::Handle<pat::JetCollection> jets_;
  edm::Handle<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<pat::MuonCollection> muons_;
  edm::Handle<pat::METCollection>  reclusteredMETs_;
  edm::Handle<edm::View<reco::PFMET> >     pfMET_ ;
  edm::EDGetTokenT<pat::JetCollection> prunedjetInputToken_;
  edm::EDGetTokenT<pat::JetCollection> softdropjetInputToken_;
  edm::EDGetTokenT<pat::JetCollection> fatjetInputToken_;
  edm::EDGetTokenT<pat::JetCollection> puppijetInputToken_;

  edm::EDGetTokenT<pat::METCollection>  metInputToken_;
  edm::EDGetTokenT<pat::METCollection>  reclusteredmetInputToken_;
  std::vector<edm::EDGetTokenT<pat::METCollection>> mettokens;
  edm::Handle<pat::JetCollection> prunedjets_;
  edm::Handle<pat::JetCollection> softdropjets_;
  edm::Handle<pat::JetCollection> puppijets_;

  std::vector<edm::EDGetTokenT<pat::JetCollection>> jetTokens;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<pat::METCollection> reclusteredmetToken_;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
  edm::EDGetTokenT<pat::JetCollection> prunedjetToken_;
  edm::EDGetTokenT<pat::JetCollection> softdropjetToken_;
  edm::EDGetTokenT<pat::JetCollection> puppijetToken_;
  edm::Handle<pat::JetCollection> fatjets_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  std::vector<std::string> jetCorrLabel_;
  std::vector<std::string> jecAK4Labels;
  std::vector<std::string> jecAK8Labels;
  bool doCorrOnTheFly_;
  edm::EDGetTokenT<edm::TriggerResults> 		     noiseFilterToken_;
  edm::Handle< edm::TriggerResults> 			     noiseFilterBits_;
  std::string HBHENoiseFilter_Selector_;
  std::string HBHENoiseIsoFilter_Selector_;
  std::string GlobalHaloNoiseFilter_Selector_;
  std::string ECALDeadCellNoiseFilter_Selector_;
  std::string GoodVtxNoiseFilter_Selector_;
  std::string EEBadScNoiseFilter_Selector_;
  edm::EDGetTokenT<bool> badMuonNoiseFilter_Selector_;
  edm::EDGetTokenT<bool> badChargedHadronNoiseFilter_Selector_;
    // DeepAK8
          deepntuples::FatJetNN* fatjetNN_ = nullptr;
          deepntuples::FatJetNN* decorrNN_ = nullptr;
  // ----------member data ---------------------------
  TTree* outTree_;
  TTree* outTreew_;
  double MW_;
  int nmetmatch, nmetno;
  int nevent, run, ls;
  int nVtx;
  int numCands;
  int nLooseEle, nLooseMu;
  int nVetoEle, nVetoMu;
  int qnumber, qnumber2, qnumber3;
  int channel, lep;
  double vbfeta, vbfmjj;
  int      vbftag;
  int nj1, nj2;
  double met, metPhi;
  double jetAK8puppi_ptJEC, jetAK8puppi_eta, jetAK8puppi_phi;
  double jetAK8puppi_tau1,  jetAK8puppi_tau2, jetAK8puppi_tau3, jetAK8puppi_tau4;
  double jetAK8puppi_tau21, jetAK8puppi_tau42;
  double jetAK8puppi_sd, jetAK8puppi_sdJEC, jetAK8puppi_sdcorr;
  double jetAK8puppi_ptJEC_2, jetAK8puppi_eta_2, jetAK8puppi_phi_2;
  double jetAK8puppi_tau1_2,  jetAK8puppi_tau2_2, jetAK8puppi_tau3_2, jetAK8puppi_tau4_2;    
  double jetAK8puppi_tau21_2, jetAK8puppi_tau42_2;
  double jetAK8puppi_sd_2, jetAK8puppi_sdJEC_2, jetAK8puppi_sdcorr_2;
  double jetAK8puppi_ptJEC_3, jetAK8puppi_eta_3, jetAK8puppi_phi_3;
  double jetAK8puppi_tau1_3,  jetAK8puppi_tau2_3, jetAK8puppi_tau3_3, jetAK8puppi_tau4_3;    
  double jetAK8puppi_tau21_3, jetAK8puppi_tau42_3;
  double jetAK8puppi_sd_3, jetAK8puppi_sdJEC_3, jetAK8puppi_sdcorr_3;
  double jetAK8puppi_ptJEC_4, jetAK8puppi_eta_4, jetAK8puppi_phi_4;
  double jetAK8puppi_tau1_4,  jetAK8puppi_tau2_4, jetAK8puppi_tau3_4, jetAK8puppi_tau4_4;
  double jetAK8puppi_tau21_4, jetAK8puppi_tau42_4;
  double jetAK8puppi_sd_4, jetAK8puppi_sdJEC_4, jetAK8puppi_sdcorr_4;
  double jetAK8puppi_ptJEC_5, jetAK8puppi_eta_5, jetAK8puppi_phi_5;
  double jetAK8puppi_tau1_5,  jetAK8puppi_tau2_5, jetAK8puppi_tau3_5, jetAK8puppi_tau4_5;
  double jetAK8puppi_tau21_5, jetAK8puppi_tau42_5;
  double jetAK8puppi_sd_5, jetAK8puppi_sdJEC_5, jetAK8puppi_sdcorr_5;
  double jetAK8puppi_ptJEC_6, jetAK8puppi_eta_6, jetAK8puppi_phi_6;
  double jetAK8puppi_tau1_6,  jetAK8puppi_tau2_6, jetAK8puppi_tau3_6, jetAK8puppi_tau4_6;
  double jetAK8puppi_tau21_6, jetAK8puppi_tau42_6;
  double jetAK8puppi_sd_6, jetAK8puppi_sdJEC_6, jetAK8puppi_sdcorr_6;
  double jetAK8puppi_ptJEC_7, jetAK8puppi_eta_7, jetAK8puppi_phi_7;
  double jetAK8puppi_tau1_7,  jetAK8puppi_tau2_7, jetAK8puppi_tau3_7, jetAK8puppi_tau4_7;
  double jetAK8puppi_tau21_7, jetAK8puppi_tau42_7;
  double jetAK8puppi_sd_7, jetAK8puppi_sdJEC_7, jetAK8puppi_sdcorr_7;
  double jetAK8puppi_ptJEC_8, jetAK8puppi_eta_8, jetAK8puppi_phi_8;
  double jetAK8puppi_tau1_8,  jetAK8puppi_tau2_8, jetAK8puppi_tau3_8, jetAK8puppi_tau4_8;
  double jetAK8puppi_tau21_8, jetAK8puppi_tau42_8;
  double jetAK8puppi_sd_8, jetAK8puppi_sdJEC_8, jetAK8puppi_sdcorr_8;
  double triggerWeight, lumiWeight, pileupWeight;
  double ptq11, etaq11, phiq11, massq11;
  double ptq21, etaq21, phiq21, massq21;
  double ptq31, etaq31, phiq31, massq31;
  double ptq12, etaq12, phiq12, massq12;
  double ptq22, etaq22, phiq22, massq22;
  double ptq32, etaq32, phiq32, massq32;
  double delPhijetmet, delPhijetmet_2, delPhijetmet_3, delPhijetmet_4, delPhijetmet_5, delPhijetmet_6, delPhijetmet_7, delPhijetmet_8;
  double pt_graviton,pt_graviton1;
  double candMasspuppiJEC;
  double massww[3];
  double pweight[211];
    TLorentzVector ak8sj11,ak8sj21,ak8sj31,ak8sj41,ak8sj51;
    TLorentzVector ak8sj12,ak8sj22,ak8sj32,ak8sj42,ak8sj52;
    TLorentzVector ak8sj13,ak8sj23,ak8sj33,ak8sj43,ak8sj53;
    TLorentzVector ak8sj14,ak8sj24,ak8sj34,ak8sj44,ak8sj54;
    TLorentzVector ak8sj15,ak8sj25,ak8sj35,ak8sj45,ak8sj55;

    double jetAK8puppi_dnnTop, jetAK8puppi_dnnW,jetAK8puppi_dnnH4q,jetAK8puppi_dnnTop_2, jetAK8puppi_dnnW_2,jetAK8puppi_dnnH4q_2,jetAK8puppi_dnnTop_3, jetAK8puppi_dnnW_3,jetAK8puppi_dnnH4q_3,jetAK8puppi_dnnTop_4, jetAK8puppi_dnnW_4,jetAK8puppi_dnnH4q_4; //DeepAK8
    double jetAK8puppi_dnnZ,jetAK8puppi_dnnZbb,jetAK8puppi_dnnHbb,jetAK8puppi_dnnZ_2,jetAK8puppi_dnnZbb_2,jetAK8puppi_dnnHbb_2,jetAK8puppi_dnnZ_3,jetAK8puppi_dnnZbb_3,jetAK8puppi_dnnHbb_3,jetAK8puppi_dnnZ_4,jetAK8puppi_dnnZbb_4,jetAK8puppi_dnnHbb_4;
    double jetAK8puppi_dnnDecorrTop, jetAK8puppi_dnnDecorrW,jetAK8puppi_dnnDecorrH4q,jetAK8puppi_dnnDecorrTop_2, jetAK8puppi_dnnDecorrW_2, jetAK8puppi_dnnDecorrH4q_2,jetAK8puppi_dnnDecorrTop_3, jetAK8puppi_dnnDecorrW_3, jetAK8puppi_dnnDecorrH4q_3,jetAK8puppi_dnnDecorrTop_4, jetAK8puppi_dnnDecorrW_4,jetAK8puppi_dnnDecorrH4q_4; //Decorrelated DeepAK8
    double jetAK8puppi_dnnDecorrZ,jetAK8puppi_dnnDecorrZbb,jetAK8puppi_dnnDecorrHbb,jetAK8puppi_dnnDecorrZ_2,jetAK8puppi_dnnDecorrZbb_2,jetAK8puppi_dnnDecorrHbb_2,jetAK8puppi_dnnDecorrZ_3,jetAK8puppi_dnnDecorrZbb_3,jetAK8puppi_dnnDecorrHbb_3,jetAK8puppi_dnnDecorrZ_4,jetAK8puppi_dnnDecorrZbb_4,jetAK8puppi_dnnDecorrHbb_4;
    double jetAK8puppi_dnnDecorrbb,jetAK8puppi_dnnDecorrcc,jetAK8puppi_dnnDecorrbbnog,jetAK8puppi_dnnDecorrccnog,jetAK8puppi_dnnDecorrbb_2,jetAK8puppi_dnnDecorrcc_2,jetAK8puppi_dnnDecorrbbnog_2,jetAK8puppi_dnnDecorrccnog_2,jetAK8puppi_dnnDecorrbb_3,jetAK8puppi_dnnDecorrcc_3,jetAK8puppi_dnnDecorrbbnog_3,jetAK8puppi_dnnDecorrccnog_3,jetAK8puppi_dnnDecorrbb_4,jetAK8puppi_dnnDecorrcc_4,jetAK8puppi_dnnDecorrbbnog_4,jetAK8puppi_dnnDecorrccnog_4;
    double jetAK8puppi_dnnqcd,jetAK8puppi_dnntop,jetAK8puppi_dnnw,jetAK8puppi_dnnz,jetAK8puppi_dnnzbb,jetAK8puppi_dnnhbb,jetAK8puppi_dnnh4q,jetAK8puppi_dnnqcd_2,jetAK8puppi_dnntop_2,jetAK8puppi_dnnw_2,jetAK8puppi_dnnz_2,jetAK8puppi_dnnzbb_2,jetAK8puppi_dnnhbb_2,jetAK8puppi_dnnh4q_2,jetAK8puppi_dnnqcd_3,jetAK8puppi_dnntop_3,jetAK8puppi_dnnw_3,jetAK8puppi_dnnz_3,jetAK8puppi_dnnzbb_3,jetAK8puppi_dnnhbb_3,jetAK8puppi_dnnh4q_3,jetAK8puppi_dnnqcd_4,jetAK8puppi_dnntop_4,jetAK8puppi_dnnw_4,jetAK8puppi_dnnz_4,jetAK8puppi_dnnzbb_4,jetAK8puppi_dnnhbb_4,jetAK8puppi_dnnh4q_4;
    double jetAK8puppi_dnnDecorrqcd,jetAK8puppi_dnnDecorrtop,jetAK8puppi_dnnDecorrw,jetAK8puppi_dnnDecorrz,jetAK8puppi_dnnDecorrzbb,jetAK8puppi_dnnDecorrhbb,jetAK8puppi_dnnDecorrh4q,jetAK8puppi_dnnDecorrqcd_2,jetAK8puppi_dnnDecorrtop_2,jetAK8puppi_dnnDecorrw_2,jetAK8puppi_dnnDecorrz_2,jetAK8puppi_dnnDecorrzbb_2,jetAK8puppi_dnnDecorrhbb_2,jetAK8puppi_dnnDecorrh4q_2,jetAK8puppi_dnnDecorrqcd_3,jetAK8puppi_dnnDecorrtop_3,jetAK8puppi_dnnDecorrw_3,jetAK8puppi_dnnDecorrz_3,jetAK8puppi_dnnDecorrzbb_3,jetAK8puppi_dnnDecorrhbb_3,jetAK8puppi_dnnDecorrh4q_3,jetAK8puppi_dnnDecorrqcd_4,jetAK8puppi_dnnDecorrtop_4,jetAK8puppi_dnnDecorrw_4,jetAK8puppi_dnnDecorrz_4,jetAK8puppi_dnnDecorrzbb_4,jetAK8puppi_dnnDecorrhbb_4,jetAK8puppi_dnnDecorrh4q_4;

    double ptgenwl[5],etagenwl[5],phigenwl[5],massgenwl[5],taggenwl[5],taggenwmother[5];
    double genw_q1_pt[5],genw_q1_eta[5],genw_q1_phi[5],genw_q1_e[5],genw_q1_pdg[5];
    double genw_q2_pt[5],genw_q2_eta[5],genw_q2_phi[5],genw_q2_e[5],genw_q2_pdg[5];
    double ptgenzl[5],etagenzl[5],phigenzl[5],massgenzl[5],taggenzl[5];
    double ptgengl[15],etagengl[15],phigengl[15],egengl[15];
    double ptgenwf[5],etagenwf[5],phigenwf[5],massgenwf[5];
    double ptgenzf[5],etagenzf[5],phigenzf[5],massgenzf[5];
    double ptgengf[15],etagengf[15],phigengf[15],egengf[15];
    double gent_b_pt,gent_b_phi,gent_b_eta,gent_b_mass;
    double genantit_b_pt,genantit_b_phi,genantit_b_eta,genantit_b_mass;
    double gent_w_pt,gent_w_phi,gent_w_eta,gent_w_mass;
    double genantit_w_pt,genantit_w_phi,genantit_w_eta,genantit_w_mass;
    double gent_w_q1_pt,gent_w_q1_phi,gent_w_q1_eta,gent_w_q1_e,gent_w_q1_pdg;
    double genantit_w_q1_pt,genantit_w_q1_phi,genantit_w_q1_eta,genantit_w_q1_e,genantit_w_q1_pdg;
    double gent_w_q2_pt,gent_w_q2_phi,gent_w_q2_eta,gent_w_q2_e,gent_w_q2_pdg;
    double genantit_w_q2_pt,genantit_w_q2_phi,genantit_w_q2_eta,genantit_w_q2_e,genantit_w_q2_pdg;
    double ptgenq1l[5],etagenq1l[5],phigenq1l[5],egenq1l[5];
    double ptgenq1f[5],etagenq1f[5],phigenq1f[5],egenq1f[5];
    double ptgenq2l[5],etagenq2l[5],phigenq2l[5],egenq2l[5];
    double ptgenq2f[5],etagenq2f[5],phigenq2f[5],egenq2f[5];
    double ptgenq3l[5],etagenq3l[5],phigenq3l[5],egenq3l[5];
    double ptgenq3f[5],etagenq3f[5],phigenq3f[5],egenq3f[5];
    double ptgenq4l[5],etagenq4l[5],phigenq4l[5],egenq4l[5];
    double ptgenq4f[5],etagenq4f[5],phigenq4f[5],egenq4f[5];
    double ptgenq5l[5],etagenq5l[5],phigenq5l[5],egenq5l[5];
    double ptgenq5f[5],etagenq5f[5],phigenq5f[5],egenq5f[5];
    double mothergenq1f[5],mothergenq2f[5],mothergenq3f[5],mothergenq4f[5],mothergenq5f[5];
    double gent_w_tag,genantit_w_tag,mothergengf[15],mmothergengf[15],mmothergenq1f[5],mmothergenq2f[5],mmothergenq3f[5],mmothergenq4f[5],mmothergenq5f[5];

  double theWeight;
  double  nump=0;
  double  numm=0;
  double  npT, npIT;
  int     nBX;
  bool isHEEP;
  double iso, isoCut, et, trackIso;
  //Gen Level
  double gen_gra_m, gen_gra_pt, gen_gra_eta, gen_gra_phi;
  double gen_rad_m, gen_rad_pt, gen_rad_eta, gen_rad_phi;
  double gen_ele_pt, gen_ele_eta, gen_ele_phi, gen_ele_e;
  double gen_mu_pt, gen_mu_eta, gen_mu_phi, gen_mu_e;
  double gen_ele_pt_2, gen_ele_eta_2, gen_ele_phi_2, gen_ele_e_2;
  double gen_mu_pt_2, gen_mu_eta_2, gen_mu_phi_2, gen_mu_e_2;
  double gen_ele_pt_3, gen_ele_eta_3, gen_ele_phi_3, gen_ele_e_3;
  double gen_mu_pt_3, gen_mu_eta_3, gen_mu_phi_3, gen_mu_e_3;
  double gen_tau_pt, gen_tau_eta, gen_tau_phi, gen_tau_e;
  double gen_tau_pt_2, gen_tau_eta_2, gen_tau_phi_2, gen_tau_e_2;
  double gen_tau_pt_3, gen_tau_eta_3, gen_tau_phi_3, gen_tau_e_3;
  double gentop_pt, gentop_eta, gentop_phi, gentop_mass;
  double genantitop_pt, genantitop_eta, genantitop_phi, genantitop_mass;
  double ptGenVlep, etaGenVlep, phiGenVlep, massGenVlep;
  double ptGenVlep_2, etaGenVlep_2, phiGenVlep_2, massGenVlep_2;
  double ptGenVlep_3, etaGenVlep_3, phiGenVlep_3, massGenVlep_3;
  double ptGenVhad, etaGenVhad, phiGenVhad, massGenVhad;
  double ptGenV_2, etaGenV_2, phiGenV_2, massGenV_2;
  double ptGenV_3, etaGenV_3, phiGenV_3, massGenV_3;
  int status_1,status_2, status_3;

  double useless;

//  JEC
//  double corr_AK8, corr_AK81[3];
 //// double corr_AK8puppi[3],corr_AK8puppiSD[3];
//  double jetAK8_pt,jetAK8_mass,jetAK8_jec,jetAK8_e,jetAK8_eta,jetAK8_phi;
// double jetAK8_pt1[3], jetAK8_mass1[3], jetAK8_SF_mass1[3], jetAK8_SF_mass2[3], jetAK8_jec1[3],jetAK8_eta1[3];
//  double jetAK8puppi_pt1[4], jetAK8puppi_mass1[4], jetAK8puppi_eta1[4], jetAK8puppi_jec1[4], jetAK8puppiSD_jec1[4];
  double corr_AK8, corr_AK81[8];
  double corr_AK8puppi[8],corr_AK8puppiSD[8];
  double jetAK8_pt,jetAK8_mass,jetAK8_jec,jetAK8_e,jetAK8_eta,jetAK8_phi;
  double jetAK8_pt1[8], jetAK8_mass1[8], jetAK8_SF_mass1[8], jetAK8_SF_mass2[8], jetAK8_jec1[8],jetAK8_eta1[8];
  double jetAK8puppi_pt1[8], jetAK8puppi_mass1[8], jetAK8puppi_eta1[8], jetAK8puppi_jec1[8], jetAK8puppiSD_jec1[8];

  double corr;
  double METraw_et, METraw_phi, METraw_sumEt;
  double MET_et, MET_phi, MET_sumEt, MET_corrPx, MET_corrPy;
  // AK4 Jets
  int ak4jet_hf[8],ak4jet_pf[8],ak4jet_hf_2[8],ak4jet_pf_2[8];
  double ak4jet_pt[8],ak4jet_pt_uncorr[8],ak4jet_eta[8],ak4jet_phi[8],ak4jet_e[8], ak4jet_dr[8]; 
  double ak4jet_csv[8],ak4jet_icsv[8], ak4jet_IDLoose[8], ak4jet_IDTight[8];
  double ak4jet_deepcsvudsg[8],ak4jet_deepcsvb[8],ak4jet_deepcsvc[8],ak4jet_deepcsvbb[8],ak4jet_deepcsvcc[8];
    


  void setDummyValues();

  /// Parameters to steer the treeDumper
  int originalNEvents_;
  double crossSectionPb_;
  double targetLumiInvPb_;
  std::string EDBRChannel_;
  bool isGen_;
  bool isJEC_;
  bool RunOnSig_,RunOnMC_;
  std::vector<JetCorrectorParameters> vPar;
  std::map<std::string,double>  TypeICorrMap_;
  edm::InputTag mets_;

  //High Level Trigger
  HLTConfigProvider hltConfig;
  edm::EDGetTokenT<edm::TriggerResults> hltToken_;
  std::vector<std::string> muPaths1_, muPaths2_, muPaths3_, muPaths4_, muPaths5_, muPaths6_, muPaths7_, muPaths8_, muPaths9_;
  std::vector<std::string> muPaths1, muPaths2, muPaths3, muPaths4, muPaths5, muPaths6, muPaths7, muPaths8, muPaths9;
  std::vector<std::string> el1_, el2_, el3_, mu1_, mu2_, mu3_, mu4_;
  std::vector<std::string> el1, el2, el3, mu1, mu2, mu3, mu4;
  int  HLT_Mu1, HLT_Mu2, HLT_Mu3, HLT_Mu4, HLT_Mu5, HLT_Mu6, HLT_Mu7, HLT_Mu8, HLT_Mu9, HLT_Mu10, HLT_Mu11, HLT_Mu12, HLT_Mu13, HLT_Mu14, HLT_Mu15, HLT_Mu16;
  bool IDLoose, IDTight, IDLoose_2, IDTight_2, IDLoose_3, IDTight_3, IDLoose_4, IDTight_4;
          
// filter
  bool passFilter_HBHE_                   ;
  bool passFilter_HBHEIso_                ;
  bool passFilter_GlobalHalo_             ;
  bool passFilter_ECALDeadCell_           ;
  bool passFilter_GoodVtx_                ;
  bool passFilter_EEBadSc_                ;
  bool passFilter_badMuon_                ;
  bool passFilter_badChargedHadron_       ;

  edm::EDGetTokenT<edm::View<reco::Candidate>> leptonicVSrc_;
  edm::EDGetTokenT<edm::View<pat::Jet>> hadronicVSrc_;
//  edm::EDGetTokenT<pat::JetCollection> hadronicVSoftDropSrc_;
  edm::EDGetTokenT<edm::View<pat::Jet>> ak4jetsSrc_;
  edm::EDGetTokenT<edm::View<pat::Electron> > looseelectronToken_ ;
  edm::EDGetTokenT<edm::View<pat::Muon>> loosemuonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron> > vetoelectronToken_ ;
  edm::EDGetTokenT<edm::View<pat::Muon>> vetomuonToken_;
  edm::EDGetTokenT<edm::View<pat::Muon>> t1muSrc_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> gravitonSrc_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> metSrc_;
  edm::EDGetTokenT<GenEventInfoProduct> GenToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genSrc_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUToken_;
  edm::EDGetTokenT<edm::View<pat::Jet>> jetsAK8Label_;
  edm::EDGetTokenT<LHEEventProduct> LheToken_;
  edm::EDGetTokenT<LHERunInfoProduct> LhestrToken_;
};

//
// constructors and destructor
//
EDBRTreeMaker::EDBRTreeMaker(const edm::ParameterSet& iConfig):
  hltToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltToken"))),
  muPaths1_(iConfig.getParameter<std::vector<std::string>>("muPaths1")),
  muPaths2_(iConfig.getParameter<std::vector<std::string>>("muPaths2")),
  muPaths3_(iConfig.getParameter<std::vector<std::string>>("muPaths3")),
  muPaths4_(iConfig.getParameter<std::vector<std::string>>("muPaths4")),
  muPaths5_(iConfig.getParameter<std::vector<std::string>>("muPaths5")),
  muPaths6_(iConfig.getParameter<std::vector<std::string>>("muPaths6")),
  muPaths7_(iConfig.getParameter<std::vector<std::string>>("muPaths7")),
  muPaths8_(iConfig.getParameter<std::vector<std::string>>("muPaths8")),
  muPaths9_(iConfig.getParameter<std::vector<std::string>>("muPaths9")),//  noiseFilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilter")))
  el1_(iConfig.getParameter<std::vector<std::string>>("el1")),
  el2_(iConfig.getParameter<std::vector<std::string>>("el2")),
  el3_(iConfig.getParameter<std::vector<std::string>>("el3")),
  mu1_(iConfig.getParameter<std::vector<std::string>>("mu1")),
  mu2_(iConfig.getParameter<std::vector<std::string>>("mu2")),
  mu3_(iConfig.getParameter<std::vector<std::string>>("mu3")),
  mu4_(iConfig.getParameter<std::vector<std::string>>("mu4"))
{
  LheToken_        = consumes<LHEEventProduct> (iConfig.getParameter<edm::InputTag>( "lhe") );
  LhestrToken_=consumes<LHERunInfoProduct,edm::InRun> (iConfig.getParameter<edm::InputTag>( "lhe") );
  originalNEvents_ = iConfig.getParameter<int>("originalNEvents");
  crossSectionPb_  = iConfig.getParameter<double>("crossSectionPb");
  targetLumiInvPb_ = iConfig.getParameter<double>("targetLumiInvPb");
  EDBRChannel_     = iConfig.getParameter<std::string>("EDBRChannel");
  isGen_           = iConfig.getParameter<bool>("isGen");
  isJEC_           = iConfig.getParameter<bool>("isJEC");
  RunOnSig_        = iConfig.getParameter<bool>("RunOnSig");
  RunOnMC_           = iConfig.getParameter<bool>("RunOnMC");
  leptonicVSrc_=consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>( "leptonicVSrc") ) ;
  jetsAK8Label_      = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>( "ak8JetSrc") ) ;
  looseelectronToken_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("looseElectronSrc"))) ;
  loosemuonToken_    = (consumes<edm::View<pat::Muon> > (iConfig.getParameter<edm::InputTag>("looseMuonSrc")));
  vetoelectronToken_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("vetoElectronSrc"))) ;
  vetomuonToken_    = (consumes<edm::View<pat::Muon> > (iConfig.getParameter<edm::InputTag>("vetoMuonSrc")));
  ak4jetsSrc_      = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>( "ak4jetsSrc") ) ;
  hadronicVSrc_ = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("hadronicVSrc") ) ;
  //hadronicVSoftDropSrc_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("hadronicVSoftDropSrc") ) ; //80X

  jetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  puppijetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("puppijets"));
  fatjetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"));
  prunedjetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("prunedjets"));
  softdropjetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("softdropjets"));
  rhoToken_  = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  GenToken_=consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>( "generator") ) ;
  genSrc_      = consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>( "genSrc") ) ;
  PUToken_=consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileup") ) ;
  metSrc_      = consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>( "metSrc") ) ;
  gravitonSrc_      = consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>( "gravitonSrc") ) ;
  metToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  t1muSrc_      = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>( "t1muSrc") ) ;
  noiseFilterToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilter"));
  HBHENoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_HBHENoiseFilter");
  HBHENoiseIsoFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_HBHENoiseIsoFilter");
  GlobalHaloNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_GlobalTightHaloFilter");
  ECALDeadCellNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter");
  GoodVtxNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_goodVertices");
  EEBadScNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_eeBadScFilter");
  badMuonNoiseFilter_Selector_ = consumes<bool>(iConfig.getParameter<edm::InputTag> ("noiseFilterSelection_badMuon"));
  badChargedHadronNoiseFilter_Selector_  = consumes<bool>(iConfig.getParameter<edm::InputTag> ("noiseFilterSelection_badChargedHadron"));

  std::string jecpath = iConfig.getParameter<std::string>("jecpath");
  std::string tmpString;
  std::vector<std::string> tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8chsPayloadNames");
  for( unsigned int v = 0; v < tmpVec.size(); ++v ){
     tmpString = jecpath + tmpVec[v];
     jecAK8Labels.push_back(tmpString);
  }
  std::vector<std::string> jecAK8LabelsGroomed;
  tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8chsPayloadNamesGroomed");
  for( unsigned int v = 0; v < tmpVec.size(); ++v ){
     tmpString = jecpath + tmpVec[v];
     jecAK8LabelsGroomed.push_back(tmpString);
  }

  std::vector<std::string> jecAK8Labelspuppi;
  tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8puppiPayloadNames");
  for( unsigned int v = 0; v < tmpVec.size(); ++v ){
     tmpString = jecpath + tmpVec[v];
     jecAK8Labelspuppi.push_back(tmpString);
  }

  std::vector<std::string> jecAK8LabelspuppiGroomed;
  tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8puppiPayloadNamesGroomed");
  for( unsigned int v = 0; v < tmpVec.size(); ++v ){
     tmpString = jecpath + tmpVec[v];
     jecAK8LabelspuppiGroomed.push_back(tmpString);
  }



  std::vector<std::string> jecAK4Labels;
  tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK4chsPayloadNames");
  for( unsigned int v = 0; v < tmpVec.size(); ++v ){
     tmpString = jecpath + tmpVec[v];
     jecAK4Labels.push_back(tmpString);
  }

  /*=======================================================================================*/
  MW_=80.385;
  nmetmatch = 0;
  nmetno = 0;
  mettokens.push_back( metToken_ );
  mettokens.push_back( reclusteredmetToken_ );
  jetTokens.push_back( jetToken_ );
  jetTokens.push_back( fatjetToken_         );
  jetTokens.push_back( prunedjetToken_      );
  jetTokens.push_back( softdropjetToken_    );
  jetTokens.push_back( puppijetToken_      );

// add 3 up
 

  metInputToken_ = mettokens[0]; 
  reclusteredmetInputToken_ = mettokens[1];

  jetCorrLabel_ = jecAK4Labels;
  offsetCorrLabel_.push_back(jetCorrLabel_[0]);
 
  doCorrOnTheFly_ = false;
  if( jecAK4Labels.size() != 0 && jecAK8Labels.size() != 0 ){

     jecAK4PayloadNames_ = jecAK4Labels;
     //jecAK4PayloadNames_.pop_back();

     jecAK8PayloadNames_ = jecAK8Labels;
     //jecAK8PayloadNames_.pop_back();

     jecAK8PayloadNamesGroomed_ = jecAK8LabelsGroomed;
     //jecAK8PayloadNamesGroomed_.pop_back();

     jecAK8puppiPayloadNames_ = jecAK8Labelspuppi;
     jecAK8puppiPayloadNamesGroomed_ = jecAK8LabelspuppiGroomed;


  fatjetInputToken_ = jetTokens[1];
  prunedjetInputToken_ = jetTokens[2];
  softdropjetInputToken_ = jetTokens[3];
  puppijetInputToken_ = jetTokens[4];
// add 3 up
     initJetCorrFactors();

     doCorrOnTheFly_ = true;

  }

 
    if(EDBRChannel_ == "VZ_CHANNEL")
        channel=VZ_CHANNEL;
    else if(EDBRChannel_ == "VW_CHANNEL")
        channel=VW_CHANNEL;
    else if(EDBRChannel_ == "VH_CHANNEL")
        channel=VH_CHANNEL;
    else {
        cms::Exception ex("InvalidConfiguration");
        ex << "Unknown channel " << EDBRChannel_<< ". Please check EDBRTreeMaker.cc for allowed values.";
    throw ex;
    }
  
    // DeepAK8
    // initialize the FatJetNN class in the constructor
    auto cc = consumesCollector();
    fatjetNN_ = new deepntuples::FatJetNN(iConfig, cc, 0.8); // jetR=0.8
    // load json for input variable transformation
    fatjetNN_->load_json(edm::FileInPath("NNKit/data/ak8/full/preprocessing.json").fullPath());
    // load DNN model and parameter files
    fatjetNN_->load_model(edm::FileInPath("NNKit/data/ak8/full/resnet-symbol.json").fullPath(),                  edm::FileInPath("NNKit/data/ak8/full/resnet.params").fullPath());
    // -----------
    
    // Decorrelated DeepAK8
    decorrNN_ = new FatJetNN(iConfig, cc, 0.8);
    decorrNN_->load_json(edm::FileInPath("NNKit/data/ak8/decorrelated/preprocessing.json").fullPath());
    decorrNN_->load_model(edm::FileInPath("NNKit/data/ak8/decorrelated/resnet-symbol.json").fullPath(),
    edm::FileInPath("NNKit/data/ak8/decorrelated/resnet.params").fullPath());
    // -----------
    //now do what ever initialization is needed
    edm::Service<TFileService> fs;

    outTree_ = fs->make<TTree>("EDBRCandidates","EDBR Candidates");
    outTreew_ = fs->make<TTree>("EDBRCandidatesw","EDBR Candidates");
    outTree_->Branch("vbfeta"    ,&vbfeta   ,"vbfeta/D"   );
    outTree_->Branch("vbfmjj"    ,&vbfmjj   ,"vbfmjj/D"   );
    outTree_->Branch("vbftag"    ,&vbftag   ,"vbftag/I"   );
    outTree_->Branch("ak8sj11"    ,&ak8sj11 ) ;
    outTree_->Branch("ak8sj21"    ,&ak8sj21  );
    outTree_->Branch("ak8sj31"    ,&ak8sj31   );
  //  outTree_->Branch("ak8sj41"    ,&ak8sj41  );
  //  outTree_->Branch("ak8sj51"    ,&ak8sj51  );
    outTree_->Branch("ak8sj12"    ,&ak8sj12 ) ;
    outTree_->Branch("ak8sj22"    ,&ak8sj22  );
    outTree_->Branch("ak8sj32"    ,&ak8sj32   );
 //   outTree_->Branch("ak8sj42"    ,&ak8sj42  );
 //   outTree_->Branch("ak8sj52"    ,&ak8sj52  );   
    outTree_->Branch("ak8sj13"    ,&ak8sj13 ) ;
    outTree_->Branch("ak8sj23"    ,&ak8sj23  );
    outTree_->Branch("ak8sj33"    ,&ak8sj33   );
 //   outTree_->Branch("ak8sj43"    ,&ak8sj43  );
 //   outTree_->Branch("ak8sj53"    ,&ak8sj53  );   
    outTree_->Branch("ak8sj14"    ,&ak8sj14 ) ;
    outTree_->Branch("ak8sj24"    ,&ak8sj24  );
    outTree_->Branch("ak8sj34"    ,&ak8sj34   );
 //   outTree_->Branch("ak8sj44"    ,&ak8sj44  );
 //   outTree_->Branch("ak8sj54"    ,&ak8sj54  );
    outTree_->Branch("ak8sj15"    ,&ak8sj14 ) ;
    outTree_->Branch("ak8sj25"    ,&ak8sj24  );
    outTree_->Branch("ak8sj35"    ,&ak8sj34   );
// DeepAK8

    outTree_->Branch("jetAK8puppi_dnnTop"         ,&jetAK8puppi_dnnTop       ,"jetAK8puppi_dnnTop/D"         );
    outTree_->Branch("jetAK8puppi_dnnW"           ,&jetAK8puppi_dnnW         ,"jetAK8puppi_dnnW/D"           );
    outTree_->Branch("jetAK8puppi_dnnH4q"         ,&jetAK8puppi_dnnH4q       ,"jetAK8puppi_dnnH4q/D"         );
    outTree_->Branch("jetAK8puppi_dnnTop_2"         ,&jetAK8puppi_dnnTop_2       ,"jetAK8puppi_dnnTop_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnW_2"           ,&jetAK8puppi_dnnW_2         ,"jetAK8puppi_dnnW_2/D"           );
    outTree_->Branch("jetAK8puppi_dnnH4q_2"         ,&jetAK8puppi_dnnH4q_2       ,"jetAK8puppi_dnnH4q_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnTop_3"         ,&jetAK8puppi_dnnTop_3       ,"jetAK8puppi_dnnTop_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnW_3"           ,&jetAK8puppi_dnnW_3         ,"jetAK8puppi_dnnW_3/D"           );
    outTree_->Branch("jetAK8puppi_dnnH4q_3"         ,&jetAK8puppi_dnnH4q_3       ,"jetAK8puppi_dnnH4q_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnTop_4"         ,&jetAK8puppi_dnnTop_4       ,"jetAK8puppi_dnnTop_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnW_4"           ,&jetAK8puppi_dnnW_4         ,"jetAK8puppi_dnnW_4/D"           );
    outTree_->Branch("jetAK8puppi_dnnH4q_4"         ,&jetAK8puppi_dnnH4q_4       ,"jetAK8puppi_dnnH4q_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnZ"         ,&jetAK8puppi_dnnZ       ,"jetAK8puppi_dnnZ/D"         );
    outTree_->Branch("jetAK8puppi_dnnZbb"         ,&jetAK8puppi_dnnZbb       ,"jetAK8puppi_dnnZbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnHbb"         ,&jetAK8puppi_dnnHbb       ,"jetAK8puppi_dnnHbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnZ_2"         ,&jetAK8puppi_dnnZ_2       ,"jetAK8puppi_dnnZ_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnZbb_2"         ,&jetAK8puppi_dnnZbb_2       ,"jetAK8puppi_dnnZbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnHbb_2"         ,&jetAK8puppi_dnnHbb_2       ,"jetAK8puppi_dnnHbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnZ_3"         ,&jetAK8puppi_dnnZ_3       ,"jetAK8puppi_dnnZ_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnZbb_3"         ,&jetAK8puppi_dnnZbb_3       ,"jetAK8puppi_dnnZbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnHbb_3"         ,&jetAK8puppi_dnnHbb_3       ,"jetAK8puppi_dnnHbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnZ_4"         ,&jetAK8puppi_dnnZ_4       ,"jetAK8puppi_dnnZ_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnZbb_4"         ,&jetAK8puppi_dnnZbb_4       ,"jetAK8puppi_dnnZbb_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnHbb_4"         ,&jetAK8puppi_dnnHbb_4       ,"jetAK8puppi_dnnHbb_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnqcd"         ,&jetAK8puppi_dnnqcd       ,"jetAK8puppi_dnnqcd/D"         );
    outTree_->Branch("jetAK8puppi_dnntop"         ,&jetAK8puppi_dnntop       ,"jetAK8puppi_dnntop/D"         );
    outTree_->Branch("jetAK8puppi_dnnw"         ,&jetAK8puppi_dnnw       ,"jetAK8puppi_dnnw/D"         );
    outTree_->Branch("jetAK8puppi_dnnz"         ,&jetAK8puppi_dnnz       ,"jetAK8puppi_dnnz/D"         );
    outTree_->Branch("jetAK8puppi_dnnzbb"         ,&jetAK8puppi_dnnzbb       ,"jetAK8puppi_dnnzbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnhbb"         ,&jetAK8puppi_dnnhbb       ,"jetAK8puppi_dnnhbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnh4q"         ,&jetAK8puppi_dnnh4q       ,"jetAK8puppi_dnnh4q/D"         );
    outTree_->Branch("jetAK8puppi_dnnqcd_2"         ,&jetAK8puppi_dnnqcd_2       ,"jetAK8puppi_dnnqcd_2/D"         );
    outTree_->Branch("jetAK8puppi_dnntop_2"         ,&jetAK8puppi_dnntop_2       ,"jetAK8puppi_dnntop_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnw_2"         ,&jetAK8puppi_dnnw_2       ,"jetAK8puppi_dnnw_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnz_2"         ,&jetAK8puppi_dnnz_2       ,"jetAK8puppi_dnnz_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnzbb_2"         ,&jetAK8puppi_dnnzbb_2       ,"jetAK8puppi_dnnzbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnhbb_2"         ,&jetAK8puppi_dnnhbb_2       ,"jetAK8puppi_dnnhbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnh4q_2"         ,&jetAK8puppi_dnnh4q_2       ,"jetAK8puppi_dnnh4q_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnqcd_3"         ,&jetAK8puppi_dnnqcd_3       ,"jetAK8puppi_dnnqcd_3/D"         );
    outTree_->Branch("jetAK8puppi_dnntop_3"         ,&jetAK8puppi_dnntop_3       ,"jetAK8puppi_dnntop_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnw_3"         ,&jetAK8puppi_dnnw_3       ,"jetAK8puppi_dnnw_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnz_3"         ,&jetAK8puppi_dnnz_3       ,"jetAK8puppi_dnnz_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnzbb_3"         ,&jetAK8puppi_dnnzbb_3       ,"jetAK8puppi_dnnzbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnhbb_3"         ,&jetAK8puppi_dnnhbb_3       ,"jetAK8puppi_dnnhbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnh4q_3"         ,&jetAK8puppi_dnnh4q_3       ,"jetAK8puppi_dnnh4q_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnqcd_4"         ,&jetAK8puppi_dnnqcd_4       ,"jetAK8puppi_dnnqcd_4/D"         );
    outTree_->Branch("jetAK8puppi_dnntop_4"         ,&jetAK8puppi_dnntop_4       ,"jetAK8puppi_dnntop_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnw_4"         ,&jetAK8puppi_dnnw_4       ,"jetAK8puppi_dnnw_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnz_4"         ,&jetAK8puppi_dnnz_4       ,"jetAK8puppi_dnnz_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnzbb_4"         ,&jetAK8puppi_dnnzbb_4       ,"jetAK8puppi_dnnzbb_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnhbb_4"         ,&jetAK8puppi_dnnhbb_4       ,"jetAK8puppi_dnnhbb_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnh4q_4"         ,&jetAK8puppi_dnnh4q_4       ,"jetAK8puppi_dnnh4q_4/D"         );

//Decorrelated DeepAK8
    outTree_->Branch("jetAK8puppi_dnnDecorrTop"         ,&jetAK8puppi_dnnDecorrTop       ,"jetAK8puppi_dnnDecorrTop/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrW"           ,&jetAK8puppi_dnnDecorrW         ,"jetAK8puppi_dnnDecorrW/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrH4q"         ,&jetAK8puppi_dnnDecorrH4q       ,"jetAK8puppi_dnnDecorrH4q/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrTop_2"         ,&jetAK8puppi_dnnDecorrTop_2       ,"jetAK8puppi_dnnDecorrTop_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrW_2"           ,&jetAK8puppi_dnnDecorrW_2         ,"jetAK8puppi_dnnDecorrW_2/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrH4q_2"         ,&jetAK8puppi_dnnDecorrH4q_2       ,"jetAK8puppi_dnnDecorrH4q_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrTop_3"         ,&jetAK8puppi_dnnDecorrTop_3       ,"jetAK8puppi_dnnDecorrTop_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrW_3"           ,&jetAK8puppi_dnnDecorrW_3         ,"jetAK8puppi_dnnDecorrW_3/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrH4q_3"         ,&jetAK8puppi_dnnDecorrH4q_3       ,"jetAK8puppi_dnnDecorrH4q_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrTop_4"         ,&jetAK8puppi_dnnDecorrTop_4       ,"jetAK8puppi_dnnDecorrTop_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrW_4"           ,&jetAK8puppi_dnnDecorrW_4         ,"jetAK8puppi_dnnDecorrW_4/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrH4q_4"         ,&jetAK8puppi_dnnDecorrH4q_4       ,"jetAK8puppi_dnnDecorrH4q_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZ"         ,&jetAK8puppi_dnnDecorrZ       ,"jetAK8puppi_dnnDecorrZ/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZbb"         ,&jetAK8puppi_dnnDecorrZbb       ,"jetAK8puppi_dnnDecorrZbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrHbb"         ,&jetAK8puppi_dnnDecorrHbb       ,"jetAK8puppi_dnnDecorrHbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZ_2"         ,&jetAK8puppi_dnnDecorrZ_2       ,"jetAK8puppi_dnnDecorrZ_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZbb_2"         ,&jetAK8puppi_dnnDecorrZbb_2       ,"jetAK8puppi_dnnDecorrZbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrHbb_2"         ,&jetAK8puppi_dnnDecorrHbb_2       ,"jetAK8puppi_dnnDecorrHbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZ_3"         ,&jetAK8puppi_dnnDecorrZ_3       ,"jetAK8puppi_dnnDecorrZ_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZbb_3"         ,&jetAK8puppi_dnnDecorrZbb_3       ,"jetAK8puppi_dnnDecorrZbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrHbb_3"         ,&jetAK8puppi_dnnDecorrHbb_3       ,"jetAK8puppi_dnnDecorrHbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZ_4"         ,&jetAK8puppi_dnnDecorrZ_4       ,"jetAK8puppi_dnnDecorrZ_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZbb_4"         ,&jetAK8puppi_dnnDecorrZbb_4       ,"jetAK8puppi_dnnDecorrZbb_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrHbb_4"         ,&jetAK8puppi_dnnDecorrHbb_4       ,"jetAK8puppi_dnnDecorrHbb_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbb"         ,&jetAK8puppi_dnnDecorrbb       ,"jetAK8puppi_dnnDecorrbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrcc"         ,&jetAK8puppi_dnnDecorrcc       ,"jetAK8puppi_dnnDecorrcc/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbbnog"         ,&jetAK8puppi_dnnDecorrbbnog       ,"jetAK8puppi_dnnDecorrbbnog/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrccnog"         ,&jetAK8puppi_dnnDecorrccnog       ,"jetAK8puppi_dnnDecorrccnog/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbb_2"         ,&jetAK8puppi_dnnDecorrbb_2       ,"jetAK8puppi_dnnDecorrbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrcc_2"         ,&jetAK8puppi_dnnDecorrcc_2       ,"jetAK8puppi_dnnDecorrcc_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbbnog_2"         ,&jetAK8puppi_dnnDecorrbbnog_2       ,"jetAK8puppi_dnnDecorrbbnog_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrccnog_2"         ,&jetAK8puppi_dnnDecorrccnog_2       ,"jetAK8puppi_dnnDecorrccnog_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbb_3"         ,&jetAK8puppi_dnnDecorrbb_3       ,"jetAK8puppi_dnnDecorrbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrcc_3"         ,&jetAK8puppi_dnnDecorrcc_3       ,"jetAK8puppi_dnnDecorrcc_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbbnog_3"         ,&jetAK8puppi_dnnDecorrbbnog_3       ,"jetAK8puppi_dnnDecorrbbnog_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrccnog_3"         ,&jetAK8puppi_dnnDecorrccnog_3       ,"jetAK8puppi_dnnDecorrccnog_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbb_4"         ,&jetAK8puppi_dnnDecorrbb_4       ,"jetAK8puppi_dnnDecorrbb_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrcc_4"         ,&jetAK8puppi_dnnDecorrcc_4       ,"jetAK8puppi_dnnDecorrcc_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbbnog_4"         ,&jetAK8puppi_dnnDecorrbbnog_4       ,"jetAK8puppi_dnnDecorrbbnog_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrccnog_4"         ,&jetAK8puppi_dnnDecorrccnog_4       ,"jetAK8puppi_dnnDecorrccnog_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrqcd"         ,&jetAK8puppi_dnnDecorrqcd       ,"jetAK8puppi_dnnDecorrqcd/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrtop"         ,&jetAK8puppi_dnnDecorrtop       ,"jetAK8puppi_dnnDecorrtop/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrw"         ,&jetAK8puppi_dnnDecorrw       ,"jetAK8puppi_dnnDecorrw/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrz"         ,&jetAK8puppi_dnnDecorrz       ,"jetAK8puppi_dnnDecorrz/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrzbb"         ,&jetAK8puppi_dnnDecorrzbb       ,"jetAK8puppi_dnnDecorrzbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrhbb"         ,&jetAK8puppi_dnnDecorrhbb       ,"jetAK8puppi_dnnDecorrhbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrh4q"         ,&jetAK8puppi_dnnDecorrh4q       ,"jetAK8puppi_dnnDecorrh4q/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrqcd_2"         ,&jetAK8puppi_dnnDecorrqcd_2       ,"jetAK8puppi_dnnDecorrqcd_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrtop_2"         ,&jetAK8puppi_dnnDecorrtop_2       ,"jetAK8puppi_dnnDecorrtop_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrw_2"         ,&jetAK8puppi_dnnDecorrw_2       ,"jetAK8puppi_dnnDecorrw_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrz_2"         ,&jetAK8puppi_dnnDecorrz_2       ,"jetAK8puppi_dnnDecorrz_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrzbb_2"         ,&jetAK8puppi_dnnDecorrzbb_2       ,"jetAK8puppi_dnnDecorrzbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrhbb_2"         ,&jetAK8puppi_dnnDecorrhbb_2       ,"jetAK8puppi_dnnDecorrhbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrh4q_2"         ,&jetAK8puppi_dnnDecorrh4q_2       ,"jetAK8puppi_dnnDecorrh4q_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrqcd_3"         ,&jetAK8puppi_dnnDecorrqcd_3       ,"jetAK8puppi_dnnDecorrqcd_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrtop_3"         ,&jetAK8puppi_dnnDecorrtop_3       ,"jetAK8puppi_dnnDecorrtop_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrw_3"         ,&jetAK8puppi_dnnDecorrw_3       ,"jetAK8puppi_dnnDecorrw_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrz_3"         ,&jetAK8puppi_dnnDecorrz_3       ,"jetAK8puppi_dnnDecorrz_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrzbb_3"         ,&jetAK8puppi_dnnDecorrzbb_3       ,"jetAK8puppi_dnnDecorrzbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrhbb_3"         ,&jetAK8puppi_dnnDecorrhbb_3       ,"jetAK8puppi_dnnDecorrhbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrh4q_3"         ,&jetAK8puppi_dnnDecorrh4q_3       ,"jetAK8puppi_dnnDecorrh4q_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrqcd_4"         ,&jetAK8puppi_dnnDecorrqcd_4       ,"jetAK8puppi_dnnDecorrqcd_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrtop_4"         ,&jetAK8puppi_dnnDecorrtop_4       ,"jetAK8puppi_dnnDecorrtop_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrw_4"         ,&jetAK8puppi_dnnDecorrw_4       ,"jetAK8puppi_dnnDecorrw_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrz_4"         ,&jetAK8puppi_dnnDecorrz_4       ,"jetAK8puppi_dnnDecorrz_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrzbb_4"         ,&jetAK8puppi_dnnDecorrzbb_4       ,"jetAK8puppi_dnnDecorrzbb_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrhbb_4"         ,&jetAK8puppi_dnnDecorrhbb_4       ,"jetAK8puppi_dnnDecorrhbb_4/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrh4q_4"         ,&jetAK8puppi_dnnDecorrh4q_4       ,"jetAK8puppi_dnnDecorrh4q_4/D"         );

    outTree_->Branch("nj1"    ,&nj1   ,"nj1/I"   );
    outTree_->Branch("nj2"    ,&nj2   ,"nj2/I"   ); 
/// Basic event quantities
  outTree_->Branch("run"             ,&run            ,"run/I");//
  outTree_->Branch("ls"              ,&ls             ,"ls/I"             );//Synch
  outTree_->Branch("event"           ,&nevent         ,"event/I"          );
  outTree_->Branch("nVtx"            ,&nVtx           ,"nVtx/I"           );
//  outTree_->Branch("nLooseEle"       ,&nLooseEle      ,"nLooseEle/I");//
//  outTree_->Branch("nLooseMu"        ,&nLooseMu       ,"nLooseMu/I");//
//  outTree_->Branch("nVetoEle"       ,&nVetoEle      ,"nVetoEle/I");//
//  outTree_->Branch("nVetoMu"        ,&nVetoMu       ,"nVetoMu/I");//
  outTree_->Branch("jetAK8puppi_ptJEC"          ,&jetAK8puppi_ptJEC         ,"jetAK8puppi_ptJEC/D"      );
  outTree_->Branch("jetAK8puppi_eta"          ,&jetAK8puppi_eta         ,"jetAK8puppi_eta/D"         );
  outTree_->Branch("jetAK8puppi_phi"          ,&jetAK8puppi_phi         ,"jetAK8puppi_phi/D"         );
  outTree_->Branch("jetAK8puppi_tau1"          ,&jetAK8puppi_tau1         ,"jetAK8puppi_tau1/D"         );
  outTree_->Branch("jetAK8puppi_tau2"          ,&jetAK8puppi_tau2         ,"jetAK8puppi_tau2/D"         );
  outTree_->Branch("jetAK8puppi_tau3"          ,&jetAK8puppi_tau3         ,"jetAK8puppi_tau3/D"         );
  outTree_->Branch("jetAK8puppi_tau21"          ,&jetAK8puppi_tau21         ,"jetAK8puppi_tau21/D"         );
  outTree_->Branch("jetAK8puppi_tau4"          ,&jetAK8puppi_tau4         ,"jetAK8puppi_tau4/D"         );
  outTree_->Branch("jetAK8puppi_tau42"          ,&jetAK8puppi_tau42         ,"jetAK8puppi_tau42/D"         );
  outTree_->Branch("jetAK8puppi_sd"          ,&jetAK8puppi_sd         ,"jetAK8puppi_sd/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC"          ,&jetAK8puppi_sdJEC         ,"jetAK8puppi_sdJEC/D"         );
  outTree_->Branch("jetAK8puppi_sdcorr"          ,&jetAK8puppi_sdcorr         ,"jetAK8puppi_sdcorr/D"         );

  outTree_->Branch("jetAK8puppi_ptJEC_2"          ,&jetAK8puppi_ptJEC_2         ,"jetAK8puppi_ptJEC_2/D"    );
  outTree_->Branch("jetAK8puppi_eta_2"          ,&jetAK8puppi_eta_2         ,"jetAK8puppi_eta_2/D"         );
  outTree_->Branch("jetAK8puppi_phi_2"          ,&jetAK8puppi_phi_2         ,"jetAK8puppi_phi_2/D"         );
  outTree_->Branch("jetAK8puppi_tau1_2"          ,&jetAK8puppi_tau1_2         ,"jetAK8puppi_tau1_2/D"         );
  outTree_->Branch("jetAK8puppi_tau2_2"          ,&jetAK8puppi_tau2_2         ,"jetAK8puppi_tau2_2/D"         );
  outTree_->Branch("jetAK8puppi_tau3_2"          ,&jetAK8puppi_tau3_2         ,"jetAK8puppi_tau3_2/D"         );
  outTree_->Branch("jetAK8puppi_tau21_2"          ,&jetAK8puppi_tau21_2         ,"jetAK8puppi_tau21_2/D"    );
  outTree_->Branch("jetAK8puppi_tau4_2"          ,&jetAK8puppi_tau4_2         ,"jetAK8puppi_tau4_2/D"       );
  outTree_->Branch("jetAK8puppi_tau42_2"          ,&jetAK8puppi_tau42_2         ,"jetAK8puppi_tau42_2/D"  );
  outTree_->Branch("jetAK8puppi_sd_2"          ,&jetAK8puppi_sd_2         ,"jetAK8puppi_sd_2/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_2"          ,&jetAK8puppi_sdJEC_2         ,"jetAK8puppi_sdJEC_2/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_2"          ,&jetAK8puppi_sdcorr_2         ,"jetAK8puppi_sdcorr_2/D"  );

  outTree_->Branch("jetAK8puppi_ptJEC_3"          ,&jetAK8puppi_ptJEC_3         ,"jetAK8puppi_ptJEC_3/D"    );
  outTree_->Branch("jetAK8puppi_eta_3"          ,&jetAK8puppi_eta_3         ,"jetAK8puppi_eta_3/D"         );
  outTree_->Branch("jetAK8puppi_phi_3"          ,&jetAK8puppi_phi_3         ,"jetAK8puppi_phi_3/D"         );
  outTree_->Branch("jetAK8puppi_tau1_3"          ,&jetAK8puppi_tau1_3         ,"jetAK8puppi_tau1_3/D"         );
  outTree_->Branch("jetAK8puppi_tau2_3"          ,&jetAK8puppi_tau2_3         ,"jetAK8puppi_tau2_3/D"         );
  outTree_->Branch("jetAK8puppi_tau3_3"          ,&jetAK8puppi_tau3_3         ,"jetAK8puppi_tau3_3/D"         );
  outTree_->Branch("jetAK8puppi_tau21_3"          ,&jetAK8puppi_tau21_3         ,"jetAK8puppi_tau21_3/D"    );
  outTree_->Branch("jetAK8puppi_tau4_3"          ,&jetAK8puppi_tau4_3         ,"jetAK8puppi_tau4_3/D"       );
  outTree_->Branch("jetAK8puppi_tau42_3"          ,&jetAK8puppi_tau42_3         ,"jetAK8puppi_tau42_3/D"  );
  outTree_->Branch("jetAK8puppi_sd_3"          ,&jetAK8puppi_sd_3         ,"jetAK8puppi_sd_3/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_3"          ,&jetAK8puppi_sdJEC_3         ,"jetAK8puppi_sdJEC_3/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_3"          ,&jetAK8puppi_sdcorr_3         ,"jetAK8puppi_sdcorr_3/D"  );

  outTree_->Branch("jetAK8puppi_ptJEC_4"          ,&jetAK8puppi_ptJEC_4         ,"jetAK8puppi_ptJEC_4/D"    );
  outTree_->Branch("jetAK8puppi_eta_4"          ,&jetAK8puppi_eta_4         ,"jetAK8puppi_eta_4/D"         );
  outTree_->Branch("jetAK8puppi_phi_4"          ,&jetAK8puppi_phi_4         ,"jetAK8puppi_phi_4/D"         );
  outTree_->Branch("jetAK8puppi_tau1_4"          ,&jetAK8puppi_tau1_4         ,"jetAK8puppi_tau1_4/D"         );
  outTree_->Branch("jetAK8puppi_tau2_4"          ,&jetAK8puppi_tau2_4         ,"jetAK8puppi_tau2_4/D"         );
  outTree_->Branch("jetAK8puppi_tau3_4"          ,&jetAK8puppi_tau3_4         ,"jetAK8puppi_tau3_4/D"         );
  outTree_->Branch("jetAK8puppi_tau21_4"          ,&jetAK8puppi_tau21_4         ,"jetAK8puppi_tau21_4/D"    );
  outTree_->Branch("jetAK8puppi_tau4_4"          ,&jetAK8puppi_tau4_4         ,"jetAK8puppi_tau4_4/D"       );
  outTree_->Branch("jetAK8puppi_tau42_4"          ,&jetAK8puppi_tau42_4         ,"jetAK8puppi_tau42_4/D"  );
  outTree_->Branch("jetAK8puppi_sd_4"          ,&jetAK8puppi_sd_4         ,"jetAK8puppi_sd_4/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_4"          ,&jetAK8puppi_sdJEC_4         ,"jetAK8puppi_sdJEC_4/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_4"          ,&jetAK8puppi_sdcorr_4         ,"jetAK8puppi_sdcorr_4/D"  );

/*
  outTree_->Branch("jetAK8puppi_ptJEC_5"          ,&jetAK8puppi_ptJEC_5         ,"jetAK8puppi_ptJEC_5/D"    );
  outTree_->Branch("jetAK8puppi_eta_5"          ,&jetAK8puppi_eta_5         ,"jetAK8puppi_eta_5/D"         );
  outTree_->Branch("jetAK8puppi_phi_5"          ,&jetAK8puppi_phi_5         ,"jetAK8puppi_phi_5/D"         );
  outTree_->Branch("jetAK8puppi_tau1_5"          ,&jetAK8puppi_tau1_5         ,"jetAK8puppi_tau1_5/D"         );
  outTree_->Branch("jetAK8puppi_tau2_5"          ,&jetAK8puppi_tau2_5         ,"jetAK8puppi_tau2_5/D"         );
  outTree_->Branch("jetAK8puppi_tau3_5"          ,&jetAK8puppi_tau3_5         ,"jetAK8puppi_tau3_5/D"         );
  outTree_->Branch("jetAK8puppi_tau21_5"          ,&jetAK8puppi_tau21_5         ,"jetAK8puppi_tau21_5/D"    );
  outTree_->Branch("jetAK8puppi_tau4_5"          ,&jetAK8puppi_tau4_5         ,"jetAK8puppi_tau4_5/D"       );
  outTree_->Branch("jetAK8puppi_tau42_5"          ,&jetAK8puppi_tau42_5         ,"jetAK8puppi_tau42_5/D"  );
  outTree_->Branch("jetAK8puppi_sd_5"          ,&jetAK8puppi_sd_5         ,"jetAK8puppi_sd_5/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_5"          ,&jetAK8puppi_sdJEC_5         ,"jetAK8puppi_sdJEC_5/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_5"          ,&jetAK8puppi_sdcorr_5         ,"jetAK8puppi_sdcorr_5/D"  );


  outTree_->Branch("jetAK8puppi_ptJEC_6"          ,&jetAK8puppi_ptJEC_6         ,"jetAK8puppi_ptJEC_6/D"    );
  outTree_->Branch("jetAK8puppi_eta_6"          ,&jetAK8puppi_eta_6         ,"jetAK8puppi_eta_6/D"         );
  outTree_->Branch("jetAK8puppi_phi_6"          ,&jetAK8puppi_phi_6         ,"jetAK8puppi_phi_6/D"         );
  outTree_->Branch("jetAK8puppi_tau1_6"          ,&jetAK8puppi_tau1_6         ,"jetAK8puppi_tau1_6/D"         );
  outTree_->Branch("jetAK8puppi_tau2_6"          ,&jetAK8puppi_tau2_6         ,"jetAK8puppi_tau2_6/D"         );
  outTree_->Branch("jetAK8puppi_tau3_6"          ,&jetAK8puppi_tau3_6         ,"jetAK8puppi_tau3_6/D"         );
  outTree_->Branch("jetAK8puppi_tau21_6"          ,&jetAK8puppi_tau21_6         ,"jetAK8puppi_tau21_6/D"    );
  outTree_->Branch("jetAK8puppi_tau4_6"          ,&jetAK8puppi_tau4_6         ,"jetAK8puppi_tau4_6/D"       );
  outTree_->Branch("jetAK8puppi_tau42_6"          ,&jetAK8puppi_tau42_6         ,"jetAK8puppi_tau42_6/D"  );
  outTree_->Branch("jetAK8puppi_sd_6"          ,&jetAK8puppi_sd_6         ,"jetAK8puppi_sd_6/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_6"          ,&jetAK8puppi_sdJEC_6         ,"jetAK8puppi_sdJEC_6/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_6"          ,&jetAK8puppi_sdcorr_6         ,"jetAK8puppi_sdcorr_6/D"  );


  outTree_->Branch("jetAK8puppi_ptJEC_7"          ,&jetAK8puppi_ptJEC_7         ,"jetAK8puppi_ptJEC_7/D"    );
  outTree_->Branch("jetAK8puppi_eta_7"          ,&jetAK8puppi_eta_7         ,"jetAK8puppi_eta_7/D"         );
  outTree_->Branch("jetAK8puppi_phi_7"          ,&jetAK8puppi_phi_7         ,"jetAK8puppi_phi_7/D"         );
  outTree_->Branch("jetAK8puppi_tau1_7"          ,&jetAK8puppi_tau1_7         ,"jetAK8puppi_tau1_7/D"         );
  outTree_->Branch("jetAK8puppi_tau2_7"          ,&jetAK8puppi_tau2_7         ,"jetAK8puppi_tau2_7/D"         );
  outTree_->Branch("jetAK8puppi_tau3_7"          ,&jetAK8puppi_tau3_7         ,"jetAK8puppi_tau3_7/D"         );
  outTree_->Branch("jetAK8puppi_tau21_7"          ,&jetAK8puppi_tau21_7         ,"jetAK8puppi_tau21_7/D"    );
  outTree_->Branch("jetAK8puppi_tau4_7"          ,&jetAK8puppi_tau4_7         ,"jetAK8puppi_tau4_7/D"       );
  outTree_->Branch("jetAK8puppi_tau42_7"          ,&jetAK8puppi_tau42_7         ,"jetAK8puppi_tau42_7/D"  );
  outTree_->Branch("jetAK8puppi_sd_7"          ,&jetAK8puppi_sd_7         ,"jetAK8puppi_sd_7/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_7"          ,&jetAK8puppi_sdJEC_7         ,"jetAK8puppi_sdJEC_7/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_7"          ,&jetAK8puppi_sdcorr_7         ,"jetAK8puppi_sdcorr_7/D"  );


  outTree_->Branch("jetAK8puppi_ptJEC_8"          ,&jetAK8puppi_ptJEC_8         ,"jetAK8puppi_ptJEC_8/D"    );
  outTree_->Branch("jetAK8puppi_eta_8"          ,&jetAK8puppi_eta_8         ,"jetAK8puppi_eta_8/D"         );
  outTree_->Branch("jetAK8puppi_phi_8"          ,&jetAK8puppi_phi_8         ,"jetAK8puppi_phi_8/D"         );
  outTree_->Branch("jetAK8puppi_tau1_8"          ,&jetAK8puppi_tau1_8         ,"jetAK8puppi_tau1_8/D"         );
  outTree_->Branch("jetAK8puppi_tau2_8"          ,&jetAK8puppi_tau2_8         ,"jetAK8puppi_tau2_8/D"         );
  outTree_->Branch("jetAK8puppi_tau3_8"          ,&jetAK8puppi_tau3_8         ,"jetAK8puppi_tau3_8/D"         );
  outTree_->Branch("jetAK8puppi_tau21_8"          ,&jetAK8puppi_tau21_8         ,"jetAK8puppi_tau21_8/D"    );
  outTree_->Branch("jetAK8puppi_tau4_8"          ,&jetAK8puppi_tau4_8         ,"jetAK8puppi_tau4_8/D"       );
  outTree_->Branch("jetAK8puppi_tau42_8"          ,&jetAK8puppi_tau42_8         ,"jetAK8puppi_tau42_8/D"  );
  outTree_->Branch("jetAK8puppi_sd_8"          ,&jetAK8puppi_sd_8         ,"jetAK8puppi_sd_8/D"         );
  outTree_->Branch("jetAK8puppi_sdJEC_8"          ,&jetAK8puppi_sdJEC_8         ,"jetAK8puppi_sdJEC_8/D"   );
  outTree_->Branch("jetAK8puppi_sdcorr_8"          ,&jetAK8puppi_sdcorr_8         ,"jetAK8puppi_sdcorr_8/D"  );
*/
  /// Other quantities
  outTree_->Branch("theWeight", &theWeight, "theWeight/D");  
  outTreew_->Branch("theWeight", &theWeight, "theWeight/D");
  outTree_->Branch("nump", &nump, "nump/D");  
  outTree_->Branch("numm", &numm, "numm/D");  
  outTree_->Branch("npT"           ,&npT         ,"npT/D"          );
  outTree_->Branch("npIT"           ,&npIT         ,"npIT/D"          );
  outTree_->Branch("nBX"           ,&nBX         ,"nBX/I"          );
  outTree_->Branch("triggerWeight"   ,&triggerWeight  ,"triggerWeight/D"  );
  outTree_->Branch("lumiWeight"      ,&lumiWeight     ,"lumiWeight/D"     );
  outTree_->Branch("pileupWeight"    ,&pileupWeight   ,"pileupWeight/D"   );
/*
  outTree_->Branch("delPhijetmet"    ,&delPhijetmet   ,"delPhijetmet/D"   );
  outTree_->Branch("delPhijetmet_2"    ,&delPhijetmet_2   ,"delPhijetmet_2/D"   );
  outTree_->Branch("delPhijetmet_3"    ,&delPhijetmet_3   ,"delPhijetmet_3/D"   );
  outTree_->Branch("delPhijetmet_4"    ,&delPhijetmet_4   ,"delPhijetmet_4/D"   );
  outTree_->Branch("delPhijetmet_5"    ,&delPhijetmet_5   ,"delPhijetmet_5/D"   );
  outTree_->Branch("delPhijetmet_6"    ,&delPhijetmet_6   ,"delPhijetmet_6/D"   );
  outTree_->Branch("delPhijetmet_7"    ,&delPhijetmet_7   ,"delPhijetmet_7/D"   );
  outTree_->Branch("delPhijetmet_8"    ,&delPhijetmet_8   ,"delPhijetmet_8/D"   );
*/
//after JEC varible
  outTree_->Branch("met"             ,&met            ,"met/D"            );
  outTree_->Branch("metPhi"          ,&metPhi         ,"metPhi/D"         );
  outTree_->Branch("METraw_et",&METraw_et,"METraw_et/D");
  outTree_->Branch("METraw_phi",&METraw_phi,"METraw_phi/D");
  outTree_->Branch("METraw_sumEt",&METraw_sumEt,"METraw_sumEt/D");
  outTree_->Branch("MET_et",&MET_et,"MET_et/D");
  outTree_->Branch("MET_phi",&MET_phi,"MET_phi/D");
  outTree_->Branch("MET_sumEt",&MET_sumEt,"MET_sumEt/D");

  outTree_->Branch("jetAK8_pt",&jetAK8_pt,"jetAK8_pt/D");
  outTree_->Branch("jetAK8_mass",&jetAK8_mass,"jetAK8_mass/D");
  outTree_->Branch("jetAK8_jec",&jetAK8_jec,"jetAK8_jec/D");
  outTree_->Branch("jetAK8_pt1",&jetAK8_pt1,"jetAK8_pt1[3]/D");
  outTree_->Branch("jetAK8_eta1",&jetAK8_eta1,"jetAK8_eta1[3]/D");
  outTree_->Branch("jetAK8_mass1",&jetAK8_mass1,"jetAK8_mass1[3]/D");
  outTree_->Branch("jetAK8_SF_mass1",&jetAK8_SF_mass1,"jetAK8_SF_mass1[3]/D");
  outTree_->Branch("jetAK8_SF_mass2",&jetAK8_SF_mass2,"jetAK8_SF_mass2[3]/D");
  outTree_->Branch("jetAK8_jec1",&jetAK8_jec1,"jetAK8_jec1[3]/D");
  outTree_->Branch("jetAK8_eta",&jetAK8_eta,"jetAK8_eta/D");
  outTree_->Branch("jetAK8_phi",&jetAK8_phi,"jetAK8_phi/D");

  outTree_->Branch("candMasspuppiJEC",&candMasspuppiJEC,"candMasspuppiJEC/D");
  outTree_->Branch("massww",&massww,"massww[3]/D");

 
  ///HLT bits
  outTree_->Branch("HLT_Mu1"   ,&HLT_Mu1  ,"HLT_Mu1/I"  );
  outTree_->Branch("HLT_Mu2"   ,&HLT_Mu2  ,"HLT_Mu2/I"  );
  outTree_->Branch("HLT_Mu3"   ,&HLT_Mu3  ,"HLT_Mu3/I"  );
  outTree_->Branch("HLT_Mu4"   ,&HLT_Mu4  ,"HLT_Mu4/I"  );
  outTree_->Branch("HLT_Mu5"   ,&HLT_Mu5  ,"HLT_Mu5/I"  );
  outTree_->Branch("HLT_Mu6"   ,&HLT_Mu6  ,"HLT_Mu6/I"  );
  outTree_->Branch("HLT_Mu7"   ,&HLT_Mu7  ,"HLT_Mu7/I"  );
  outTree_->Branch("HLT_Mu8"   ,&HLT_Mu8  ,"HLT_Mu8/I"  );
  outTree_->Branch("HLT_Mu9"   ,&HLT_Mu9  ,"HLT_Mu9/I"  );
  outTree_->Branch("HLT_Mu10"   ,&HLT_Mu10  ,"HLT_Mu10/I"  );
  outTree_->Branch("HLT_Mu11"   ,&HLT_Mu11  ,"HLT_Mu11/I"  );
  outTree_->Branch("HLT_Mu12"   ,&HLT_Mu12  ,"HLT_Mu12/I"  );
  outTree_->Branch("HLT_Mu13"   ,&HLT_Mu13  ,"HLT_Mu13/I"  );
  outTree_->Branch("HLT_Mu14"   ,&HLT_Mu14  ,"HLT_Mu14/I"  );
  outTree_->Branch("HLT_Mu15"   ,&HLT_Mu15  ,"HLT_Mu15/I"  );
  outTree_->Branch("HLT_Mu16"   ,&HLT_Mu16  ,"HLT_Mu16/I"  );
// filter
  outTree_->Branch("passFilter_HBHE"                 ,&passFilter_HBHE_                ,"passFilter_HBHE_/O");
  outTree_->Branch("passFilter_HBHEIso"                 ,&passFilter_HBHEIso_                ,"passFilter_HBHEIso_/O");
  outTree_->Branch("passFilter_GlobalHalo"              ,&passFilter_GlobalHalo_             ,"passFilter_GlobalHalo_/O");
  outTree_->Branch("passFilter_ECALDeadCell"         ,&passFilter_ECALDeadCell_        ,"passFilter_ECALDeadCell_/O");
  outTree_->Branch("passFilter_GoodVtx"              ,&passFilter_GoodVtx_             ,"passFilter_GoodVtx_/O");
  outTree_->Branch("passFilter_EEBadSc"              ,&passFilter_EEBadSc_             ,"passFilter_EEBadSc_/O");
  outTree_->Branch("passFilter_badMuon"                 ,&passFilter_badMuon_                ,"passFilter_badMuon_/O");
  outTree_->Branch("passFilter_badChargedHadron"                 ,&passFilter_badChargedHadron_                ,"passFilter_badChargedHadron_/O");

  /// AK4 Jets Info
  outTree_->Branch("ak4jet_hf"        , ak4jet_hf       ,"ak4jet_hf[8]/I"       );
  outTree_->Branch("ak4jet_pf"        , ak4jet_pf       ,"ak4jet_pf[8]/I"       );
  outTree_->Branch("ak4jet_pt"        , ak4jet_pt       ,"ak4jet_pt[8]/D"       );
  outTree_->Branch("ak4jet_pt_uncorr"        , ak4jet_pt_uncorr       ,"ak4jet_pt_uncorr[8]/D"       );
  outTree_->Branch("ak4jet_eta"        , ak4jet_eta       ,"ak4jet_eta[8]/D"       );
  outTree_->Branch("ak4jet_phi"        , ak4jet_phi       ,"ak4jet_phi[8]/D"       );
  outTree_->Branch("ak4jet_e"        , ak4jet_e       ,"ak4jet_e[8]/D"       );
  outTree_->Branch("ak4jet_dr"        , ak4jet_dr       ,"ak4jet_dr[8]/D"       );
  outTree_->Branch("ak4jet_csv"        , ak4jet_csv       ,"ak4jet_csv[8]/D"       );
  outTree_->Branch("ak4jet_icsv"        , ak4jet_icsv       ,"ak4jet_icsv[8]/D"       );
  outTree_->Branch("ak4jet_IDLoose"        , ak4jet_IDLoose       ,"ak4jet_IDLoose[8]/D"       );
  outTree_->Branch("ak4jet_IDTight"        , ak4jet_IDTight       ,"ak4jet_IDTight[8]/D"       );
   outTree_->Branch("ak4jet_deepcsvudsg"        , ak4jet_deepcsvudsg       ,"ak4jet_deepcsvudsg[8]/D"       );
    outTree_->Branch("ak4jet_deepcsvb"        , ak4jet_deepcsvb       ,"ak4jet_deepcsvb[8]/D"       );
    outTree_->Branch("ak4jet_deepcsvc"        , ak4jet_deepcsvc       ,"ak4jet_deepcsvc[8]/D"       );
    outTree_->Branch("ak4jet_deepcsvbb"        , ak4jet_deepcsvbb       ,"ak4jet_deepcsvbb[8]/D"       );
    outTree_->Branch("ak4jet_deepcsvcc"        , ak4jet_deepcsvcc       ,"ak4jet_deepcsvcc[8]/D"       );

  /// Gen Level quantities
  /*
  outTree_->Branch("ptgenwl"        ,&ptgenwl       ,"ptgenwl[5]/D"       );
  outTree_->Branch("etagenwl"        ,&etagenwl       ,"etagenwl[5]/D"       );
  outTree_->Branch("phigenwl"        ,&phigenwl       ,"phigenwl[5]/D"       );
  outTree_->Branch("massgenwl"        ,&massgenwl       ,"massgenwl[5]/D"       );
  outTree_->Branch("taggenwl"        ,&taggenwl       ,"taggenwl[5]/D"       );
  outTree_->Branch("taggenwmother"   ,taggenwmother   ,"taggenwmother[5]/D"   );
  outTree_->Branch("genw_q1_pt"           ,genw_q1_pt         ,"genw_q1_pt[5]/D"          );
  outTree_->Branch("genw_q1_phi"           ,genw_q1_phi         ,"genw_q1_phi[5]/D"          );
  outTree_->Branch("genw_q1_eta"           ,genw_q1_eta         ,"genw_q1_eta[5]/D"          );
  outTree_->Branch("genw_q1_e"           ,genw_q1_e         ,"genw_q1_e[5]/D"          );
  outTree_->Branch("genw_q1_pdg"           ,genw_q1_pdg         ,"genw_q1_pdg[5]/D"          );
  outTree_->Branch("genw_q2_pt"           ,genw_q2_pt         ,"genw_q2_pt[5]/D"          );
  outTree_->Branch("genw_q2_phi"           ,genw_q2_phi         ,"genw_q2_phi[5]/D"          );
  outTree_->Branch("genw_q2_eta"           ,genw_q2_eta         ,"genw_q2_eta[5]/D"          );
  outTree_->Branch("genw_q2_e"           ,genw_q2_e         ,"genw_q2_e[5]/D"          );
  outTree_->Branch("genw_q2_pdg"           ,genw_q2_pdg         ,"genw_q2_pdg[5]/D"          );
  outTree_->Branch("ptgenzl"        ,&ptgenzl       ,"ptgenzl[5]/D"       );
  outTree_->Branch("etagenzl"        ,&etagenzl       ,"etagenzl[5]/D"       );
  outTree_->Branch("phigenzl"        ,&phigenzl       ,"phigenzl[5]/D"       );
  outTree_->Branch("massgenzl"        ,&massgenzl       ,"massgenzl[5]/D"       );
  outTree_->Branch("taggenzl"        ,&taggenzl       ,"taggenzl[5]/D"       );
  outTree_->Branch("ptgengl"        ,&ptgengl       ,"ptgengl[15]/D"       );
  outTree_->Branch("etagengl"        ,&etagengl       ,"etagengl[15]/D"       );
  outTree_->Branch("phigengl"        ,&phigengl       ,"phigengl[15]/D"       );
  outTree_->Branch("egengl"        ,&egengl       ,"egengl[15]/D"       );
  outTree_->Branch("ptgenq1l"        ,&ptgenq1l       ,"ptgenq1l[5]/D"       );
  outTree_->Branch("etagenq1l"        ,&etagenq1l       ,"etagenq1l[5]/D"       );
  outTree_->Branch("phigenq1l"        ,&phigenq1l       ,"phigenq1l[5]/D"       );
  outTree_->Branch("egenq1l"        ,&egenq1l       ,"egenq1l[5]/D"       );
  outTree_->Branch("ptgenq2l"        ,&ptgenq2l       ,"ptgenq2l[5]/D"       );
  outTree_->Branch("etagenq2l"        ,&etagenq2l       ,"etagenq2l[5]/D"       );
  outTree_->Branch("phigenq2l"        ,&phigenq2l       ,"phigenq2l[5]/D"       );
  outTree_->Branch("egenq2l"        ,&egenq2l       ,"egenq2l[5]/D"       );
  outTree_->Branch("ptgenq3l"        ,&ptgenq3l       ,"ptgenq3l[5]/D"       );
  outTree_->Branch("etagenq3l"        ,&etagenq3l       ,"etagenq3l[5]/D"       );
  outTree_->Branch("phigenq3l"        ,&phigenq3l       ,"phigenq3l[5]/D"       );
  outTree_->Branch("egenq3l"        ,&egenq3l       ,"egenq3l[5]/D"       );
  outTree_->Branch("ptgenq4l"        ,&ptgenq4l       ,"ptgenq4l[5]/D"       );
  outTree_->Branch("etagenq4l"        ,&etagenq4l       ,"etagenq4l[5]/D"       );
  outTree_->Branch("phigenq4l"        ,&phigenq4l       ,"phigenq4l[5]/D"       );
  outTree_->Branch("egenq4l"        ,&egenq4l       ,"egenq4l[5]/D"       );
  outTree_->Branch("ptgenq5l"        ,&ptgenq5l       ,"ptgenq5l[5]/D"       );
  outTree_->Branch("etagenq5l"        ,&etagenq5l       ,"etagenq5l[5]/D"       );
  outTree_->Branch("phigenq5l"        ,&phigenq5l       ,"phigenq5l[5]/D"       );
  outTree_->Branch("egenq5l"        ,&egenq5l       ,"egenq5l[5]/D"       );
  outTree_->Branch("ptgenwf"        ,&ptgenwf       ,"ptgenwf[5]/D"       );
  outTree_->Branch("etagenwf"        ,&etagenwf       ,"etagenwf[5]/D"       );
  outTree_->Branch("phigenwf"        ,&phigenwf       ,"phigenwf[5]/D"       );
  outTree_->Branch("massgenwf"        ,&massgenwf       ,"massgenwf[5]/D"       );
  outTree_->Branch("ptgenzf"        ,&ptgenzf       ,"ptgenzf[5]/D"       );
  outTree_->Branch("etagenzf"        ,&etagenzf       ,"etagenzf[5]/D"       );
  outTree_->Branch("phigenzf"        ,&phigenzf       ,"phigenzf[5]/D"       );
  outTree_->Branch("massgenzf"        ,&massgenzf       ,"massgenzf[5]/D"       );
  outTree_->Branch("ptgengf"        ,&ptgengf       ,"ptgengf[5]/D"       );
  outTree_->Branch("etagengf"        ,&etagengf       ,"etagengf[5]/D"       );
  outTree_->Branch("phigengf"        ,&phigengf       ,"phigengf[5]/D"       );
  outTree_->Branch("egengf"        ,&egengf       ,"egengf[5]/D"       );
  outTree_->Branch("ptgenq1f"        ,&ptgenq1f       ,"ptgenq1f[5]/D"       );
  outTree_->Branch("etagenq1f"        ,&etagenq1f       ,"etagenq1f[5]/D"       );
  outTree_->Branch("phigenq1f"        ,&phigenq1f       ,"phigenq1f[5]/D"       );
  outTree_->Branch("egenq1f"        ,&egenq1f       ,"egenq1f[5]/D"       );
  outTree_->Branch("ptgenq2f"        ,&ptgenq2f       ,"ptgenq2f[5]/D"       );
  outTree_->Branch("etagenq2f"        ,&etagenq2f       ,"etagenq2f[5]/D"       );
  outTree_->Branch("phigenq2f"        ,&phigenq2f       ,"phigenq2f[5]/D"       );
  outTree_->Branch("egenq2f"        ,&egenq2f       ,"egenq2f[5]/D"       );
  outTree_->Branch("ptgenq3f"        ,&ptgenq3f       ,"ptgenq3f[5]/D"       );
  outTree_->Branch("etagenq3f"        ,&etagenq3f       ,"etagenq3f[5]/D"       );
  outTree_->Branch("phigenq3f"        ,&phigenq3f       ,"phigenq3f[5]/D"       );
  outTree_->Branch("egenq3f"        ,&egenq3f       ,"egenq3f[5]/D"       );
  outTree_->Branch("ptgenq4f"        ,&ptgenq4f       ,"ptgenq4f[5]/D"       );
  outTree_->Branch("etagenq4f"        ,&etagenq4f       ,"etagenq4f[5]/D"       );
  outTree_->Branch("phigenq4f"        ,&phigenq4f       ,"phigenq4f[5]/D"       );
  outTree_->Branch("egenq4f"        ,&egenq4f       ,"egenq4f[5]/D"       );
  outTree_->Branch("ptgenq5f"        ,&ptgenq5f       ,"ptgenq5f[5]/D"       );
  outTree_->Branch("etagenq5f"        ,&etagenq5f       ,"etagenq5f[5]/D"       );
  outTree_->Branch("phigenq5f"        ,&phigenq5f       ,"phigenq5f[5]/D"       );
  outTree_->Branch("egenq5f"        ,&egenq5f       ,"egenq5f[5]/D"       );
  */
    if (RunOnMC_){
        outTree_->Branch("gent_b_pt"           ,&gent_b_pt         ,"gent_b_pt/D"          );
        outTree_->Branch("gent_b_eta"           ,&gent_b_eta         ,"gent_b_eta/D"          );
        outTree_->Branch("gent_b_phi"           ,&gent_b_phi         ,"gent_b_phi/D"          );
        outTree_->Branch("gent_b_mass"           ,&gent_b_mass         ,"gent_b_mass/D"          );
        outTree_->Branch("genantit_b_pt"           ,&genantit_b_pt         ,"genantit_b_pt/D"          );
        outTree_->Branch("genantit_b_eta"           ,&genantit_b_eta         ,"genantit_b_eta/D"          );
        outTree_->Branch("genantit_b_phi"           ,&genantit_b_phi         ,"genantit_b_phi/D"          );
        outTree_->Branch("genantit_b_mass"           ,&genantit_b_mass         ,"genantit_b_mass/D"          );
        outTree_->Branch("gent_w_pt"           ,&gent_w_pt         ,"gent_w_pt/D"          );
        outTree_->Branch("gent_w_eta"           ,&gent_w_eta         ,"gent_w_eta/D"          );
        outTree_->Branch("gent_w_phi"           ,&gent_w_phi         ,"gent_w_phi/D"          );
        outTree_->Branch("gent_w_mass"           ,&gent_w_mass         ,"gent_w_mass/D"          );
        outTree_->Branch("genantit_w_pt"           ,&genantit_w_pt         ,"genantit_w_pt/D"          );
        outTree_->Branch("genantit_w_eta"           ,&genantit_w_eta         ,"genantit_w_eta/D"          );
        outTree_->Branch("genantit_w_phi"           ,&genantit_w_phi         ,"genantit_w_phi/D"          );
        outTree_->Branch("genantit_w_mass"           ,&genantit_w_mass         ,"genantit_w_mass/D"          );
        outTree_->Branch("gent_w_tag"           ,&gent_w_tag         ,"gent_w_tag/D"          );
        outTree_->Branch("gent_w_q1_pt"           ,&gent_w_q1_pt         ,"gent_w_q1_pt/D"          );
        outTree_->Branch("gent_w_q1_eta"           ,&gent_w_q1_eta         ,"gent_w_q1_eta/D"          );
        outTree_->Branch("gent_w_q1_phi"           ,&gent_w_q1_phi         ,"gent_w_q1_phi/D"          );
        outTree_->Branch("gent_w_q1_e"           ,&gent_w_q1_e         ,"gent_w_q1_e/D"          );
        outTree_->Branch("gent_w_q1_pdg"           ,&gent_w_q1_pdg         ,"gent_w_q1_pdg/D"          );
        outTree_->Branch("gent_w_q2_pt"           ,&gent_w_q2_pt         ,"gent_w_q2_pt/D"          );
        outTree_->Branch("gent_w_q2_eta"           ,&gent_w_q2_eta         ,"gent_w_q2_eta/D"          );
        outTree_->Branch("gent_w_q2_phi"           ,&gent_w_q2_phi         ,"gent_w_q2_phi/D"          );
        outTree_->Branch("gent_w_q2_e"           ,&gent_w_q2_e         ,"gent_w_q2_e/D"          );
        outTree_->Branch("gent_w_q2_pdg"           ,&gent_w_q2_pdg         ,"gent_w_q2_pdg/D"          );
        outTree_->Branch("genantit_w_tag"           ,&genantit_w_tag         ,"genantit_w_tag/D"          );
        outTree_->Branch("genantit_w_q1_pt"           ,&genantit_w_q1_pt         ,"genantit_w_q1_pt/D"          );
        outTree_->Branch("genantit_w_q1_eta"           ,&genantit_w_q1_eta         ,"genantit_w_q1_eta/D"          );
        outTree_->Branch("genantit_w_q1_phi"           ,&genantit_w_q1_phi         ,"genantit_w_q1_phi/D"          );
        outTree_->Branch("genantit_w_q1_e"           ,&genantit_w_q1_e         ,"genantit_w_q1_e/D"          );
        outTree_->Branch("genantit_w_q1_pdg"           ,&genantit_w_q1_pdg         ,"genantit_w_q1_pdg/D"          );
        outTree_->Branch("genantit_w_q2_pt"           ,&genantit_w_q2_pt         ,"genantit_w_q2_pt/D"          );
        outTree_->Branch("genantit_w_q2_eta"           ,&genantit_w_q2_eta         ,"genantit_w_q2_eta/D"          );
        outTree_->Branch("genantit_w_q2_phi"           ,&genantit_w_q2_phi         ,"genantit_w_q2_phi/D"          );
        outTree_->Branch("genantit_w_q2_e"           ,&genantit_w_q2_e         ,"gent_w_q2_e/D"          );
        outTree_->Branch("genantit_w_q2_pdg"           ,&genantit_w_q2_pdg         ,"genantit_w_q2_pdg/D"          );

        outTree_->Branch("ptgenwl"           ,ptgenwl         ,"ptgenwl[5]/D"          );
        outTree_->Branch("etagenwl"           ,etagenwl         ,"etagenwl[5]/D"          );
        outTree_->Branch("phigenwl"           ,phigenwl       ,"phigenwl[5]/D"          );
        outTree_->Branch("massgenwl"           ,massgenwl         ,"massgenwl[5]/D"          );
        outTree_->Branch("taggenwl"           ,taggenwl         ,"taggenwl[5]/D"          );
        outTree_->Branch("taggenwmother"           ,taggenwmother         ,"taggenwmother[5]/D"          );
        outTree_->Branch("genw_q1_pt"           ,genw_q1_pt         ,"genw_q1_pt[5]/D"          );
        outTree_->Branch("genw_q1_phi"           ,genw_q1_phi         ,"genw_q1_phi[5]/D"          );
        outTree_->Branch("genw_q1_eta"           ,genw_q1_eta         ,"genw_q1_eta[5]/D"          );
        outTree_->Branch("genw_q1_e"           ,genw_q1_e         ,"genw_q1_e[5]/D"          );
        outTree_->Branch("genw_q1_pdg"           ,genw_q1_pdg         ,"genw_q1_pdg[5]/D"          );
        outTree_->Branch("genw_q2_pt"           ,genw_q2_pt         ,"genw_q2_pt[5]/D"          );
        outTree_->Branch("genw_q2_phi"           ,genw_q2_phi         ,"genw_q2_phi[5]/D"          );
        outTree_->Branch("genw_q2_eta"           ,genw_q2_eta         ,"genw_q2_eta[5]/D"          );
        outTree_->Branch("genw_q2_e"           ,genw_q2_e         ,"genw_q2_e[5]/D"          );
        outTree_->Branch("genw_q2_pdg"           ,genw_q2_pdg         ,"genw_q2_pdg[5]/D"          );

        outTree_->Branch("ptgenzl"           ,ptgenzl         ,"ptgenzl[5]/D"          );
        outTree_->Branch("etagenzl"           ,etagenzl         ,"etagenzl[5]/D"          );
        outTree_->Branch("phigenzl"           ,phigenzl       ,"phigenzl[5]/D"          );
        outTree_->Branch("massgenzl"           ,massgenzl         ,"massgenzl[5]/D"          );
        outTree_->Branch("taggenzl"           ,taggenzl         ,"taggenzl[5]/D"          );
        outTree_->Branch("ptgengl"           ,ptgengl         ,"ptgengl[15]/D"          );
        outTree_->Branch("etagengl"           ,etagengl         ,"etagengl[15]/D"          );
        outTree_->Branch("phigengl"           ,phigengl       ,"phigengl[15]/D"          );
        outTree_->Branch("egengl"           ,egengl         ,"egengl[15]/D"          );
        outTree_->Branch("ptgenwf"           ,ptgenwf         ,"ptgenwf[5]/D"          );
        outTree_->Branch("etagenwf"           ,etagenwf         ,"etagenwf[5]/D"          );
        outTree_->Branch("phigenwf"           ,phigenwf       ,"phigenwf[5]/D"          );
        outTree_->Branch("massgenwf"           ,massgenwf         ,"massgenwf[5]/D"          );
        outTree_->Branch("ptgenzf"           ,ptgenzf         ,"ptgenzf[5]/D"          );
        outTree_->Branch("etagenzf"           ,etagenzf         ,"etagenzf[5]/D"          );
        outTree_->Branch("phigenzf"           ,phigenzf       ,"phigenzf[5]/D"          );
        outTree_->Branch("massgenzf"           ,massgenzf         ,"massgenzf[5]/D"          );
        outTree_->Branch("ptgengf"           ,ptgengf         ,"ptgengf[15]/D"          );
        outTree_->Branch("etagengf"           ,etagengf         ,"etagengf[15]/D"          );
        outTree_->Branch("phigengf"           ,phigengf       ,"phigengf[15]/D"          );
        outTree_->Branch("egengf"           ,egengf         ,"egengf[15]/D"          );

        outTree_->Branch("gent_b_pt"           ,&gent_b_pt         ,"gent_b_pt/D"          );
        outTree_->Branch("gent_b_eta"           ,&gent_b_eta         ,"gent_b_eta/D"          );
        outTree_->Branch("gent_b_phi"           ,&gent_b_phi         ,"gent_b_phi/D"          );
        outTree_->Branch("gent_b_mass"           ,&gent_b_mass         ,"gent_b_mass/D"          );
        outTree_->Branch("genantit_b_pt"           ,&genantit_b_pt         ,"genantit_b_pt/D"          );
        outTree_->Branch("genantit_b_eta"           ,&genantit_b_eta         ,"genantit_b_eta/D"          );
        outTree_->Branch("genantit_b_phi"           ,&genantit_b_phi         ,"genantit_b_phi/D"          );
        outTree_->Branch("genantit_b_mass"           ,&genantit_b_mass         ,"genantit_b_mass/D"          );
        outTree_->Branch("gent_w_pt"           ,&gent_w_pt         ,"gent_w_pt/D"          );
        outTree_->Branch("gent_w_eta"           ,&gent_w_eta         ,"gent_w_eta/D"          );
        outTree_->Branch("gent_w_phi"           ,&gent_w_phi         ,"gent_w_phi/D"          );
        outTree_->Branch("gent_w_mass"           ,&gent_w_mass         ,"gent_w_mass/D"          );
        outTree_->Branch("genantit_w_pt"           ,&genantit_w_pt         ,"genantit_w_pt/D"          );
        outTree_->Branch("genantit_w_eta"           ,&genantit_w_eta         ,"genantit_w_eta/D"          );
        outTree_->Branch("genantit_w_phi"           ,&genantit_w_phi         ,"genantit_w_phi/D"          );
        outTree_->Branch("genantit_w_mass"           ,&genantit_w_mass         ,"genantit_w_mass/D"          );
        outTree_->Branch("gent_w_tag"           ,&gent_w_tag         ,"gent_w_tag/D"          );
        outTree_->Branch("gent_w_q1_pt"           ,&gent_w_q1_pt         ,"gent_w_q1_pt/D"          );
        outTree_->Branch("gent_w_q1_eta"           ,&gent_w_q1_eta         ,"gent_w_q1_eta/D"          );
        outTree_->Branch("gent_w_q1_phi"           ,&gent_w_q1_phi         ,"gent_w_q1_phi/D"          );
        outTree_->Branch("gent_w_q1_e"           ,&gent_w_q1_e         ,"gent_w_q1_e/D"          );
        outTree_->Branch("gent_w_q1_pdg"           ,&gent_w_q1_pdg         ,"gent_w_q1_pdg/D"          );
        outTree_->Branch("gent_w_q2_pt"           ,&gent_w_q2_pt         ,"gent_w_q2_pt/D"          );
        outTree_->Branch("gent_w_q2_eta"           ,&gent_w_q2_eta         ,"gent_w_q2_eta/D"          );
        outTree_->Branch("gent_w_q2_phi"           ,&gent_w_q2_phi         ,"gent_w_q2_phi/D"          );
        outTree_->Branch("gent_w_q2_e"           ,&gent_w_q2_e         ,"gent_w_q2_e/D"          );
        outTree_->Branch("gent_w_q2_pdg"           ,&gent_w_q2_pdg         ,"gent_w_q2_pdg/D"          );
        outTree_->Branch("genantit_w_tag"           ,&genantit_w_tag         ,"genantit_w_tag/D"          );
        outTree_->Branch("genantit_w_q1_pt"           ,&genantit_w_q1_pt         ,"genantit_w_q1_pt/D"          );
        outTree_->Branch("genantit_w_q1_eta"           ,&genantit_w_q1_eta         ,"genantit_w_q1_eta/D"          );
        outTree_->Branch("genantit_w_q1_phi"           ,&genantit_w_q1_phi         ,"genantit_w_q1_phi/D"          );
        outTree_->Branch("genantit_w_q1_e"           ,&genantit_w_q1_e         ,"genantit_w_q1_e/D"          );
        outTree_->Branch("genantit_w_q1_pdg"           ,&genantit_w_q1_pdg         ,"genantit_w_q1_pdg/D"          );
        outTree_->Branch("genantit_w_q2_pt"           ,&genantit_w_q2_pt         ,"genantit_w_q2_pt/D"          );
        outTree_->Branch("genantit_w_q2_eta"           ,&genantit_w_q2_eta         ,"genantit_w_q2_eta/D"          );
        outTree_->Branch("genantit_w_q2_phi"           ,&genantit_w_q2_phi         ,"genantit_w_q2_phi/D"          );
        outTree_->Branch("genantit_w_q2_e"           ,&genantit_w_q2_e         ,"gent_w_q2_e/D"          );
        outTree_->Branch("genantit_w_q2_pdg"           ,&genantit_w_q2_pdg         ,"genantit_w_q2_pdg/D"          );

        outTree_->Branch("ptgenq1l"           ,ptgenq1l         ,"ptgenq1l[5]/D"          );
        outTree_->Branch("etagenq1l"           ,etagenq1l         ,"etagenq1l[5]/D"          );
        outTree_->Branch("phigenq1l"           ,phigenq1l       ,"phigenq1l[5]/D"          );
        outTree_->Branch("egenq1l"           ,egenq1l         ,"egenq1l[5]/D"          );
        outTree_->Branch("ptgenq1f"           ,ptgenq1f         ,"ptgenq1f[5]/D"          );
        outTree_->Branch("etagenq1f"           ,etagenq1f         ,"etagenq1f[5]/D"          );
        outTree_->Branch("phigenq1f"           ,phigenq1f       ,"phigenq1f[5]/D"          );
        outTree_->Branch("egenq1f"           ,egenq1f         ,"egenq1f[5]/D"          );
        outTree_->Branch("ptgenq2l"           ,ptgenq2l         ,"ptgenq2l[5]/D"          );
        outTree_->Branch("etagenq2l"           ,etagenq2l         ,"etagenq2l[5]/D"          );
        outTree_->Branch("phigenq2l"           ,phigenq2l       ,"phigenq2l[5]/D"          );
        outTree_->Branch("egenq2l"           ,egenq2l         ,"egenq2l[5]/D"          );
        outTree_->Branch("ptgenq2f"           ,ptgenq2f         ,"ptgenq2f[5]/D"          );
        outTree_->Branch("etagenq2f"           ,etagenq2f         ,"etagenq2f[5]/D"          );
        outTree_->Branch("phigenq2f"           ,phigenq2f       ,"phigenq2f[5]/D"          );
        outTree_->Branch("egenq2f"           ,egenq2f         ,"egenq2f[5]/D"          );
        outTree_->Branch("ptgenq3l"           ,ptgenq3l         ,"ptgenq3l[5]/D"          );
        outTree_->Branch("etagenq3l"           ,etagenq3l         ,"etagenq3l[5]/D"          );
        outTree_->Branch("phigenq3l"           ,phigenq3l       ,"phigenq3l[5]/D"          );
        outTree_->Branch("egenq3l"           ,egenq3l         ,"egenq3l[5]/D"          );
        outTree_->Branch("ptgenq3f"           ,ptgenq3f         ,"ptgenq3f[5]/D"          );
        outTree_->Branch("etagenq3f"           ,etagenq3f         ,"etagenq3f[5]/D"          );
        outTree_->Branch("phigenq3f"           ,phigenq3f       ,"phigenq3f[5]/D"          );
        outTree_->Branch("egenq3f"           ,egenq3f         ,"egenq3f[5]/D"          );
        outTree_->Branch("ptgenq4l"           ,ptgenq4l         ,"ptgenq4l[5]/D"          );
        outTree_->Branch("etagenq4l"           ,etagenq4l         ,"etagenq4l[5]/D"          );
        outTree_->Branch("phigenq4l"           ,phigenq4l       ,"phigenq4l[5]/D"          );
        outTree_->Branch("egenq4l"           ,egenq4l         ,"egenq4l[5]/D"          );
        outTree_->Branch("ptgenq4f"           ,ptgenq4f         ,"ptgenq4f[5]/D"          );
        outTree_->Branch("etagenq4f"           ,etagenq4f         ,"etagenq4f[5]/D"          );
        outTree_->Branch("phigenq4f"           ,phigenq4f       ,"phigenq4f[5]/D"          );
        outTree_->Branch("egenq4f"           ,egenq4f         ,"egenq4f[5]/D"          );
        outTree_->Branch("ptgenq5l"           ,ptgenq5l         ,"ptgenq5l[5]/D"          );
        outTree_->Branch("etagenq5l"           ,etagenq5l         ,"etagenq5l[5]/D"          );
        outTree_->Branch("phigenq5l"           ,phigenq5l       ,"phigenq5l[5]/D"          );
        outTree_->Branch("egenq5l"           ,egenq5l         ,"egenq5l[5]/D"          );
        outTree_->Branch("ptgenq5f"           ,ptgenq5f         ,"ptgenq5f[5]/D"          );
        outTree_->Branch("etagenq5f"           ,etagenq5f         ,"etagenq5f[5]/D"          );
        outTree_->Branch("phigenq5f"           ,phigenq5f       ,"phigenq5f[5]/D"          );
        outTree_->Branch("egenq5f"           ,egenq5f         ,"egenq5f[5]/D"          );
        outTree_->Branch("mothergenq1f"           ,mothergenq1f         ,"mothergenq1f[5]/D"          );
        outTree_->Branch("mothergenq2f"           ,mothergenq2f         ,"mothergenq2f[5]/D"          );
        outTree_->Branch("mothergenq3f"           ,mothergenq3f         ,"mothergenq3f[5]/D"          );
        outTree_->Branch("mothergenq4f"           ,mothergenq4f         ,"mothergenq4f[5]/D"          );
        outTree_->Branch("mothergenq5f"           ,mothergenq5f         ,"mothergenq5f[5]/D"          );

        outTree_->Branch("mothergengf"           ,mothergengf         ,"mothergengf[15]/D"          );
        outTree_->Branch("mmothergengf"           ,mmothergengf         ,"mmothergengf[15]/D"          );

        outTree_->Branch("mmothergenq1f"           ,mmothergenq1f         ,"mmothergenq1f[5]/D"          );
        outTree_->Branch("mmothergenq2f"           ,mmothergenq2f         ,"mmothergenq2f[5]/D"          );
        outTree_->Branch("mmothergenq3f"           ,mmothergenq3f         ,"mmothergenq3f[5]/D"          );
        outTree_->Branch("mmothergenq4f"           ,mmothergenq4f         ,"mmothergenq4f[5]/D"          );
        outTree_->Branch("mmothergenq5f"           ,mmothergenq5f         ,"mmothergenq5f[5]/D"          );
  outTree_->Branch("gen_gra_m"        ,&gen_gra_m       ,"gen_gra_m/D"       );
  outTree_->Branch("gen_gra_pt"        ,&gen_gra_pt       ,"gen_gra_pt/D"       );
  outTree_->Branch("gen_gra_eta"        ,&gen_gra_eta       ,"gen_gra_eta/D"       );
  outTree_->Branch("gen_gra_phi"        ,&gen_gra_phi       ,"gen_gra_phi/D"       );
  outTree_->Branch("gen_rad_m"        ,&gen_rad_m       ,"gen_rad_m/D"       );
  outTree_->Branch("gen_rad_pt"        ,&gen_rad_pt       ,"gen_rad_pt/D"       );
  outTree_->Branch("gen_rad_eta"        ,&gen_rad_eta       ,"gen_rad_eta/D"       );
  outTree_->Branch("gen_rad_phi"        ,&gen_rad_phi       ,"gen_rad_phi/D"       );
/*  outTree_->Branch("gen_ele_pt"        ,&gen_ele_pt       ,"gen_ele_pt/D"       );
  outTree_->Branch("gen_ele_eta"        ,&gen_ele_eta       ,"gen_ele_eta/D"       );
  outTree_->Branch("gen_ele_phi"        ,&gen_ele_phi       ,"gen_ele_phi/D"       );
  outTree_->Branch("gen_ele_e"        ,&gen_ele_e       ,"gen_ele_e/D"       );
  outTree_->Branch("gen_mu_pt"        ,&gen_mu_pt       ,"gen_mu_pt/D"       );
  outTree_->Branch("gen_mu_eta"        ,&gen_mu_eta       ,"gen_mu_eta/D"       );
  outTree_->Branch("gen_mu_phi"        ,&gen_mu_phi       ,"gen_mu_phi/D"       );
  outTree_->Branch("gen_mu_e"        ,&gen_mu_e       ,"gen_mu_e/D"       );
  outTree_->Branch("gen_ele_pt_2"        ,&gen_ele_pt_2       ,"gen_ele_pt_2/D"       );
  outTree_->Branch("gen_ele_eta_2"        ,&gen_ele_eta_2       ,"gen_ele_eta_2/D"       );
  outTree_->Branch("gen_ele_phi_2"        ,&gen_ele_phi_2       ,"gen_ele_phi_2/D"       );
  outTree_->Branch("gen_ele_e_2"        ,&gen_ele_e_2       ,"gen_ele_e_2/D"       );
  outTree_->Branch("gen_mu_pt_2"        ,&gen_mu_pt_2       ,"gen_mu_pt_2/D"       );
  outTree_->Branch("gen_mu_eta_2"        ,&gen_mu_eta_2       ,"gen_mu_eta_2/D"       );
  outTree_->Branch("gen_mu_phi_2"        ,&gen_mu_phi_2       ,"gen_mu_phi_2/D"       );
  outTree_->Branch("gen_mu_e_2"        ,&gen_mu_e_2       ,"gen_mu_e_2/D"       );
  outTree_->Branch("gen_ele_pt_3"        ,&gen_ele_pt_3       ,"gen_ele_pt_3/D"       );
  outTree_->Branch("gen_ele_eta_3"        ,&gen_ele_eta_3       ,"gen_ele_eta_3/D"       );
  outTree_->Branch("gen_ele_phi_3"        ,&gen_ele_phi_3       ,"gen_ele_phi_3/D"       );
  outTree_->Branch("gen_ele_e_3"        ,&gen_ele_e_3       ,"gen_ele_e_3/D"       );
  outTree_->Branch("gen_mu_pt_3"        ,&gen_mu_pt_3       ,"gen_mu_pt_3/D"       );
  outTree_->Branch("gen_mu_eta_3"        ,&gen_mu_eta_3       ,"gen_mu_eta_3/D"       );
  outTree_->Branch("gen_mu_phi_3"        ,&gen_mu_phi_3       ,"gen_mu_phi_3/D"       );
  outTree_->Branch("gen_mu_e_3"        ,&gen_mu_e_3       ,"gen_mu_e_3/D"       );
*/ 
  outTree_->Branch("gen_tau_pt"        ,&gen_tau_pt       ,"gen_tau_pt/D"       );
  outTree_->Branch("gen_tau_eta"        ,&gen_tau_eta       ,"gen_tau_eta/D"       );
  outTree_->Branch("gen_tau_phi"        ,&gen_tau_phi       ,"gen_tau_phi/D"       );
  outTree_->Branch("gen_tau_e"        ,&gen_tau_e       ,"gen_tau_e/D"       );
  outTree_->Branch("gen_tau_pt_2"        ,&gen_tau_pt_2       ,"gen_tau_pt_2/D"       );
  outTree_->Branch("gen_tau_eta_2"        ,&gen_tau_eta_2       ,"gen_tau_eta_2/D"       );
  outTree_->Branch("gen_tau_phi_2"        ,&gen_tau_phi_2       ,"gen_tau_phi_2/D"       );
  outTree_->Branch("gen_tau_e_2"        ,&gen_tau_e_2       ,"gen_tau_e_2/D"       );
  outTree_->Branch("gen_tau_pt_3"        ,&gen_tau_pt_3       ,"gen_tau_pt_3/D"       );
  outTree_->Branch("gen_tau_eta_3"        ,&gen_tau_eta_3       ,"gen_tau_eta_3/D"       );
  outTree_->Branch("gen_tau_phi_3"        ,&gen_tau_phi_3       ,"gen_tau_phi_3/D"       );
  outTree_->Branch("gen_tau_e_3"        ,&gen_tau_e_3       ,"gen_tau_e_3/D"       );
 
  outTree_->Branch("gentop_pt"        ,&gentop_pt       ,"gentop_pt/D"       );
  outTree_->Branch("gentop_eta"        ,&gentop_eta       ,"gentop_eta/D"       );
  outTree_->Branch("gentop_phi"        ,&gentop_phi       ,"gentop_phi/D"       );
  outTree_->Branch("gentop_mass"        ,&gentop_mass       ,"gentop_mass/D"       );
  outTree_->Branch("genantitop_pt"        ,&genantitop_pt       ,"genantitop_pt/D"       );
  outTree_->Branch("genantitop_eta"        ,&genantitop_eta       ,"genantitop_eta/D"       );
  outTree_->Branch("genantitop_phi"        ,&genantitop_phi       ,"genantitop_phi/D"       );
  outTree_->Branch("genantitop_mass"        ,&genantitop_mass       ,"genantitop_mass/D"       );
  //outTree_->Branch("ptGenVlep"        ,&ptGenVlep       ,"ptGenVlep/D"       );
//  outTree_->Branch("etaGenVlep"        ,&etaGenVlep       ,"etaGenVlep/D"       );
//  outTree_->Branch("phiGenVlep"        ,&phiGenVlep       ,"phiGenVlep/D"       );
//  outTree_->Branch("massGenVlep"        ,&massGenVlep       ,"massGenVlep/D"       );
  outTree_->Branch("ptq11"        ,&ptq11       ,"ptq11/D"       );
  outTree_->Branch("etaq11"        ,&etaq11       ,"etaq11/D"       );
  outTree_->Branch("phiq11"        ,&phiq11       ,"phiq11/D"       );
  outTree_->Branch("massq11"        ,&massq11       ,"massq11/D"       );
  outTree_->Branch("ptq21"        ,&ptq21       ,"ptq21/D"       );
  outTree_->Branch("etaq21"        ,&etaq21       ,"etaq21/D"       );
  outTree_->Branch("phiq21"        ,&phiq21       ,"phiq21/D"       );
  outTree_->Branch("massq21"        ,&massq21       ,"massq21/D"       );
  outTree_->Branch("ptq31"        ,&ptq31       ,"ptq31/D"       );
  outTree_->Branch("etaq31"        ,&etaq31       ,"etaq31/D"       );
  outTree_->Branch("phiq31"        ,&phiq31       ,"phiq31/D"       );
  outTree_->Branch("massq31"        ,&massq31       ,"massq31/D"       );

  outTree_->Branch("ptq12"        ,&ptq12       ,"ptq12/D"       );
  outTree_->Branch("etaq12"        ,&etaq12       ,"etaq12/D"       );
  outTree_->Branch("phiq12"        ,&phiq12       ,"phiq12/D"       );
  outTree_->Branch("massq12"        ,&massq12       ,"massq12/D"       );
  outTree_->Branch("ptq22"        ,&ptq22       ,"ptq22/D"       );
  outTree_->Branch("etaq22"        ,&etaq22       ,"etaq22/D"       );
  outTree_->Branch("phiq22"        ,&phiq22       ,"phiq22/D"       );
  outTree_->Branch("massq22"        ,&massq22       ,"massq22/D"       );
  outTree_->Branch("ptq32"        ,&ptq32       ,"ptq32/D"       );
  outTree_->Branch("etaq32"        ,&etaq32       ,"etaq32/D"       );
  outTree_->Branch("phiq32"        ,&phiq32       ,"phiq32/D"       );
  outTree_->Branch("massq32"        ,&massq32       ,"massq32/D"       );

  outTree_->Branch("ptGenVhad"        ,&ptGenVhad       ,"ptGenVhad/D"       );
  outTree_->Branch("etaGenVhad"        ,&etaGenVhad       ,"etaGenVhad/D"       );
  outTree_->Branch("phiGenVhad"        ,&phiGenVhad       ,"phiGenVhad/D"       );
  outTree_->Branch("massGenVhad"        ,&massGenVhad       ,"massGenVhad/D"       );
  outTree_->Branch("ptGenV_2"        ,&ptGenV_2       ,"ptGenV_2/D"       );
  outTree_->Branch("etaGenV_2"        ,&etaGenV_2       ,"etaGenV_2/D"       );
  outTree_->Branch("phiGenV_2"        ,&phiGenV_2       ,"phiGenV_2/D"       );
  outTree_->Branch("massGenV_2"        ,&massGenV_2       ,"massGenV_2/D"       );
  outTree_->Branch("ptGenV_3"        ,&ptGenV_3       ,"ptGenV_3/D"       );
  outTree_->Branch("etaGenV_3"        ,&etaGenV_3       ,"etaGenV_3/D"       );
  outTree_->Branch("phiGenV_3"        ,&phiGenV_3       ,"phiGenV_3/D"       );
  outTree_->Branch("massGenV_3"        ,&massGenV_3       ,"massGenV_3/D"       );

  outTree_->Branch("ptGenVlep"        ,&ptGenVlep       ,"ptGenVlep/D"       );
  outTree_->Branch("etaGenVlep"        ,&etaGenVlep       ,"etaGenVlep/D"       );
  outTree_->Branch("phiGenVlep"        ,&phiGenVlep       ,"phiGenVlep/D"       );
  outTree_->Branch("massGenVlep"        ,&massGenVlep       ,"massGenVlep/D"       );
  outTree_->Branch("ptGenVlep_2"        ,&ptGenVlep_2       ,"ptGenVlep_2/D"       );
  outTree_->Branch("etaGenVlep_2"        ,&etaGenVlep_2       ,"etaGenVlep_2/D"       );
  outTree_->Branch("phiGenVlep_2"        ,&phiGenVlep_2       ,"phiGenVlep_2/D"       );
  outTree_->Branch("massGenVlep_2"        ,&massGenVlep_2       ,"massGenVlep_2/D"       );
  outTree_->Branch("ptGenVlep_3"        ,&ptGenVlep_3       ,"ptGenVlep_3/D"       );
  outTree_->Branch("etaGenVlep_3"        ,&etaGenVlep_3       ,"etaGenVlep_3/D"       );
  outTree_->Branch("phiGenVlep_3"        ,&phiGenVlep_3       ,"phiGenVlep_3/D"       );
  outTree_->Branch("massGenVlep_3"        ,&massGenVlep_3       ,"massGenVlep_3/D"       );
  
  outTree_->Branch("status_1"           ,&status_1         ,"status_1/I"          );
  outTree_->Branch("status_2"           ,&status_2         ,"status_2/I"          );
  outTree_->Branch("status_3"           ,&status_3         ,"status_3/I"          );
  outTree_->Branch("pweight" ,pweight ,"pweight[211]/D" );
  }
}
/*
const reco::Candidate*  EDBRTreeMaker::findLastParticle(const reco::Candidate *particle,int IDpdg){
    int iw=0;
    int pidw=0;
    const reco::Candidate* pw=particle;
    for(int ii=0;particle->daughter(ii)!=NULL;ii++){
        if(abs(particle->daughter(ii)->pdgId())>pidw) {
            iw=ii;
            pidw=abs(particle->daughter(ii)->pdgId());
        }
    }
    if( abs(pidw) == IDpdg ){
        
        pw = particle->daughter(iw);
        return (findLastParticle(pw,IDpdg));
    }
    return pw;
    }
*/

const reco::Candidate*  EDBRTreeMaker::findLastParticle(const reco::Candidate *particle,int IDpdg){
    int iw=0;
    int pidw=0;
    const reco::Candidate* pw=particle;
    for(int ii=0;particle->daughter(ii)!=NULL;ii++){
        if(abs(particle->daughter(ii)->pdgId())== IDpdg) {
            iw=ii;
            pidw=abs(particle->daughter(ii)->pdgId());
        }
    }
    if( abs(pidw) == IDpdg ){
        pw = particle->daughter(iw);
        return (findLastParticle(pw,IDpdg));
    }
    return pw;
}


const reco::Candidate*  EDBRTreeMaker::findFirstParticle(const reco::Candidate *particle,int IDpdg){
    if (particle->mother(0)!=NULL){
        if(abs(particle->mother(0)->pdgId()) == IDpdg )
        return (findFirstParticle(particle->mother(0),IDpdg));
    }
    return particle;
}




EDBRTreeMaker::~EDBRTreeMaker()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


bool
EDBRTreeMaker::looseJetID( const pat::Jet& j ) {
// refer to https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
        if(j.pt()>170.){
	double NHF = j.neutralHadronEnergyFraction();
	double NEMF = j.neutralEmEnergyFraction();
	double CHF = j.chargedHadronEnergyFraction();
	double CEMF = j.chargedEmEnergyFraction();
	int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
	int NumNeutralParticle =j.neutralMultiplicity();
	int CHM = j.chargedMultiplicity(); 
	double eta = j.eta();
          return ((  (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7 ) || (NHF<0.98 && NEMF>0.01 && NumNeutralParticle>2 && abs(eta)>2.7 && abs(eta)<=3.0 ) || (NEMF<0.90 && NumNeutralParticle>10 && abs(eta)>3.0) ) ;
}
else{
return (1);
    }
}

bool
EDBRTreeMaker::tightJetID( const pat::Jet& j ) {
// refer to https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
        if(j.pt()>170.){
        double NHF = j.neutralHadronEnergyFraction();
        double NEMF = j.neutralEmEnergyFraction();
        double CHF = j.chargedHadronEnergyFraction();
        double CEMF = j.chargedEmEnergyFraction();
        int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
        int NumNeutralParticle =j.neutralMultiplicity();
        int CHM = j.chargedMultiplicity();
        double eta = j.eta();
        return ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 ) || (NEMF<0.90 && NumNeutralParticle>10 && abs(eta)>3.0 )  ;
}
else{
return (1);
    }
}

float
EDBRTreeMaker::dEtaInSeed( const pat::Electron*  ele ){
  return ele->superCluster().isNonnull() && ele->superCluster()->seed().isNonnull() ?
    ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta() : std::numeric_limits<float>::max();
}

void EDBRTreeMaker::initJetCorrFactors( void ){
   std::vector<JetCorrectorParameters> vPar;
   for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(), payloadEnd = jecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
   }
  // Make the FactorizedJetCorrector
  jecAK8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNamesGroomed_.begin(), payloadEnd = jecAK8PayloadNamesGroomed_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  // Make the FactorizedJetCorrector
  jecAK8Groomed_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNamesGroomed_.begin(), payloadEnd = jecAK8PayloadNamesGroomed_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  jecAK8GroomedSD_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8puppiPayloadNames_.begin(), payloadEnd = jecAK8puppiPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
 jecAK8puppi_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8puppiPayloadNamesGroomed_.begin(), payloadEnd = jecAK8puppiPayloadNamesGroomed_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
 jecAK8puppiGroomed_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4PayloadNames_.begin(), payloadEnd = jecAK4PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  // Make the FactorizedJetCorrector
  jecAK4_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = offsetCorrLabel_.begin(), payloadEnd = offsetCorrLabel_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  jecOffset_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
}


double EDBRTreeMaker::getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
   double jetCorrFactor = 1.;
   if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
      jecAK4_->setJetEta( rawJetP4.eta() );
      jecAK4_->setJetPt ( rawJetP4.pt() );
      jecAK4_->setJetE  ( rawJetP4.energy() );
      jecAK4_->setJetPhi( rawJetP4.phi()    );
      jecAK4_->setJetA  ( jet.jetArea() );
      jecAK4_->setRho   ( *(rho_.product()) );
      jecAK4_->setNPV   ( nVtx );
      jetCorrFactor = jecAK4_->getCorrection();
   }
   reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
   corrJetP4 *= jetCorrFactor;
   return jetCorrFactor;
}

double EDBRTreeMaker::getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
   double jetCorrFactor = 1.;
   if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
      jecOffset_->setJetEta( rawJetP4.eta()     );
      jecOffset_->setJetPt ( rawJetP4.pt()      );
      jecOffset_->setJetE  ( rawJetP4.energy()  );
      jecOffset_->setJetPhi( rawJetP4.phi()     );
      jecOffset_->setJetA  ( jet.jetArea()      );
      jecOffset_->setRho   ( *(rho_.product())  );
      jecOffset_->setNPV   ( nVtx  );
      jetCorrFactor = jecOffset_->getCorrection();
   }
   reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
   corrJetP4 *= jetCorrFactor;
   return jetCorrFactor;
}

//-------------------------------------------------------------------------------------------------------------------------------------//
//
// member functions
//
void EDBRTreeMaker::addTypeICorr( edm::Event const & event ){
   TypeICorrMap_.clear();

   event.getByToken(jetToken_      , jets_    );
   event.getByToken(rhoToken_      , rho_     );
   edm::Handle<reco::VertexCollection> vertices_;
   event.getByToken(vtxToken_, vertices_);
   edm::Handle<edm::View<pat::Muon>> muons_;
   event.getByToken(t1muSrc_,muons_);
   bool skipEM_                    = true;
   double skipEMfractionThreshold_ = 0.9;
   bool skipMuons_                 = true;
   std::string skipMuonSelection_string = "isGlobalMuon | isStandAloneMuon";
   StringCutObjectSelector<reco::Candidate>* skipMuonSelection_ = new StringCutObjectSelector<reco::Candidate>(skipMuonSelection_string,true);
   double jetCorrEtaMax_           = 9.9;
   double type1JetPtThreshold_     = 15.0; //10.0;
   double corrEx    = 0;
   double corrEy    = 0;
   double corrSumEt = 0;
   for (const pat::Jet &jet : *jets_) {
     double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
     if ( skipEM_ && emEnergyFraction > skipEMfractionThreshold_ ) continue;
     reco::Candidate::LorentzVector rawJetP4 = jet.correctedP4(0);
     double corr = getJEC(rawJetP4, jet, jetCorrEtaMax_, jetCorrLabel_);
     if ( skipMuons_ ) {
          const std::vector<reco::CandidatePtr> & cands = jet.daughterPtrVector();
          for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();
             cand != cands.end(); ++cand ) {
     	     const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
     	     const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
               if ( mu != 0 && (*skipMuonSelection_)(*mu) ) {
                  reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
                  rawJetP4 -= muonP4;
               }
           }
      }
     reco::Candidate::LorentzVector corrJetP4 = corr*rawJetP4;
     if ( corrJetP4.pt() > type1JetPtThreshold_ ) {
                 reco::Candidate::LorentzVector tmpP4 = jet.correctedP4(0);
                 corr = getJECOffset(tmpP4, jet, jetCorrEtaMax_, offsetCorrLabel_);
                 reco::Candidate::LorentzVector rawJetP4offsetCorr = corr*rawJetP4;
                 corrEx    -= (corrJetP4.px() - rawJetP4offsetCorr.px());
                 corrEy    -= (corrJetP4.py() - rawJetP4offsetCorr.py());
                 corrSumEt += (corrJetP4.Et() - rawJetP4offsetCorr.Et());
      }
  }
 TypeICorrMap_["corrEx"]    = corrEx;
 TypeICorrMap_["corrEy"]    = corrEy;
 TypeICorrMap_["corrSumEt"] = corrSumEt;
}

//-------------------------------------------------------------------------------------------------------------------------------------//
math::XYZTLorentzVector
EDBRTreeMaker::getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType){
    double leppt = lep.Pt();
    double lepphi = lep.Phi();
    double lepeta = lep.Eta();
    double lepenergy = lep.Energy();
    double metpt = MetPt;
    double metphi = MetPhi;
    double  px = metpt*cos(metphi);
    double  py = metpt*sin(metphi);
    double  pz = 0;
    double  pxl= leppt*cos(lepphi);
    double  pyl= leppt*sin(lepphi);
    double  pzl= leppt*sinh(lepeta);
    double  El = lepenergy;
    double  a = pow(MW_,2) + pow(px+pxl,2) + pow(py+pyl,2) - px*px - py*py - El*El + pzl*pzl;
    double  b = 2.*pzl;
    double  A = b*b -4.*El*El;
    double  B = 2.*a*b;
    double  C = a*a-4.*(px*px+py*py)*El*El;
    double M_mu =  0;
    int type=2; // use the small abs real root
    a = MW_*MW_ - M_mu*M_mu + 2.0*pxl*px + 2.0*pyl*py;
    A = 4.0*(El*El - pzl*pzl);
    B = -4.0*a*pzl;
    C = 4.0*El*El*(px*px + py*py) - a*a;
    double tmproot = B*B - 4.0*A*C;
    if (tmproot<0) {
        pz = - B/(2*A); // take real part of complex roots
    }
    else {
        double tmpsol1 = (-B + sqrt(tmproot))/(2.0*A);
        double tmpsol2 = (-B - sqrt(tmproot))/(2.0*A);
        if (type == 0 ) {
            if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
            else { pz = tmpsol1; }
            if ( abs(pz) > 300. ) {
                if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
                else { pz = tmpsol2; }
            }
        }
        if (type == 1 ) {
            // two real roots, pick the one closest to pz of muon
            if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
            else {pz = tmpsol1; }
        }
        if (type == 2 ) {
            // pick the most central root.
            if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
            else { pz = tmpsol2; }
        }
        /*if (type == 3 ) {
         // pick the largest value of the cosine
         TVector3 p3w, p3mu;
         p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol1);
         p3mu.SetXYZ(pxl, pyl, pzl );
         double sinthcm1 = 2.*(p3mu.Perp(p3w))/MW_;
         p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol2);
         double sinthcm2 = 2.*(p3mu.Perp(p3w))/MW_;
         double costhcm1 = sqrt(1. - sinthcm1*sinthcm1);
         double costhcm2 = sqrt(1. - sinthcm2*sinthcm2);
         if ( costhcm1 > costhcm2 ) { pz = tmpsol1; otherSol_ = tmpsol2; }
         else { pz = tmpsol2;otherSol_ = tmpsol1; }
         }*///end of type3
    }//endl of if real root
    //dont correct pt neutrino
    math::XYZTLorentzVector outP4(px,py,pz,sqrt(px*px+py*py+pz*pz));
    return outP4;
}//end neutrinoP4


// ------------ method called for each event  ------------
void
EDBRTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // DeepAK8
         fatjetNN_->readEvent(iEvent, iSetup);
            
                // Decorrelated DeepAK8
         decorrNN_->readEvent(iEvent, iSetup);
    
   using namespace edm;
   setDummyValues(); //Initalize variables with dummy values
   nevent = iEvent.eventAuxiliary().event();
   run    = iEvent.eventAuxiliary().run();
   ls     = iEvent.eventAuxiliary().luminosityBlock();
   Handle<TriggerResults> trigRes;
   iEvent.getByToken(hltToken_, trigRes);
 

   int mtemp1=0;
   for (size_t i=0; i<muPaths1.size();i++) {
      mtemp1 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths1[i]));
      if(HLT_Mu1<mtemp1) HLT_Mu1=mtemp1;
   }
   int mtemp2=0;
   for (size_t i=0; i<muPaths2.size();i++) {
      mtemp2 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths2[i]));
      if(HLT_Mu2<mtemp2) HLT_Mu2=mtemp2;
   }
   int mtemp3=0;
   for (size_t i=0; i<muPaths3.size();i++) {
      mtemp3 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths3[i]));
      if(HLT_Mu3<mtemp3) HLT_Mu3=mtemp3;
   }
   int mtemp4=0;
   for (size_t i=0; i<muPaths4.size();i++) {
      mtemp4 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths4[i]));
      if(HLT_Mu4<mtemp4) HLT_Mu4=mtemp4;
   }
   int mtemp5=0;
   for (size_t i=0; i<muPaths5.size();i++) {
      mtemp5 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths5[i]));
      if(HLT_Mu5<mtemp5) HLT_Mu5=mtemp5;
   }
   int mtemp6=0;
   for (size_t i=0; i<muPaths6.size();i++) {
      mtemp6 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths6[i]));
      if(HLT_Mu6<mtemp6) HLT_Mu6=mtemp6;
   }
   int mtemp7=0;
   for (size_t i=0; i<muPaths7.size();i++) {
      mtemp7 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths7[i]));
      if(HLT_Mu7<mtemp7) HLT_Mu7=mtemp7;
   }
   int mtemp8=0;
   for (size_t i=0; i<muPaths8.size();i++) {
      mtemp8 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths8[i]));
      if(HLT_Mu8<mtemp8) HLT_Mu8=mtemp8;
   }
   int mtemp9=0;
   for (size_t i=0; i<muPaths9.size();i++) {
      mtemp9 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths9[i]));
      if(HLT_Mu9<mtemp9) HLT_Mu9=mtemp9;
   }
   int mtemp10=0;
   for (size_t i=0; i<el1.size();i++) {
      mtemp10 = (int)trigRes->accept(hltConfig.triggerIndex(el1[i]));
      if(HLT_Mu10<mtemp10) HLT_Mu10=mtemp10;
   }
   int mtemp11=0;
   for (size_t i=0; i<el2.size();i++) {
      mtemp11 = (int)trigRes->accept(hltConfig.triggerIndex(el2[i]));
      if(HLT_Mu11<mtemp11) HLT_Mu11=mtemp11;
   }

   int mtemp12=0;
   for (size_t i=0; i<el3.size();i++) {
      mtemp12 = (int)trigRes->accept(hltConfig.triggerIndex(el3[i]));
      if(HLT_Mu12<mtemp12) HLT_Mu12=mtemp12;
   }

   int mtemp13=0;
   for (size_t i=0; i<mu1.size();i++) {
      mtemp13 = (int)trigRes->accept(hltConfig.triggerIndex(mu1[i]));
      if(HLT_Mu13<mtemp13) HLT_Mu13=mtemp13;
   }

   int mtemp14=0;
   for (size_t i=0; i<mu2.size();i++) {
      mtemp14 = (int)trigRes->accept(hltConfig.triggerIndex(mu2[i]));
      if(HLT_Mu14<mtemp14) HLT_Mu14=mtemp14;
   }  

   int mtemp15=0;
   for (size_t i=0; i<mu3.size();i++) {
      mtemp15 = (int)trigRes->accept(hltConfig.triggerIndex(mu3[i]));
      if(HLT_Mu15<mtemp15) HLT_Mu15=mtemp15;
   }  

   int mtemp16=0;
   for (size_t i=0; i<mu4.size();i++) {
      mtemp16 = (int)trigRes->accept(hltConfig.triggerIndex(mu4[i]));
      if(HLT_Mu16<mtemp16) HLT_Mu16=mtemp16;
   }  

   edm::Handle<edm::View<pat::Jet> > hadronicVs;
   iEvent.getByToken(hadronicVSrc_, hadronicVs);
   //edm::Handle<pat::JetCollection> hadronicVSoftDrop;
   //iEvent.getByToken(hadronicVSoftDropSrc_, hadronicVSoftDrop);
   edm::Handle<edm::View<reco::Candidate> > leptonicVs;
   iEvent.getByToken(leptonicVSrc_, leptonicVs);
   edm::Handle<double> rho;

   iEvent.getByToken(rhoToken_      , rho     );
   double fastJetRho = *(rho.product());
   useless = fastJetRho;
   edm::Handle<edm::View<pat::Jet> > ak4jets;
   iEvent.getByToken(ak4jetsSrc_, ak4jets);
   edm::Handle<edm::View<reco::Candidate> > gravitons;

//   iEvent.getByLabel(gravitonSrc_.c_str(), gravitons);
   iEvent.getByToken(gravitonSrc_, gravitons);
   edm::Handle<edm::View<reco::Candidate> > metHandle;
   iEvent.getByToken(metSrc_, metHandle);
   edm::Handle<edm::View<pat::Muon>> loosemus;
   iEvent.getByToken(loosemuonToken_,loosemus);

   edm::Handle<edm::View<pat::Muon>> vetomus;
   iEvent.getByToken(vetomuonToken_,vetomus);
   edm::Handle<edm::View<pat::Electron>> looseels;

   iEvent.getByToken(looseelectronToken_, looseels);
   edm::Handle<edm::View<pat::Electron>> vetoels;

   iEvent.getByToken(vetoelectronToken_, vetoels);
   edm::Handle<edm::View<reco::GenParticle> > genParticles;//define genParticle
   iEvent.getByToken(genSrc_, genParticles);
   if (RunOnSig_||RunOnMC_){
     // edm::Handle<LHEEventProduct> wgtsource;
     // iEvent.getByLabel("externalLHEProducer", wgtsource);
     // iEvent.getByLabel("source", wgtsource);
     edm::Handle<LHEEventProduct> wgtsource;
     iEvent.getByToken(LheToken_, wgtsource);

//std::cout<<"weight number "<<wgtsource->weights().size()<<std::endl;
     for ( int i=0; i<211; ++i) {
       pweight[i]= -99;//
     //  pweight[i]= wgtsource->weights()[i].wgt/wgtsource->originalXWGTUP();
////cout<<wgtsource->weights()[i].id<<" "<<pweight[i]<<endl;
                              }
  //   for ( int i=9; i<110; ++i) {
    //    pweight[i]= wgtsource->weights()[i+101].wgt/wgtsource->originalXWGTUP();
////cout<<wgtsource->weights()[i].id<<" "<<pweight[i]<<endl;
    //                            }
     edm::Handle<GenEventInfoProduct> genEvtInfo;
     iEvent.getByToken(GenToken_,genEvtInfo);
     theWeight = genEvtInfo->weight();
     if(theWeight>0) nump = nump+1;
     if(theWeight<0) numm = numm+1;

     edm::Handle<std::vector<PileupSummaryInfo>>  PupInfo;
     iEvent.getByToken(PUToken_, PupInfo);
     std::vector<PileupSummaryInfo>::const_iterator PVI;
     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         nBX = PVI->getBunchCrossing();
         if(nBX == 0) { // "0" is the in-time crossing, negative values are the early crossings, positive are late
            npT = PVI->getTrueNumInteractions();
            npIT = PVI->getPU_NumInteractions();
         }
      }
   }
//filter
    iEvent.getByToken(noiseFilterToken_, noiseFilterBits_);
    const edm::TriggerNames &names = iEvent.triggerNames(*noiseFilterBits_);
    for (unsigned int i = 0, n = noiseFilterBits_->size(); i < n; ++i) {
      if (names.triggerName(i) == HBHENoiseFilter_Selector_)
        passFilter_HBHE_ = noiseFilterBits_->accept(i); // TO BE USED
      if (names.triggerName(i) == HBHENoiseIsoFilter_Selector_)
        passFilter_HBHEIso_ = noiseFilterBits_->accept(i); // TO BE USED
      if (names.triggerName(i) == GlobalHaloNoiseFilter_Selector_)
        passFilter_GlobalHalo_ = noiseFilterBits_->accept(i); // TO BE USED
      if (names.triggerName(i) == ECALDeadCellNoiseFilter_Selector_)
        passFilter_ECALDeadCell_ = noiseFilterBits_->accept(i); // under scrutiny
      if (names.triggerName(i) == GoodVtxNoiseFilter_Selector_)
        passFilter_GoodVtx_ = noiseFilterBits_->accept(i); // TO BE USED
      if (names.triggerName(i) == EEBadScNoiseFilter_Selector_)
        passFilter_EEBadSc_ = noiseFilterBits_->accept(i); // under scrutiny  
    }
     edm::Handle<bool> badMuonResultHandle;
     edm::Handle<bool> badChargedHadronResultHandle;
     iEvent.getByToken(badMuonNoiseFilter_Selector_, badMuonResultHandle);
     iEvent.getByToken(badChargedHadronNoiseFilter_Selector_, badChargedHadronResultHandle);
     passFilter_badMuon_ = *badMuonResultHandle;
     passFilter_badChargedHadron_ = *badChargedHadronResultHandle;

   numCands = gravitons->size();
// ************************* Gen Level Information******************//
   if(RunOnMC_)
   {//MC Info
       Int_t havegra=0;
	  for(size_t ik=0; ik<genParticles->size();ik++)
	{// loop on gen
            const reco::Candidate* ptop0 = &(*genParticles)[ik];
            const reco::Candidate* ptop=findLastParticle(ptop0,6);
                if(ptop0->pdgId()== 6 && gentop_pt==-99) {
                    gentop_pt = ptop->pt();
                    gentop_eta = ptop->eta();
                    gentop_phi = ptop->phi();
                    gentop_mass = ptop->mass();
                    for(int i=0;ptop->daughter(i)!=NULL;i++){
                        if(abs(ptop->daughter(i)->pdgId())==24){
                            gent_w_pt=ptop->daughter(i)->pt();
                            gent_w_eta=ptop->daughter(i)->eta();
                            gent_w_phi=ptop->daughter(i)->phi();
                            gent_w_mass=ptop->daughter(i)->mass();
                            const reco::Candidate* ptw0 = ptop->daughter(i);
                            const reco::Candidate* ptw= findLastParticle(ptw0,24);
                            if(ptw->daughter(0)!=NULL)
                            {
                                if( abs(ptw->daughter(0)->pdgId())<=6 ){
                                    gent_w_tag=4;
                                    gent_w_q1_pt=ptw->daughter(0)->pt();
                                    gent_w_q1_eta=ptw->daughter(0)->eta();
                                    gent_w_q1_phi=ptw->daughter(0)->phi();
                                    gent_w_q1_e=ptw->daughter(0)->energy();
                                    gent_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    gent_w_q2_pt=ptw->daughter(1)->pt();
                                    gent_w_q2_eta=ptw->daughter(1)->eta();
                                    gent_w_q2_phi=ptw->daughter(1)->phi();
                                    gent_w_q2_e=ptw->daughter(1)->energy();
                                    gent_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                                //11-15      ->       11-16
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12 ) gent_w_tag=1;
                                if( abs(ptw->daughter(0)->pdgId())==12 ||abs(ptw->daughter(0)->pdgId())==13 ) gent_w_tag=2;
                                if( abs(ptw->daughter(0)->pdgId())==14 ||abs(ptw->daughter(0)->pdgId())==15 ) gent_w_tag=3;
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12||abs(ptw->daughter(0)->pdgId())==13 ||abs(ptw->daughter(0)->pdgId())==14||abs(ptw->daughter(0)->pdgId())==15 ||abs(ptw->daughter(0)->pdgId())==16)
                                {
                                    gent_w_q1_pt=ptw->daughter(0)->pt();
                                    gent_w_q1_eta=ptw->daughter(0)->eta();
                                    gent_w_q1_phi=ptw->daughter(0)->phi();
                                    gent_w_q1_e=ptw->daughter(0)->energy();
                                    gent_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    gent_w_q2_pt=ptw->daughter(1)->pt();
                                    gent_w_q2_eta=ptw->daughter(1)->eta();
                                    gent_w_q2_phi=ptw->daughter(1)->phi();
                                    gent_w_q2_e=ptw->daughter(1)->energy();
                                    gent_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                        }
                        }
                        if(abs(ptop->daughter(i)->pdgId())==5){
                            gent_b_pt=ptop->daughter(i)->pt();
                            gent_b_eta=ptop->daughter(i)->eta();
                            gent_b_phi=ptop->daughter(i)->phi();
                            gent_b_mass=ptop->daughter(i)->mass();
                        }
                }
                }
                if(ptop0->pdgId()== -6 && genantitop_pt==-99) {
                    genantitop_pt = ptop->pt();
                    genantitop_eta = ptop->eta();
                    genantitop_phi = ptop->phi();
                    genantitop_mass = ptop->mass();
                    for(int i=0;ptop->daughter(i)!=NULL;i++){
                        if(abs(ptop->daughter(i)->pdgId())==24){
                            genantit_w_pt=ptop->daughter(i)->pt();
                            genantit_w_eta=ptop->daughter(i)->eta();
                            genantit_w_phi=ptop->daughter(i)->phi();
                            genantit_w_mass=ptop->daughter(i)->mass();
                            const reco::Candidate* ptw0 = ptop->daughter(i);
                            const reco::Candidate* ptw=findLastParticle(ptw0,24);
                            if(ptw->daughter(0)!=NULL)
                            {
                                if( abs(ptw->daughter(0)->pdgId())<=6 ){
                                    genantit_w_tag=4;
                                    genantit_w_q1_pt=ptw->daughter(0)->pt();
                                    genantit_w_q1_eta=ptw->daughter(0)->eta();
                                    genantit_w_q1_phi=ptw->daughter(0)->phi();
                                    genantit_w_q1_e=ptw->daughter(0)->energy();
                                    genantit_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    genantit_w_q2_pt=ptw->daughter(1)->pt();
                                    genantit_w_q2_eta=ptw->daughter(1)->eta();
                                    genantit_w_q2_phi=ptw->daughter(1)->phi();
                                    genantit_w_q2_e=ptw->daughter(1)->energy();
                                    genantit_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12 ) genantit_w_tag=1;
                                if( abs(ptw->daughter(0)->pdgId())==12 ||abs(ptw->daughter(0)->pdgId())==13 ) genantit_w_tag=2;
                                if( abs(ptw->daughter(0)->pdgId())==14 ||abs(ptw->daughter(0)->pdgId())==15 ) genantit_w_tag=3;
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12||abs(ptw->daughter(0)->pdgId())==13 ||abs(ptw->daughter(0)->pdgId())==14||abs(ptw->daughter(0)->pdgId())==15 ||abs(ptw->daughter(0)->pdgId())==16)
                                {
                                    genantit_w_q1_pt=ptw->daughter(0)->pt();
                                    genantit_w_q1_eta=ptw->daughter(0)->eta();
                                    genantit_w_q1_phi=ptw->daughter(0)->phi();
                                    genantit_w_q1_e=ptw->daughter(0)->energy();
                                    genantit_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    genantit_w_q2_pt=ptw->daughter(1)->pt();
                                    genantit_w_q2_eta=ptw->daughter(1)->eta();
                                    genantit_w_q2_phi=ptw->daughter(1)->phi();
                                    genantit_w_q2_e=ptw->daughter(1)->energy();
                                    genantit_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                            }
                        }
                        if(abs(ptop->daughter(i)->pdgId())==5){
                            genantit_b_pt=ptop->daughter(i)->pt();
                            genantit_b_eta=ptop->daughter(i)->eta();
                            genantit_b_phi=ptop->daughter(i)->phi();
                            genantit_b_mass=ptop->daughter(i)->mass();
                        }
                    }
                }

		if( abs((*genParticles)[ik].pdgId())==9000024 || abs((*genParticles)[ik].pdgId())==6 )
		{//if Wkk
                havegra=1;
                const reco::Candidate* pwkk0 = &(*genParticles)[ik];
                const reco::Candidate* pwkk=findLastParticle(pwkk0,9000024);
                   gen_gra_m=pwkk->mass();
		   gen_gra_pt=pwkk->pt();
		   gen_gra_eta=pwkk->eta();
                   gen_gra_phi=pwkk->phi();
                     qnumber=1;
                for(int i=0;pwkk->daughter(i)!=NULL;i++)
                   {//loop on Wkk daughter
                    if(abs(pwkk->daughter(i)->pdgId())==24)
                       {//if w
                         const reco::Candidate* pw0 = pwkk->daughter(i);
                           const reco::Candidate* pw= findLastParticle(pw0,24);  
                         if(pw->daughter(0)!=NULL)
                         {//loop on w daughter
                            const reco::Candidate* pl = pw->daughter(0);
                            if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)|| (abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16))
                            {//beign of lep-w
                                ptGenVlep = pw->pt();////////////something wrong?
                                etaGenVlep = pw->eta();
                                phiGenVlep = pw->phi();
                                massGenVlep = pw->mass();
                                status_1=0;
                            for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())==11)
                               {
                                 gen_ele_pt=pl->pt();
                                 gen_ele_eta=pl->eta();
                                 gen_ele_phi=pl->phi();
                                 gen_ele_e=pl->energy();
                                   status_1=1;
                               }
                               if(abs(pl->pdgId())==13)
                               {
                                 gen_mu_pt=pl->pt();
                                 gen_mu_eta=pl->eta();
                                 gen_mu_phi=pl->phi();
                                 gen_mu_e=pl->energy();
                                   status_1=2;
                               }
                                if(abs(pl->pdgId())==15)
                                    {
                            gen_tau_pt=pl->pt();
                            gen_tau_eta=pl->eta();
                            gen_tau_phi=pl->phi();
                            gen_tau_e=pl->energy();
                                        status_1=3;
                                    }
                                }
                             }//end of if lep-w

                             for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())<6)
                               {
                              if(qnumber==1)
                                 {
                                 ptq11 = pl->pt();
                                 etaq11 = pl->eta();
                                 phiq11 = pl->phi();
                                 massq11 = pl->mass();
                                 }
                                 if(qnumber==2)
                                 {
                                 ptq12 = pl->pt();
                                 etaq12 = pl->eta();
                                 phiq12 = pl->phi();
                                 massq12 = pl->mass();
                                 }
                                 qnumber=qnumber+1;
                               }                                          
                                                                     }
                         //    pl = pw->daughter(0);    // @@@@@@@@@@  delete

        		     if(abs(pl->pdgId())<6) 
                             {
			         ptGenVhad = pw->pt();
			         etaGenVhad = pw->eta();
			         phiGenVhad = pw->phi();
			         massGenVhad = pw->mass();
                                 status_1=4;
                             }
        		     if(abs(pl->pdgId())==24) 
                             {
                                 status_1=5;
                             }
                     //if(status_1<0) cout<<"pw->daughter(0)  "<<pl->pdgId()<<endl;
			   }//end of loop on w daughter
                       }//end of if w				 
                   }//end of loop on Wkk daughter
 		}//end of if Wkk
        if( havegra==0&&abs((*genParticles)[ik].pdgId())==24 )
        {//if W
             qnumber=1;
            const reco::Candidate* pw0 = &(*genParticles)[ik];
            const reco::Candidate* pw=findFirstParticle(pw0,24);
            if(pw->mother(0)->pdgId()!=9000025){
                const reco::Candidate* pw1= findLastParticle(pw0,24);
                if(pw1->daughter(0)!=NULL)
                {//loop on w daughter
                const reco::Candidate* pl = pw1->daughter(0);
                if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)|| (abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16))
                {//beign of lep-w
                    ptGenVlep = pw1->pt();
                    etaGenVlep = pw1->eta();
                    phiGenVlep = pw1->phi();
                    massGenVlep = pw1->mass();
                    status_1=0;
                    for(int ii=0;pw1->daughter(ii)!=NULL;ii++){
                        const reco::Candidate* pl = pw1->daughter(ii);
                        if(abs(pl->pdgId())==11)
                        {
                            gen_ele_pt=pl->pt();
                            gen_ele_eta=pl->eta();
                            gen_ele_phi=pl->phi();
                            gen_ele_e=pl->energy();
                            status_1=1;
                        }
                        if(abs(pl->pdgId())==13)
                        {
                            gen_mu_pt=pl->pt();
                            gen_mu_eta=pl->eta();
                            gen_mu_phi=pl->phi();
                            gen_mu_e=pl->energy();
                            status_1=2;
                        }
                        if(abs(pl->pdgId())==15)
                        {
                            gen_tau_pt=pl->pt();
                            gen_tau_eta=pl->eta();
                            gen_tau_phi=pl->phi();
                            gen_tau_e=pl->energy();
                            status_1=3;
                        }
                    }
                }//end of if lep-w

                             for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())<6)
                               {
                              if(qnumber==1)
                                 {
                                 ptq11 = pl->pt();
                                 etaq11 = pl->eta();
                                 phiq11 = pl->phi();
                                 massq11 = pl->mass();
                                 }
                                 if(qnumber==2)
                                 {
                                 ptq12 = pl->pt();
                                 etaq12 = pl->eta();
                                 phiq12 = pl->phi();
                                 massq12 = pl->mass();
                                 }
                                 qnumber=qnumber+1;
                               }
                                                                     }
                             pl = pw->daughter(0);


                if(abs(pl->pdgId())<6)
                {
                    ptGenVhad = pw1->pt();
                    etaGenVhad = pw1->eta();
                    phiGenVhad = pw1->phi();
                    massGenVhad = pw1->mass();
                    status_1=4;
                }
                if(abs(pl->pdgId())==24)
                {
                    status_1=5;
                }
                }
            }
        }
		if( abs((*genParticles)[ik].pdgId())==9000025 )
		{//if Radion
                   gen_rad_m=(*genParticles)[ik].mass();
		   gen_rad_pt=(*genParticles)[ik].pt();
	           gen_rad_eta=(*genParticles)[ik].eta();
                   gen_rad_phi=(*genParticles)[ik].phi();
                   qnumber2=1;
                   qnumber3=1;
                   for(int i=0;(*genParticles)[ik].daughter(i)!=NULL;i++)
                   {//loop on Radion daughter
                                 
                      if(((*genParticles)[ik].daughter(i)->pdgId())==24)
                       {//if w-
                         const reco::Candidate* pw0 = (*genParticles)[ik].daughter(i);
                           const reco::Candidate* pw= findLastParticle(pw0,24);
                         if(pw->daughter(0)!=NULL)
                         {//loop on w daughter
                            const reco::Candidate* pl = pw->daughter(0);
                            if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)||(abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16))
                            {//beign of lep-w
                               ptGenVlep_2 = pw->pt();
                               etaGenVlep_2 = pw->eta();
                               phiGenVlep_2 = pw->phi();
                               massGenVlep_2 = pw->mass();
                               status_2=0;
                            for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                               const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())==11)
                               {
                                 gen_ele_pt_2=pl->pt();
                                 gen_ele_eta_2=pl->eta();
                                 gen_ele_phi_2=pl->phi();
                                 gen_ele_e_2=pl->energy();
                                 status_2=1;  
                               }
                               if(abs(pl->pdgId())==13)
                               {
                                 gen_mu_pt_2=pl->pt();
                                 gen_mu_eta_2=pl->eta();
                                 gen_mu_phi_2=pl->phi();
                                 gen_mu_e_2=pl->energy();
                                 status_2=2;  
                               }
				if(abs(pl->pdgId())==15)
                               {
                                 gen_tau_pt_2=pl->pt();
                                 gen_tau_eta_2=pl->eta();
                                 gen_tau_phi_2=pl->phi();
                                 gen_tau_e_2=pl->energy();
                                 status_2=3;
                               }
                            }
                             }//end of if lep-w

                             for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())<6)
                               {
                                if(qnumber2==1)
                                 {
                                 ptq21 = pl->pt();
                                 etaq21 = pl->eta();
                                 phiq21 = pl->phi();
                                 massq21 = pl->mass();
                                 }
                                 if(qnumber2==2)
                                 {
                                 ptq22 = pl->pt();
                                 etaq22 = pl->eta();
                                 phiq22 = pl->phi();
                                 massq22 = pl->mass();
                                 }
                                 qnumber2=qnumber2+1;                               
                               }
                                                                     }
                             pl = pw->daughter(0);

        		     if(abs(pl->pdgId())<6) 
                             {
                                 ptGenV_2 = pw->pt();
                                 etaGenV_2 = pw->eta();
                                 phiGenV_2 = pw->phi();
                                 massGenV_2 = pw->mass();
                                 status_2=4;  
                             }
        		     if(abs(pl->pdgId())==24) 
                             {
                                 status_2=5;  
                             }
			   }//end of loop on w daughter
                       }//end of if w-	
                      if(((*genParticles)[ik].daughter(i)->pdgId())==-24)
                       {//if w+
                         const reco::Candidate* pw0 = (*genParticles)[ik].daughter(i);
                           const reco::Candidate* pw= findLastParticle(pw0,24);
			//cout<<((*genParticles)[ik].daughter(i)->pdgId())<<endl;
                         if(pw->daughter(0)!=NULL)
                         {//loop on w daughter
                            const reco::Candidate* pl = pw->daughter(0);
			//cout<<(pl->pdgId())<<endl;
                            if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)|| (abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16))
                            {//beign of lep-w
                               ptGenVlep_3 = pw->pt();
                               etaGenVlep_3 = pw->eta();
                               phiGenVlep_3 = pw->phi();
                               massGenVlep_3 = pw->mass();
                               status_3=0;
                            for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                              const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())==11)
                               {
                                 gen_ele_pt_3=pl->pt();
                                 gen_ele_eta_3=pl->eta();
                                 gen_ele_phi_3=pl->phi();
                                 gen_ele_e_3=pl->energy();
                                 status_3=1;  
                               }
                               if(abs(pl->pdgId())==13)
                               {
                                 gen_mu_pt_3=pl->pt();
                                 gen_mu_eta_3=pl->eta();
                                 gen_mu_phi_3=pl->phi();
                                 gen_mu_e_3=pl->energy();
                                 status_3=2;  
                               }
                               if(abs(pl->pdgId())==15)
                               {
                                 gen_tau_pt_3=pl->pt();
                                 gen_tau_eta_3=pl->eta();
                                 gen_tau_phi_3=pl->phi();
                                 gen_tau_e_3=pl->energy();
                                 status_3=3;  
                               }
                            }
                             }//end of if lep-w

                             for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                               if(abs(pl->pdgId())<6)
                               {
                                if(qnumber3==1)
                                 {
                                 ptq31 = pl->pt();
                                 etaq31 = pl->eta();
                                 phiq31 = pl->phi();
                                 massq31 = pl->mass();
                                 }
                                 if(qnumber3==2)
                                 {
                                 ptq32 = pl->pt();
                                 etaq32 = pl->eta();
                                 phiq32 = pl->phi();
                                 massq32 = pl->mass();
                                 }
                                 qnumber3=qnumber3+1;
                               }
                                                                      }
                             pl = pw->daughter(0);


        		     if(abs(pl->pdgId())<6) 
                             {
                                 ptGenV_3 = pw->pt();
                                 etaGenV_3 = pw->eta();
                                 phiGenV_3 = pw->phi();
                                 massGenV_3 = pw->mass();
                                 status_3=4;  
                             }
        		     if(abs(pl->pdgId())==24) 
                             {
                                 status_3=5;  
                             }
			   }//end of loop on w daughter
                       }//end of if w+	
			 
                   }//end of loop on Radion daughter
 		}//end of if Radion

        }//end of loop on gen

    //w and top info
        for( auto p=genParticles->begin(); p!= genParticles->end(); ++p)
        {}
        int igenw=0;
        int sizew=5;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
            if(abs((*genParticles)[ik].pdgId())==24)
                {
                    const reco::Candidate* pwtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pwtmp=findLastParticle(pwtmp1,24);
                    int woverlap=0;
                    for (int ia=0;ia<igenw;ia++){
                        if(pwtmp->pt()==ptgenwl[ia]) woverlap=1;
                    }
                    const reco::Candidate* pwtmp2=findFirstParticle(pwtmp1,24);
                    if(pwtmp->pt()>50&&igenw<sizew&&woverlap==0){
                    ptgenwl[igenw] = pwtmp->pt();
                    etagenwl[igenw] = pwtmp->eta();
                    phigenwl[igenw] = pwtmp->phi();
                    massgenwl[igenw] = pwtmp->mass();
                    ptgenwf[igenw] = pwtmp2->pt();
                    etagenwf[igenw] = pwtmp2->eta();
                    phigenwf[igenw] = pwtmp2->phi();
                    massgenwf[igenw] = pwtmp2->mass();
                        taggenwmother[igenw]=pwtmp2->mother(0)->pdgId();
                    if(pwtmp->daughter(0)!=NULL)//loop on w daughter
                    {
                        const reco::Candidate* pltmp = pwtmp->daughter(0);
                         if( (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ){
                                taggenwl[igenw]=1;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ){
                            taggenwl[igenw]=2;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16) ){
                            taggenwl[igenw]=3;                    }//end of w daugter loop
                        if(abs(pltmp->pdgId())<6 ) {
                            taggenwl[igenw]=4;
                            genw_q1_pt[igenw]=pwtmp->daughter(0)->pt();
                            genw_q1_eta[igenw]=pwtmp->daughter(0)->eta();
                            genw_q1_phi[igenw]=pwtmp->daughter(0)->phi();
                            genw_q1_e[igenw]=pwtmp->daughter(0)->energy();
                            genw_q1_pdg[igenw]=pwtmp->daughter(0)->pdgId();
                            genw_q2_pt[igenw]=pwtmp->daughter(1)->pt();
                            genw_q2_eta[igenw]=pwtmp->daughter(1)->eta();
                            genw_q2_phi[igenw]=pwtmp->daughter(1)->phi();
                            genw_q2_e[igenw]=pwtmp->daughter(1)->energy();
                            genw_q2_pdg[igenw]=pwtmp->daughter(1)->pdgId();
                            }
                        if( (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ||(abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ||(abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16)){
                            genw_q1_pt[igenw]=pwtmp->daughter(0)->pt();
                            genw_q1_eta[igenw]=pwtmp->daughter(0)->eta();
                            genw_q1_phi[igenw]=pwtmp->daughter(0)->phi();
                            genw_q1_e[igenw]=pwtmp->daughter(0)->energy();
                            genw_q1_pdg[igenw]=pwtmp->daughter(0)->pdgId();
                            genw_q2_pt[igenw]=pwtmp->daughter(1)->pt();
                            genw_q2_eta[igenw]=pwtmp->daughter(1)->eta();
                            genw_q2_phi[igenw]=pwtmp->daughter(1)->phi();
                            genw_q2_e[igenw]=pwtmp->daughter(1)->energy();
                            genw_q2_pdg[igenw]=pwtmp->daughter(1)->pdgId();

                        }
                    }
                    igenw+=1;
                    }
                }//end of if w

        }

        int igenz=0;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
            if(abs((*genParticles)[ik].pdgId())==23)
            {
                const reco::Candidate* pztmp1 = &(*genParticles)[ik];
                const reco::Candidate* pztmp=findLastParticle(pztmp1,23);
                int zoverlap=0;
                for (int ia=0;ia<igenz;ia++){
                    if(pztmp->pt()==ptgenzl[ia]) zoverlap=1;}
                if(pztmp->pt()>50&&igenz<sizew&&zoverlap==0){
                    ptgenzl[igenz] = pztmp->pt();
                    etagenzl[igenz] = pztmp->eta();
                    phigenzl[igenz] = pztmp->phi();
                    massgenzl[igenz] = pztmp->mass();
                    const reco::Candidate* pztmp2=findFirstParticle(pztmp1,23);
                    ptgenzf[igenz] = pztmp2->pt();
                    etagenzf[igenz] = pztmp2->eta();
                    phigenzf[igenz] = pztmp2->phi();
                    massgenzf[igenz] = pztmp2->mass();
                    if(pztmp->daughter(0)!=NULL)//loop on w daughter
                    {
                        const reco::Candidate* pltmp = pztmp->daughter(0);
                        if( (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ){
                            taggenzl[igenz]=1;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ){
                            taggenzl[igenz]=2;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16) ){
                            taggenzl[igenz]=3;                    }//end of w daugter loop
                        if(abs(pltmp->pdgId())<6 ) {
                            taggenzl[igenz]=4;}
                    }
                    igenz+=1;
                }
            }//end of if w
        }

        int igeng=0;
        int sizeg=15;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
            if(abs((*genParticles)[ik].pdgId())==21)
            {
                const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                const reco::Candidate* pgtmp=findLastParticle(pgtmp1,21);
                int goverlap=0;
                for (int ia=0;ia<igeng;ia++){
                    if(pgtmp->pt()==ptgengl[ia]) goverlap=1;}
                if(pgtmp->pt()>50&&igeng<sizeg&&goverlap==0){
                    ptgengl[igeng] = pgtmp->pt();
                    etagengl[igeng] = pgtmp->eta();
                    phigengl[igeng] = pgtmp->phi();
                    egengl[igeng] = pgtmp->energy();
                    const reco::Candidate* pgtmp2=findFirstParticle(pgtmp,21);
                    ptgengf[igeng] = pgtmp2->pt();
                    etagengf[igeng] = pgtmp2->eta();
                    phigengf[igeng] = pgtmp2->phi();
                    egengf[igeng] = pgtmp2->energy();
                    mothergengf[igeng] = pgtmp2->mother(0)->pdgId();
                    const reco::Candidate* pgtmp3=findLastParticle(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                    if (pgtmp3->mother(0)!=NULL) mmothergengf[igeng] = pgtmp3->mother(0)->pdgId();
                    igeng+=1;
                }
            }//end of if w

        }

        int igenq1=0;
        int sizeq1=5;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
                if(abs((*genParticles)[ik].pdgId())==1)
                {
                    const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pgtmp=findLastParticle(pgtmp1,1);
                    int goverlap=0;
                    for (int ia=0;ia<igenq1;ia++){
                        if(pgtmp->pt()==ptgenq1l[ia]) goverlap=1;}
                    const reco::Candidate* pgtmp2=findFirstParticle(pgtmp,1);
                    if(pgtmp->pt()>50&&igenq1<sizeq1&&goverlap==0&&abs(pgtmp2->mother(0)->pdgId())!=24){
                        ptgenq1l[igenq1] = pgtmp->pt();
                        etagenq1l[igenq1] = pgtmp->eta();
                        phigenq1l[igenq1] = pgtmp->phi();
                        egenq1l[igenq1] = pgtmp->energy();
                        ptgenq1f[igenq1] = pgtmp2->pt();
                        etagenq1f[igenq1] = pgtmp2->eta();
                        phigenq1f[igenq1] = pgtmp2->phi();
                        egenq1f[igenq1] = pgtmp2->energy();
                        mothergenq1f[igenq1] = pgtmp2->mother(0)->pdgId();
                        const reco::Candidate* pgtmp3=findLastParticle(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                        if (pgtmp3->mother(0)!=NULL) mmothergenq1f[igenq1] = pgtmp3->mother(0)->pdgId();
                        igenq1+=1;
                    }
            }
        }

            int igenq2=0;
            int sizeq2=5;
            for(size_t ik=0; ik<genParticles->size();ik++)
            {
                if(abs((*genParticles)[ik].pdgId())==2)
                {
                    const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pgtmp=findLastParticle(pgtmp1,2);
                    int goverlap=0;
                    for (int ia=0;ia<igenq2;ia++){
                        if(pgtmp->pt()==ptgenq2l[ia]) goverlap=1;}
                    const reco::Candidate* pgtmp2=findFirstParticle(pgtmp1,2);
                    if(pgtmp->pt()>50&&igenq2<sizeq2&&goverlap==0&&abs(pgtmp2->mother(0)->pdgId())!=24){
                        ptgenq2l[igenq2] = pgtmp->pt();
                        etagenq2l[igenq2] = pgtmp->eta();
                        phigenq2l[igenq2] = pgtmp->phi();
                        egenq2l[igenq2] = pgtmp->energy();
                        ptgenq2f[igenq2] = pgtmp2->pt();
                        etagenq2f[igenq2] = pgtmp2->eta();
                        phigenq2f[igenq2] = pgtmp2->phi();
                        egenq2f[igenq2] = pgtmp2->energy();
                        mothergenq2f[igenq2] = pgtmp2->mother(0)->pdgId();
                        const reco::Candidate* pgtmp3=findLastParticle(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                        if (pgtmp3->mother(0)!=NULL) mmothergenq2f[igenq2] = pgtmp3->mother(0)->pdgId();
                        igenq2+=1;
                    }
                }
            }

            int igenq3=0;
            int sizeq3=5;
            for(size_t ik=0; ik<genParticles->size();ik++)
            {
                if(abs((*genParticles)[ik].pdgId())==3)
                {
                    const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pgtmp=findLastParticle(pgtmp1,3);
                    int goverlap=0;
                    for (int ia=0;ia<igenq3;ia++){
                        if(pgtmp->pt()==ptgenq3l[ia]) goverlap=1;}
                    const reco::Candidate* pgtmp2=findFirstParticle(pgtmp1,3);
                    if(pgtmp->pt()>50&&igenq3<sizeq3&&goverlap==0&&abs(pgtmp2->mother(0)->pdgId())!=24){
                        ptgenq3l[igenq3] = pgtmp->pt();
                        etagenq3l[igenq3] = pgtmp->eta();
                        phigenq3l[igenq3] = pgtmp->phi();
                        egenq3l[igenq3] = pgtmp->energy();
                        ptgenq3f[igenq3] = pgtmp2->pt();
                        etagenq3f[igenq3] = pgtmp2->eta();
                        phigenq3f[igenq3] = pgtmp2->phi();
                        egenq3f[igenq3] = pgtmp2->energy();
                        mothergenq3f[igenq3] = pgtmp2->mother(0)->pdgId();
                        const reco::Candidate* pgtmp3=findLastParticle(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                        if (pgtmp3->mother(0)!=NULL) mmothergenq3f[igenq3] = pgtmp3->mother(0)->pdgId();
                        igenq3+=1;
                    }
                }
            }

            int igenq4=0;
            int sizeq4=5;
            for(size_t ik=0; ik<genParticles->size();ik++)
            {
                if(abs((*genParticles)[ik].pdgId())==4)
                {
                    const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pgtmp=findLastParticle(pgtmp1,4);
                    int goverlap=0;
                    for (int ia=0;ia<igenq4;ia++){
                        if(pgtmp->pt()==ptgenq4l[ia]) goverlap=1;}
                    const reco::Candidate* pgtmp2=findFirstParticle(pgtmp1,4);
                    if(pgtmp->pt()>50&&igenq4<sizeq4&&goverlap==0&&abs(pgtmp2->mother(0)->pdgId())!=24){
                        ptgenq4l[igenq4] = pgtmp->pt();
                        etagenq4l[igenq4] = pgtmp->eta();
                        phigenq4l[igenq4] = pgtmp->phi();
                        egenq4l[igenq4] = pgtmp->energy();
                        ptgenq4f[igenq4] = pgtmp2->pt();
                        etagenq4f[igenq4] = pgtmp2->eta();
                        phigenq4f[igenq4] = pgtmp2->phi();
                        egenq4f[igenq4] = pgtmp2->energy();
                        mothergenq4f[igenq4] = pgtmp2->mother(0)->pdgId();
                        const reco::Candidate* pgtmp3=findLastParticle(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                        if (pgtmp3->mother(0)!=NULL) mmothergenq4f[igenq4] = pgtmp3->mother(0)->pdgId();
                        igenq4+=1;
                    }
                }
            }

            int igenq5=0;
            int sizeq5=5;
            for(size_t ik=0; ik<genParticles->size();ik++)
            {
                if(abs((*genParticles)[ik].pdgId())==5)
                {
                    const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pgtmp=findLastParticle(pgtmp1,5);
                    int goverlap=0;
                    for (int ia=0;ia<igenq5;ia++){
                        if(pgtmp->pt()==ptgenq5l[ia]) goverlap=1;}
                    const reco::Candidate* pgtmp2=findFirstParticle(pgtmp1,5);
                    if(pgtmp->pt()>50&&igenq5<sizeq5&&goverlap==0&&abs(pgtmp2->mother(0)->pdgId())!=24&&abs(pgtmp2->mother(0)->pdgId())!=6){
                        ptgenq5l[igenq5] = pgtmp->pt();
                        etagenq5l[igenq5] = pgtmp->eta();
                        phigenq5l[igenq5] = pgtmp->phi();
                        egenq5l[igenq5] = pgtmp->energy();
                        ptgenq5f[igenq5] = pgtmp2->pt();
                        etagenq5f[igenq5] = pgtmp2->eta();
                        phigenq5f[igenq5] = pgtmp2->phi();
                        egenq5f[igenq5] = pgtmp2->energy();
                        mothergenq5f[igenq5] = pgtmp2->mother(0)->pdgId();
                        const reco::Candidate* pgtmp3=findLastParticle(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                        if (pgtmp3->mother(0)!=NULL) mmothergenq5f[igenq5] = pgtmp3->mother(0)->pdgId();
                        igenq5+=1;
                    }
                }
            }
     //   cout<<"nng   "<<gentop_pt<<"  "<<genantitop_pt<<"  "<<igenw<<"  "<<igenq1<<"  "<<igenq2<<"  "<<igenq3<<"  "<<igenq4<<"  "<<igenq5<<"  "<<igeng<<"  "<<endl;
      
   }//end of MC Info
// *************************End of Gen Level Information******************//

    if(hadronicVs->size()!= 0){
       const reco::Candidate& metCand = metHandle->at(0);
       nLooseMu = loosemus->size();
       nLooseEle = looseels->size();
       nVetoMu = vetomus->size();
       nVetoEle = vetoels->size();
       edm::Handle<reco::VertexCollection> vertices;
       iEvent.getByToken(vtxToken_, vertices);
     //  cout<<"lllllllll"<<nLooseMu<<" "<<nLooseEle<<" "<<nVetoMu<<" "<<nVetoEle<<endl;
       if (nVetoMu != 0) return;
       if (hadronicVs->size()<2) return;
       if (vertices->empty()) return; // skip the event if no PV found

       nVtx = vertices->size();
       reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
       for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx) {
               // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
               // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
               if (  /// !vtx->isFake() &&
                     !(vtx->chi2()==0 && vtx->ndof()==0) 
	             &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	             && fabs(vtx->position().Z())<=24.0) {
                     firstGoodVertex = vtx;
                     break;
                    }           
       }
       if ( firstGoodVertex==vertices->end() ) return; // skip event if there are no good PVs


           // ************************* MET ********************** //
              iEvent.getByToken(metInputToken_ , METs_ );
                addTypeICorr(iEvent);
                for (const pat::MET &met : *METs_) {
                        const float rawPt	 = met.uncorPt();
		    const float rawPhi   = met.uncorPhi();
		    const float rawSumEt = met.uncorSumEt();
                        TVector2 rawMET_;
                        rawMET_.SetMagPhi (rawPt, rawPhi );
                        Double_t rawPx = rawMET_.Px();
                        Double_t rawPy = rawMET_.Py();
                        Double_t rawEt = std::hypot(rawPx,rawPy);
            	    METraw_et = rawEt;
        	   	    METraw_phi = rawPhi;
        	        	    METraw_sumEt = rawSumEt;
                        double pxcorr = rawPx+TypeICorrMap_["corrEx"];
                        double pycorr = rawPy+TypeICorrMap_["corrEy"];
                        double et     = std::hypot(pxcorr,pycorr);
                        double sumEtcorr = rawSumEt+TypeICorrMap_["corrSumEt"];
                        TLorentzVector corrmet; corrmet.SetPxPyPzE(pxcorr,pycorr,0.,et);
            	    useless = sumEtcorr;
            	    useless = rawEt;
            	    MET_et = et;
            	    MET_phi = corrmet.Phi();
            	    MET_sumEt = sumEtcorr;
            	    MET_corrPx = TypeICorrMap_["corrEx"];
            	    MET_corrPy = TypeICorrMap_["corrEy"]; 
                }
           // ***************************************************************** //  
 
       /// For the time being, set these to 1
        triggerWeight=1.0;
        pileupWeight=1.0;
        double targetEvents = targetLumiInvPb_*crossSectionPb_;
        lumiWeight = targetEvents/originalNEvents_;
        met          = metCand.pt();
        metPhi       = metCand.phi();

////////////////////////lep ID  ////////////////////////////////////
   /*     if( leptonicV.daughter(0)->isMuon()||leptonicV.daughter(1)->isMuon()){

                       const pat::Muon *mu1 = abs(leptonicV.daughter(0)->pdgId())==13 ?
                                                  (pat::Muon*)leptonicV.daughter(0):
                                                  (pat::Muon*)leptonicV.daughter(1);
                isHighPt = mu1->isHighPtMuon(vertices->at(0));
                trackIso = mu1->trackIso();
                muchaiso=mu1->pfIsolationR04().sumChargedHadronPt;
                muneuiso=mu1->pfIsolationR04().sumNeutralHadronEt;
                muphoiso=mu1->pfIsolationR04().sumPhotonEt;
                muPU=mu1->pfIsolationR04().sumPUPt;
                    muisolation = (muchaiso+ std::max(0.0,muneuiso+muphoiso-0.5*muPU))/mu1->pt();

}*/
                              if(leptonicVs->size()!= 0){
  
        const reco::Candidate& leptonicV = leptonicVs->at(0);
 
        if( leptonicV.daughter(0)->isElectron()||leptonicV.daughter(1)->isElectron() ) {
                       const pat::Electron *el1 = leptonicV.daughter(0)->isElectron() ?
                                                  (pat::Electron*)leptonicV.daughter(0):
                                                  (pat::Electron*)leptonicV.daughter(1);
                double etaSC1         = el1->superCluster()->eta();
                double d01            = (-1)*el1->gsfTrack()->dxy(firstGoodVertex->position());
                isHEEP = false;
                et = el1->energy()!=0. ? el1->et()/el1->energy()*el1->caloEnergy() : 0.;
                if( et > 35. ) {
                     if( fabs(etaSC1) < 1.4442 ){
                        iso = el1->dr03EcalRecHitSumEt() + el1->dr03HcalDepth1TowerSumEt();
                        isoCut = 2 + 0.03*et + 0.28*fastJetRho;
                        if( el1->ecalDriven() == 1 && dEtaInSeed( el1 ) < 0.004 && el1->deltaPhiSuperClusterTrackAtVtx() < 0.06 &&
                         el1->hadronicOverEm() < (1./el1->superCluster()->energy()+0.05) &&
                         (el1->full5x5_e2x5Max()/el1->full5x5_e5x5() > 0.94 || el1->full5x5_e1x5()/el1->full5x5_e5x5() > 0.83) &&
                         el1->dr03TkSumPt() < 5. && el1->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 && iso < isoCut && fabs(d01) < 0.02 ) isHEEP = true;
                     }
                     if( fabs(etaSC1) > 1.566 && fabs(etaSC1) < 2.5 ){
                        iso = el1->dr03EcalRecHitSumEt() + el1->dr03HcalDepth1TowerSumEt();
                        if( et <= 50 )
                                isoCut = 2.5 + 0.28*fastJetRho;
                        else
                                isoCut = 2.5+0.03*(et-50.) + 0.28*fastJetRho;
                        if( el1->ecalDriven() == 1 && dEtaInSeed( el1 ) < 0.006 && el1->deltaPhiSuperClusterTrackAtVtx() < 0.06 &&
                         el1->hadronicOverEm() < (5./el1->superCluster()->energy()+0.05) && el1->full5x5_sigmaIetaIeta() < 0.03 &&
                         el1->dr03TkSumPt() < 5. && el1->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&//numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
                         iso < isoCut && fabs(d01) < 0.05 ) isHEEP = true;
                     }
        }
        }
                                     }
                                else{
                              isHEEP = false;
                                    }
//cout<<"isHEEP:"<<isHEEP<<endl;
////////////////////////JEC for AK8/////////////////////////////////
        reco::Candidate::LorentzVector uncorrPrunedJet;
        bool doPuppi  = iEvent.getByToken(puppijetInputToken_, puppijets_ );
        if( doPuppi ){//1
         for(size_t ij=0; ij<puppijets_->size();ij++){
           corr_AK8puppi[ij] = 1;
           corr_AK8puppiSD[ij] = 1;
           const pat::Jet& hadronicVa = puppijets_->at(ij);
           reco::Candidate::LorentzVector uncorrJet;
           if(not isJEC_) doCorrOnTheFly_ = false;
           if( doCorrOnTheFly_ ){
              uncorrJet = hadronicVa.correctedP4(0);
              jecAK8puppi_->setJetEta( uncorrJet.eta()          );
              jecAK8puppi_->setJetPt ( uncorrJet.pt()           );
              jecAK8puppi_->setJetE  ( uncorrJet.energy()       );
              jecAK8puppi_->setRho   (fastJetRho);
              jecAK8puppi_->setNPV   (nVtx);
              jecAK8puppi_->setJetA  (hadronicVa.jetArea());
              corr_AK8puppi[ij] = jecAK8puppi_->getCorrection();
              jecAK8puppiGroomed_->setJetEta( uncorrJet.eta()          );
              jecAK8puppiGroomed_->setJetPt ( uncorrJet.pt()           );
              jecAK8puppiGroomed_->setJetE  ( uncorrJet.energy()       );
              jecAK8puppiGroomed_->setRho   (fastJetRho);
              jecAK8puppiGroomed_->setNPV   (nVtx);
              jecAK8puppiGroomed_->setJetA  (hadronicVa.jetArea());
              corr_AK8puppiSD[ij] = jecAK8puppiGroomed_->getCorrection();
           }
           else{uncorrJet = hadronicVa.p4();}
           if(ij<8){
              jetAK8puppi_pt1[ij] = corr_AK8puppi[ij]*uncorrJet.pt();
              jetAK8puppi_mass1[ij] = corr_AK8puppi[ij]*uncorrJet.mass();
              jetAK8puppi_eta1[ij] = uncorrJet.eta();
              jetAK8puppi_jec1[ij] = corr_AK8puppi[ij];
              jetAK8puppiSD_jec1[ij] = corr_AK8puppiSD[ij];
           }
         }
         int usenumber3 = -1; double pt_larger=0;
         int numvhad = puppijets_->size();
         for( int inum = 0; inum< numvhad; inum++){
           const pat::Jet& Vpuppi = puppijets_->at(inum);
      //cout<<looseJetID(Vpuppi)<<endl;
      //     continue; 
           if(looseJetID(Vpuppi)<1) continue;     
           if(jetAK8puppi_pt1[inum] > pt_larger && fabs(jetAK8puppi_eta1[inum])<2.4 && inum<8) {pt_larger = jetAK8puppi_pt1[inum]; usenumber3 = inum; continue;}
        }
       if (usenumber3>-1) {//2
        const pat::Jet& hadronicVpuppi = puppijets_->at(usenumber3);
                // Decorrelated DeepAK8
                JetHelper jet_helper(&hadronicVpuppi);
                //jet_helper.setSubjets(*hadronicVSoftDrop, 0.8); // jetR=0.8, only for 80X 
                const auto& nnpreds = fatjetNN_->predict(jet_helper);
                FatJetNNHelper nn(nnpreds);
                jetAK8puppi_dnnTop       = nn.get_binarized_score_top();
                jetAK8puppi_dnnW         = nn.get_binarized_score_w();
                jetAK8puppi_dnnH4q       = nn.get_binarized_score_h4q();
                jetAK8puppi_dnnZ         = nn.get_binarized_score_z();
                jetAK8puppi_dnnZbb       = nn.get_binarized_score_zbb();
                jetAK8puppi_dnnHbb       = nn.get_binarized_score_hbb();
                jetAK8puppi_dnnqcd       = nn.get_raw_score_qcd();
                jetAK8puppi_dnntop       = nn.get_raw_score_top();
                jetAK8puppi_dnnw         = nn.get_raw_score_w();
                jetAK8puppi_dnnz         = nn.get_raw_score_z();
                jetAK8puppi_dnnzbb         = nn.get_raw_score_zbb();
                jetAK8puppi_dnnhbb         = nn.get_raw_score_hbb();
                jetAK8puppi_dnnh4q         = nn.get_raw_score_h4q();
                // Decorrelated DeepAK8
                const auto& mdpreds = decorrNN_->predict(jet_helper);
                FatJetNNHelper md(mdpreds);
                jetAK8puppi_dnnDecorrTop       = md.get_binarized_score_top();
                jetAK8puppi_dnnDecorrW         = md.get_binarized_score_w();
                jetAK8puppi_dnnDecorrH4q       = md.get_binarized_score_h4q();
                jetAK8puppi_dnnDecorrZ         = md.get_binarized_score_z();
                jetAK8puppi_dnnDecorrZbb       = md.get_binarized_score_zbb();
                jetAK8puppi_dnnDecorrHbb       = md.get_binarized_score_hbb();
                jetAK8puppi_dnnDecorrbb        = md.get_flavor_score_bb();
                jetAK8puppi_dnnDecorrcc        = md.get_flavor_score_cc();
                jetAK8puppi_dnnDecorrbbnog        = md.get_flavor_score_bb_no_gluon();
                jetAK8puppi_dnnDecorrbbnog        = md.get_flavor_score_cc_no_gluon();
                jetAK8puppi_dnnDecorrqcd       = md.get_raw_score_qcd();
                jetAK8puppi_dnnDecorrtop       = md.get_raw_score_top();
                jetAK8puppi_dnnDecorrw         = md.get_raw_score_w();
                jetAK8puppi_dnnDecorrz         = md.get_raw_score_z();
                jetAK8puppi_dnnDecorrzbb         = md.get_raw_score_zbb();
                jetAK8puppi_dnnDecorrhbb         = md.get_raw_score_hbb();
                jetAK8puppi_dnnDecorrh4q         = md.get_raw_score_h4q();
                IDLoose = looseJetID(hadronicVpuppi);
                IDTight = tightJetID(hadronicVpuppi);
                jetAK8puppi_ptJEC       = jetAK8puppi_pt1[usenumber3]; // unpruned corrected jet pt
                jetAK8puppi_eta     = jetAK8puppi_eta1[usenumber3]; // unpruned (w/o jec) jet eta
                jetAK8puppi_phi      = hadronicVpuppi.phi(); // unpruned (w/o jec) jet phi
                jetAK8puppi_tau1         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau1");
                jetAK8puppi_tau2         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau2");
                jetAK8puppi_tau3         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau21        = jetAK8puppi_tau2/jetAK8puppi_tau1;
                jetAK8puppi_tau4         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau4");
                jetAK8puppi_tau42        = jetAK8puppi_tau4/jetAK8puppi_tau2;
                jetAK8puppi_sd       =  hadronicVpuppi.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass
                jetAK8puppi_sdJEC  =corr_AK8puppiSD[usenumber3]*jetAK8puppi_sd;
                Double_t gencorrect=1.0;
                Double_t recocorrect_0eta1p3=1.0;
                Double_t recocorrect_1p3eta2p5=1.0;
                gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC*0.08,-1.2);
                recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC+3.449e-07*pow(jetAK8puppi_ptJEC,2)-2.681e-10*pow(jetAK8puppi_ptJEC,3)+8.674e-14*pow(jetAK8puppi_ptJEC,4)-1.001e-17*pow(jetAK8puppi_ptJEC,5);
                recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC+8.37e-07*pow(jetAK8puppi_ptJEC,2)-5.204e-10*pow(jetAK8puppi_ptJEC,3)+1.454e-13*pow(jetAK8puppi_ptJEC,4)-1.504e-17*pow(jetAK8puppi_ptJEC,5);
                if (fabs(jetAK8puppi_eta)<=1.3){jetAK8puppi_sdcorr=jetAK8puppi_sd*gencorrect*recocorrect_0eta1p3;}
                else if (fabs(jetAK8puppi_eta)<2.5 && fabs(jetAK8puppi_eta)>1.3){jetAK8puppi_sdcorr=jetAK8puppi_sd*gencorrect*recocorrect_1p3eta2p5;}
           int usenumber2 = -1; double pt_larger2=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger2 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber3 && inum<8) {pt_larger2 = jetAK8puppi_pt1[inum]; usenumber2 = inum; continue;}
           }
           
          if(usenumber2>-1)  {
               const pat::Jet& hadronicVpuppi_2 = puppijets_->at(usenumber2);
                // DeepAK8
                JetHelper jet_helper_2(&hadronicVpuppi_2);
                // jet_helper_2.setSubjets(*hadronicVSoftDrop, 0.8); // jetR=0.8
                const auto& nnpreds_2 = fatjetNN_->predict(jet_helper_2);
                FatJetNNHelper nn_2(nnpreds_2);
                jetAK8puppi_dnnTop_2       = nn_2.get_binarized_score_top();
                jetAK8puppi_dnnW_2         = nn_2.get_binarized_score_w();
                jetAK8puppi_dnnH4q_2       = nn_2.get_binarized_score_h4q();
                jetAK8puppi_dnnZ_2         = nn_2.get_binarized_score_z();
                jetAK8puppi_dnnZbb_2       = nn_2.get_binarized_score_zbb();
                jetAK8puppi_dnnHbb_2       = nn_2.get_binarized_score_hbb();
                jetAK8puppi_dnnqcd_2       = nn_2.get_raw_score_qcd();
                jetAK8puppi_dnntop_2       = nn_2.get_raw_score_top();
                jetAK8puppi_dnnw_2         = nn_2.get_raw_score_w();
                jetAK8puppi_dnnz_2         = nn_2.get_raw_score_z();
                jetAK8puppi_dnnzbb_2         = nn_2.get_raw_score_zbb();
                jetAK8puppi_dnnhbb_2         = nn_2.get_raw_score_hbb();
                jetAK8puppi_dnnh4q_2         = nn_2.get_raw_score_h4q();
                // Decorrelated DeepAK8
                const auto& mdpreds_2 = decorrNN_->predict(jet_helper_2);
                FatJetNNHelper md_2(mdpreds_2);
                jetAK8puppi_dnnDecorrTop_2       = md_2.get_binarized_score_top();
                jetAK8puppi_dnnDecorrW_2         = md_2.get_binarized_score_w();
                jetAK8puppi_dnnDecorrH4q_2       = md_2.get_binarized_score_h4q();
                jetAK8puppi_dnnDecorrZ_2         = md_2.get_binarized_score_z();
                jetAK8puppi_dnnDecorrZbb_2       = md_2.get_binarized_score_zbb();
                jetAK8puppi_dnnDecorrHbb_2       = md_2.get_binarized_score_hbb();
                jetAK8puppi_dnnDecorrbb_2        = md_2.get_flavor_score_bb();
                jetAK8puppi_dnnDecorrcc_2        = md_2.get_flavor_score_cc();
                jetAK8puppi_dnnDecorrbbnog_2        = md_2.get_flavor_score_bb_no_gluon();
                jetAK8puppi_dnnDecorrbbnog_2        = md_2.get_flavor_score_cc_no_gluon();
                jetAK8puppi_dnnDecorrqcd_2       = md_2.get_raw_score_qcd();
                jetAK8puppi_dnnDecorrtop_2       = md_2.get_raw_score_top();
                jetAK8puppi_dnnDecorrw_2         = md_2.get_raw_score_w();
                jetAK8puppi_dnnDecorrz_2         = md_2.get_raw_score_z();
                jetAK8puppi_dnnDecorrzbb_2         = md_2.get_raw_score_zbb();
                jetAK8puppi_dnnDecorrhbb_2         = md_2.get_raw_score_hbb();
                jetAK8puppi_dnnDecorrh4q_2         = md_2.get_raw_score_h4q();

                IDLoose_2 = looseJetID(hadronicVpuppi_2);
                IDTight_2 = tightJetID(hadronicVpuppi_2);
               jetAK8puppi_ptJEC_2       = jetAK8puppi_pt1[usenumber2]; // unpruned corrected jet pt
               jetAK8puppi_eta_2     = jetAK8puppi_eta1[usenumber2]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_2      = hadronicVpuppi_2.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau1");
               jetAK8puppi_tau2_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau2");
               jetAK8puppi_tau3_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau3");
               jetAK8puppi_tau21_2        = jetAK8puppi_tau2_2/jetAK8puppi_tau1_2;
               jetAK8puppi_tau4_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau4");
               jetAK8puppi_tau42_2        = jetAK8puppi_tau4_2/jetAK8puppi_tau2_2;
               
               jetAK8puppi_sd_2       =  hadronicVpuppi_2.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_2  =corr_AK8puppiSD[usenumber2]*jetAK8puppi_sd_2;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_2*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_2+3.449e-07*pow(jetAK8puppi_ptJEC_2,2)-2.681e-10*pow(jetAK8puppi_ptJEC_2,3)+8.674e-14*pow(jetAK8puppi_ptJEC_2,4)-1.001e-17*pow(jetAK8puppi_ptJEC_2,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_2+8.37e-07*pow(jetAK8puppi_ptJEC_2,2)-5.204e-10*pow(jetAK8puppi_ptJEC_2,3)+1.454e-13*pow(jetAK8puppi_ptJEC_2,4)-1.504e-17*pow(jetAK8puppi_ptJEC_2,5);
               if (fabs(jetAK8puppi_eta_2)<=1.3){jetAK8puppi_sdcorr_2=jetAK8puppi_sd_2*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_2)<2.5 && fabs(jetAK8puppi_eta_2)>1.3){jetAK8puppi_sdcorr_2=jetAK8puppi_sd_2*gencorrect*recocorrect_1p3eta2p5;}
           }
           


           int usenumber1 = -1; double pt_larger3=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger3 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber3 && inum != usenumber2 && inum<8) {pt_larger3 = jetAK8puppi_pt1[inum]; usenumber1 = inum; continue;}
           }
           
           
          if(usenumber1>-1)  {
               const pat::Jet& hadronicVpuppi_3 = puppijets_->at(usenumber1);
                // DeepAK8
                JetHelper jet_helper_3(&hadronicVpuppi_3);
                // jet_helper_3.setSubjets(*hadronicVSoftDrop, 0.8); // jetR=0.8
                const auto& nnpreds_3 = fatjetNN_->predict(jet_helper_3);
                FatJetNNHelper nn_3(nnpreds_3);
                jetAK8puppi_dnnTop_3       = nn_3.get_binarized_score_top();
                jetAK8puppi_dnnW_3         = nn_3.get_binarized_score_w();
                jetAK8puppi_dnnH4q_3       = nn_3.get_binarized_score_h4q();
                jetAK8puppi_dnnZ_3         = nn_3.get_binarized_score_z();
                jetAK8puppi_dnnZbb_3       = nn_3.get_binarized_score_zbb();
                jetAK8puppi_dnnHbb_3       = nn_3.get_binarized_score_hbb();
                jetAK8puppi_dnnqcd_3       = nn_3.get_raw_score_qcd();
                jetAK8puppi_dnntop_3       = nn_3.get_raw_score_top();
                jetAK8puppi_dnnw_3         = nn_3.get_raw_score_w();
                jetAK8puppi_dnnz_3         = nn_3.get_raw_score_z();
                jetAK8puppi_dnnzbb_3         = nn_3.get_raw_score_zbb();
                jetAK8puppi_dnnhbb_3         = nn_3.get_raw_score_hbb();
                jetAK8puppi_dnnh4q_3         = nn_3.get_raw_score_h4q();
                // Decorrelated DeepAK8
                const auto& mdpreds_3 = decorrNN_->predict(jet_helper_3);
                FatJetNNHelper md_3(mdpreds_3);
                jetAK8puppi_dnnDecorrTop_3       = md_3.get_binarized_score_top();
                jetAK8puppi_dnnDecorrW_3         = md_3.get_binarized_score_w();
                jetAK8puppi_dnnDecorrH4q_3       = md_3.get_binarized_score_h4q();
                jetAK8puppi_dnnDecorrZ_3         = md_3.get_binarized_score_z();
                jetAK8puppi_dnnDecorrZbb_3       = md_3.get_binarized_score_zbb();
                jetAK8puppi_dnnDecorrHbb_3       = md_3.get_binarized_score_hbb();
                jetAK8puppi_dnnDecorrbb_3        = md_3.get_flavor_score_bb();
                jetAK8puppi_dnnDecorrcc_3        = md_3.get_flavor_score_cc();
                jetAK8puppi_dnnDecorrbbnog_3        = md_3.get_flavor_score_bb_no_gluon();
                jetAK8puppi_dnnDecorrbbnog_3        = md_3.get_flavor_score_cc_no_gluon();
                jetAK8puppi_dnnDecorrqcd_3       = md_3.get_raw_score_qcd();
                jetAK8puppi_dnnDecorrtop_3       = md_3.get_raw_score_top();
                jetAK8puppi_dnnDecorrw_3         = md_3.get_raw_score_w();
                jetAK8puppi_dnnDecorrz_3         = md_3.get_raw_score_z();
                jetAK8puppi_dnnDecorrzbb_3         = md_3.get_raw_score_zbb();
                jetAK8puppi_dnnDecorrhbb_3         = md_3.get_raw_score_hbb();
                jetAK8puppi_dnnDecorrh4q_3         = md_3.get_raw_score_h4q();

                IDLoose_3 = looseJetID(hadronicVpuppi_3);
                IDTight_3 = tightJetID(hadronicVpuppi_3);
               jetAK8puppi_ptJEC_3       = jetAK8puppi_pt1[usenumber1]; // unpruned corrected jet pt
               jetAK8puppi_eta_3     = jetAK8puppi_eta1[usenumber1]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_3      = hadronicVpuppi_3.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau1");
               jetAK8puppi_tau2_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau2");
               jetAK8puppi_tau3_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau3");
               jetAK8puppi_tau21_3        = jetAK8puppi_tau2_3/jetAK8puppi_tau1_3;
               jetAK8puppi_tau4_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau4");
               jetAK8puppi_tau42_3        = jetAK8puppi_tau4_3/jetAK8puppi_tau2_3;
               
               jetAK8puppi_sd_3       =  hadronicVpuppi_3.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_3  =corr_AK8puppiSD[usenumber1]*jetAK8puppi_sd_3;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_3*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_3+3.449e-07*pow(jetAK8puppi_ptJEC_3,2)-2.681e-10*pow(jetAK8puppi_ptJEC_3,3)+8.674e-14*pow(jetAK8puppi_ptJEC_3,4)-1.001e-17*pow(jetAK8puppi_ptJEC_3,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_3+8.37e-07*pow(jetAK8puppi_ptJEC_3,2)-5.204e-10*pow(jetAK8puppi_ptJEC_3,3)+1.454e-13*pow(jetAK8puppi_ptJEC_3,4)-1.504e-17*pow(jetAK8puppi_ptJEC_3,5);
               if (fabs(jetAK8puppi_eta_3)<=1.3){jetAK8puppi_sdcorr_3=jetAK8puppi_sd_3*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_3)<2.5 && fabs(jetAK8puppi_eta_3)>1.3){jetAK8puppi_sdcorr_3=jetAK8puppi_sd_3*gencorrect*recocorrect_1p3eta2p5;}
           }
           

           int usenumber4 = -1; double pt_larger4=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger4 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber3 && inum != usenumber2 && inum != usenumber1 && inum<8) {pt_larger4 = jetAK8puppi_pt1[inum]; usenumber4 = inum; continue;}
           }


          if(usenumber4>-1)  {
                const pat::Jet& hadronicVpuppi_4 = puppijets_->at(usenumber4);
                // DeepAK8
                JetHelper jet_helper_4(&hadronicVpuppi_4);
                //  jet_helper_4.setSubjets(*hadronicVSoftDrop, 0.8); // jetR=0.8
                const auto& nnpreds_4 = fatjetNN_->predict(jet_helper_4);
                FatJetNNHelper nn_4(nnpreds_4);
                jetAK8puppi_dnnTop_4       = nn_4.get_binarized_score_top();
                jetAK8puppi_dnnW_4         = nn_4.get_binarized_score_w();
                jetAK8puppi_dnnH4q_4       = nn_4.get_binarized_score_h4q();
                jetAK8puppi_dnnZ_4         = nn_4.get_binarized_score_z();
                jetAK8puppi_dnnZbb_4       = nn_4.get_binarized_score_zbb();
                jetAK8puppi_dnnHbb_4       = nn_4.get_binarized_score_hbb();
                jetAK8puppi_dnnqcd_4       = nn_4.get_raw_score_qcd();
                jetAK8puppi_dnntop_4       = nn_4.get_raw_score_top();
                jetAK8puppi_dnnw_4         = nn_4.get_raw_score_w();
                jetAK8puppi_dnnz_4         = nn_4.get_raw_score_z();
                jetAK8puppi_dnnzbb_4         = nn_4.get_raw_score_zbb();
                jetAK8puppi_dnnhbb_4         = nn_4.get_raw_score_hbb();
                jetAK8puppi_dnnh4q_4         = nn_4.get_raw_score_h4q();
                // Decorrelated DeepAK8
                const auto& mdpreds_4 = decorrNN_->predict(jet_helper_4);
                FatJetNNHelper md_4(mdpreds_4);
                jetAK8puppi_dnnDecorrTop_4       = md_4.get_binarized_score_top();
                jetAK8puppi_dnnDecorrW_4         = md_4.get_binarized_score_w();
                jetAK8puppi_dnnDecorrH4q_4       = md_4.get_binarized_score_h4q();
                jetAK8puppi_dnnDecorrZ_4         = md_4.get_binarized_score_z();
                jetAK8puppi_dnnDecorrZbb_4       = md_4.get_binarized_score_zbb();
                jetAK8puppi_dnnDecorrHbb_4       = md_4.get_binarized_score_hbb();
                jetAK8puppi_dnnDecorrbb_4        = md_4.get_flavor_score_bb();
                jetAK8puppi_dnnDecorrcc_4        = md_4.get_flavor_score_cc();
                jetAK8puppi_dnnDecorrbbnog_4        = md_4.get_flavor_score_bb_no_gluon();
                jetAK8puppi_dnnDecorrbbnog_4        = md_4.get_flavor_score_cc_no_gluon();
                jetAK8puppi_dnnDecorrqcd_4       = md_4.get_raw_score_qcd();
                jetAK8puppi_dnnDecorrtop_4       = md_4.get_raw_score_top();
                jetAK8puppi_dnnDecorrw_4         = md_4.get_raw_score_w();
                jetAK8puppi_dnnDecorrz_4         = md_4.get_raw_score_z();
                jetAK8puppi_dnnDecorrzbb_4         = md_4.get_raw_score_zbb();
                jetAK8puppi_dnnDecorrhbb_4         = md_4.get_raw_score_hbb();
                jetAK8puppi_dnnDecorrh4q_4         = md_4.get_raw_score_h4q();

                IDLoose_4 = looseJetID(hadronicVpuppi_4);
                IDTight_4 = tightJetID(hadronicVpuppi_4);
               jetAK8puppi_ptJEC_4       = jetAK8puppi_pt1[usenumber4]; // unpruned corrected jet pt
               jetAK8puppi_eta_4     = jetAK8puppi_eta1[usenumber4]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_4      = hadronicVpuppi_4.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_4         = hadronicVpuppi_4.userFloat("NjettinessAK8Puppi:tau1");
               jetAK8puppi_tau2_4         = hadronicVpuppi_4.userFloat("NjettinessAK8Puppi:tau2");
               jetAK8puppi_tau3_4         = hadronicVpuppi_4.userFloat("NjettinessAK8Puppi:tau3");
               jetAK8puppi_tau21_4        = jetAK8puppi_tau2_4/jetAK8puppi_tau1_4;
               jetAK8puppi_tau4_4         = hadronicVpuppi_4.userFloat("NjettinessAK8Puppi:tau4");
               jetAK8puppi_tau42_4        = jetAK8puppi_tau4_4/jetAK8puppi_tau2_4;

               jetAK8puppi_sd_4       =  hadronicVpuppi_4.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_4  =corr_AK8puppiSD[usenumber4]*jetAK8puppi_sd_4;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_4*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_4+3.449e-07*pow(jetAK8puppi_ptJEC_4,2)-2.681e-10*pow(jetAK8puppi_ptJEC_4,3)+8.674e-14*pow(jetAK8puppi_ptJEC_4,4)-1.001e-17*pow(jetAK8puppi_ptJEC_4,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_4+8.37e-07*pow(jetAK8puppi_ptJEC_4,2)-5.204e-10*pow(jetAK8puppi_ptJEC_4,3)+1.454e-13*pow(jetAK8puppi_ptJEC_4,4)-1.504e-17*pow(jetAK8puppi_ptJEC_4,5);
               if (fabs(jetAK8puppi_eta_4)<=1.3){jetAK8puppi_sdcorr_4=jetAK8puppi_sd_4*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_4)<2.5 && fabs(jetAK8puppi_eta_4)>1.3){jetAK8puppi_sdcorr_4=jetAK8puppi_sd_4*gencorrect*recocorrect_1p3eta2p5;}
           }


           int usenumber5 = -1; double pt_larger5=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger5 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber4 && inum != usenumber3 && inum != usenumber2 && inum != usenumber1 && inum<8) {pt_larger5 = jetAK8puppi_pt1[inum]; usenumber5 = inum; continue;}
           }


          if(usenumber5>-1)  {
               const pat::Jet& hadronicVpuppi_5 = puppijets_->at(usenumber5);
               jetAK8puppi_ptJEC_5       = jetAK8puppi_pt1[usenumber5]; // unpruned corrected jet pt
               jetAK8puppi_eta_5     = jetAK8puppi_eta1[usenumber5]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_5      = hadronicVpuppi_5.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_5         = hadronicVpuppi_5.userFloat("NjettinessAK8Puppi:tau1");
               jetAK8puppi_tau2_5         = hadronicVpuppi_5.userFloat("NjettinessAK8Puppi:tau2");
               jetAK8puppi_tau3_5         = hadronicVpuppi_5.userFloat("NjettinessAK8Puppi:tau3");
               jetAK8puppi_tau21_5        = jetAK8puppi_tau2_5/jetAK8puppi_tau1_5;
               jetAK8puppi_tau4_5         = hadronicVpuppi_5.userFloat("NjettinessAK8Puppi:tau4");
               jetAK8puppi_tau42_5        = jetAK8puppi_tau4_5/jetAK8puppi_tau2_5;

               jetAK8puppi_sd_5       =  hadronicVpuppi_5.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_5  =corr_AK8puppiSD[usenumber5]*jetAK8puppi_sd_5;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_5*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_5+3.449e-07*pow(jetAK8puppi_ptJEC_5,2)-2.681e-10*pow(jetAK8puppi_ptJEC_5,3)+8.674e-14*pow(jetAK8puppi_ptJEC_5,4)-1.001e-17*pow(jetAK8puppi_ptJEC_5,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_5+8.37e-07*pow(jetAK8puppi_ptJEC_5,2)-5.204e-10*pow(jetAK8puppi_ptJEC_5,3)+1.454e-13*pow(jetAK8puppi_ptJEC_5,4)-1.504e-17*pow(jetAK8puppi_ptJEC_5,5);
               if (fabs(jetAK8puppi_eta_5)<=1.3){jetAK8puppi_sdcorr_5=jetAK8puppi_sd_5*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_5)<2.5 && fabs(jetAK8puppi_eta_5)>1.3){jetAK8puppi_sdcorr_5=jetAK8puppi_sd_5*gencorrect*recocorrect_1p3eta2p5;}
           }



           int usenumber6 = -1; double pt_larger6=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger6 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber5 && inum != usenumber4 && inum != usenumber3 && inum != usenumber2 && inum != usenumber1 && inum<8) {pt_larger6 = jetAK8puppi_pt1[inum]; usenumber6 = inum; continue;}
           }


          if(usenumber6>-1)  {
               const pat::Jet& hadronicVpuppi_6 = puppijets_->at(usenumber6);
               jetAK8puppi_ptJEC_6       = jetAK8puppi_pt1[usenumber6]; // unpruned corrected jet pt
               jetAK8puppi_eta_6     = jetAK8puppi_eta1[usenumber6]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_6      = hadronicVpuppi_6.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_6         = hadronicVpuppi_6.userFloat("NjettinessAK8Puppi:tau1");
               jetAK8puppi_tau2_6         = hadronicVpuppi_6.userFloat("NjettinessAK8Puppi:tau2");
               jetAK8puppi_tau3_6         = hadronicVpuppi_6.userFloat("NjettinessAK8Puppi:tau3");
               jetAK8puppi_tau21_6        = jetAK8puppi_tau2_6/jetAK8puppi_tau1_6;
               jetAK8puppi_tau4_6         = hadronicVpuppi_6.userFloat("NjettinessAK8Puppi:tau4");
               jetAK8puppi_tau42_6        = jetAK8puppi_tau4_6/jetAK8puppi_tau2_6;

               jetAK8puppi_sd_6       =  hadronicVpuppi_6.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_6  =corr_AK8puppiSD[usenumber6]*jetAK8puppi_sd_6;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_6*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_6+3.449e-07*pow(jetAK8puppi_ptJEC_6,2)-2.681e-10*pow(jetAK8puppi_ptJEC_6,3)+8.674e-14*pow(jetAK8puppi_ptJEC_6,4)-1.001e-17*pow(jetAK8puppi_ptJEC_6,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_6+8.37e-07*pow(jetAK8puppi_ptJEC_6,2)-5.204e-10*pow(jetAK8puppi_ptJEC_6,3)+1.454e-13*pow(jetAK8puppi_ptJEC_6,4)-1.504e-17*pow(jetAK8puppi_ptJEC_6,5);
               if (fabs(jetAK8puppi_eta_6)<=1.3){jetAK8puppi_sdcorr_6=jetAK8puppi_sd_6*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_6)<2.5 && fabs(jetAK8puppi_eta_6)>1.3){jetAK8puppi_sdcorr_6=jetAK8puppi_sd_6*gencorrect*recocorrect_1p3eta2p5;}
           }


           int usenumber7 = -1; double pt_larger7=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger7 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber6 && inum != usenumber5 && inum != usenumber4 && inum != usenumber3 && inum != usenumber2 && inum != usenumber1 && inum<8) {pt_larger7 = jetAK8puppi_pt1[inum]; usenumber7 = inum; continue;}
           }


          if(usenumber7>-1)  {
               const pat::Jet& hadronicVpuppi_7 = puppijets_->at(usenumber7);
               jetAK8puppi_ptJEC_7       = jetAK8puppi_pt1[usenumber7]; // unpruned corrected jet pt
               jetAK8puppi_eta_7     = jetAK8puppi_eta1[usenumber7]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_7      = hadronicVpuppi_7.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_7         = hadronicVpuppi_7.userFloat("NjettinessAK8Puppi:tau1");
               jetAK8puppi_tau2_7         = hadronicVpuppi_7.userFloat("NjettinessAK8Puppi:tau2");
               jetAK8puppi_tau3_7         = hadronicVpuppi_7.userFloat("NjettinessAK8Puppi:tau3");
               jetAK8puppi_tau21_7        = jetAK8puppi_tau2_7/jetAK8puppi_tau1_7;
               jetAK8puppi_tau4_7         = hadronicVpuppi_7.userFloat("NjettinessAK8Puppi:tau4");
               jetAK8puppi_tau42_7        = jetAK8puppi_tau4_7/jetAK8puppi_tau2_7;

               jetAK8puppi_sd_7       =  hadronicVpuppi_7.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_7  =corr_AK8puppiSD[usenumber7]*jetAK8puppi_sd_7;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_7*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_7+3.449e-07*pow(jetAK8puppi_ptJEC_7,2)-2.681e-10*pow(jetAK8puppi_ptJEC_7,3)+8.674e-14*pow(jetAK8puppi_ptJEC_7,4)-1.001e-17*pow(jetAK8puppi_ptJEC_7,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_7+8.37e-07*pow(jetAK8puppi_ptJEC_7,2)-5.204e-10*pow(jetAK8puppi_ptJEC_7,3)+1.454e-13*pow(jetAK8puppi_ptJEC_7,4)-1.504e-17*pow(jetAK8puppi_ptJEC_7,5);
               if (fabs(jetAK8puppi_eta_7)<=1.3){jetAK8puppi_sdcorr_7=jetAK8puppi_sd_7*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_7)<2.5 && fabs(jetAK8puppi_eta_7)>1.3){jetAK8puppi_sdcorr_7=jetAK8puppi_sd_7*gencorrect*recocorrect_1p3eta2p5;}
           }


           int usenumber8 = -1; double pt_larger8=0;
           for( int inum = 0; inum< numvhad; inum++){
               const pat::Jet& Vpuppi = puppijets_->at(inum);
               if(looseJetID(Vpuppi)<1) continue;
               if(jetAK8puppi_pt1[inum] > pt_larger8 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber7 && inum != usenumber6 && inum != usenumber5 && inum != usenumber4 && inum != usenumber3 && inum != usenumber2 && inum != usenumber1 && inum<8) {pt_larger8 = jetAK8puppi_pt1[inum]; usenumber8 = inum; continue;}
           }


          if(usenumber8>-1)  {
               const pat::Jet& hadronicVpuppi_8 = puppijets_->at(usenumber8);
               jetAK8puppi_ptJEC_8       = jetAK8puppi_pt1[usenumber8]; // unpruned corrected jet pt
               jetAK8puppi_eta_8     = jetAK8puppi_eta1[usenumber8]; // unpruned (w/o jec) jet eta
               jetAK8puppi_phi_8      = hadronicVpuppi_8.phi(); // unpruned (w/o jec) jet phi
               jetAK8puppi_tau1_8         = hadronicVpuppi_8.userFloat("NjettinessAK8Puppi:tau1");
               jetAK8puppi_tau2_8         = hadronicVpuppi_8.userFloat("NjettinessAK8Puppi:tau2");
               jetAK8puppi_tau3_8         = hadronicVpuppi_8.userFloat("NjettinessAK8Puppi:tau3");
               jetAK8puppi_tau21_8        = jetAK8puppi_tau2_8/jetAK8puppi_tau1_8;
               jetAK8puppi_tau4_8         = hadronicVpuppi_8.userFloat("NjettinessAK8Puppi:tau4");
               jetAK8puppi_tau42_8        = jetAK8puppi_tau4_8/jetAK8puppi_tau2_8;

               jetAK8puppi_sd_8       =  hadronicVpuppi_8.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass
               jetAK8puppi_sdJEC_8  =corr_AK8puppiSD[usenumber8]*jetAK8puppi_sd_8;
               Double_t gencorrect=1.0;
               Double_t recocorrect_0eta1p3=1.0;
               Double_t recocorrect_1p3eta2p5=1.0;
               gencorrect=1.006-1.062*pow(jetAK8puppi_ptJEC_8*0.08,-1.2);
               recocorrect_0eta1p3=1.093-1.501e-04*jetAK8puppi_ptJEC_8+3.449e-07*pow(jetAK8puppi_ptJEC_8,2)-2.681e-10*pow(jetAK8puppi_ptJEC_8,3)+8.674e-14*pow(jetAK8puppi_ptJEC_8,4)-1.001e-17*pow(jetAK8puppi_ptJEC_8,5);
               recocorrect_1p3eta2p5=1.272-5.72e-04*jetAK8puppi_ptJEC_8+8.37e-07*pow(jetAK8puppi_ptJEC_8,2)-5.204e-10*pow(jetAK8puppi_ptJEC_8,3)+1.454e-13*pow(jetAK8puppi_ptJEC_8,4)-1.504e-17*pow(jetAK8puppi_ptJEC_8,5);
               if (fabs(jetAK8puppi_eta_8)<=1.3){jetAK8puppi_sdcorr_8=jetAK8puppi_sd_8*gencorrect*recocorrect_0eta1p3;}
               else if (fabs(jetAK8puppi_eta_8)<2.5 && fabs(jetAK8puppi_eta_8)>1.3){jetAK8puppi_sdcorr_8=jetAK8puppi_sd_8*gencorrect*recocorrect_1p3eta2p5;}
           }
	int nak4 = 0;
        double tj1=-10.0, tj2=-10.0; 
        for (size_t ik=0; ik<ak4jets->size();ik++)
         {//3
            double corr = 1;
            reco::Candidate::LorentzVector uncorrJet;
             if( doCorrOnTheFly_ ){
	    uncorrJet = (*ak4jets)[ik].correctedP4(0);
              jecAK4_->setJetEta( uncorrJet.eta() );
              jecAK4_->setJetPt ( uncorrJet.pt() );
              jecAK4_->setJetE ( uncorrJet.energy() );
              jecAK4_->setRho ( fastJetRho );
              jecAK4_->setNPV ( vertices->size() );
              jecAK4_->setJetA ( (*ak4jets)[ik].jetArea() );
              corr = jecAK4_->getCorrection();
	    } else {uncorrJet = (*ak4jets)[ik].p4();}
            if( (corr*uncorrJet.pt())>20 && (fabs((*ak4jets)[ik].eta()) < 5.0) && looseJetID((*ak4jets)[ik])>0 && nak4<8){
                ak4jet_hf[nak4]=(*ak4jets)[ik].hadronFlavour();
                ak4jet_pf[nak4]=(*ak4jets)[ik].partonFlavour();
                ak4jet_pt[nak4] =  corr*uncorrJet.pt();
                ak4jet_pt_uncorr[nak4] =  uncorrJet.pt();  
                ak4jet_eta[nak4] = (*ak4jets)[ik].eta();
                ak4jet_phi[nak4] = (*ak4jets)[ik].phi();
                ak4jet_e[nak4] =   corr*uncorrJet.energy();
                ak4jet_csv[nak4] = (*ak4jets)[ik].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
                ak4jet_icsv[nak4] = (*ak4jets)[ik].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");   
                ak4jet_IDLoose[nak4] = looseJetID((*ak4jets)[ik]); 
                ak4jet_IDTight[nak4] = tightJetID((*ak4jets)[ik]);
             //   cout<<"Run:"<<RunOnSig_<<endl;
                if(!RunOnSig_){
                        ak4jet_deepcsvudsg[nak4] = (*ak4jets)[ik].bDiscriminator("deepFlavourJetTags:probudsg");
                        ak4jet_deepcsvb[nak4] = (*ak4jets)[ik].bDiscriminator("deepFlavourJetTags:probb");
                        ak4jet_deepcsvc[nak4] = (*ak4jets)[ik].bDiscriminator("deepFlavourJetTags:probc");
                        ak4jet_deepcsvbb[nak4] = (*ak4jets)[ik].bDiscriminator("deepFlavourJetTags:probbb");
                        ak4jet_deepcsvcc[nak4] = (*ak4jets)[ik].bDiscriminator("deepFlavourJetTags:probcc");}
                else{ 
                    ak4jet_deepcsvudsg[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probudsg");
                    //cout<<(*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probudsg")<<"   "<<(*ak4jets)[ik].bDiscriminator("deepFlavourJetTags:probudsg")<<endl;
                    ak4jet_deepcsvb[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probb");
                    ak4jet_deepcsvc[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probc");
                    ak4jet_deepcsvbb[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probbb");
                    ak4jet_deepcsvcc[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probcc");
                    }
                    ak4jet_IDLoose[nak4] = looseJetID((*ak4jets)[ik]);
                    ak4jet_IDTight[nak4] = tightJetID((*ak4jets)[ik]);
                    if(ak4jet_pt[nak4]>tj1 ) {
                        if(tj1>tj2) {tj2=tj1; nj2=nj1;}
                        tj1=ak4jet_pt[nak4]; nj1=nak4;
                    }
                    else if(ak4jet_pt[nak4]>tj2){
                        tj2=ak4jet_pt[nak4]; nj2=nak4;}
                nak4 = nak4 + 1;
            }
         }//3
            if(nj1>-1 && nj2>-1 && ak4jet_pt[nj1]>30. && ak4jet_pt[nj2]>30.) {
                vbfeta=fabs(ak4jet_eta[nj1]-ak4jet_eta[nj2]);
                TLorentzVector vbfj1, vbfj2;
                vbfj1.SetPtEtaPhiE(ak4jet_pt[nj1], ak4jet_eta[nj1], ak4jet_phi[nj1], ak4jet_e[nj1]);
                vbfj2.SetPtEtaPhiE(ak4jet_pt[nj2], ak4jet_eta[nj2], ak4jet_phi[nj2], ak4jet_e[nj2]);
                vbfmjj=(vbfj1+vbfj2).Mag();
            }
            if(vbfeta>4.0 && vbfmjj>400) {vbftag=1;}

                TLorentzVector ghadronicVpuppi, gravitonpuppiJEC,ghadronicVpuppi_2, gravitonpuppiJEC_2, ghadronicVpuppi_3, gravitonpuppiJEC_3;
                ghadronicVpuppi.SetPtEtaPhiM(jetAK8puppi_ptJEC, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sdJEC);
                ghadronicVpuppi_2.SetPtEtaPhiM(jetAK8puppi_ptJEC_2, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sdJEC_2);
                ghadronicVpuppi_3.SetPtEtaPhiM(jetAK8puppi_ptJEC_3, jetAK8puppi_eta_3, jetAK8puppi_phi_3, jetAK8puppi_sdJEC_3);
	        gravitonpuppiJEC =  ghadronicVpuppi+ ghadronicVpuppi_2+ ghadronicVpuppi_3;
                candMasspuppiJEC     = gravitonpuppiJEC.Mag();
                delPhijetmet = deltaPhi(jetAK8puppi_phi, MET_phi);
                delPhijetmet_2 = deltaPhi(jetAK8puppi_phi_2, MET_phi);
                delPhijetmet_3 = deltaPhi(jetAK8puppi_phi_3, MET_phi);
           
     /*      bool IDLoose = looseJetID(hadronicVpuppi);
           bool IDTight = tightJetID(hadronicVpuppi);
           bool IDLoose_2 = looseJetID(hadronicVpuppi_2);
           bool IDTight_2 = tightJetID(hadronicVpuppi_2);
           bool IDLoose_3 = looseJetID(hadronicVpuppi_3);
           bool IDTight_3 = tightJetID(hadronicVpuppi_3);
           bool IDLoose_4 = looseJetID(hadronicVpuppi_4);
           bool IDTight_4 = tightJetID(hadronicVpuppi_4);
*/
           TLorentzVector lvw[3];
           lvw[0] = ghadronicVpuppi;
           lvw[1] = ghadronicVpuppi_2;
           lvw[2] = ghadronicVpuppi_3;
           Double_t Wpt[3];
           Wpt[0]=jetAK8puppi_ptJEC;
           Wpt[1]=jetAK8puppi_ptJEC_2;
           Wpt[2]=jetAK8puppi_ptJEC_3;

           Int_t *indexx=new Int_t[3];
           TMath::Sort(3,Wpt,indexx,1);
           //cout<<Wpt[indexx[0]]<<"   "<<Wpt[indexx[1]]<<"   "<<Wpt[indexx[2]]<<"   "<<endl;
           massww[0] = (lvw[0]+lvw[1]).Mag();
           massww[1] = (lvw[0]+lvw[2]).Mag();
           massww[2] = (lvw[1]+lvw[2]).Mag();
           edm::Handle<edm::View<pat::Jet> > jetsAK8;
           
            iEvent.getByToken(jetsAK8Label_, jetsAK8);
            
            edm::View<pat::Jet>::const_iterator beginAK8 = jetsAK8->begin();
            edm::View<pat::Jet>::const_iterator endAK8 = jetsAK8->end();
            edm::View<pat::Jet>::const_iterator ijetAK8 = beginAK8;
           edm::View<pat::Jet>::const_iterator ijetAK8_j1,ijetAK8_j2,ijetAK8_j3;
           double drak8jetmatch1=10000.,drak8jetmatch2=10000.,drak8jetmatch3=10000.;
         
            for(ijetAK8 = beginAK8; ijetAK8 != endAK8; ++ijetAK8 ) {
                    if(ijetAK8->pt()>0){
                    double tmpdrak8jet1=deltaR(ijetAK8->eta(),ijetAK8->phi(),jetAK8puppi_eta,jetAK8puppi_phi);                  
                    double tmpdrak8jet2=deltaR(ijetAK8->eta(),ijetAK8->phi(),jetAK8puppi_eta_2,jetAK8puppi_phi_2);
                    double tmpdrak8jet3=deltaR(ijetAK8->eta(),ijetAK8->phi(),jetAK8puppi_eta_3,jetAK8puppi_phi_3);
                    if (tmpdrak8jet1<drak8jetmatch1) {drak8jetmatch1=tmpdrak8jet1; ijetAK8_j1=ijetAK8;}
                    if (tmpdrak8jet2<drak8jetmatch2) {drak8jetmatch2=tmpdrak8jet2; ijetAK8_j2=ijetAK8;}
                    if (tmpdrak8jet3<drak8jetmatch3) {drak8jetmatch3=tmpdrak8jet3; ijetAK8_j3=ijetAK8;}
                                       }  
            }
           auto const & sdSubjetsPuppi_1 = ijetAK8_j1->subjets("SoftDropPuppi");
           Int_t nsj1=0;
           for ( auto const & puppiSDSJ_1 : sdSubjetsPuppi_1 ) {
                      if (nsj1==0)
                           ak8sj11.SetPtEtaPhiM(puppiSDSJ_1->correctedP4(0).pt(),puppiSDSJ_1->correctedP4(0).eta(),puppiSDSJ_1->correctedP4(0).phi(),puppiSDSJ_1->correctedP4(0).mass());
                       if (nsj1==1)
                           ak8sj12.SetPtEtaPhiM(puppiSDSJ_1->correctedP4(0).pt(),puppiSDSJ_1->correctedP4(0).eta(),puppiSDSJ_1->correctedP4(0).phi(),puppiSDSJ_1->correctedP4(0).mass());
                       if (nsj1==2)
                           ak8sj13.SetPtEtaPhiM(puppiSDSJ_1->correctedP4(0).pt(),puppiSDSJ_1->correctedP4(0).eta(),puppiSDSJ_1->correctedP4(0).phi(),puppiSDSJ_1->correctedP4(0).mass());
                       if (nsj1==3)
                           ak8sj14.SetPtEtaPhiM(puppiSDSJ_1->correctedP4(0).pt(),puppiSDSJ_1->correctedP4(0).eta(),puppiSDSJ_1->correctedP4(0).phi(),puppiSDSJ_1->correctedP4(0).mass());
                       if (nsj1==4)
                           ak8sj15.SetPtEtaPhiM(puppiSDSJ_1->correctedP4(0).pt(),puppiSDSJ_1->correctedP4(0).eta(),puppiSDSJ_1->correctedP4(0).phi(),puppiSDSJ_1->correctedP4(0).mass());
                   nsj1++;
               }

           auto const & sdSubjetsPuppi_2 = ijetAK8_j2->subjets("SoftDropPuppi");
           Int_t nsj2=0;
           for ( auto const & puppiSDSJ_2 : sdSubjetsPuppi_2 ) {
                      if (nsj2==0)
                           ak8sj21.SetPtEtaPhiM(puppiSDSJ_2->correctedP4(0).pt(),puppiSDSJ_2->correctedP4(0).eta(),puppiSDSJ_2->correctedP4(0).phi(),puppiSDSJ_2->correctedP4(0).mass());
                       if (nsj2==1)
                           ak8sj22.SetPtEtaPhiM(puppiSDSJ_2->correctedP4(0).pt(),puppiSDSJ_2->correctedP4(0).eta(),puppiSDSJ_2->correctedP4(0).phi(),puppiSDSJ_2->correctedP4(0).mass());
                       if (nsj2==2)
                           ak8sj23.SetPtEtaPhiM(puppiSDSJ_2->correctedP4(0).pt(),puppiSDSJ_2->correctedP4(0).eta(),puppiSDSJ_2->correctedP4(0).phi(),puppiSDSJ_2->correctedP4(0).mass());
                       if (nsj2==3)
                           ak8sj24.SetPtEtaPhiM(puppiSDSJ_2->correctedP4(0).pt(),puppiSDSJ_2->correctedP4(0).eta(),puppiSDSJ_2->correctedP4(0).phi(),puppiSDSJ_2->correctedP4(0).mass());
                       if (nsj2==4)
                           ak8sj25.SetPtEtaPhiM(puppiSDSJ_2->correctedP4(0).pt(),puppiSDSJ_2->correctedP4(0).eta(),puppiSDSJ_2->correctedP4(0).phi(),puppiSDSJ_2->correctedP4(0).mass());
                   nsj2++;
               }
 
           auto const & sdSubjetsPuppi_3 = ijetAK8_j3->subjets("SoftDropPuppi");
           Int_t nsj3=0;
           for ( auto const & puppiSDSJ_3 : sdSubjetsPuppi_3 ) {
                      if (nsj3==0)
                           ak8sj31.SetPtEtaPhiM(puppiSDSJ_3->correctedP4(0).pt(),puppiSDSJ_3->correctedP4(0).eta(),puppiSDSJ_3->correctedP4(0).phi(),puppiSDSJ_3->correctedP4(0).mass());
                       if (nsj3==1)
                           ak8sj32.SetPtEtaPhiM(puppiSDSJ_3->correctedP4(0).pt(),puppiSDSJ_3->correctedP4(0).eta(),puppiSDSJ_3->correctedP4(0).phi(),puppiSDSJ_3->correctedP4(0).mass());
                       if (nsj3==2)
                           ak8sj33.SetPtEtaPhiM(puppiSDSJ_3->correctedP4(0).pt(),puppiSDSJ_3->correctedP4(0).eta(),puppiSDSJ_3->correctedP4(0).phi(),puppiSDSJ_3->correctedP4(0).mass());
                       if (nsj3==3)
                           ak8sj34.SetPtEtaPhiM(puppiSDSJ_3->correctedP4(0).pt(),puppiSDSJ_3->correctedP4(0).eta(),puppiSDSJ_3->correctedP4(0).phi(),puppiSDSJ_3->correctedP4(0).mass());
                       if (nsj3==4)
                           ak8sj35.SetPtEtaPhiM(puppiSDSJ_3->correctedP4(0).pt(),puppiSDSJ_3->correctedP4(0).eta(),puppiSDSJ_3->correctedP4(0).phi(),puppiSDSJ_3->correctedP4(0).mass());
                   nsj3++;
               }
         }//2
        }//1
Int_t Nj8=(jetAK8puppi_ptJEC>0)+(jetAK8puppi_ptJEC_2>0)+(jetAK8puppi_ptJEC_3>0)+(jetAK8puppi_ptJEC_4>0)+(jetAK8puppi_ptJEC_5>0)+(jetAK8puppi_ptJEC_6>0)+(jetAK8puppi_ptJEC_7>0)+(jetAK8puppi_ptJEC_8>0);

TLorentzVector AK81,AK82,AK83,AK84;
AK81.SetPtEtaPhiM(0,-99,-99,-99);
AK82.SetPtEtaPhiM(0,-99,-99,-99);
AK83.SetPtEtaPhiM(0,-99,-99,-99);
AK84.SetPtEtaPhiM(0,-99,-99,-99);
Double_t MJJ=-99;
Double_t MJJJ=-99;
Double_t MJJJJ=-99;
if(Nj8==2){
AK81.SetPtEtaPhiM(jetAK8puppi_ptJEC,jetAK8puppi_eta,jetAK8puppi_phi,jetAK8puppi_sdcorr);
AK82.SetPtEtaPhiM(jetAK8puppi_ptJEC_2,jetAK8puppi_eta_2,jetAK8puppi_phi_2,jetAK8puppi_sdcorr_2);
MJJ=(AK81+AK82).M();
          }
if(Nj8==3){
MJJJ=candMasspuppiJEC;
          }
if(Nj8==4){
AK81.SetPtEtaPhiM(jetAK8puppi_ptJEC,jetAK8puppi_eta,jetAK8puppi_phi,jetAK8puppi_sdcorr);
AK82.SetPtEtaPhiM(jetAK8puppi_ptJEC_2,jetAK8puppi_eta_2,jetAK8puppi_phi_2,jetAK8puppi_sdcorr_2);
AK83.SetPtEtaPhiM(jetAK8puppi_ptJEC_3,jetAK8puppi_eta_3,jetAK8puppi_phi_3,jetAK8puppi_sdcorr_3);
AK84.SetPtEtaPhiM(jetAK8puppi_ptJEC_4,jetAK8puppi_eta_4,jetAK8puppi_phi_4,jetAK8puppi_sdcorr_4);
MJJJJ=(AK81+AK82+AK83+AK84).M();
          }

//cout<<passFilter_EEBadSc_<<"  "<<passFilter_badMuon_<<" "<<passFilter_badChargedHadron_<<endl;
//cout<<"IDLoose:"<<IDLoose<<endl;
//cout<<"IDLoose2:"<<IDLoose_2<<endl;
         //cout<<"isHEEP:"<<isHEEP<<endl;
         //cout<<"nVetoEle:"<<nVetoEle<<endl;
         //cout<<"nVetoMu:"<<nVetoMu<<endl;
         //nVetoEle=isHEEP*nVetoEle;//Consider HEEP7 condition
         if(RunOnMC_){
              passFilter_EEBadSc_=true;//Only used for Data
                     }
         if(nVetoEle==0&&nVetoMu==0&&((Nj8==2&&MJJ>-1000&&IDLoose>0&&IDLoose_2>0)||(Nj8==3&&MJJJ>-1000&&IDLoose>0&&IDLoose_2>0&&IDLoose_3>0&&fabs(jetAK8puppi_eta_3)<2.4)||(Nj8==4&&MJJJJ>-1000&&IDLoose>0&&IDLoose_2>0&&IDLoose_3>0&&IDLoose_4>0&&fabs(jetAK8puppi_eta_4)<2.4)) && jetAK8puppi_ptJEC>400 && jetAK8puppi_ptJEC_2>200 && jetAK8puppi_sdcorr>40 && jetAK8puppi_sdcorr_2>40&& (HLT_Mu1>0 || HLT_Mu2>0 || HLT_Mu3>0 || HLT_Mu4>0 || HLT_Mu5>0 || HLT_Mu6>0 || HLT_Mu7>0 || HLT_Mu8>0 || HLT_Mu9>0) && fabs(jetAK8puppi_eta)<2.4&& fabs(jetAK8puppi_eta_2)<2.4 && passFilter_HBHE_>0 && passFilter_GlobalHalo_>0 && passFilter_HBHEIso_>0 && passFilter_ECALDeadCell_>0 && passFilter_GoodVtx_>0 && passFilter_badMuon_>0 && passFilter_EEBadSc_>0){
              outTree_->Fill();
                     }
       outTreew_->Fill();
   }
   else {
       outTreew_->Fill();
        }
}
//-------------------------------------------------------------------------------------------------------------------------------------//


void EDBRTreeMaker::setDummyValues() {
     npT=-1.;
     npIT=-1.;
     nBX=-1;
    ak8sj11.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj21.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj31.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj12.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj22.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj32.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj13.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj23.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj33.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj14.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj24.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj34.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj15.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj25.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj35.SetPtEtaPhiM(0,-99,-99,-99);

    for(int p=0;p<211;p++){
                pweight[p]=-99; 
                          }
    for(Int_t ii=0;ii<8;ii++){
        ak4jet_hf[ii] = -99;
        ak4jet_pf[ii] = -99;
        ak4jet_pt[ii] = -99;
        ak4jet_pt_uncorr[ii] = -99;
        ak4jet_eta[ii] = -99;
        ak4jet_phi[ii] = -99;
        ak4jet_e[ii] = -99;
        ak4jet_dr[ii] = -99;
        ak4jet_csv[ii] = -99;
        ak4jet_icsv[ii] = -99;
        ak4jet_deepcsvudsg[ii] = -99;
        ak4jet_deepcsvb[ii] = -99;
        ak4jet_deepcsvc[ii] = -99;
        ak4jet_deepcsvbb[ii] = -99;
        ak4jet_deepcsvcc[ii] = -99;
        ak4jet_IDLoose[ii] = -99;
        ak4jet_IDTight[ii] = -99;
    }
     vbfeta=-10.;
     vbfmjj=-10.;
     vbftag=0;
     nLooseEle      =-99;
     nLooseMu       =-99;
     nVetoEle      =-99;
     nVetoMu       =-99;
     nVtx           = -99;
     triggerWeight  = -99;
     pileupWeight   = -99;
     lumiWeight     = -99;
 
    for(int i=0;i<5;i++){
        ptgenwl[i]=-99;etagenwl[i]=-99;phigenwl[i]=-99;massgenwl[i]=-99;taggenwl[i]=-99;taggenwmother[i]=-99;
        genw_q1_pt[i]=-99;genw_q1_phi[i]=-99;genw_q1_eta[i]=-99;genw_q1_e[i]=-99;genw_q1_pdg[i]=-99;
        genw_q2_pt[i]=-99;genw_q2_phi[i]=-99;genw_q2_eta[i]=-99;genw_q2_e[i]=-99;genw_q2_pdg[i]=-99;
        ptgenzl[i]=-99;etagenzl[i]=-99;phigenzl[i]=-99;massgenzl[i]=-99;taggenzl[i]=-99;
        ptgenwf[i]=-99;etagenwf[i]=-99;phigenwf[i]=-99;massgenwf[i]=-99;
        ptgenzf[i]=-99;etagenzf[i]=-99;phigenzf[i]=-99;massgenzf[i]=-99;
    }
    for(int i=0;i<15;i++){
        ptgengl[i]=-99;etagengl[i]=-99;phigengl[i]=-99;egengl[i]=-99;
        ptgengf[i]=-99;etagengf[i]=-99;phigengf[i]=-99;egengf[i]=-99;
        mothergengf[i]=-99;mmothergengf[i]=-99;
    }
    for(int i=0;i<5;i++){
        ptgenq1l[i]=-99;etagenq1l[i]=-99;phigenq1l[i]=-99;egenq1l[i]=-99;
        ptgenq1f[i]=-99;etagenq1f[i]=-99;phigenq1f[i]=-99;egenq1f[i]=-99;
        ptgenq2l[i]=-99;etagenq2l[i]=-99;phigenq2l[i]=-99;egenq2l[i]=-99;
        ptgenq2f[i]=-99;etagenq2f[i]=-99;phigenq2f[i]=-99;egenq2f[i]=-99;
        ptgenq3l[i]=-99;etagenq3l[i]=-99;phigenq3l[i]=-99;egenq3l[i]=-99;
        ptgenq3f[i]=-99;etagenq3f[i]=-99;phigenq3f[i]=-99;egenq3f[i]=-99;
        ptgenq4l[i]=-99;etagenq4l[i]=-99;phigenq4l[i]=-99;egenq4l[i]=-99;
        ptgenq4f[i]=-99;etagenq4f[i]=-99;phigenq4f[i]=-99;egenq4f[i]=-99;
        ptgenq5l[i]=-99;etagenq5l[i]=-99;phigenq5l[i]=-99;egenq5l[i]=-99;
        ptgenq5f[i]=-99;etagenq5f[i]=-99;phigenq5f[i]=-99;egenq5f[i]=-99;
        mmothergenq1f[i]=-99;mmothergenq2f[i]=-99;mmothergenq3f[i]=-99;mmothergenq4f[i]=-99;mmothergenq5f[i]=-99;

    }
    gent_b_pt=-99;gent_b_phi=-99;gent_b_eta=-99;gent_b_mass=-99;
    genantit_b_pt=-99;genantit_b_phi=-99;genantit_b_eta=-99;genantit_b_mass=-99;
    gent_w_pt=-99;gent_w_phi=-99;gent_w_eta=-99;gent_w_mass=-99;
    genantit_w_pt=-99;genantit_w_phi=-99;genantit_w_eta=-99;genantit_w_mass=-99;
    gent_w_q1_pt=-99;gent_w_q1_phi=-99;gent_w_q1_eta=-99;gent_w_q1_e=-99;gent_w_q1_pdg=-99;
    genantit_w_q1_pt=-99;genantit_w_q1_phi=-99;genantit_w_q1_eta=-99;genantit_w_q1_e=-99;genantit_w_q1_pdg=-99;
    gent_w_q2_pt=-99;gent_w_q2_phi=-99;gent_w_q2_eta=-99;gent_w_q2_e=-99;gent_w_q2_pdg=-99;
    genantit_w_q2_pt=-99;genantit_w_q2_phi=-99;genantit_w_q2_eta=-99;genantit_w_q2_e=-99;genantit_w_q2_pdg=-99;
    gent_w_tag=-99;genantit_w_tag=-99;

    // DeepAK8
    jetAK8puppi_dnnTop        = -99;
    jetAK8puppi_dnnW          = -99;
    jetAK8puppi_dnnH4q        = -99;
    jetAK8puppi_dnnTop_2        = -99;
    jetAK8puppi_dnnW_2          = -99;
    jetAK8puppi_dnnH4q_2        = -99;
    jetAK8puppi_dnnTop_3        = -99;
    jetAK8puppi_dnnW_3          = -99;
    jetAK8puppi_dnnH4q_3        = -99;
    jetAK8puppi_dnnTop_4        = -99;
    jetAK8puppi_dnnW_4          = -99;
    jetAK8puppi_dnnH4q_4        = -99;
    jetAK8puppi_dnnZ        = -99;
    jetAK8puppi_dnnZbb        = -99;
    jetAK8puppi_dnnHbb        = -99;
    jetAK8puppi_dnnZ_2        = -99;
    jetAK8puppi_dnnZbb_2        = -99;
    jetAK8puppi_dnnHbb_2        = -99;
    jetAK8puppi_dnnZ_3        = -99;
    jetAK8puppi_dnnZbb_3        = -99;
    jetAK8puppi_dnnHbb_3        = -99;
    jetAK8puppi_dnnZ_4        = -99;
    jetAK8puppi_dnnZbb_4        = -99;
    jetAK8puppi_dnnHbb_4        = -99;
    jetAK8puppi_dnnqcd        = -99;
    jetAK8puppi_dnntop        = -99;
    jetAK8puppi_dnnw          = -99;
    jetAK8puppi_dnnz          = -99;
    jetAK8puppi_dnnzbb        = -99;
    jetAK8puppi_dnnhbb        = -99;
    jetAK8puppi_dnnh4q        = -99;
    jetAK8puppi_dnnqcd_2        = -99;
    jetAK8puppi_dnntop_2        = -99;
    jetAK8puppi_dnnw_2          = -99;
    jetAK8puppi_dnnz_2          = -99;
    jetAK8puppi_dnnzbb_2        = -99;
    jetAK8puppi_dnnhbb_2        = -99;
    jetAK8puppi_dnnh4q_2        = -99;
    jetAK8puppi_dnnqcd_3        = -99;
    jetAK8puppi_dnntop_3        = -99;
    jetAK8puppi_dnnw_3          = -99;
    jetAK8puppi_dnnz_3          = -99;
    jetAK8puppi_dnnzbb_3        = -99;
    jetAK8puppi_dnnhbb_3        = -99;
    jetAK8puppi_dnnh4q_3        = -99;
    jetAK8puppi_dnnqcd_4        = -99;
    jetAK8puppi_dnntop_4        = -99;
    jetAK8puppi_dnnw_4          = -99;
    jetAK8puppi_dnnz_4          = -99;
    jetAK8puppi_dnnzbb_4        = -99;
    jetAK8puppi_dnnhbb_4        = -99;
    jetAK8puppi_dnnh4q_4        = -99;
    // Decorrelated DeepAK8
    jetAK8puppi_dnnDecorrTop        = -99;
    jetAK8puppi_dnnDecorrW          = -99;
    jetAK8puppi_dnnDecorrH4q        = -99;
    jetAK8puppi_dnnDecorrTop_2        = -99;
    jetAK8puppi_dnnDecorrW_2          = -99;
    jetAK8puppi_dnnDecorrH4q_2        = -99;
    jetAK8puppi_dnnDecorrTop_3        = -99;
    jetAK8puppi_dnnDecorrW_3          = -99;
    jetAK8puppi_dnnDecorrH4q_3        = -99;
    jetAK8puppi_dnnDecorrTop_4        = -99;
    jetAK8puppi_dnnDecorrW_4          = -99;
    jetAK8puppi_dnnDecorrH4q_4        = -99;
    jetAK8puppi_dnnDecorrZ        = -99;
    jetAK8puppi_dnnDecorrZbb        = -99;
    jetAK8puppi_dnnDecorrHbb        = -99;
    jetAK8puppi_dnnDecorrZ_2        = -99;
    jetAK8puppi_dnnDecorrZbb_2        = -99;
    jetAK8puppi_dnnDecorrHbb_2        = -99;
    jetAK8puppi_dnnDecorrZ_3        = -99;
    jetAK8puppi_dnnDecorrZbb_3        = -99;
    jetAK8puppi_dnnDecorrHbb_3        = -99;
    jetAK8puppi_dnnDecorrZ_4        = -99;
    jetAK8puppi_dnnDecorrZbb_4        = -99;
    jetAK8puppi_dnnDecorrHbb_4        = -99;
    jetAK8puppi_dnnDecorrbb        = -99;
    jetAK8puppi_dnnDecorrcc        = -99;
    jetAK8puppi_dnnDecorrbbnog        = -99;
    jetAK8puppi_dnnDecorrccnog        = -99;
    jetAK8puppi_dnnDecorrbb_2        = -99;
    jetAK8puppi_dnnDecorrcc_2        = -99;
    jetAK8puppi_dnnDecorrbbnog_2        = -99;
    jetAK8puppi_dnnDecorrccnog_2        = -99;
    jetAK8puppi_dnnDecorrbb_3        = -99;
    jetAK8puppi_dnnDecorrcc_3        = -99;
    jetAK8puppi_dnnDecorrbbnog_3        = -99;
    jetAK8puppi_dnnDecorrccnog_3        = -99;
    jetAK8puppi_dnnDecorrbb_4        = -99;
    jetAK8puppi_dnnDecorrcc_4        = -99;
    jetAK8puppi_dnnDecorrbbnog_4        = -99;
    jetAK8puppi_dnnDecorrccnog_4        = -99;
    jetAK8puppi_dnnDecorrqcd        = -99;
    jetAK8puppi_dnnDecorrtop        = -99;
    jetAK8puppi_dnnDecorrw          = -99;
    jetAK8puppi_dnnDecorrz          = -99;
    jetAK8puppi_dnnDecorrzbb        = -99;
    jetAK8puppi_dnnDecorrhbb        = -99;
    jetAK8puppi_dnnDecorrh4q        = -99;
    jetAK8puppi_dnnDecorrqcd_2        = -99;
    jetAK8puppi_dnnDecorrtop_2        = -99;
    jetAK8puppi_dnnDecorrw_2          = -99;
    jetAK8puppi_dnnDecorrz_2          = -99;
    jetAK8puppi_dnnDecorrzbb_2        = -99;
    jetAK8puppi_dnnDecorrhbb_2        = -99;
    jetAK8puppi_dnnDecorrh4q_2        = -99;
    jetAK8puppi_dnnDecorrqcd_3        = -99;
    jetAK8puppi_dnnDecorrtop_3        = -99;
    jetAK8puppi_dnnDecorrw_3          = -99;
    jetAK8puppi_dnnDecorrz_3          = -99;
    jetAK8puppi_dnnDecorrzbb_3        = -99;
    jetAK8puppi_dnnDecorrhbb_3        = -99;
    jetAK8puppi_dnnDecorrh4q_3        = -99;
    jetAK8puppi_dnnDecorrqcd_4        = -99;
    jetAK8puppi_dnnDecorrtop_4        = -99;
    jetAK8puppi_dnnDecorrw_4          = -99;
    jetAK8puppi_dnnDecorrz_4          = -99;
    jetAK8puppi_dnnDecorrzbb_4        = -99;
    jetAK8puppi_dnnDecorrhbb_4        = -99;
    jetAK8puppi_dnnDecorrh4q_4        = -99;

     jetAK8puppi_ptJEC         = -99;
     jetAK8puppi_eta         = -99;
     jetAK8puppi_phi         = -99;
     jetAK8puppi_tau1         = -99;
     jetAK8puppi_tau2         = -99;
     jetAK8puppi_tau3         = -99;
     jetAK8puppi_tau21         = -99;
     jetAK8puppi_tau4         = -99;
     jetAK8puppi_tau42         = -99;
     jetAK8puppi_sd         = -99;
     jetAK8puppi_sdJEC         = -99; 
     jetAK8puppi_sdcorr         = -99;

     jetAK8puppi_ptJEC_2         = -99;
     jetAK8puppi_eta_2         = -99;
     jetAK8puppi_phi_2         = -99;
     jetAK8puppi_tau1_2         = -99;
     jetAK8puppi_tau2_2         = -99;
     jetAK8puppi_tau3_2         = -99;
     jetAK8puppi_tau21_2         = -99;
     jetAK8puppi_tau4_2         = -99;
     jetAK8puppi_tau42_2         = -99;
     jetAK8puppi_sd_2         = -99;
     jetAK8puppi_sdJEC_2         = -99;
     jetAK8puppi_sdcorr_2         = -99;
 
     jetAK8puppi_ptJEC_3         = -99;
     jetAK8puppi_eta_3         = -99;
     jetAK8puppi_phi_3         = -99;
     jetAK8puppi_tau1_3         = -99;
     jetAK8puppi_tau2_3         = -99;
     jetAK8puppi_tau3_3         = -99;
     jetAK8puppi_tau21_3         = -99;
     jetAK8puppi_tau4_3         = -99;
     jetAK8puppi_tau42_3         = -99;
     jetAK8puppi_sd_3         = -99;
     jetAK8puppi_sdJEC_3         = -99;
     jetAK8puppi_sdcorr_3         = -99;

     jetAK8puppi_ptJEC_4         = -99;
     jetAK8puppi_eta_4         = -99;
     jetAK8puppi_phi_4         = -99;
     jetAK8puppi_tau1_4         = -99;
     jetAK8puppi_tau2_4         = -99;
     jetAK8puppi_tau3_4         = -99;
     jetAK8puppi_tau21_4         = -99;
     jetAK8puppi_tau4_4         = -99;
     jetAK8puppi_tau42_4         = -99;
     jetAK8puppi_sd_4         = -99;
     jetAK8puppi_sdJEC_4         = -99;
     jetAK8puppi_sdcorr_4         = -99;

     jetAK8puppi_ptJEC_5         = -99;
     jetAK8puppi_eta_5         = -99;
     jetAK8puppi_phi_5         = -99;
     jetAK8puppi_tau1_5         = -99;
     jetAK8puppi_tau2_5         = -99;
     jetAK8puppi_tau3_5         = -99;
     jetAK8puppi_tau21_5         = -99;
     jetAK8puppi_tau4_5         = -99;
     jetAK8puppi_tau42_5         = -99;
     jetAK8puppi_sd_5         = -99;
     jetAK8puppi_sdJEC_5         = -99;
     jetAK8puppi_sdcorr_5         = -99;

     jetAK8puppi_ptJEC_6         = -99;
     jetAK8puppi_eta_6         = -99;
     jetAK8puppi_phi_6         = -99;
     jetAK8puppi_tau1_6         = -99;
     jetAK8puppi_tau2_6         = -99;
     jetAK8puppi_tau3_6         = -99;
     jetAK8puppi_tau21_6         = -99;
     jetAK8puppi_tau4_6         = -99;
     jetAK8puppi_tau42_6         = -99;
     jetAK8puppi_sd_6         = -99;
     jetAK8puppi_sdJEC_6         = -99;
     jetAK8puppi_sdcorr_6         = -99;

     jetAK8puppi_ptJEC_7         = -99;
     jetAK8puppi_eta_7         = -99;
     jetAK8puppi_phi_7         = -99;
     jetAK8puppi_tau1_7         = -99;
     jetAK8puppi_tau2_7         = -99;
     jetAK8puppi_tau3_7         = -99;
     jetAK8puppi_tau21_7         = -99;
     jetAK8puppi_tau4_7         = -99;
     jetAK8puppi_tau42_7         = -99;
     jetAK8puppi_sd_7         = -99;
     jetAK8puppi_sdJEC_7         = -99;
     jetAK8puppi_sdcorr_7         = -99;

     jetAK8puppi_ptJEC_8         = -99;
     jetAK8puppi_eta_8         = -99;
     jetAK8puppi_phi_8         = -99;
     jetAK8puppi_tau1_8         = -99;
     jetAK8puppi_tau2_8         = -99;
     jetAK8puppi_tau3_8         = -99;
     jetAK8puppi_tau21_8         = -99;
     jetAK8puppi_tau4_8         = -99;
     jetAK8puppi_tau42_8         = -99;
     jetAK8puppi_sd_8         = -99;
     jetAK8puppi_sdJEC_8         = -99;
     jetAK8puppi_sdcorr_8         = -99;

     met            = -99;
     metPhi         = -99;
     delPhijetmet =  -99;
     delPhijetmet_2 =  -99;
     delPhijetmet_3 =  -99;

     gen_gra_m      = -99;
     gen_gra_pt     = -99;
     gen_gra_eta     = -99;
     gen_gra_phi     = -99;
     gen_rad_m      = -99;
     gen_rad_pt     = -99;
     gen_rad_eta     = -99;
     gen_rad_phi     = -99;
     gen_ele_pt     = -99;
     gen_ele_eta    = -99;
     gen_ele_phi    = -99;
     gen_ele_e      = -99;
     gen_mu_pt     = -99;
     gen_mu_eta    = -99;
     gen_mu_phi    = -99;
     gen_mu_e      = -99;
     gen_ele_pt_2     = -99;
     gen_ele_eta_2    = -99;
     gen_ele_phi_2    = -99;
     gen_ele_e_2      = -99;
     gen_mu_pt_2     = -99;
     gen_mu_eta_2    = -99;
     gen_mu_phi_2    = -99;
     gen_mu_e_2      = -99;
     gen_ele_pt_3     = -99;
     gen_ele_eta_3    = -99;
     gen_ele_phi_3    = -99;
     gen_ele_e_3      = -99;
     gen_mu_pt_3     = -99;
     gen_mu_eta_3    = -99;
     gen_mu_phi_3    = -99;
     gen_mu_e_3      = -99;

     gen_tau_pt     = -99;
     gen_tau_eta    = -99;
     gen_tau_phi    = -99;
     gen_tau_e      = -99;
     gen_tau_pt_2     = -99;
     gen_tau_eta_2    = -99;
     gen_tau_phi_2    = -99;
     gen_tau_e_2      = -99;
     gen_tau_pt_3     = -99;
     gen_tau_eta_3    = -99;
     gen_tau_phi_3    = -99;
     gen_tau_e_3      = -99;

     gentop_pt  = -99;
     gentop_eta  = -99;
     gentop_phi  = -99;
     gentop_mass  = -99;
     genantitop_pt  = -99;
     genantitop_eta  = -99;
     genantitop_phi  = -99;
     genantitop_mass  = -99;
     ptGenVlep      = -99;
     etaGenVlep      = -99;
     phiGenVlep      = -99;
     massGenVlep      = -99;
     ptGenVlep_2      = -99;
     etaGenVlep_2      = -99;
     phiGenVlep_2      = -99;
     massGenVlep_2      = -99;
     ptGenVlep_3      = -99;
     etaGenVlep_3      = -99;
     phiGenVlep_3      = -99;
     massGenVlep_3      = -99;     
     ptGenV_2      = -99;
     etaGenV_2      = -99;
     phiGenV_2      = -99;
     massGenV_2      = -99;
     ptGenV_3      = -99;
     etaGenV_3      = -99;
     phiGenV_3      = -99;
     massGenV_3      = -99;
     ptGenVhad      = -99;
     etaGenVhad      = -99;
     phiGenVhad      = -99;
     massGenVhad      = -99;
     ptq11 =-99;
     etaq11 =-99;
     phiq11 =-99;
     massq11 =-99;
     ptq12 =-99;
     etaq12 =-99;
     phiq12 =-99;
     massq12 =-99; 
     ptq21 =-99;
     etaq21 =-99;
     phiq21 =-99;
     massq21 =-99;
     ptq22 =-99;
     etaq22 =-99;
     phiq22 =-99;
     massq22 =-99;
     ptq31 =-99;
     etaq31 =-99;
     phiq31 =-99;
     massq31 =-99;
     ptq32 =-99;
     etaq32 =-99;
     phiq32 =-99;
     massq32 =-99;
     status_1       =  -1;
     status_2       =  -1;
     status_3       =  -1;

     for(Int_t ii=0;ii<8;ii++){
        ak4jet_hf[ii] = -99;
        ak4jet_pf[ii] = -99;
        ak4jet_pt[ii] = -99;
        ak4jet_pt_uncorr[ii] = -99;
        ak4jet_eta[ii] = -99;
        ak4jet_phi[ii] = -99;
        ak4jet_e[ii] = -99;
        ak4jet_dr[ii] = -99;
        ak4jet_csv[ii] = -99;
        ak4jet_icsv[ii] = -99;
        ak4jet_IDLoose[ii] = -99;
        ak4jet_IDTight[ii] = -99;
    }
 

     jetAK8_mass = -99;
     jetAK8_pt = -99;
     jetAK8_jec = -99;
     jetAK8_mass1[0] = -99;
     jetAK8_mass1[1] = -99;
     jetAK8_mass1[2] = -99;
     jetAK8_SF_mass1[0] = -99;
     jetAK8_SF_mass1[1] = -99;
     jetAK8_SF_mass1[2] = -99;
     jetAK8_SF_mass2[0] = -99;
     jetAK8_SF_mass2[1] = -99;
     jetAK8_SF_mass2[2] = -99;
     jetAK8_pt1[0] = -99;
     jetAK8_pt1[1] = -99;
     jetAK8_pt1[2] = -99;
     jetAK8_jec1[0] = -99;
     jetAK8_jec1[1] = -99;
     jetAK8_jec1[2] = -99;
     corr_AK81[0] = -99;
     corr_AK81[1] = -99;
     corr_AK81[2] = -99;
     jetAK8_eta = -99;
     jetAK8_eta1[0] = -99;
     jetAK8_eta1[1] = -99;
     jetAK8_eta1[2] = -99;
     jetAK8_phi = -99;

     METraw_et = -99;
     METraw_phi = -99;
     METraw_sumEt = -99;
     MET_et = -99;
     MET_phi = -99;
     MET_sumEt = -99;
     MET_corrPx = -99;
     MET_corrPy = -99;

     candMasspuppiJEC     =  -99;
     massww[0] = -99;
     massww[1] = -99;
     massww[2] = -99;

     HLT_Mu1=-99;
     HLT_Mu2=-99;
     HLT_Mu3=-99;
     HLT_Mu4=-99;
     HLT_Mu5=-99;
     HLT_Mu6=-99;
     HLT_Mu7=-99;
     HLT_Mu8=-99;
     HLT_Mu9=-99;
     HLT_Mu10=-99;
     HLT_Mu11=-99;
     HLT_Mu12=-99;
     HLT_Mu13=-99;
     HLT_Mu14=-99;
     HLT_Mu15=-99;
     HLT_Mu16=-99;

     theWeight = -99;
     //nump = 0;
     //numm = 0;
     passFilter_HBHE_                  = false;
     passFilter_HBHEIso_               = false;
     passFilter_GlobalHalo_            = false;
     passFilter_ECALDeadCell_          = false;
     passFilter_GoodVtx_               = false;
     passFilter_EEBadSc_               = false;
     passFilter_badMuon_               = false;
     passFilter_badChargedHadron_      = false;

}

// ------------ method called once each job just before starting event loop  ------------
void 
EDBRTreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void EDBRTreeMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  muPaths1.clear();
  muPaths2.clear();
  muPaths3.clear();
  muPaths4.clear();
  muPaths5.clear();
  muPaths6.clear();
  muPaths7.clear();
  muPaths8.clear();
  muPaths9.clear();
  el1.clear();
  el2.clear();
  el3.clear();
  mu1.clear();
  mu2.clear();
  mu3.clear();
  mu4.clear();

  std::cout<<"-----begin-----"<<std::endl;
   bool changed;
   if ( !hltConfig.init(iRun, iSetup, "HLT", changed) ) {
        edm::LogError("HltAnalysis") << "Initialization of HLTConfigProvider failed!!";
       return;
      }
   for (size_t i = 0; i < muPaths1_.size(); i++) {
         std::vector<std::string> foundPaths1 = hltConfig.matched( hltConfig.triggerNames(), muPaths1_[i] );
         while ( !foundPaths1.empty() ){
               muPaths1.push_back( foundPaths1.back() );
               foundPaths1.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-1 Information **************\n";
   for (size_t i=0; i < muPaths1.size(); i++) std::cout << "\n HLT paths-1:   " << i<<"  "<<muPaths1[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths2_.size(); i++) {
         std::vector<std::string> foundPaths2 = hltConfig.matched( hltConfig.triggerNames(), muPaths2_[i] );
         while ( !foundPaths2.empty() ){
               muPaths2.push_back( foundPaths2.back() );
               foundPaths2.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-2 Information **************\n";
   for (size_t i=0; i < muPaths2.size(); i++) std::cout << "\n Muon paths-2:   " << i<<"  "<<muPaths2[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths3_.size(); i++) {
         std::vector<std::string> foundPaths3 = hltConfig.matched( hltConfig.triggerNames(), muPaths3_[i] );
         while ( !foundPaths3.empty() ){
               muPaths3.push_back( foundPaths3.back() );
               foundPaths3.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-3 Information **************\n";
   for (size_t i=0; i < muPaths3.size(); i++) std::cout << "\n Muon paths-3:   " << i<<"  "<<muPaths3[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths4_.size(); i++) {
         std::vector<std::string> foundPaths4 = hltConfig.matched( hltConfig.triggerNames(), muPaths4_[i] );
         while ( !foundPaths4.empty() ){
               muPaths4.push_back( foundPaths4.back() );
               foundPaths4.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-4 Information **************\n";
   for (size_t i=0; i < muPaths4.size(); i++) std::cout << "\n Muon paths-4:   " << i<<"  "<<muPaths4[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths5_.size(); i++) {
         std::vector<std::string> foundPaths5 = hltConfig.matched( hltConfig.triggerNames(), muPaths5_[i] );
         while ( !foundPaths5.empty() ){
               muPaths5.push_back( foundPaths5.back() );
               foundPaths5.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-5 Information **************\n";
   for (size_t i=0; i < muPaths5.size(); i++) std::cout << "\n Muon paths-5:   " << i<<"  "<<muPaths5[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths6_.size(); i++) {
         std::vector<std::string> foundPaths6 = hltConfig.matched( hltConfig.triggerNames(), muPaths6_[i] );
         while ( !foundPaths6.empty() ){
               muPaths6.push_back( foundPaths6.back() );
               foundPaths6.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-6 Information **************\n";
   for (size_t i=0; i < muPaths6.size(); i++) std::cout << "\n Muon paths-6:   " << i<<"  "<<muPaths6[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths7_.size(); i++) {
         std::vector<std::string> foundPaths7 = hltConfig.matched( hltConfig.triggerNames(), muPaths7_[i] );
         while ( !foundPaths7.empty() ){
               muPaths7.push_back( foundPaths7.back() );
               foundPaths7.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-7 Information **************\n";
   for (size_t i=0; i < muPaths7.size(); i++) std::cout << "\n Muon paths-7:   " << i<<"  "<<muPaths7[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths8_.size(); i++) {
         std::vector<std::string> foundPaths8 = hltConfig.matched( hltConfig.triggerNames(), muPaths8_[i] );
         while ( !foundPaths8.empty() ){
               muPaths8.push_back( foundPaths8.back() );
               foundPaths8.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-8 Information **************\n";
   for (size_t i=0; i < muPaths8.size(); i++) std::cout << "\n Muon paths-8:   " << i<<"  "<<muPaths8[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < muPaths9_.size(); i++) {
         std::vector<std::string> foundPaths9 = hltConfig.matched( hltConfig.triggerNames(), muPaths9_[i] );
         while ( !foundPaths9.empty() ){
               muPaths9.push_back( foundPaths9.back() );
               foundPaths9.pop_back();
                                      }
                                                }

   std::cout<<"\n************** HLT-9 Information **************\n";
   for (size_t i=0; i < muPaths9.size(); i++) std::cout << "\n Muon paths-9:   " << i<<"  "<<muPaths9[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < el1_.size(); i++) {
         std::vector<std::string> foundPaths10 = hltConfig.matched( hltConfig.triggerNames(), el1_[i] );
         while ( !foundPaths10.empty() ){
               el1.push_back( foundPaths10.back() );
               foundPaths10.pop_back();
                                      }         
                                                }
   std::cout<<"\n************** HLT-10 Information **************\n";
   for (size_t i=0; i < el1.size(); i++) std::cout << "\n HLT paths-10:   " << i<<"  "<<el1[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < el2_.size(); i++) {
         std::vector<std::string> foundPaths11 = hltConfig.matched( hltConfig.triggerNames(), el2_[i] );
         while ( !foundPaths11.empty() ){
               el2.push_back( foundPaths11.back() );
               foundPaths11.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-11 Information **************\n";
   for (size_t i=0; i < el2.size(); i++) std::cout << "\n HLT paths-11:   " << i<<"  "<<el2[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < el3_.size(); i++) {
         std::vector<std::string> foundPaths12 = hltConfig.matched( hltConfig.triggerNames(), el3_[i] );
         while ( !foundPaths12.empty() ){
               el3.push_back( foundPaths12.back() );
               foundPaths12.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-12 Information **************\n";
   for (size_t i=0; i < el3.size(); i++) std::cout << "\n HLT paths-12:   " << i<<"  "<<el3[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < mu1_.size(); i++) {
         std::vector<std::string> foundPaths13 = hltConfig.matched( hltConfig.triggerNames(), mu1_[i] );
         while ( !foundPaths13.empty() ){
               mu1.push_back( foundPaths13.back() );
               foundPaths13.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-13 Information **************\n";
   for (size_t i=0; i < mu1.size(); i++) std::cout << "\n HLT paths-13:   " << i<<"  "<<mu1[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";


   for (size_t i = 0; i < mu2_.size(); i++) {
         std::vector<std::string> foundPaths14 = hltConfig.matched( hltConfig.triggerNames(), mu2_[i] );
         while ( !foundPaths14.empty() ){
               mu2.push_back( foundPaths14.back() );
               foundPaths14.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-14 Information **************\n";
   for (size_t i=0; i < mu2.size(); i++) std::cout << "\n HLT paths-14:   " << i<<"  "<<mu2[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < mu3_.size(); i++) {
         std::vector<std::string> foundPaths15 = hltConfig.matched( hltConfig.triggerNames(), mu3_[i] );
         while ( !foundPaths15.empty() ){
               mu3.push_back( foundPaths15.back() );
               foundPaths15.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-15 Information **************\n";
   for (size_t i=0; i < mu3.size(); i++) std::cout << "\n HLT paths-15:   " << i<<"  "<<mu3[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";

   for (size_t i = 0; i < mu4_.size(); i++) {
         std::vector<std::string> foundPaths16 = hltConfig.matched( hltConfig.triggerNames(), mu4_[i] );
         while ( !foundPaths16.empty() ){
               mu4.push_back( foundPaths16.back() );
               foundPaths16.pop_back();
                                      }
                                                }
   std::cout<<"\n************** HLT-16 Information **************\n";
   for (size_t i=0; i < mu4.size(); i++) std::cout << "\n HLT paths-16:   " << i<<"  "<<mu4[i].c_str() <<"\t"<< std::endl;
   std::cout<<"\n*********************************************\n\n";



}

void EDBRTreeMaker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
std::cout << "EDBRTreeMaker endJob()... endRun" << std::endl;
}


void
EDBRTreeMaker::endJob() {
  std::cout << "EDBRTreeMaker endJob()..." << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(EDBRTreeMaker);
