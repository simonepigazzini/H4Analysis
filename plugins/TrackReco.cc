#include "TrackReco.h"
#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

double TrackReco::Track::chi2(const double* par)
{
  if (par) //set the actual fit parameters
    {
      // std::cout << "X " << par[0] << " Y " << par[1] << " alpha " << par[2] << " beta " << par[3] << endl;
      position_(0)=par[0];
      position_(1)=par[1];
      if (fitAngle_)
	{
	  angle_(0)=par[2];
	  angle_(1)=par[3];
	}
    }

   double chi2=0;
   int i=0;
   for (auto& hit: TrackReco::Track::hits_)
     {
       //calculate the residual
       Measurement_t pos;
       MeasurementErrorMatrix_t posErrInverse;
       hit.globalPosition(pos);
       hit.globalPositionErrorInverse(posErrInverse);
       Measurement_t delta= pos-this->statusAt((this->hodo_.layers_[hit.layer_].position_(2))); //calculation of the residual
       chi2+=ROOT::Math::Dot(delta,posErrInverse*delta); // delta^T * (V^-1) * delta
       //       cout << i << " " << hit.layer_ << ":" << pos << ":" << delta << ":" << chi2 << ":" << this->hodo_.layers_[hit.layer_].position_(2) << endl;
       ++i;
     }
   return chi2;
}

bool TrackReco::Track::fitTrack()
{
  //---setup minimization
  int nFreeParamters = fitAngle_? 4:2;
  ROOT::Math::Functor chi2(this, &TrackReco::Track::chi2, nFreeParamters);
  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetMaxIterations(1000000);
  minimizer->SetTolerance(1e-6);
  minimizer->SetPrintLevel(0);
  minimizer->SetFunction(chi2);
  minimizer->SetLimitedVariable(0, "X", 0., 1E-6, -30, 30);
  minimizer->SetLimitedVariable(1, "Y", 0., 1E-6, -30, 30);
  if (fitAngle_)
    {
      minimizer->SetLimitedVariable(2, "alpha",0, 1E-6,  -0.1, 0.1);
      minimizer->SetLimitedVariable(3, "beta", 0, 1E-6,  -0.1, 0.1);
    }

  //---fit
  minimizer->Minimize();
  position_(0) = minimizer->X()[0];
  position_(1) = minimizer->X()[1];
  angle_(0) = fitAngle_ ? minimizer->X()[2] : 0; 
  angle_(1) = fitAngle_ ? minimizer->X()[3] : 0;
  
  delete minimizer;        

  return true;
}

//----------Begin-------------------------------------------------------------------------
bool TrackReco::Begin(CfgManager& opts, uint64* index)
{
    //---create a position tree
//    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
//        opts.GetOpt<bool>(instanceName_+".storeTree") : true;


//  Now build the hodoscope
    GlobalCoord_t layer0Pos(0,0,0);
    GlobalCoord_t layer1Pos(0.,0.,0);
    GlobalCoord_t layer2Pos(0.,0.,500);
    GlobalCoord_t layer3Pos(0,0,1000);
    GlobalCoord_t layer4Pos(0.,0.,1000);

    TrackLayer layer_0(layer0Pos);
    TrackLayer layer_1(layer1Pos);
    TrackLayer layer_2(layer2Pos);
    TrackLayer layer_3(layer3Pos);
    TrackLayer layer_4(layer4Pos);

    hodo_.addLayer(layer_0);
    hodo_.addLayer(layer_1);
    hodo_.addLayer(layer_2);
    hodo_.addLayer(layer_3);
    hodo_.addLayer(layer_4);

    std::cout << "BUILT HODO" << std::endl;
    int i=0;
    for (auto& layer : hodo_.layers_)
      {
	cout << "LAYER "<< i << ":" << layer.position_ << std::endl;
	++i;
      }
    return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool TrackReco::ProcessEvent(H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{

    tracks_.clear();

    //Get track measurements
    TrackMeasurement aHit_0(0.5,0.,hodo_,0);
    aHit_0.setVarianceX(0.04);
    aHit_0.setVarianceY(9999.);
    aHit_0.calculateInverseVariance();
   TrackMeasurement aHit_1(0.,0.5,hodo_,1);
    aHit_1.setVarianceX(9999.);
    aHit_1.setVarianceY(0.04);
    aHit_1.calculateInverseVariance();
    TrackMeasurement aHit_2(0.3,0.3,hodo_,2);
    aHit_2.setVarianceX(0.0004);
    aHit_2.setVarianceY(0.0004);
    aHit_2.calculateInverseVariance();
    TrackMeasurement aHit_3(0.5,0.,hodo_,3);
    aHit_3.setVarianceX(0.04);
    aHit_3.setVarianceY(9999.);
    aHit_3.calculateInverseVariance();
   TrackMeasurement aHit_4(0.,0.5,hodo_,4);
    aHit_4.setVarianceX(9999.);
    aHit_4.setVarianceY(0.04);
    aHit_4.calculateInverseVariance();

    Track aTrack(hodo_);
    aTrack.addMeasurement(aHit_0);
    aTrack.addMeasurement(aHit_1);
    aTrack.addMeasurement(aHit_2);
    aTrack.addMeasurement(aHit_3);
    aTrack.addMeasurement(aHit_4);

    //fit track
    aTrack.fitAngle_ = true;
    aTrack.fitTrack();

    std::cout << aTrack.position_(0) << "," << aTrack.position_(1) ;
    if (aTrack.fitAngle_)
      std::cout << ":" << aTrack.angle_(0) << "," << aTrack.angle_(1) ;
    std::cout << std::endl;

    std::cout << "CHI2 " << aTrack.chi2() << std::endl;
    //---fill output tree
    return true;
}
