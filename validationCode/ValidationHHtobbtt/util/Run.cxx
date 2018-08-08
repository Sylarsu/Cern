#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include "SampleHandler/DiskListEOS.h"
#include "EventLoopGrid/PrunDriver.h"
#include <TSystem.h>

#include "ValidationHHtobbtt/MyEventLoop.h"

int main( int argc, char* argv[] ) {

   // Take the submit directory from the input if provided:
   int runGrid = 0;
   std::string submitDir = "submitDir";
   std::string jobId = "test.01";
   std::string user(gSystem->GetFromPipe("whoami"));
   double nFiles = 0;

   if( argc > 1 ) submitDir = argv[ 1 ];
   if( argc > 2 ) runGrid = std::stoi( argv[ 2 ]);
   if( argc > 3 ) jobId = argv[ 3 ]; 
   if( argc > 4 ) nFiles = (double)std::stoi( argv[ 4 ] );

   // Set up the job for xAOD access:
   xAOD::Init().ignore();

   // Construct the samples to run on:
   SH::SampleHandler sh;

   if( runGrid ) {

      std::string jobname=""+user+"."+jobId+"";

      //SH::scanDQ2 (sh, "mc15_13TeV:mc15_13TeV.301054.PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYtautau_3000M3500.recon.AOD.e3649_s2576_s2132_r6545/");
      //SH::scanDQ2 (sh, "mc15_13TeV:mc15_13TeV.361108.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau.recon.AOD.e3601_s2576_s2132_r6545/");

      // Set the name of the input TTree. It's always "CollectionTree"
      // for xAOD files.
      sh.setMetaString( "nc_tree", "CollectionTree" );

      // Print what we found:
      sh.print();

      // Create an EventLoop job:
      EL::Job job;
      job.sampleHandler( sh );

      // Add our analysis to the job:
      MyEventLoop* alg = new MyEventLoop();

      EL::PrunDriver driver;

      //const char* calibFilePath = gSystem->ExpandPathName("$ROOTCOREBIN/data/TESCalibration/EnergyCalibrationLC2012_retuned.root");
      //alg->m_filePath=calibFilePath;
      //alg->setCalibrationFilePath(calibFilePath);
      
      job.useXAOD();
      job.algsAdd( alg );

      driver.options()->setString("nc_outputSampleName", ("user."+jobname+".%in:name[2]%.ValidationHHtobbtt").c_str());
      if(nFiles != 0) driver.options()->setDouble("nc_nFiles", nFiles);
      driver.options()->setDouble("nc_mergeOutput", 1);
      driver.options()->setString(EL::Job::optGridNGBPerJob, "MAX");
      driver.options()->setDouble(EL::Job::optGridMaxNFilesPerJob, 50);

      driver.submit(job, "uniqueJobDirectory"); 

   }
   else {

     //const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/c/carquin/public/bbtautau/validation/samplesHadHad/");
     const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/c/carquin/public/bbtautau/validation/2HDMNLO/");
     //const char* inputFilePath = gSystem->ExpandPathName ("/tmp/carquin/");
     //const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/c/carquin/private/AtlasPhysicsMadGraphProd/sudir6/");
     //const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/d/dimicco/public/");
     //const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/c/carquin/public/bbtautau/validation/");
     SH::DiskListLocal list (inputFilePath);
     //SH::scanDir (sh, list, "DAOD_TRUTH1.341937*.root");
     //SH::scanDir (sh, list, "DAOD_TRUTH1.341938*.root");
     //SH::scanDir (sh, list, "Validation_*");
     
     //SH::scanDir (sh, list, "DAOD_TRUTH0.test_303351_hh.pool.root");
     //SH::scanDir (sh, list, "DAOD_TRUTH0.test_303*_lh.pool.root");
     //SH::scanDir (sh,list);
     //SH::scanDir (sh, list, "Validation_303*");
     SH::scanDir (sh, list, "DAOD_TRUTH1.Xhh_m*_ttbb_hh.pool.root*");
     //SH::scanDir (sh, list, "DAOD_TRUTH0.test.pool.root"); // specifying one particular file for testing

     //SH::DiskListEOS list ("/eos/atlas/user/c/carquin/TauCP", "root://eosatlas//eos/atlas/user/c/carquin/TauCP");
     //SH::scanDir (sh, list, "*023391.pool.root.*");

     //SH::DiskListEOS list ("/eos/atlas/user/m/mbecking/trigger/AOD/", "root://eosatlas//eos/atlas/user/m/mbecking/trigger/AOD/");
     //SH::scanDir (sh, list, "*05064171*.pool.root.*");
     //SH::scanDir (sh, list, "*05064171*.pool.root.*");

     // Set the name of the input TTree. It's always "CollectionTree"
     // for xAOD files.
     sh.setMetaString( "nc_tree", "CollectionTree" );

     // Print what we found:
     sh.print();

     // Create an EventLoop job:
     EL::Job job;
     job.sampleHandler( sh );

     // Add our analysis to the job:
     MyEventLoop* alg = new MyEventLoop();

     //const char* calibFilePath = gSystem->ExpandPathName("$ROOTCOREBIN/data/TESCalibration/EnergyCalibrationLC2012_retuned.root");

     //alg->m_filePath=calibFilePath;
     //alg->setCalibrationFilePath(calibFilePath);

     job.useXAOD();
     job.algsAdd( alg );
     // By default run only over 1000 events locally
     //job.options()->setDouble (EL::Job::optMaxEvents, 1000);

     // Run the job using the local/direct driver:
     EL::DirectDriver driver;

     //job.options()->setDouble (EL::Job::optMaxEvents, 10);

     driver.submit( job, submitDir );

   }

   return 0;
}
