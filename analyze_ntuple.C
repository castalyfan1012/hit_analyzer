#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <cmath>

void analyze_ntuple()
{
  ////////////////////////////////////
  /// Load input files into RDataFrame
  /// (This is the ROOT equivalent of
  ///  a pandas dataframe)
  ////////////////////////////////////
  std::vector<std::string> filenames;
  const char* fileListPath = "filelist_xrootd_small.txt";
  //const char* fileListPath = "filelist_small.txt";

  std::ifstream fileList(fileListPath);
  
  std::string line;
  while (std::getline(fileList, line)) {
    if (!line.empty()) {
      filenames.push_back(line);
    }
  }
  fileList.close();
  ROOT::RDataFrame rdf("caloskim/TrackCaloSkim", filenames);

  ////////////////////////////////////
  /// Define new columns
  /// - Track length based on trk.start and trk.end
  /// - Number of hits per wire based on trk.hits0.h.wire size
  ////////////////////////////////////
  auto rdf_with_length = rdf.Define("track_length_cm", [](float start_x, float start_y, float start_z,
                                                          float end_x, float end_y, float end_z) {
    float dx = end_x - start_x;
    float dy = end_y - start_y;
    float dz = end_z - start_z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }, {"trk.start.x", "trk.start.y", "trk.start.z", "trk.end.x", "trk.end.y", "trk.end.z"})
  .Define("num_hits_per_wire", [](ROOT::RVec<unsigned short> wires) {
    return static_cast<float>(wires.size());
  }, {"trk.hits0.h.wire"});

  ////////////////////////////////////
  /// Print contents of RDataFrame
  /// --> only showing a few columns
  ///     (track end points)
  ////////////////////////////////////
  rdf_with_length.Display({
    /* "track_length_cm", */
    /* "num_hits_per_wire", */
    /* "trk.length", */
    "trk.meta.evt",
    "trk.hits0.h.plane",
    "trk.hits0.pitch",
    "trk.hits0.h.time",
    "trk.end.x",
    "trk.hits0.h.sp.x",
    "trk.hits0.tp.x",
    "trk.hits0.h.wire",
    "trk.id",
    "trk.hits0.h.id",
    "trk.hits0.h.tpc",
    "trk.hits0.ontraj",
    "trk.end.x",
  }, 11)->Print();

  ////////////////////////////////////
  /// Do some analysis
  /// --> Example: Loop over 
  ///     collection plane hit wires
  ////////////////////////////////////
  rdf_with_length.Foreach([&](ROOT::RVec<unsigned short> col_wires) {
      for(int i = 0; i < col_wires.size(); i++)
      { 
        // Do something here
        //std::cout << col_wires[i] << std::endl;
      }
    }, {"trk.hits2.h.wire"});
}