#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TH1F.h>
#include <TFile.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

void dead_wires()
{
  // Load input files into RDataFrame
  std::vector<std::string> filenames;
  std::ifstream fileList("filelist_xrootd_small.txt");
  std::string line;
  while (std::getline(fileList, line)) {
    if (!line.empty()) {
      filenames.push_back(line);
    }
  }
  fileList.close();
  ROOT::RDataFrame rdf("caloskim/TrackCaloSkim", filenames);

  // Define track length column and filter
  auto rdf_filtered = rdf.Define("track_length_cm", [](float start_x, float start_y, float start_z,
                                                       float end_x, float end_y, float end_z) {
    float dx = end_x - start_x;
    float dy = end_y - start_y;
    float dz = end_z - start_z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }, {"trk.start.x", "trk.start.y", "trk.start.z", "trk.end.x", "trk.end.y", "trk.end.z"})
  .Filter([](float length) { return length > 50.0; }, {"track_length_cm"});

  // Count total entries for progress tracking
  auto n_entries = rdf_filtered.Count();
  Long64_t total_entries = *n_entries;
  Long64_t processed_entries = 0;
  Long64_t progress_step = total_entries / 10; // Update progress every 10%

  // Determine min and max wire numbers for each plane
  auto min_wire0 = rdf_filtered.Min<ROOT::RVec<unsigned short>>("trk.hits0.h.wire");
  auto max_wire0 = rdf_filtered.Max<ROOT::RVec<unsigned short>>("trk.hits0.h.wire");
  auto min_wire1 = rdf_filtered.Min<ROOT::RVec<unsigned short>>("trk.hits1.h.wire");
  auto max_wire1 = rdf_filtered.Max<ROOT::RVec<unsigned short>>("trk.hits1.h.wire");
  auto min_wire2 = rdf_filtered.Min<ROOT::RVec<unsigned short>>("trk.hits2.h.wire");
  auto max_wire2 = rdf_filtered.Max<ROOT::RVec<unsigned short>>("trk.hits2.h.wire");
  int min_wires0 = *min_wire0;
  int max_wires0 = *max_wire0 + 1; // Add 1 for histogram binning
  int min_wires1 = *min_wire1;
  int max_wires1 = *max_wire1 + 1;
  int min_wires2 = *min_wire2;
  int max_wires2 = *max_wire2 + 1;

  // --- Separate TPC Results ---
  // Define histograms for each TPC and plane with individual wire ranges
  TH1F* h_wires_tpc0_plane0 = new TH1F("h_wires_tpc0_plane0", "Hit Wires TPC 0 Plane 0;Wire Number;Entries", max_wires0 - min_wires0, min_wires0, max_wires0);
  TH1F* h_wires_tpc0_plane1 = new TH1F("h_wires_tpc0_plane1", "Hit Wires TPC 0 Plane 1;Wire Number;Entries", max_wires1 - min_wires1, min_wires1, max_wires1);
  TH1F* h_wires_tpc0_plane2 = new TH1F("h_wires_tpc0_plane2", "Hit Wires TPC 0 Plane 2;Wire Number;Entries", max_wires2 - min_wires2, min_wires2, max_wires2);
  TH1F* h_wires_tpc1_plane0 = new TH1F("h_wires_tpc1_plane0", "Hit Wires TPC 1 Plane 0;Wire Number;Entries", max_wires0 - min_wires0, min_wires0, max_wires0);
  TH1F* h_wires_tpc1_plane1 = new TH1F("h_wires_tpc1_plane1", "Hit Wires TPC 1 Plane 1;Wire Number;Entries", max_wires1 - min_wires1, min_wires1, max_wires1);
  TH1F* h_wires_tpc1_plane2 = new TH1F("h_wires_tpc1_plane2", "Hit Wires TPC 1 Plane 2;Wire Number;Entries", max_wires2 - min_wires2, min_wires2, max_wires2);

  // --- Combined TPC Results ---
  // Define histograms for combined TPCs (3 planes)
  TH1F* h_wires_plane0 = new TH1F("h_wires_plane0", "Hit Wires Plane 0 (Combined TPC);Wire Number;Entries", max_wires0 - min_wires0, min_wires0, max_wires0);
  TH1F* h_wires_plane1 = new TH1F("h_wires_plane1", "Hit Wires Plane 1 (Combined TPC);Wire Number;Entries", max_wires1 - min_wires1, min_wires1, max_wires1);
  TH1F* h_wires_plane2 = new TH1F("h_wires_plane2", "Hit Wires Plane 2 (Combined TPC);Wire Number;Entries", max_wires2 - min_wires2, min_wires2, max_wires2);

  // Fill histograms for both separate and combined TPC results
  rdf_filtered.Foreach([&](ROOT::RVec<unsigned short> wires0, ROOT::RVec<unsigned short> tpcs0, ROOT::RVec<unsigned short> planes0,
                           ROOT::RVec<unsigned short> wires1, ROOT::RVec<unsigned short> tpcs1, ROOT::RVec<unsigned short> planes1,
                           ROOT::RVec<unsigned short> wires2, ROOT::RVec<unsigned short> tpcs2, ROOT::RVec<unsigned short> planes2) {
    // Process plane 0
    for (size_t i = 0; i < wires0.size() && i < tpcs0.size() && i < planes0.size(); ++i) {
      if (planes0[i] == 0 && wires0[i] >= min_wires0 && wires0[i] < max_wires0) {
        if (tpcs0[i] == 0) h_wires_tpc0_plane0->Fill(wires0[i]);
        else if (tpcs0[i] == 1) h_wires_tpc1_plane0->Fill(wires0[i]);
        h_wires_plane0->Fill(wires0[i]); // Combined TPC
      }
    }
    // Process plane 1
    for (size_t i = 0; i < wires1.size() && i < tpcs1.size() && i < planes1.size(); ++i) {
      if (planes1[i] == 1 && wires1[i] >= min_wires1 && wires1[i] < max_wires1) {
        if (tpcs1[i] == 0) h_wires_tpc0_plane1->Fill(wires1[i]);
        else if (tpcs1[i] == 1) h_wires_tpc1_plane1->Fill(wires1[i]);
        h_wires_plane1->Fill(wires1[i]); // Combined TPC
      }
    }
    // Process plane 2
    for (size_t i = 0; i < wires2.size() && i < tpcs2.size() && i < planes2.size(); ++i) {
      if (planes2[i] == 2 && wires2[i] >= min_wires2 && wires2[i] < max_wires2) {
        if (tpcs2[i] == 0) h_wires_tpc0_plane2->Fill(wires2[i]);
        else if (tpcs2[i] == 1) h_wires_tpc1_plane2->Fill(wires2[i]);
        h_wires_plane2->Fill(wires2[i]); // Combined TPC
      }
    }
    // Progress update
    processed_entries++;
    if (progress_step > 0 && processed_entries % progress_step == 0) {
      std::cout << "Processed " << (processed_entries * 100 / total_entries) << "% of entries (" 
                << processed_entries << "/" << total_entries << ")\n";
    }
  }, {"trk.hits0.h.wire", "trk.hits0.h.tpc", "trk.hits0.h.plane",
      "trk.hits1.h.wire", "trk.hits1.h.tpc", "trk.hits1.h.plane",
      "trk.hits2.h.wire", "trk.hits2.h.tpc", "trk.hits2.h.plane"});

  // Final progress update
  std::cout << "Processed 100% of entries (" << processed_entries << "/" << total_entries << ")\n";

  // --- Separate TPC: Identify dead channels (bins with 0 entries) ---
  std::vector<int> dead_channels_tpc0_plane0, dead_channels_tpc0_plane1, dead_channels_tpc0_plane2;
  std::vector<int> dead_channels_tpc1_plane0, dead_channels_tpc1_plane1, dead_channels_tpc1_plane2;
  for (int i = 1; i <= h_wires_tpc0_plane0->GetNbinsX(); ++i) {
    int wire_num = min_wires0 + (i - 1); // Convert bin to wire number
    if (h_wires_tpc0_plane0->GetBinContent(i) == 0 && wire_num >= min_wires0 && wire_num < max_wires0) dead_channels_tpc0_plane0.push_back(wire_num);
    if (h_wires_tpc0_plane1->GetBinContent(i) == 0 && wire_num >= min_wires1 && wire_num < max_wires1) dead_channels_tpc0_plane1.push_back(wire_num);
    if (h_wires_tpc0_plane2->GetBinContent(i) == 0 && wire_num >= min_wires2 && wire_num < max_wires2) dead_channels_tpc0_plane2.push_back(wire_num);
    if (h_wires_tpc1_plane0->GetBinContent(i) == 0 && wire_num >= min_wires0 && wire_num < max_wires0) dead_channels_tpc1_plane0.push_back(wire_num);
    if (h_wires_tpc1_plane1->GetBinContent(i) == 0 && wire_num >= min_wires1 && wire_num < max_wires1) dead_channels_tpc1_plane1.push_back(wire_num);
    if (h_wires_tpc1_plane2->GetBinContent(i) == 0 && wire_num >= min_wires2 && wire_num < max_wires2) dead_channels_tpc1_plane2.push_back(wire_num);
  }

  // Save dead channels to CSV (Separate TPC)
  std::ofstream dead_file_tpc("dead_channels.csv");
  if (dead_file_tpc.is_open()) {
    dead_file_tpc << "Wire,Plane,TPC\n";
    for (auto ch : dead_channels_tpc0_plane0) dead_file_tpc << ch << ",0,0\n";
    for (auto ch : dead_channels_tpc0_plane1) dead_file_tpc << ch << ",1,0\n";
    for (auto ch : dead_channels_tpc0_plane2) dead_file_tpc << ch << ",2,0\n";
    for (auto ch : dead_channels_tpc1_plane0) dead_file_tpc << ch << ",0,1\n";
    for (auto ch : dead_channels_tpc1_plane1) dead_file_tpc << ch << ",1,1\n";
    for (auto ch : dead_channels_tpc1_plane2) dead_file_tpc << ch << ",2,1\n";
    dead_file_tpc.close();
  } else {
    std::cout << "Error: Could not open dead_channels.csv for writing\n";
  }


  // Save histograms to ROOT file
  TFile* out_file = new TFile("hit_wires.root", "RECREATE");
  h_wires_tpc0_plane0->Write();
  h_wires_tpc0_plane1->Write();
  h_wires_tpc0_plane2->Write();
  h_wires_tpc1_plane0->Write();
  h_wires_tpc1_plane1->Write();
  h_wires_tpc1_plane2->Write();
  h_wires_plane0->Write();
  h_wires_plane1->Write();
  h_wires_plane2->Write();
  out_file->Close();
}