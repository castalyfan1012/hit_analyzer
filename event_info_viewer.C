#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cmath>

void display_event_info(const std::string& file_path) {
    // Open output file
    std::ofstream out_file("event_info.txt");
    if (!out_file.is_open()) {
        std::cout << "Error: Cannot open event_info.txt for writing" << std::endl;
        return;
    }

    // Load dead channels from CSV
    std::set<std::tuple<unsigned short, unsigned short, unsigned short>> dead_channels;
    std::ifstream dead_csv("dead_channels.csv");
    bool has_dead_channels = false;
    
    if (dead_csv.is_open()) {
        has_dead_channels = true;
        std::string header, line;
        std::getline(dead_csv, header);
        while (std::getline(dead_csv, line)) {
            std::stringstream ss(line);
            std::string wire_str, plane_str, tpc_str;
            std::getline(ss, wire_str, ',');
            std::getline(ss, plane_str, ',');
            std::getline(ss, tpc_str, ',');
            try {
                dead_channels.emplace(std::stoi(wire_str), std::stoi(plane_str), std::stoi(tpc_str));
            } catch (...) {
                out_file << "Warning: Invalid line in dead_channels.csv: " << line << std::endl;
            }
        }
        dead_csv.close();
        out_file << "Loaded " << dead_channels.size() << " dead channels for efficiency calculation" << std::endl;
    } else {
        out_file << "Note: dead_channels.csv not found - efficiency calculation will be skipped" << std::endl;
    }

    // Define required columns
    std::vector<std::string> required_columns = {
        "trk.meta.evt", "trk.meta.run", "trk.meta.subrun", "trk.id", "trk.length",
        "trk.start.x", "trk.start.y", "trk.start.z", "trk.end.x", "trk.end.y", "trk.end.z",
        "trk.hits0.h.wire", "trk.hits0.pitch", "trk.hits0.h.tpc", "trk.hits0.ontraj",
        "trk.hits0.h.time", "trk.hits0.h.id",
        "trk.hits1.h.wire", "trk.hits1.pitch", "trk.hits1.h.tpc", "trk.hits1.ontraj",
        "trk.hits1.h.time", "trk.hits1.h.id",
        "trk.hits2.h.wire", "trk.hits2.pitch", "trk.hits2.h.tpc", "trk.hits2.ontraj",
        "trk.hits2.h.time", "trk.hits2.h.id"
    };

    // Validate columns in ROOT file
    TFile file(file_path.c_str(), "READ");
    if (!file.IsOpen()) {
        out_file << "Error: Cannot open ROOT file: " << file_path << std::endl;
        out_file.close();
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(file.Get("caloskim/TrackCaloSkim"));
    if (!tree) {
        out_file << "Error: Cannot find tree 'caloskim/TrackCaloSkim' in file: " << file_path << std::endl;
        file.Close();
        out_file.close();
        return;
    }

    std::vector<std::string> missing_columns;
    for (const auto& col : required_columns) {
        if (!tree->GetBranch(col.c_str())) {
            missing_columns.push_back(col);
        }
    }
    file.Close();

    if (!missing_columns.empty()) {
        out_file << "Error: The following required columns are missing in the ROOT file:" << std::endl;
        for (const auto& col : missing_columns) {
            out_file << "  - " << col << std::endl;
        }
        out_file << "Please check the ROOT file schema and update the column names in the script." << std::endl;
        out_file.close();
        return;
    }

    // Load the ROOT file
    ROOT::RDataFrame rdf("caloskim/TrackCaloSkim", {file_path});
    
    // Helper function to check for holes > 10 consecutive wires
    auto has_large_holes = [](const std::vector<unsigned short>& sorted_wires) {
        if (sorted_wires.size() < 2) return false;
        for (size_t i = 1; i < sorted_wires.size(); ++i) {
            if (sorted_wires[i] - sorted_wires[i-1] > 11) return true;
        }
        return false;
    };

    // Function to calculate efficiency for a plane
    auto calculate_plane_efficiency = [&](const ROOT::RVec<unsigned short>& wires, 
                                         const ROOT::RVec<float>& pitches,
                                         const ROOT::RVec<unsigned short>& tpcs,
                                         const ROOT::RVec<bool>& ontraj,
                                         int plane) -> std::pair<float, float> {
        if (!has_dead_channels || wires.empty()) return {-1.0f, -1.0f};
        
        // Filter valid hits (on-trajectory, valid pitch, not dead)
        std::vector<unsigned short> valid_wires;
        std::vector<float> valid_pitches;
        std::vector<unsigned short> valid_tpcs;
        
        for (size_t i = 0; i < wires.size() && i < pitches.size() && i < tpcs.size() && i < ontraj.size(); ++i) {
            if (ontraj[i] && pitches[i] > 0.0f && 
                dead_channels.find({wires[i], static_cast<unsigned short>(plane), tpcs[i]}) == dead_channels.end()) {
                valid_wires.push_back(wires[i]);
                valid_pitches.push_back(pitches[i]);
                valid_tpcs.push_back(tpcs[i]);
            }
        }
        
        if (valid_wires.size() < 25) return {-1.0f, -1.0f};
        
        std::set<unsigned short> unique_wires(valid_wires.begin(), valid_wires.end());
        std::vector<unsigned short> sorted_wires(unique_wires.begin(), unique_wires.end());
        std::sort(sorted_wires.begin(), sorted_wires.end());
        
        if (has_large_holes(sorted_wires)) return {-1.0f, -1.0f};
        
        // Calculate non-dead wires per TPC
        std::set<unsigned short> unique_tpcs(valid_tpcs.begin(), valid_tpcs.end());
        int n_non_dead_wires = 0;
        std::set<unsigned short> unique_hit_wires;
        
        for (unsigned short tpc_id : unique_tpcs) {
            unsigned short min_wire = std::numeric_limits<unsigned short>::max();
            unsigned short max_wire = 0;
            for (size_t i = 0; i < valid_wires.size(); ++i) {
                if (valid_tpcs[i] == tpc_id) {
                    min_wire = std::min(min_wire, valid_wires[i]);
                    max_wire = std::max(max_wire, valid_wires[i]);
                    unique_hit_wires.insert(valid_wires[i]);
                }
            }
            if (min_wire <= max_wire) {
                for (unsigned short wire = min_wire; wire <= max_wire; ++wire) {
                    if (dead_channels.find({wire, static_cast<unsigned short>(plane), tpc_id}) == dead_channels.end()) {
                        n_non_dead_wires++;
                    }
                }
            }
        }
        
        if (n_non_dead_wires < 25) return {-1.0f, -1.0f};
        
        // Calculate efficiency and average pitch
        float efficiency = n_non_dead_wires > 0 ? std::min(static_cast<float>(unique_hit_wires.size()) / n_non_dead_wires, 1.0f) : 0.0f;
        
        float avg_pitch = 0.0f;
        if (!valid_pitches.empty()) {
            float sum_pitch = 0.0f;
            int n_valid_pitches = 0;
            for (float p : valid_pitches) {
                if (p > 0) {
                    sum_pitch += p;
                    n_valid_pitches++;
                }
            }
            avg_pitch = n_valid_pitches > 0 ? sum_pitch / n_valid_pitches : 0.0f;
        }
        
        return {efficiency, avg_pitch};
    };

    // Display detailed information for each track
    out_file << "\n" << std::string(80, '=') << "\n";
    out_file << "DETAILED TRACK INFORMATION\n";
    out_file << "File: " << file_path << "\n";
    out_file << std::string(80, '=') << "\n";

    int track_count = 0;
    
    rdf.Foreach([&](int evt, int run, int subrun, int trk_id, float trk_length,
                   float start_x, float start_y, float start_z,
                   float end_x, float end_y, float end_z,
                   ROOT::RVec<unsigned short> wires0, ROOT::RVec<float> pitches0, 
                   ROOT::RVec<unsigned short> tpcs0, ROOT::RVec<bool> ontraj0,
                   ROOT::RVec<float> times0, ROOT::RVec<int> hit_ids0,
                   ROOT::RVec<unsigned short> wires1, ROOT::RVec<float> pitches1,
                   ROOT::RVec<unsigned short> tpcs1, ROOT::RVec<bool> ontraj1,
                   ROOT::RVec<float> times1, ROOT::RVec<int> hit_ids1,
                   ROOT::RVec<unsigned short> wires2, ROOT::RVec<float> pitches2,
                   ROOT::RVec<unsigned short> tpcs2, ROOT::RVec<bool> ontraj2,
                   ROOT::RVec<float> times2, ROOT::RVec<int> hit_ids2) {
        
        track_count++;
        
        out_file << "\n" << std::string(60, '-') << "\n";
        out_file << "TRACK #" << track_count << "\n";
        out_file << std::string(60, '-') << "\n";
        
        out_file << std::fixed << std::setprecision(3);
        out_file << "Event Info:\n";
        out_file << "  Run: " << run << ", Subrun: " << subrun << ", Event: " << evt << "\n";
        out_file << "  Track ID: " << trk_id << "\n";
        out_file << "  Track Length: " << std::setprecision(2) << trk_length << " cm\n";
        
        out_file << "\nTrack Endpoints:\n";
        out_file << "  Start: (" << start_x << ", " << start_y << ", " << start_z << ") cm\n";
        out_file << "  End:   (" << end_x << ", " << end_y << ", " << end_z << ") cm\n";
        
        float dx = end_x - start_x;
        float dy = end_y - start_y; 
        float dz = end_z - start_z;
        float length = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (length > 0) {
            out_file << "  Direction: (" << dx/length << ", " << dy/length << ", " << dz/length << ")\n";
        }
        
        std::vector<std::tuple<int, ROOT::RVec<unsigned short>, ROOT::RVec<float>, 
                              ROOT::RVec<unsigned short>, ROOT::RVec<bool>, 
                              ROOT::RVec<float>, ROOT::RVec<int>>> plane_data = {
            {0, wires0, pitches0, tpcs0, ontraj0, times0, hit_ids0},
            {1, wires1, pitches1, tpcs1, ontraj1, times1, hit_ids1},
            {2, wires2, pitches2, tpcs2, ontraj2, times2, hit_ids2}
        };
        
        for (const auto& [plane, wires, pitches, tpcs, ontraj, times, hit_ids] : plane_data) {
            out_file << "\n--- Plane " << plane << " ---\n";
            out_file << "Total hits: " << wires.size() << "\n";
            
            if (wires.empty()) {
                out_file << "No hits on this plane\n";
                continue;
            }
            
            int on_traj_count = 0;
            for (bool ot : ontraj) {
                if (ot) on_traj_count++;
            }
            out_file << "On-trajectory hits: " << on_traj_count << "\n";
            
            auto [min_wire, max_wire] = std::minmax_element(wires.begin(), wires.end());
            std::set<unsigned short> unique_wires(wires.begin(), wires.end());
            out_file << "Wire range: " << *min_wire << " - " << *max_wire 
                     << " (span: " << (*max_wire - *min_wire + 1) << ", unique: " << unique_wires.size() << ")\n";
            
            if (!tpcs.empty()) {
                std::set<unsigned short> unique_tpcs(tpcs.begin(), tpcs.end());
                out_file << "TPCs: ";
                for (auto tpc : unique_tpcs) {
                    out_file << tpc << " ";
                }
                out_file << "\n";
            }
            
            std::vector<float> valid_pitches;
            for (size_t i = 0; i < pitches.size(); ++i) {
                if (ontraj[i] && pitches[i] > 0) {
                    valid_pitches.push_back(pitches[i]);
                }
            }
            
            if (!valid_pitches.empty()) {
                auto [min_pitch, max_pitch] = std::minmax_element(valid_pitches.begin(), valid_pitches.end());
                float sum_pitch = 0.0f;
                for (float p : valid_pitches) sum_pitch += p;
                float avg_pitch = sum_pitch / valid_pitches.size();
                
                out_file << std::fixed << std::setprecision(6);
                out_file << "Pitch: min=" << *min_pitch << ", max=" << *max_pitch 
                         << ", avg=" << avg_pitch << " cm (" << valid_pitches.size() << " valid)\n";
                
                auto [efficiency, eff_avg_pitch] = calculate_plane_efficiency(wires, pitches, tpcs, ontraj, plane);
                if (efficiency >= 0) {
                    out_file << "Efficiency: " << efficiency*100 << "% (avg pitch: " << eff_avg_pitch << " cm)\n";
                }
            }
            
            if (!times.empty()) {
                auto [min_time, max_time] = std::minmax_element(times.begin(), times.end());
                out_file << "Hit time range: " << *min_time << " - " << *max_time << " μs\n";
            } else {
                out_file << "No hit times available on this plane\n";
            }
            
            std::vector<int> valid_hit_indices;
            for (int i = 0; i < wires.size(); ++i) {
                if (i < ontraj.size() && ontraj[i] && i < pitches.size() && pitches[i] > 0) {
                    valid_hit_indices.push_back(i);
                }
            }
            
            if (!valid_hit_indices.empty()) {
                out_file << "All " << valid_hit_indices.size() << " valid hits (on-trajectory, pitch > 0):\n";
                out_file << "  #    Wire  TPC  Pitch     Time      HitID\n";
                for (size_t idx = 0; idx < valid_hit_indices.size(); ++idx) {
                    int i = valid_hit_indices[idx];
                    out_file << "  " << std::setw(3) << (idx+1);
                    out_file << "  " << std::setw(4) << wires[i];
                    out_file << "  " << std::setw(3) << (tpcs.size() > i ? tpcs[i] : 0);
                    out_file << "  " << std::setw(9) << std::setprecision(6) << (pitches.size() > i ? pitches[i] : -1.0f);
                    out_file << "  " << std::setw(9) << std::setprecision(3) << (times.size() > i ? times[i] : -1.0f);
                    out_file << "  " << (hit_ids.size() > i ? hit_ids[i] : -1);
                    out_file << "\n";
                }
            } else {
                out_file << "No valid hits found on this plane\n";
            }
        }
        
    }, {"trk.meta.evt", "trk.meta.run", "trk.meta.subrun", "trk.id", "trk.length",
        "trk.start.x", "trk.start.y", "trk.start.z", "trk.end.x", "trk.end.y", "trk.end.z",
        "trk.hits0.h.wire", "trk.hits0.pitch", "trk.hits0.h.tpc", "trk.hits0.ontraj",
        "trk.hits0.h.time", "trk.hits0.h.id",
        "trk.hits1.h.wire", "trk.hits1.pitch", "trk.hits1.h.tpc", "trk.hits1.ontraj",
        "trk.hits1.h.time", "trk.hits1.h.id",
        "trk.hits2.h.wire", "trk.hits2.pitch", "trk.hits2.h.tpc", "trk.hits2.ontraj",
        "trk.hits2.h.time", "trk.hits2.h.id"});
        
    out_file << "\n" << std::string(60, '=') << "\n";
    out_file << "SUMMARY: Displayed information for " << track_count << " tracks\n";
    out_file << std::string(60, '=') << "\n";
    out_file.close();
}

void event_info_viewer() {
    std::string file_path;
    std::cout << "Enter the ROOT file path: ";
    std::getline(std::cin, file_path);
    
    if (file_path.empty()) {
        std::cout << "No file path provided. Using default example..." << std::endl;
        file_path = "/pnfs/sbn/data_add/sbn_nd/poms_production/data/MCP2025B_02/v10_06_00_02/DevSample_1e20/reco2/bnblight/60/hist_reco2_reco1_filtered_decoded-raw_filtered_data_EventBuilder5_art1_run18351_13_strmBNBLight_20250325T113833-60726fda-f325-e401-fe4f-f8b24274d492.root";
    }
    
    try {
        display_event_info(file_path);
    } catch (const std::exception& e) {
        std::cout << "Error processing file: " << e.what() << std::endl;
        std::cout << "Please check that the file path is correct and accessible." << std::endl;
    }
}

void show_event_info(const std::string& file_path) {
    display_event_info(file_path);
}