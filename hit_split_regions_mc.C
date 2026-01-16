#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>
#include <limits>
#include <map>
#include <filesystem>
#include <memory>

void hit_split_regions_mc() {
    // Create output directory
    std::filesystem::create_directories("split_regions");
    
    // Load input files into RDataFrame
    std::vector<std::string> filenames;
    std::ifstream fileList("tfilelist_xrootd.txt");
    std::string line;
    while (std::getline(fileList, line)) {
        if (!line.empty()) filenames.push_back(line);
    }
    fileList.close();

    // Statistics accumulators
    size_t total_samples = 0;
    size_t total_events = 0;
    
    // TPC geometry definitions
    struct TPCRegion {
        std::string name;
        float y_min, y_max, z_min, z_max;
        int tpc_id;
    };
    
    // Define 8 regions (4 per TPC)
    std::vector<TPCRegion> regions = {
        {"TPC0_00", -203.732f, 0.0f, -5.68434e-14f,244.7f, 0},
        {"TPC0_01", 0.0f, 203.732f, -5.68434e-14f,244.7f, 0},
        {"TPC0_10", -203.732f, 0.0f,264.7f, 500.1f, 0},
        {"TPC0_11", 0.0f, 203.732f,264.7f, 500.1f, 0},
        {"TPC1_00", -203.732f, 0.0f, -5.68434e-14f,244.7f, 1},
        {"TPC1_01", 0.0f, 203.732f, -5.68434e-14f,244.7f, 1},
        {"TPC1_10", -203.732f, 0.0f,264.7f, 500.1f, 1},
        {"TPC1_11", 0.0f, 203.732f,264.7f, 500.1f, 1}
    };
    
    // Create separate CSV files for each region
    std::vector<std::unique_ptr<std::ofstream>> csv_files;
    for (const auto& region : regions) {
        std::string filename = "split_regions/" + region.name + "_hits_mc.csv";
        auto csv_file = std::make_unique<std::ofstream>(filename);
        *csv_file << "TrackID,Plane,TPC,TrackLength,ValidHits,MinX,MaxX,MinY,MaxY,MinZ,MaxZ,AvgPitch,HitEfficiency,AnodeTPC0_Hits,Cathode_Hits,AnodeTPC1_Hits,Other_Hits\n";
        csv_files.push_back(std::move(csv_file));
    }

    // Load dead channels from CSV
    std::set<std::tuple<unsigned short, unsigned short, unsigned short>> dead_channels;
    std::ifstream dead_csv("dead_channels.csv");
    if (dead_csv.is_open()) {
        std::string header;
        std::getline(dead_csv, header); // Skip header
        while (std::getline(dead_csv, line)) {
            std::stringstream ss(line);
            std::string wire_str, plane_str, tpc_str;
            std::getline(ss, wire_str, ',');
            std::getline(ss, plane_str, ',');
            std::getline(ss, tpc_str, ',');
            try {
                dead_channels.emplace(std::stoi(wire_str), std::stoi(plane_str), std::stoi(tpc_str));
            } catch (...) {
                std::cout << "Warning: Invalid line in dead_channels.csv: " << line << std::endl;
            }
        }
        dead_csv.close();
        std::cout << "Loaded " << dead_channels.size() << " dead channels" << std::endl;
    } else {
        std::cout << "Warning: Could not open dead_channels.csv, proceeding without dead channel filtering" << std::endl;
    }

    // Helper function to classify hits by X coordinate
    auto classify_x_region = [](float x) -> std::string {
        if (x >= -202.2f && x <= -152.2f) {
            return "anode_tpc0";
        } else if (x >= -50.0f && x <= 50.0f) {
            return "cathode";
        } else if (x >= 152.2f && x <= 202.2f) {
            return "anode_tpc1";
        } else {
            return "other";
        }
    };

    // Helper function to determine which region a hit belongs to
    auto get_region_index = [&](float x, float y, float z, unsigned short tpc_id) -> int {
        for (size_t i = 0; i < regions.size(); ++i) {
            const auto& region = regions[i];
            if (region.tpc_id == tpc_id &&
                y >= region.y_min && y < region.y_max &&
                z >= region.z_min && z < region.z_max) {
                return static_cast<int>(i);
            }
        }
        return -1; // No matching region
    };

    // Helper function to check for holes > 10 consecutive wires
    auto has_large_holes = [](const std::vector<unsigned short>& sorted_wires) {
        if (sorted_wires.size() < 2) return false;
        for (size_t i = 1; i < sorted_wires.size(); ++i) {
            if (sorted_wires[i] - sorted_wires[i-1] > 11) return true;
        }
        return false;
    };

    // Process hits for each region
    auto process_hits = [&](int trk_id, float track_length, 
                           ROOT::RVec<unsigned short> wires, 
                           ROOT::RVec<float> pitches, 
                           ROOT::RVec<unsigned short> tpcs,
                           ROOT::RVec<float> x_coords,
                           ROOT::RVec<float> y_coords, 
                           ROOT::RVec<float> z_coords,
                           ROOT::RVec<bool> ontraj,
                           int plane) {
        
        if (wires.empty()) return;

        // Group hits by region
        std::map<int, std::vector<size_t>> region_hits; // region_index -> hit indices
        
        for (size_t i = 0; i < wires.size() && i < x_coords.size() && i < y_coords.size() && 
             i < z_coords.size() && i < tpcs.size() && i < ontraj.size() && i < pitches.size(); ++i) {
            
            // Skip invalid hits
            if (!ontraj[i] || pitches[i] == -1.0f) continue;
            if (std::isnan(x_coords[i]) || std::isnan(y_coords[i]) || std::isnan(z_coords[i])) continue;
            
            // Skip dead channels
            if (dead_channels.find({wires[i], static_cast<unsigned short>(plane), tpcs[i]}) != dead_channels.end()) {
                continue;
            }
            
            int region_idx = get_region_index(x_coords[i], y_coords[i], z_coords[i], tpcs[i]);
            if (region_idx >= 0) {
                region_hits[region_idx].push_back(i);
            }
        }
        
        // Process each region
        for (const auto& [region_idx, hit_indices] : region_hits) {
            if (hit_indices.size() < 10) continue; // Skip regions with too few hits
            
            // Collect unique wires and check for large holes
            std::set<unsigned short> unique_wires;
            for (size_t idx : hit_indices) {
                unique_wires.insert(wires[idx]);
            }
            
            if (unique_wires.size() < 25) continue;
            
            std::vector<unsigned short> sorted_wires(unique_wires.begin(), unique_wires.end());
            std::sort(sorted_wires.begin(), sorted_wires.end());
            if (has_large_holes(sorted_wires)) continue;
            
            // Calculate statistics for this region
            float min_x = std::numeric_limits<float>::max();
            float max_x = std::numeric_limits<float>::lowest();
            float min_y = std::numeric_limits<float>::max();
            float max_y = std::numeric_limits<float>::lowest();
            float min_z = std::numeric_limits<float>::max();
            float max_z = std::numeric_limits<float>::lowest();
            float sum_pitch = 0.0f;
            int valid_pitches = 0;
            
            // X-coordinate region counters
            int anode_tpc0_hits = 0;
            int cathode_hits = 0;
            int anode_tpc1_hits = 0;
            int other_hits = 0;
            
            for (size_t idx : hit_indices) {
                // Update coordinate ranges
                min_x = std::min(min_x, x_coords[idx]);
                max_x = std::max(max_x, x_coords[idx]);
                min_y = std::min(min_y, y_coords[idx]);
                max_y = std::max(max_y, y_coords[idx]);
                min_z = std::min(min_z, z_coords[idx]);
                max_z = std::max(max_z, z_coords[idx]);
                
                // Sum pitches for average
                if (pitches[idx] > 0) {
                    sum_pitch += pitches[idx];
                    valid_pitches++;
                }
                
                // Classify hit by X coordinate
                std::string x_region = classify_x_region(x_coords[idx]);
                if (x_region == "anode_tpc0") {
                    anode_tpc0_hits++;
                } else if (x_region == "cathode") {
                    cathode_hits++;
                } else if (x_region == "anode_tpc1") {
                    anode_tpc1_hits++;
                } else {
                    other_hits++;
                }
            }
            
            float avg_pitch = valid_pitches > 0 ? sum_pitch / valid_pitches : 0.0f;
            
            // Compute hit efficiency
            float hit_eff = 0.0f;
            if (!sorted_wires.empty()) {
                unsigned short min_wire = sorted_wires.front();
                unsigned short max_wire = sorted_wires.back();
                int num_expected = static_cast<int>(max_wire - min_wire + 1);
                int num_dead = 0;
                unsigned short tpc = static_cast<unsigned short>(regions[region_idx].tpc_id);
                unsigned short pl = static_cast<unsigned short>(plane);
                for (unsigned short w = min_wire; w <= max_wire; ++w) {
                    if (dead_channels.count({w, pl, tpc})) {
                        ++num_dead;
                    }
                }
                int num_live = num_expected - num_dead;
                if (num_live > 0) {
                    hit_eff = static_cast<float>(unique_wires.size()) / static_cast<float>(num_live);
                }
            }
            
            // Write to appropriate region CSV (including the new X-coordinate flags)
            *csv_files[region_idx] << trk_id << "," << plane << "," << regions[region_idx].tpc_id << ","
                                  << track_length << "," << hit_indices.size() << ","
                                  << min_x << "," << max_x << ","
                                  << min_y << "," << max_y << ","
                                  << min_z << "," << max_z << ","
                                  << avg_pitch << "," << hit_eff << ","
                                  << anode_tpc0_hits << "," << cathode_hits << ","
                                  << anode_tpc1_hits << "," << other_hits << "\n";
            
            total_events++;
        }
    };

    // Process each ROOT file individually
    for (size_t file_idx = 0; file_idx < filenames.size(); ++file_idx) {
        ROOT::RDataFrame rdf_file("caloskim/TrackCaloSkim", {filenames[file_idx]});

        // Filter tracks with length > 50 cm
        auto rdf_filtered = rdf_file.Filter([](float length) { return length > 50.0; }, {"trk.length"});

        // Process data for each plane
        rdf_filtered.Foreach([&](int trk_id, float track_length,
                               // Plane 0 data
                               ROOT::RVec<unsigned short> wires0, ROOT::RVec<float> pitches0, ROOT::RVec<unsigned short> tpcs0,
                               ROOT::RVec<float> x0, ROOT::RVec<float> y0, ROOT::RVec<float> z0, ROOT::RVec<bool> ontraj0,
                               // Plane 1 data
                               ROOT::RVec<unsigned short> wires1, ROOT::RVec<float> pitches1, ROOT::RVec<unsigned short> tpcs1,
                               ROOT::RVec<float> x1, ROOT::RVec<float> y1, ROOT::RVec<float> z1, ROOT::RVec<bool> ontraj1,
                               // Plane 2 data
                               ROOT::RVec<unsigned short> wires2, ROOT::RVec<float> pitches2, ROOT::RVec<unsigned short> tpcs2,
                               ROOT::RVec<float> x2, ROOT::RVec<float> y2, ROOT::RVec<float> z2, ROOT::RVec<bool> ontraj2) {
                               
            process_hits(trk_id, track_length, wires0, pitches0, tpcs0, x0, y0, z0, ontraj0, 0);
            process_hits(trk_id, track_length, wires1, pitches1, tpcs1, x1, y1, z1, ontraj1, 1);
            process_hits(trk_id, track_length, wires2, pitches2, tpcs2, x2, y2, z2, ontraj2, 2);
            
        }, {"trk.id", "trk.length",
            // Plane 0
            "trk.hits0.h.wire", "trk.hits0.pitch", "trk.hits0.h.tpc",
            "trk.hits0.h.sp.x", "trk.hits0.h.sp.y", "trk.hits0.h.sp.z", "trk.hits0.ontraj",
            // Plane 1
            "trk.hits1.h.wire", "trk.hits1.pitch", "trk.hits1.h.tpc",
            "trk.hits1.h.sp.x", "trk.hits1.h.sp.y", "trk.hits1.h.sp.z", "trk.hits1.ontraj",
            // Plane 2
            "trk.hits2.h.wire", "trk.hits2.pitch", "trk.hits2.h.tpc",
            "trk.hits2.h.sp.x", "trk.hits2.h.sp.y", "trk.hits2.h.sp.z", "trk.hits2.ontraj"});

        // Update sample count and print checkpoint
        total_samples++;
        if (total_samples % 10 == 0) {
            std::cout << "Processed " << total_samples << " samples" << std::endl;
        }
    }

    // Close all CSV files and print region statistics
    std::cout << "\n=== Processing Statistics ===\n";
    std::cout << "Total samples processed: " << total_samples << "\n";
    std::cout << "Total events recorded: " << total_events << "\n";
    
    std::cout << "\n=== Region Statistics ===\n";
    for (size_t i = 0; i < regions.size(); ++i) {
        csv_files[i]->close();
        
        // Count lines in each file to report statistics
        std::ifstream count_file("split_regions/" + regions[i].name + "_hits_mc.csv");
        size_t line_count = 0;
        std::string temp_line;
        while (std::getline(count_file, temp_line)) {
            line_count++;
        }
        line_count--; // Subtract header line
        count_file.close();
        
        std::cout << regions[i].name << " (TPC " << regions[i].tpc_id << ", y:[" 
                  << regions[i].y_min << "," << regions[i].y_max << "], z:[" 
                  << regions[i].z_min << "," << regions[i].z_max << "]): " 
                  << line_count << " entries\n";
    }
    
    std::cout << "\nCSV files saved in split_regions/ directory\n";
}