#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>
#include <limits>
#include <thread>
#include <mutex>
#include <atomic>

void hit_analyzer() {
    // ============================================================================
    // CONFIGURATION: Manually set input file lists and output CSV names here
    // ============================================================================
    std::string data_filelist = "filelist_xrootd_data.txt";        // Data file list
    std::string mc_filelist = "filelist_xrootd_mc.txt";         // MC file list

    std::string data_output_csv = "hiteff_data.csv";  // Data output CSV
    std::string mc_output_csv = "hiteff_mc.csv"; // MC output CSV
    // ============================================================================

    // Enable ROOT implicit multithreading ONCE at the start (before any threads)
    ROOT::EnableImplicitMT();

    // Load dead channels from CSV
    std::set<std::tuple<unsigned short, unsigned short, unsigned short>> dead_channels;
    std::ifstream dead_csv("dead_channels.csv");
    std::string line;
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
        std::cout << "Error: Could not open dead_channels.csv" << std::endl;
    }

    // Helper function to check for holes > 10 consecutive wires
    auto has_large_holes = [](const std::vector<unsigned short>& sorted_wires) {
        if (sorted_wires.size() < 2) return false;
        for (size_t i = 1; i < sorted_wires.size(); ++i) {
            if (sorted_wires[i] - sorted_wires[i-1] > 11) return true;
        }
        return false;
    };

    // Mutex for thread-safe console output
    std::mutex cout_mutex;

    // Lambda to process a dataset (data or MC)
    auto process_dataset = [&](const std::string& filelist_name, const std::string& output_csv_name, const std::string& dataset_type) {
        {
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << "\n=== Processing " << dataset_type << " dataset ===" << std::endl;
        }
        
        // Load input files
        std::vector<std::string> filenames;
        std::ifstream fileList(filelist_name);
        if (!fileList.is_open()) {
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << "Error: Could not open " << filelist_name << std::endl;
            return;
        }
        std::string file_line;
        while (std::getline(fileList, file_line)) {
            if (!file_line.empty()) filenames.push_back(file_line);
        }
        fileList.close();
        {
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << "[" << dataset_type << "] Loaded " << filenames.size() << " files from " << filelist_name << std::endl;
        }

        // Statistics accumulators
        std::atomic<size_t> total_samples{0};
        std::atomic<size_t> total_events{0};
        float min_pitch = std::numeric_limits<float>::max();
        float max_pitch = std::numeric_limits<float>::lowest();
        size_t min_wires = std::numeric_limits<size_t>::max();
        size_t max_wires = 0;
        size_t min_hits = std::numeric_limits<size_t>::max();
        size_t max_hits = 0;
        double total_efficiency = 0.0;
        size_t efficiency_count = 0;

        // Open CSV file for output
        std::ofstream csv_out(output_csv_name);
        csv_out << "TrackID,Plane,TPC,TrackLength,ValidHits,NonDeadWires,Efficiency,AvgPitch\n";

        // Mutex for thread-safe CSV writing and statistics updates
        std::mutex csv_mutex;

        // Calculate efficiency and average pitch
        auto calculate_efficiency = [&](int trk_id, float track_length, ROOT::RVec<unsigned short> wires, 
                                       ROOT::RVec<float> pitches, ROOT::RVec<unsigned short> tpcs, int plane) {
            if (wires.empty()) return;

            std::set<unsigned short> unique_wires(wires.begin(), wires.end());
            if (unique_wires.size() < 25) return;

            std::vector<unsigned short> sorted_wires(unique_wires.begin(), unique_wires.end());
            std::sort(sorted_wires.begin(), sorted_wires.end());
            if (has_large_holes(sorted_wires)) return;

            unsigned short min_wire = sorted_wires.front();
            unsigned short max_wire = sorted_wires.back();
            unsigned short tpc_id = tpcs.empty() ? 0 : tpcs[0];

            int n_non_dead_wires = 0;
            for (unsigned short wire = min_wire; wire <= max_wire; ++wire) {
                if (dead_channels.find({wire, static_cast<unsigned short>(plane), tpc_id}) == dead_channels.end()) {
                    n_non_dead_wires++;
                }
            }
            if (n_non_dead_wires < 25) return;

            std::set<unsigned short> unique_hit_wires;
            for (unsigned short wire : wires) {
                if (dead_channels.find({wire, static_cast<unsigned short>(plane), tpc_id}) == dead_channels.end()) {
                    unique_hit_wires.insert(wire);
                }
            }
            int n_valid_hits = unique_hit_wires.size();

            float efficiency = n_non_dead_wires > 0 ? static_cast<float>(n_valid_hits) / n_non_dead_wires : 0.0;
            float avg_pitch = 0.0;
            if (!pitches.empty()) {
                float sum_pitch = 0;
                int n_valid_pitches = 0;
                for (float p : pitches) {
                    if (p > 0) {
                        sum_pitch += p;
                        n_valid_pitches++;
                    }
                }
                avg_pitch = n_valid_pitches > 0 ? sum_pitch / n_valid_pitches : 0.0;
            }

            // Update statistics (thread-safe)
            {
                std::lock_guard<std::mutex> lock(csv_mutex);
                total_events++;
                if (avg_pitch > 0) {
                    min_pitch = std::min(min_pitch, avg_pitch);
                    max_pitch = std::max(max_pitch, avg_pitch);
                    total_efficiency += efficiency;
                    efficiency_count++;
                }
                min_wires = std::min(min_wires, unique_wires.size());
                max_wires = std::max(max_wires, unique_wires.size());
                min_hits = std::min(min_hits, static_cast<size_t>(n_valid_hits));
                max_hits = std::max(max_hits, static_cast<size_t>(n_valid_hits));

                if (avg_pitch > 0) {
                    csv_out << trk_id << "," << plane << "," << tpc_id << ","
                            << track_length << "," << n_valid_hits << "," << n_non_dead_wires << ","
                            << efficiency << "," << avg_pitch << "\n";
                }
            }
        };

        // Process each ROOT file individually to track sample progress
        for (size_t file_idx = 0; file_idx < filenames.size(); ++file_idx) {
            ROOT::RDataFrame rdf_file("caloskim/TrackCaloSkim", {filenames[file_idx]});

            // Filter tracks with length > 50 cm
            auto rdf_filtered = rdf_file.Filter([](float length) { return length > 50.0; }, {"trk.length"});

            // Define valid hits for each plane
            auto rdf_cleaned = rdf_filtered
                .Define("valid_wires0", [&](ROOT::RVec<unsigned short> wires, ROOT::RVec<unsigned short> planes, 
                                           ROOT::RVec<unsigned short> tpcs, ROOT::RVec<bool> ontraj, ROOT::RVec<float> pitches) {
                    ROOT::RVec<unsigned short> valid;
                    for (size_t i = 0; i < wires.size() && i < planes.size() && i < tpcs.size() && i < ontraj.size() && i < pitches.size(); ++i) {
                        if (ontraj[i] && pitches[i] != -1.0f && 
                            dead_channels.find({wires[i], planes[i], tpcs[i]}) == dead_channels.end()) {
                            valid.push_back(wires[i]);
                        }
                    }
                    return valid;
                }, {"trk.hits0.h.wire", "trk.hits0.h.plane", "trk.hits0.h.tpc", "trk.hits0.ontraj", "trk.hits0.pitch"})
                .Define("valid_pitches0", [&](ROOT::RVec<unsigned short> wires, ROOT::RVec<unsigned short> planes, 
                                             ROOT::RVec<unsigned short> tpcs, ROOT::RVec<bool> ontraj, ROOT::RVec<float> pitches) {
                    ROOT::RVec<float> valid;
                    for (size_t i = 0; i < wires.size() && i < planes.size() && i < tpcs.size() && i < ontraj.size() && i < pitches.size(); ++i) {
                        if (ontraj[i] && pitches[i] != -1.0f && 
                            dead_channels.find({wires[i], planes[i], tpcs[i]}) == dead_channels.end()) {
                            valid.push_back(pitches[i]);
                        }
                    }
                    return valid;
                }, {"trk.hits0.h.wire", "trk.hits0.h.plane", "trk.hits0.h.tpc", "trk.hits0.ontraj", "trk.hits0.pitch"})
                .Define("valid_tpcs0", [&](ROOT::RVec<unsigned short> wires, ROOT::RVec<unsigned short> planes, 
                                          ROOT::RVec<unsigned short> tpcs, ROOT::RVec<bool> ontraj, ROOT::RVec<float> pitches) {
                    ROOT::RVec<unsigned short> valid;
                    for (size_t i = 0; i < wires.size() && i < planes.size() && i < tpcs.size() && i < ontraj.size() && i < pitches.size(); ++i) {
                        if (ontraj[i] && pitches[i] != -1.0f && 
                            dead_channels.find({wires[i], planes[i], tpcs[i]}) == dead_channels.end()) {
                            valid.push_back(tpcs[i]);
                        }
                    }
                    return valid;
                }, {"trk.hits0.h.wire", "trk.hits0.h.plane", "trk.hits0.h.tpc", "trk.hits0.ontraj", "trk.hits0.pitch"})
                .Define("valid_wires1", [&](ROOT::RVec<unsigned short> wires, ROOT::RVec<unsigned short> planes, 
                                           ROOT::RVec<unsigned short> tpcs, ROOT::RVec<bool> ontraj, ROOT::RVec<float> pitches) {
                    ROOT::RVec<unsigned short> valid;
                    for (size_t i = 0; i < wires.size() && i < planes.size() && i < tpcs.size() && i < ontraj.size() && i < pitches.size(); ++i) {
                        if (ontraj[i] && pitches[i] != -1.0f && 
                            dead_channels.find({wires[i], planes[i], tpcs[i]}) == dead_channels.end()) {
                            valid.push_back(wires[i]);
                        }
                    }
                    return valid;
                }, {"trk.hits1.h.wire", "trk.hits1.h.plane", "trk.hits1.h.tpc", "trk.hits1.ontraj", "trk.hits1.pitch"})
                .Define("valid_pitches1", [&](ROOT::RVec<unsigned short> wires, ROOT::RVec<unsigned short> planes, 
                                             ROOT::RVec<unsigned short> tpcs, ROOT::RVec<bool> ontraj, ROOT::RVec<float> pitches) {
                    ROOT::RVec<float> valid;
                    for (size_t i = 0; i < wires.size() && i < planes.size() && i < tpcs.size() && i < ontraj.size() && i < pitches.size(); ++i) {
                        if (ontraj[i] && pitches[i] != -1.0f && 
                            dead_channels.find({wires[i], planes[i], tpcs[i]}) == dead_channels.end()) {
                            valid.push_back(pitches[i]);
                        }
                    }
                    return valid;
                }, {"trk.hits1.h.wire", "trk.hits1.h.plane", "trk.hits1.h.tpc", "trk.hits1.ontraj", "trk.hits1.pitch"})
                .Define("valid_tpcs1", [&](ROOT::RVec<unsigned short> wires, ROOT::RVec<unsigned short> planes, 
                                          ROOT::RVec<unsigned short> tpcs, ROOT::RVec<bool> ontraj, ROOT::RVec<float> pitches) {
                    ROOT::RVec<unsigned short> valid;
                    for (size_t i = 0; i < wires.size() && i < planes.size() && i < tpcs.size() && i < ontraj.size() && i < pitches.size(); ++i) {
                        if (ontraj[i] && pitches[i] != -1.0f && 
                            dead_channels.find({wires[i], planes[i], tpcs[i]}) == dead_channels.end()) {
                            valid.push_back(tpcs[i]);
                        }
                    }
                    return valid;
                }, {"trk.hits1.h.wire", "trk.hits1.h.plane", "trk.hits1.h.tpc", "trk.hits1.ontraj", "trk.hits1.pitch"})
                .Define("valid_wires2", [&](ROOT::RVec<unsigned short> wires, ROOT::RVec<unsigned short> planes, 
                                           ROOT::RVec<unsigned short> tpcs, ROOT::RVec<bool> ontraj, ROOT::RVec<float> pitches) {
                    ROOT::RVec<unsigned short> valid;
                    for (size_t i = 0; i < wires.size() && i < planes.size() && i < tpcs.size() && i < ontraj.size() && i < pitches.size(); ++i) {
                        if (ontraj[i] && pitches[i] != -1.0f && 
                            dead_channels.find({wires[i], planes[i], tpcs[i]}) == dead_channels.end()) {
                            valid.push_back(wires[i]);
                        }
                    }
                    return valid;
                }, {"trk.hits2.h.wire", "trk.hits2.h.plane", "trk.hits2.h.tpc", "trk.hits2.ontraj", "trk.hits2.pitch"})
                .Define("valid_pitches2", [&](ROOT::RVec<unsigned short> wires, ROOT::RVec<unsigned short> planes, 
                                             ROOT::RVec<unsigned short> tpcs, ROOT::RVec<bool> ontraj, ROOT::RVec<float> pitches) {
                    ROOT::RVec<float> valid;
                    for (size_t i = 0; i < wires.size() && i < planes.size() && i < tpcs.size() && i < ontraj.size() && i < pitches.size(); ++i) {
                        if (ontraj[i] && pitches[i] != -1.0f && 
                            dead_channels.find({wires[i], planes[i], tpcs[i]}) == dead_channels.end()) {
                            valid.push_back(pitches[i]);
                        }
                    }
                    return valid;
                }, {"trk.hits2.h.wire", "trk.hits2.h.plane", "trk.hits2.h.tpc", "trk.hits2.ontraj", "trk.hits2.pitch"})
                .Define("valid_tpcs2", [&](ROOT::RVec<unsigned short> wires, ROOT::RVec<unsigned short> planes, 
                                          ROOT::RVec<unsigned short> tpcs, ROOT::RVec<bool> ontraj, ROOT::RVec<float> pitches) {
                    ROOT::RVec<unsigned short> valid;
                    for (size_t i = 0; i < wires.size() && i < planes.size() && i < tpcs.size() && i < ontraj.size() && i < pitches.size(); ++i) {
                        if (ontraj[i] && pitches[i] != -1.0f && 
                            dead_channels.find({wires[i], planes[i], tpcs[i]}) == dead_channels.end()) {
                            valid.push_back(tpcs[i]);
                        }
                    }
                    return valid;
                }, {"trk.hits2.h.wire", "trk.hits2.h.plane", "trk.hits2.h.tpc", "trk.hits2.ontraj", "trk.hits2.pitch"});

            // Process data for each plane
            rdf_cleaned.Foreach([&](int trk_id, float track_length, 
                                   ROOT::RVec<unsigned short> valid_wires0, ROOT::RVec<float> valid_pitches0, ROOT::RVec<unsigned short> valid_tpcs0,
                                   ROOT::RVec<unsigned short> valid_wires1, ROOT::RVec<float> valid_pitches1, ROOT::RVec<unsigned short> valid_tpcs1,
                                   ROOT::RVec<unsigned short> valid_wires2, ROOT::RVec<float> valid_pitches2, ROOT::RVec<unsigned short> valid_tpcs2) {
                calculate_efficiency(trk_id, track_length, valid_wires0, valid_pitches0, valid_tpcs0, 0);
                calculate_efficiency(trk_id, track_length, valid_wires1, valid_pitches1, valid_tpcs1, 1);
                calculate_efficiency(trk_id, track_length, valid_wires2, valid_pitches2, valid_tpcs2, 2);
            }, {"trk.id", "trk.length", 
                "valid_wires0", "valid_pitches0", "valid_tpcs0",
                "valid_wires1", "valid_pitches1", "valid_tpcs1",
                "valid_wires2", "valid_pitches2", "valid_tpcs2"});

            // Update sample count and print checkpoint
            total_samples++;
            if (total_samples % 10 == 0) {
                std::lock_guard<std::mutex> lock(cout_mutex);
                std::cout << "[" << dataset_type << "] Processed " << total_samples << " samples" << std::endl;
            }
        }

        // Print statistics
        {
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << "\n=== Processing Statistics for " << dataset_type << " ===" << std::endl;
            std::cout << "Total samples processed: " << total_samples << "\n";
            std::cout << "Total events recorded: " << total_events << "\n";
            std::cout << "Average pitch: min = " << (min_pitch == std::numeric_limits<float>::max() ? 0 : min_pitch) 
                      << ", max = " << (max_pitch == std::numeric_limits<float>::lowest() ? 0 : max_pitch) << "\n";
            std::cout << "Wires per track: min = " << (min_wires == std::numeric_limits<size_t>::max() ? 0 : min_wires) 
                      << ", max = " << max_wires << "\n";
            std::cout << "Valid hits per track: min = " << (min_hits == std::numeric_limits<size_t>::max() ? 0 : min_hits) 
                      << ", max = " << max_hits << "\n";
            std::cout << "Average efficiency: " << (efficiency_count > 0 ? total_efficiency / efficiency_count : 0.0) << "\n";
        }

        csv_out.close();
    };

    // Process Data FIRST, then MC in parallel using threads
    std::thread data_thread([&]() { process_dataset(data_filelist, data_output_csv, "Data"); });
    std::thread mc_thread([&]() { process_dataset(mc_filelist, mc_output_csv, "MC"); });

    // Wait for both threads to complete
    data_thread.join();
    mc_thread.join();

    std::cout << "\n=== Analysis Complete ===" << std::endl;
}