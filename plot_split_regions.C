#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <map>

void plot_split_regions() {
    // Create output directory
    gSystem->mkdir("plots_split_regions", kTRUE);
    
    // Define region information
    struct RegionInfo {
        std::string name;
        std::string display_name;
        std::string coord_range;
        int tpc_id;
    };
    
    std::vector<RegionInfo> regions = {
        {"TPC0_00", "TPC 0 Region (0,0)", "y:[-203.7, 0], z:[0, 250.1]", 0},
        {"TPC0_01", "TPC 0 Region (0,1)", "y:[0, 203.7], z:[0, 250.1]", 0},
        {"TPC0_10", "TPC 0 Region (1,0)", "y:[-203.7, 0], z:[250.1, 500.1]", 0},
        {"TPC0_11", "TPC 0 Region (1,1)", "y:[0, 203.7], z:[250.1, 500.1]", 0},
        {"TPC1_00", "TPC 1 Region (0,0)", "y:[-203.7, 0], z:[0, 250.1]", 1},
        {"TPC1_01", "TPC 1 Region (0,1)", "y:[0, 203.7], z:[0, 250.1]", 1},
        {"TPC1_10", "TPC 1 Region (1,0)", "y:[-203.7, 0], z:[250.1, 500.1]", 1},
        {"TPC1_11", "TPC 1 Region (1,1)", "y:[0, 203.7], z:[250.1, 500.1]", 1}
    };
    
    // Define anode/cathode regions
    struct SpecialRegion {
        std::string name;
        std::string display_name;
    };
    
    std::vector<SpecialRegion> special_regions = {
        {"anode_tpc0", "Anode TPC0 (x:[-202.2, -152.2])"},
        {"cathode", "Cathode (x:[-50, 50])"},
        {"anode_tpc1", "Anode TPC1 (x:[152.2, 202.2])"}
    };
    
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);
    
    // Function to create and save plot
    auto create_plot = [](const std::string& canvas_name, const std::string& title, 
                         const std::vector<float>& pitches_plane0_data, const std::vector<float>& effs_plane0_data,
                         const std::vector<float>& pitches_plane1_data, const std::vector<float>& effs_plane1_data,
                         const std::vector<float>& pitches_plane2_data, const std::vector<float>& effs_plane2_data,
                         const std::vector<float>& pitches_plane0_mc, const std::vector<float>& effs_plane0_mc,
                         const std::vector<float>& pitches_plane1_mc, const std::vector<float>& effs_plane1_mc,
                         const std::vector<float>& pitches_plane2_mc, const std::vector<float>& effs_plane2_mc,
                         int total_entries_data, int total_entries_mc,
                         float min_pitch_data, float max_pitch_data, float mean_pitch_data, float mean_eff_data,
                         float min_pitch_mc, float max_pitch_mc, float mean_pitch_mc, float mean_eff_mc,
                         const std::string& output_filename) {
        
        // Create TProfile histograms - Data
        TProfile* h_eff0_data = new TProfile(("h_eff0_data_" + canvas_name).c_str(), "Hit Efficiency vs Pitch;Average Pitch [cm];Efficiency", 100, 0.295, 0.8, 0, 1);
        TProfile* h_eff1_data = new TProfile(("h_eff1_data_" + canvas_name).c_str(), "Hit Efficiency vs Pitch;Average Pitch [cm];Efficiency", 100, 0.295, 0.8, 0, 1);
        TProfile* h_eff2_data = new TProfile(("h_eff2_data_" + canvas_name).c_str(), "Hit Efficiency vs Pitch;Average Pitch [cm];Efficiency", 100, 0.295, 0.8, 0, 1);
        
        // Fill histograms - Data
        for (size_t i = 0; i < pitches_plane0_data.size(); ++i)
            h_eff0_data->Fill(pitches_plane0_data[i], effs_plane0_data[i]);
        for (size_t i = 0; i < pitches_plane1_data.size(); ++i)
            h_eff1_data->Fill(pitches_plane1_data[i], effs_plane1_data[i]);
        for (size_t i = 0; i < pitches_plane2_data.size(); ++i)
            h_eff2_data->Fill(pitches_plane2_data[i], effs_plane2_data[i]);
        
        // Create TProfile histograms - MC
        TProfile* h_eff0_mc = new TProfile(("h_eff0_mc_" + canvas_name).c_str(), "Hit Efficiency vs Pitch;Average Pitch [cm];Efficiency", 100, 0.295, 0.8, 0, 1);
        TProfile* h_eff1_mc = new TProfile(("h_eff1_mc_" + canvas_name).c_str(), "Hit Efficiency vs Pitch;Average Pitch [cm];Efficiency", 100, 0.295, 0.8, 0, 1);
        TProfile* h_eff2_mc = new TProfile(("h_eff2_mc_" + canvas_name).c_str(), "Hit Efficiency vs Pitch;Average Pitch [cm];Efficiency", 100, 0.295, 0.8, 0, 1);
        
        // Fill histograms - MC
        for (size_t i = 0; i < pitches_plane0_mc.size(); ++i)
            h_eff0_mc->Fill(pitches_plane0_mc[i], effs_plane0_mc[i]);
        for (size_t i = 0; i < pitches_plane1_mc.size(); ++i)
            h_eff1_mc->Fill(pitches_plane1_mc[i], effs_plane1_mc[i]);
        for (size_t i = 0; i < pitches_plane2_mc.size(); ++i)
            h_eff2_mc->Fill(pitches_plane2_mc[i], effs_plane2_mc[i]);
        
        // Create canvas
        TCanvas* c = new TCanvas(canvas_name.c_str(), title.c_str(), 1000, 600);
        c->SetBatch(kTRUE);
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.12);
        c->SetGrid(1, 1);
        
        // Set colors and styles - Data
        h_eff0_data->SetLineColor(kBlue);
        h_eff0_data->SetMarkerColor(kBlue);
        h_eff0_data->SetMarkerStyle(20);
        h_eff0_data->SetMarkerSize(0.5);
        h_eff0_data->GetXaxis()->SetRangeUser(0.295, 0.8);
        h_eff0_data->GetYaxis()->SetRangeUser(0.95, 1.002);
        h_eff0_data->GetXaxis()->SetTitleSize(0.05);
        h_eff0_data->GetYaxis()->SetTitleSize(0.05);
        
        // Set title
        h_eff0_data->SetTitle((title + ";Average Pitch [cm];Efficiency").c_str());
        
        // Draw first histogram
        h_eff0_data->Draw("P");
        
        h_eff1_data->SetLineColor(kRed);
        h_eff1_data->SetMarkerColor(kRed);
        h_eff1_data->SetMarkerStyle(21);
        h_eff1_data->SetMarkerSize(0.5);
        if (pitches_plane1_data.size() > 0) h_eff1_data->Draw("P SAME");
        
        h_eff2_data->SetLineColor(kGreen+2);
        h_eff2_data->SetMarkerColor(kGreen+2);
        h_eff2_data->SetMarkerStyle(22);
        h_eff2_data->SetMarkerSize(0.5);
        if (pitches_plane2_data.size() > 0) h_eff2_data->Draw("P SAME");
        
        // Set colors and styles - MC
        h_eff0_mc->SetLineColor(kCyan);
        h_eff0_mc->SetMarkerColor(kCyan);
        h_eff0_mc->SetMarkerStyle(24);
        h_eff0_mc->SetMarkerSize(0.5);
        if (pitches_plane0_mc.size() > 0) h_eff0_mc->Draw("P SAME");
        
        h_eff1_mc->SetLineColor(kMagenta);
        h_eff1_mc->SetMarkerColor(kMagenta);
        h_eff1_mc->SetMarkerStyle(25);
        h_eff1_mc->SetMarkerSize(0.5);
        if (pitches_plane1_mc.size() > 0) h_eff1_mc->Draw("P SAME");
        
        h_eff2_mc->SetLineColor(kYellow+2);
        h_eff2_mc->SetMarkerColor(kYellow+2);
        h_eff2_mc->SetMarkerStyle(26);
        h_eff2_mc->SetMarkerSize(0.5);
        if (pitches_plane2_mc.size() > 0) h_eff2_mc->Draw("P SAME");
        
        // Create legend
        TLegend* leg = new TLegend(0.55, 0.2, 0.75, 0.5);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        if (pitches_plane0_data.size() > 0) leg->AddEntry(h_eff0_data, Form("Plane 0 Data (%zu tracks)", pitches_plane0_data.size()), "lp");
        if (pitches_plane1_data.size() > 0) leg->AddEntry(h_eff1_data, Form("Plane 1 Data (%zu tracks)", pitches_plane1_data.size()), "lp");
        if (pitches_plane2_data.size() > 0) leg->AddEntry(h_eff2_data, Form("Plane 2 Data (%zu tracks)", pitches_plane2_data.size()), "lp");
        leg->AddEntry((TObject*)0, "", "");
        if (pitches_plane0_mc.size() > 0) leg->AddEntry(h_eff0_mc, Form("Plane 0 MC (%zu tracks)", pitches_plane0_mc.size()), "lp");
        if (pitches_plane1_mc.size() > 0) leg->AddEntry(h_eff1_mc, Form("Plane 1 MC (%zu tracks)", pitches_plane1_mc.size()), "lp");
        if (pitches_plane2_mc.size() > 0) leg->AddEntry(h_eff2_mc, Form("Plane 2 MC (%zu tracks)", pitches_plane2_mc.size()), "lp");
        leg->Draw();
        
        // Add statistics text
        TLatex* latex = new TLatex();
        latex->SetTextSize(0.025);
        latex->DrawLatex(0.6, 0.77, "Data Statistics:");
        latex->DrawLatex(0.6, 0.74, Form("Total Tracks: %d", total_entries_data));
        latex->DrawLatex(0.6, 0.71, Form("Pitches - Min: %.6f, Max: %.5f", min_pitch_data, max_pitch_data));
        latex->DrawLatex(0.6, 0.68, Form("Mean: %.6f", mean_pitch_data));
        latex->DrawLatex(0.6, 0.65, Form("Efficiency - Mean: %.5f", mean_eff_data));
        
        latex->DrawLatex(0.6, 0.59, "MC Statistics:");
        latex->DrawLatex(0.6, 0.56, Form("Total Tracks: %d", total_entries_mc));
        latex->DrawLatex(0.6, 0.53, Form("Pitches - Min: %.6f, Max: %.5f", min_pitch_mc, max_pitch_mc));
        latex->DrawLatex(0.6, 0.50, Form("Mean: %.6f", mean_pitch_mc));
        latex->DrawLatex(0.6, 0.47, Form("Efficiency - Mean: %.5f", mean_eff_mc));
        
        // Save plot
        c->SaveAs(output_filename.c_str());
        
        std::cout << "Saved plot: " << output_filename << " (Total tracks Data: " << total_entries_data << ", MC: " << total_entries_mc << ")" << std::endl;
        
        // Clean up
        delete h_eff0_data;
        delete h_eff1_data;
        delete h_eff2_data;
        delete h_eff0_mc;
        delete h_eff1_mc;
        delete h_eff2_mc;
        delete c;
        delete leg;
        delete latex;
    };
    
    // Process each region
    for (const auto& region : regions) {
        std::string filename_data = "split_regions/" + region.name + "_hits_data.csv";
        std::string filename_mc = "split_regions/" + region.name + "_hits_mc.csv";
        
        // Load CSV data for this region - Data
        std::vector<float> pitches_plane0_data, pitches_plane1_data, pitches_plane2_data;
        std::vector<float> effs_plane0_data, effs_plane1_data, effs_plane2_data;
        std::vector<float> all_pitches_data, all_effs_data;
        
        std::ifstream csv_file_data(filename_data);
        if (!csv_file_data.is_open()) {
            std::cout << "Warning: Cannot open " << filename_data << ", skipping..." << std::endl;
            continue;
        }
        
        std::string header;
        std::getline(csv_file_data, header); // Skip header
        std::string line;
        int total_entries_data = 0;
        
        while (std::getline(csv_file_data, line)) {
            std::stringstream ss(line);
            std::string trk_id_str, plane_str, tpc_str, length_str, hits_str, 
                       min_x_str, max_x_str, min_y_str, max_y_str, min_z_str, max_z_str, 
                       pitch_str, eff_str, anode_tpc0_hits_str, cathode_hits_str, anode_tpc1_hits_str, other_hits_str;
            
            std::getline(ss, trk_id_str, ',');
            std::getline(ss, plane_str, ',');
            std::getline(ss, tpc_str, ',');
            std::getline(ss, length_str, ',');
            std::getline(ss, hits_str, ',');
            std::getline(ss, min_x_str, ',');
            std::getline(ss, max_x_str, ',');
            std::getline(ss, min_y_str, ',');
            std::getline(ss, max_y_str, ',');
            std::getline(ss, min_z_str, ',');
            std::getline(ss, max_z_str, ',');
            std::getline(ss, pitch_str, ',');
            std::getline(ss, eff_str, ',');
            std::getline(ss, anode_tpc0_hits_str, ',');
            std::getline(ss, cathode_hits_str, ',');
            std::getline(ss, anode_tpc1_hits_str, ',');
            std::getline(ss, other_hits_str, ',');
            
            try {
                int plane = std::stoi(plane_str);
                float avg_pitch = std::stof(pitch_str);
                float efficiency = std::stof(eff_str);
                
                if (plane == 0) {
                    pitches_plane0_data.push_back(avg_pitch);
                    effs_plane0_data.push_back(efficiency);
                } else if (plane == 1) {
                    pitches_plane1_data.push_back(avg_pitch);
                    effs_plane1_data.push_back(efficiency);
                } else if (plane == 2) {
                    pitches_plane2_data.push_back(avg_pitch);
                    effs_plane2_data.push_back(efficiency);
                }
                
                all_pitches_data.push_back(avg_pitch);
                all_effs_data.push_back(efficiency);
                total_entries_data++;
                
            } catch (...) {
                std::cout << "Warning: Invalid line in " << filename_data << ": " << line << std::endl;
            }
        }
        csv_file_data.close();
        
        // Load CSV data for this region - MC
        std::vector<float> pitches_plane0_mc, pitches_plane1_mc, pitches_plane2_mc;
        std::vector<float> effs_plane0_mc, effs_plane1_mc, effs_plane2_mc;
        std::vector<float> all_pitches_mc, all_effs_mc;
        
        std::ifstream csv_file_mc(filename_mc);
        if (!csv_file_mc.is_open()) {
            std::cout << "Warning: Cannot open " << filename_mc << ", skipping..." << std::endl;
            continue;
        }
        
        std::getline(csv_file_mc, header); // Skip header
        
        int total_entries_mc = 0;
        
        while (std::getline(csv_file_mc, line)) {
            std::stringstream ss(line);
            std::string trk_id_str, plane_str, tpc_str, length_str, hits_str, 
                       min_x_str, max_x_str, min_y_str, max_y_str, min_z_str, max_z_str, 
                       pitch_str, eff_str, anode_tpc0_hits_str, cathode_hits_str, anode_tpc1_hits_str, other_hits_str;
            
            std::getline(ss, trk_id_str, ',');
            std::getline(ss, plane_str, ',');
            std::getline(ss, tpc_str, ',');
            std::getline(ss, length_str, ',');
            std::getline(ss, hits_str, ',');
            std::getline(ss, min_x_str, ',');
            std::getline(ss, max_x_str, ',');
            std::getline(ss, min_y_str, ',');
            std::getline(ss, max_y_str, ',');
            std::getline(ss, min_z_str, ',');
            std::getline(ss, max_z_str, ',');
            std::getline(ss, pitch_str, ',');
            std::getline(ss, eff_str, ',');
            std::getline(ss, anode_tpc0_hits_str, ',');
            std::getline(ss, cathode_hits_str, ',');
            std::getline(ss, anode_tpc1_hits_str, ',');
            std::getline(ss, other_hits_str, ',');
            
            try {
                int plane = std::stoi(plane_str);
                float avg_pitch = std::stof(pitch_str);
                float efficiency = std::stof(eff_str);
                
                if (plane == 0) {
                    pitches_plane0_mc.push_back(avg_pitch);
                    effs_plane0_mc.push_back(efficiency);
                } else if (plane == 1) {
                    pitches_plane1_mc.push_back(avg_pitch);
                    effs_plane1_mc.push_back(efficiency);
                } else if (plane == 2) {
                    pitches_plane2_mc.push_back(avg_pitch);
                    effs_plane2_mc.push_back(efficiency);
                }
                
                all_pitches_mc.push_back(avg_pitch);
                all_effs_mc.push_back(efficiency);
                total_entries_mc++;
                
            } catch (...) {
                std::cout << "Warning: Invalid line in " << filename_mc << ": " << line << std::endl;
            }
        }
        csv_file_mc.close();
        
        if (total_entries_data == 0 && total_entries_mc == 0) {
            std::cout << "No valid data found for region " << region.name << ", skipping..." << std::endl;
            continue;
        }
        
        // Calculate statistics - Data
        float min_pitch_data = all_pitches_data.empty() ? 0.0 : *std::min_element(all_pitches_data.begin(), all_pitches_data.end());
        float max_pitch_data = all_pitches_data.empty() ? 0.0 : *std::max_element(all_pitches_data.begin(), all_pitches_data.end());
        float mean_pitch_data = all_pitches_data.empty() ? 0.0 : std::accumulate(all_pitches_data.begin(), all_pitches_data.end(), 0.0) / all_pitches_data.size();
        float mean_eff_data = all_effs_data.empty() ? 0.0 : std::accumulate(all_effs_data.begin(), all_effs_data.end(), 0.0) / all_effs_data.size();
        
        // Calculate statistics - MC
        float min_pitch_mc = all_pitches_mc.empty() ? 0.0 : *std::min_element(all_pitches_mc.begin(), all_pitches_mc.end());
        float max_pitch_mc = all_pitches_mc.empty() ? 0.0 : *std::max_element(all_pitches_mc.begin(), all_pitches_mc.end());
        float mean_pitch_mc = all_pitches_mc.empty() ? 0.0 : std::accumulate(all_pitches_mc.begin(), all_pitches_mc.end(), 0.0) / all_pitches_mc.size();
        float mean_eff_mc = all_effs_mc.empty() ? 0.0 : std::accumulate(all_effs_mc.begin(), all_effs_mc.end(), 0.0) / all_effs_mc.size();
        
        // Create and save region plot
        create_plot("c_" + region.name, 
                    region.display_name + " - " + region.coord_range,
                    pitches_plane0_data, effs_plane0_data,
                    pitches_plane1_data, effs_plane1_data,
                    pitches_plane2_data, effs_plane2_data,
                    pitches_plane0_mc, effs_plane0_mc,
                    pitches_plane1_mc, effs_plane1_mc,
                    pitches_plane2_mc, effs_plane2_mc,
                    total_entries_data, total_entries_mc,
                    min_pitch_data, max_pitch_data, mean_pitch_data, mean_eff_data,
                    min_pitch_mc, max_pitch_mc, mean_pitch_mc, mean_eff_mc,
                    "plots_split_regions/hit_efficiency_" + region.name + ".png");
    }
    
    // Process special regions (anode_tpc0, cathode, anode_tpc1)
    for (const auto& special_region : special_regions) {
        // Initialize data vectors for special regions
        std::vector<float> pitches_plane0_data, pitches_plane1_data, pitches_plane2_data;
        std::vector<float> effs_plane0_data, effs_plane1_data, effs_plane2_data;
        std::vector<float> all_pitches_data, all_effs_data;
        std::vector<float> pitches_plane0_mc, pitches_plane1_mc, pitches_plane2_mc;
        std::vector<float> effs_plane0_mc, effs_plane1_mc, effs_plane2_mc;
        std::vector<float> all_pitches_mc, all_effs_mc;
        
        int total_entries_data = 0;
        int total_entries_mc = 0;
        
        // Process all region files to collect data for special regions
        for (const auto& region : regions) {
            std::string filename_data = "split_regions/" + region.name + "_hits_data.csv";
            std::string filename_mc = "split_regions/" + region.name + "_hits_mc.csv";
            
            // Load Data
            std::ifstream csv_file_data(filename_data);
            if (csv_file_data.is_open()) {
                std::string header;
                std::getline(csv_file_data, header); // Skip header
                std::string line;
                
                while (std::getline(csv_file_data, line)) {
                    std::stringstream ss(line);
                    std::string trk_id_str, plane_str, tpc_str, length_str, hits_str, 
                               min_x_str, max_x_str, min_y_str, max_y_str, min_z_str, max_z_str, 
                               pitch_str, eff_str, anode_tpc0_hits_str, cathode_hits_str, anode_tpc1_hits_str, other_hits_str;
                    
                    std::getline(ss, trk_id_str, ',');
                    std::getline(ss, plane_str, ',');
                    std::getline(ss, tpc_str, ',');
                    std::getline(ss, length_str, ',');
                    std::getline(ss, hits_str, ',');
                    std::getline(ss, min_x_str, ',');
                    std::getline(ss, max_x_str, ',');
                    std::getline(ss, min_y_str, ',');
                    std::getline(ss, max_y_str, ',');
                    std::getline(ss, min_z_str, ',');
                    std::getline(ss, max_z_str, ',');
                    std::getline(ss, pitch_str, ',');
                    std::getline(ss, eff_str, ',');
                    std::getline(ss, anode_tpc0_hits_str, ',');
                    std::getline(ss, cathode_hits_str, ',');
                    std::getline(ss, anode_tpc1_hits_str, ',');
                    std::getline(ss, other_hits_str, ',');
                    
                    try {
                        int plane = std::stoi(plane_str);
                        float avg_pitch = std::stof(pitch_str);
                        float efficiency = std::stof(eff_str);
                        int anode_tpc0_hits = std::stoi(anode_tpc0_hits_str);
                        int cathode_hits = std::stoi(cathode_hits_str);
                        int anode_tpc1_hits = std::stoi(anode_tpc1_hits_str);
                        
                        bool include = false;
                        if (special_region.name == "anode_tpc0" && anode_tpc0_hits > 0) include = true;
                        else if (special_region.name == "cathode" && cathode_hits > 0) include = true;
                        else if (special_region.name == "anode_tpc1" && anode_tpc1_hits > 0) include = true;
                        
                        if (include) {
                            if (plane == 0) {
                                pitches_plane0_data.push_back(avg_pitch);
                                effs_plane0_data.push_back(efficiency);
                            } else if (plane == 1) {
                                pitches_plane1_data.push_back(avg_pitch);
                                effs_plane1_data.push_back(efficiency);
                            } else if (plane == 2) {
                                pitches_plane2_data.push_back(avg_pitch);
                                effs_plane2_data.push_back(efficiency);
                            }
                            
                            all_pitches_data.push_back(avg_pitch);
                            all_effs_data.push_back(efficiency);
                            total_entries_data++;
                        }
                    } catch (...) {
                        std::cout << "Warning: Invalid line in " << filename_data << ": " << line << std::endl;
                    }
                }
                csv_file_data.close();
            }
            
            // Load MC
            std::ifstream csv_file_mc(filename_mc);
            if (csv_file_mc.is_open()) {
                std::string header;
                std::getline(csv_file_mc, header); // Skip header
                std::string line;
                
                while (std::getline(csv_file_mc, line)) {
                    std::stringstream ss(line);
                    std::string trk_id_str, plane_str, tpc_str, length_str, hits_str, 
                               min_x_str, max_x_str, min_y_str, max_y_str, min_z_str, max_z_str, 
                               pitch_str, eff_str, anode_tpc0_hits_str, cathode_hits_str, anode_tpc1_hits_str, other_hits_str;
                    
                    std::getline(ss, trk_id_str, ',');
                    std::getline(ss, plane_str, ',');
                    std::getline(ss, tpc_str, ',');
                    std::getline(ss, length_str, ',');
                    std::getline(ss, hits_str, ',');
                    std::getline(ss, min_x_str, ',');
                    std::getline(ss, max_x_str, ',');
                    std::getline(ss, min_y_str, ',');
                    std::getline(ss, max_y_str, ',');
                    std::getline(ss, min_z_str, ',');
                    std::getline(ss, max_z_str, ',');
                    std::getline(ss, pitch_str, ',');
                    std::getline(ss, eff_str, ',');
                    std::getline(ss, anode_tpc0_hits_str, ',');
                    std::getline(ss, cathode_hits_str, ',');
                    std::getline(ss, anode_tpc1_hits_str, ',');
                    std::getline(ss, other_hits_str, ',');
                    
                    try {
                        int plane = std::stoi(plane_str);
                        float avg_pitch = std::stof(pitch_str);
                        float efficiency = std::stof(eff_str);
                        int anode_tpc0_hits = std::stoi(anode_tpc0_hits_str);
                        int cathode_hits = std::stoi(cathode_hits_str);
                        int anode_tpc1_hits = std::stoi(anode_tpc1_hits_str);
                        
                        bool include = false;
                        if (special_region.name == "anode_tpc0" && anode_tpc0_hits > 0) include = true;
                        else if (special_region.name == "cathode" && cathode_hits > 0) include = true;
                        else if (special_region.name == "anode_tpc1" && anode_tpc1_hits > 0) include = true;
                        
                        if (include) {
                            if (plane == 0) {
                                pitches_plane0_mc.push_back(avg_pitch);
                                effs_plane0_mc.push_back(efficiency);
                            } else if (plane == 1) {
                                pitches_plane1_mc.push_back(avg_pitch);
                                effs_plane1_mc.push_back(efficiency);
                            } else if (plane == 2) {
                                pitches_plane2_mc.push_back(avg_pitch);
                                effs_plane2_mc.push_back(efficiency);
                            }
                            
                            all_pitches_mc.push_back(avg_pitch);
                            all_effs_mc.push_back(efficiency);
                            total_entries_mc++;
                        }
                    } catch (...) {
                        std::cout << "Warning: Invalid line in " << filename_mc << ": " << line << std::endl;
                    }
                }
                csv_file_mc.close();
            }
        }
        
        if (total_entries_data == 0 && total_entries_mc == 0) {
            std::cout << "No valid data found for special region " << special_region.name << ", skipping..." << std::endl;
            continue;
        }
        
        // Calculate statistics - Data
        float min_pitch_data = all_pitches_data.empty() ? 0.0 : *std::min_element(all_pitches_data.begin(), all_pitches_data.end());
        float max_pitch_data = all_pitches_data.empty() ? 0.0 : *std::max_element(all_pitches_data.begin(), all_pitches_data.end());
        float mean_pitch_data = all_pitches_data.empty() ? 0.0 : std::accumulate(all_pitches_data.begin(), all_pitches_data.end(), 0.0) / all_pitches_data.size();
        float mean_eff_data = all_effs_data.empty() ? 0.0 : std::accumulate(all_effs_data.begin(), all_effs_data.end(), 0.0) / all_effs_data.size();
        
        // Calculate statistics - MC
        float min_pitch_mc = all_pitches_mc.empty() ? 0.0 : *std::min_element(all_pitches_mc.begin(), all_pitches_mc.end());
        float max_pitch_mc = all_pitches_mc.empty() ? 0.0 : *std::max_element(all_pitches_mc.begin(), all_pitches_mc.end());
        float mean_pitch_mc = all_pitches_mc.empty() ? 0.0 : std::accumulate(all_pitches_mc.begin(), all_pitches_mc.end(), 0.0) / all_pitches_mc.size();
        float mean_eff_mc = all_effs_mc.empty() ? 0.0 : std::accumulate(all_effs_mc.begin(), all_effs_mc.end(), 0.0) / all_effs_mc.size();
        
        // Create and save special region plot
        create_plot("c_" + special_region.name, 
                    special_region.display_name,
                    pitches_plane0_data, effs_plane0_data,
                    pitches_plane1_data, effs_plane1_data,
                    pitches_plane2_data, effs_plane2_data,
                    pitches_plane0_mc, effs_plane0_mc,
                    pitches_plane1_mc, effs_plane1_mc,
                    pitches_plane2_mc, effs_plane2_mc,
                    total_entries_data, total_entries_mc,
                    min_pitch_data, max_pitch_data, mean_pitch_data, mean_eff_data,
                    min_pitch_mc, max_pitch_mc, mean_pitch_mc, mean_eff_mc,
                    "plots_split_regions/hit_efficiency_" + special_region.name + ".png");
    }
    
    std::cout << "\nAll region plots saved in plots_split_regions/ directory" << std::endl;
}