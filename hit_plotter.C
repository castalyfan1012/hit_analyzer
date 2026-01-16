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

void hit_plotter() {
    // ========== CONFIGURATION PARAMETERS ==========
    // Input CSV files
    const char* data_csv = "hiteff_data.csv";
    const char* mc_csv = "hiteff_mc.csv";
    
    // Output directory
    const char* output_dir = "plots_hiteff";
    
    // Plot axis ranges
    const float pitch_xmin = 0.28;
    const float pitch_xmax = 0.8;
    const float eff_ymin = 0.95;
    const float eff_ymax = 1.002;
    
    // Binned plot ranges
    const int nbins = 30;
    const float binned_xmin = 0.3;
    const float binned_xmax = 2.5;
    const float binned_ymin = 0.96;
    const float binned_ymax = 1.0;
    // ==============================================
    
    // Create output directory
    gSystem->mkdir(output_dir, kTRUE);
    
    // Load Data CSV
    std::vector<float> pitches_plane0_data, pitches_plane1_data, pitches_plane2_data;
    std::vector<float> effs_plane0_data, effs_plane1_data, effs_plane2_data;
    std::ifstream csv_data(data_csv);
    if (!csv_data.is_open()) {
        std::cout << "Error: Cannot open " << data_csv << std::endl;
        return;
    }
    
    std::string header;
    std::getline(csv_data, header);
    std::string line;
    while (std::getline(csv_data, line)) {
        std::stringstream ss(line);
        std::string trk_id_str, plane_str, tpc_str, length_str, hits_str, wires_str, eff_str, pitch_str;
        std::getline(ss, trk_id_str, ',');
        std::getline(ss, plane_str, ',');
        std::getline(ss, tpc_str, ',');
        std::getline(ss, length_str, ',');
        std::getline(ss, hits_str, ',');
        std::getline(ss, wires_str, ',');
        std::getline(ss, eff_str, ',');
        std::getline(ss, pitch_str, ',');
        
        try {
            int plane = std::stoi(plane_str);
            float efficiency = std::stof(eff_str);
            float avg_pitch = std::stof(pitch_str);
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
        } catch (...) {
            std::cout << "Warning: Invalid line in " << data_csv << ": " << line << std::endl;
        }
    }
    csv_data.close();
    
    // Load MC CSV
    std::vector<float> pitches_plane0_mc, pitches_plane1_mc, pitches_plane2_mc;
    std::vector<float> effs_plane0_mc, effs_plane1_mc, effs_plane2_mc;
    std::ifstream csv_mc(mc_csv);
    if (!csv_mc.is_open()) {
        std::cout << "Error: Cannot open " << mc_csv << std::endl;
        return;
    }
    
    std::getline(csv_mc, header);
    while (std::getline(csv_mc, line)) {
        std::stringstream ss(line);
        std::string trk_id_str, plane_str, tpc_str, length_str, hits_str, wires_str, eff_str, pitch_str;
        std::getline(ss, trk_id_str, ',');
        std::getline(ss, plane_str, ',');
        std::getline(ss, tpc_str, ',');
        std::getline(ss, length_str, ',');
        std::getline(ss, hits_str, ',');
        std::getline(ss, wires_str, ',');
        std::getline(ss, eff_str, ',');
        std::getline(ss, pitch_str, ',');
        
        try {
            int plane = std::stoi(plane_str);
            float efficiency = std::stof(eff_str);
            float avg_pitch = std::stof(pitch_str);
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
        } catch (...) {
            std::cout << "Warning: Invalid line in " << mc_csv << ": " << line << std::endl;
        }
    }
    csv_mc.close();
    
    // Combine all planes for Data
    std::vector<float> all_pitches_data, all_effs_data;
    all_pitches_data.insert(all_pitches_data.end(), pitches_plane0_data.begin(), pitches_plane0_data.end());
    all_pitches_data.insert(all_pitches_data.end(), pitches_plane1_data.begin(), pitches_plane1_data.end());
    all_pitches_data.insert(all_pitches_data.end(), pitches_plane2_data.begin(), pitches_plane2_data.end());
    all_effs_data.insert(all_effs_data.end(), effs_plane0_data.begin(), effs_plane0_data.end());
    all_effs_data.insert(all_effs_data.end(), effs_plane1_data.begin(), effs_plane1_data.end());
    all_effs_data.insert(all_effs_data.end(), effs_plane2_data.begin(), effs_plane2_data.end());
    
    // Combine all planes for MC
    std::vector<float> all_pitches_mc, all_effs_mc;
    all_pitches_mc.insert(all_pitches_mc.end(), pitches_plane0_mc.begin(), pitches_plane0_mc.end());
    all_pitches_mc.insert(all_pitches_mc.end(), pitches_plane1_mc.begin(), pitches_plane1_mc.end());
    all_pitches_mc.insert(all_pitches_mc.end(), pitches_plane2_mc.begin(), pitches_plane2_mc.end());
    all_effs_mc.insert(all_effs_mc.end(), effs_plane0_mc.begin(), effs_plane0_mc.end());
    all_effs_mc.insert(all_effs_mc.end(), effs_plane1_mc.begin(), effs_plane1_mc.end());
    all_effs_mc.insert(all_effs_mc.end(), effs_plane2_mc.begin(), effs_plane2_mc.end());
    
    // Calculate statistics
    int total_tracks_data = all_pitches_data.size();
    float min_pitch_data = total_tracks_data > 0 ? *std::min_element(all_pitches_data.begin(), all_pitches_data.end()) : 0.0;
    float max_pitch_data = total_tracks_data > 0 ? *std::max_element(all_pitches_data.begin(), all_pitches_data.end()) : 0.0;
    float mean_pitch_data = total_tracks_data > 0 ? std::accumulate(all_pitches_data.begin(), all_pitches_data.end(), 0.0) / total_tracks_data : 0.0;
    float mean_eff_data = total_tracks_data > 0 ? std::accumulate(all_effs_data.begin(), all_effs_data.end(), 0.0) / total_tracks_data : 0.0;
    
    int total_tracks_mc = all_pitches_mc.size();
    float min_pitch_mc = total_tracks_mc > 0 ? *std::min_element(all_pitches_mc.begin(), all_pitches_mc.end()) : 0.0;
    float max_pitch_mc = total_tracks_mc > 0 ? *std::max_element(all_pitches_mc.begin(), all_pitches_mc.end()) : 0.0;
    float mean_pitch_mc = total_tracks_mc > 0 ? std::accumulate(all_pitches_mc.begin(), all_pitches_mc.end(), 0.0) / total_tracks_mc : 0.0;
    float mean_eff_mc = total_tracks_mc > 0 ? std::accumulate(all_effs_mc.begin(), all_effs_mc.end(), 0.0) / total_tracks_mc : 0.0;
    
    // Create TProfile histograms for Data
    TProfile* h_eff0_data = new TProfile("h_eff0_data", "SBND TPC Hit Efficiency;Average Pitch [cm];Efficiency", 200, pitch_xmin, pitch_xmax, 0, 1);
    TProfile* h_eff1_data = new TProfile("h_eff1_data", "SBND TPC Hit Efficiency;Average Pitch [cm];Efficiency", 200, pitch_xmin, pitch_xmax, 0, 1);
    TProfile* h_eff2_data = new TProfile("h_eff2_data", "SBND TPC Hit Efficiency;Average Pitch [cm];Efficiency", 200, pitch_xmin, pitch_xmax, 0, 1);
    
    for (size_t i = 0; i < pitches_plane0_data.size(); ++i)
        h_eff0_data->Fill(pitches_plane0_data[i], effs_plane0_data[i]);
    for (size_t i = 0; i < pitches_plane1_data.size(); ++i)
        h_eff1_data->Fill(pitches_plane1_data[i], effs_plane1_data[i]);
    for (size_t i = 0; i < pitches_plane2_data.size(); ++i)
        h_eff2_data->Fill(pitches_plane2_data[i], effs_plane2_data[i]);
    
    // Create TProfile histograms for MC
    TProfile* h_eff0_mc = new TProfile("h_eff0_mc", "SBND TPC Hit Efficiency;Average Pitch [cm];Efficiency", 200, pitch_xmin, pitch_xmax, 0, 1);
    TProfile* h_eff1_mc = new TProfile("h_eff1_mc", "SBND TPC Hit Efficiency;Average Pitch [cm];Efficiency", 200, pitch_xmin, pitch_xmax, 0, 1);
    TProfile* h_eff2_mc = new TProfile("h_eff2_mc", "SBND TPC Hit Efficiency;Average Pitch [cm];Efficiency", 200, pitch_xmin, pitch_xmax, 0, 1);
    
    for (size_t i = 0; i < pitches_plane0_mc.size(); ++i)
        h_eff0_mc->Fill(pitches_plane0_mc[i], effs_plane0_mc[i]);
    for (size_t i = 0; i < pitches_plane1_mc.size(); ++i)
        h_eff1_mc->Fill(pitches_plane1_mc[i], effs_plane1_mc[i]);
    for (size_t i = 0; i < pitches_plane2_mc.size(); ++i)
        h_eff2_mc->Fill(pitches_plane2_mc[i], effs_plane2_mc[i]);
    
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);
    
    // ========== PLOT 1: Hit Efficiency vs Pitch (All Planes) ==========
    TCanvas* c1 = new TCanvas("c1", "SBND TPC Hit Efficiency", 1000, 600);
    c1->SetBatch(kTRUE);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    c1->SetGrid(1, 1);
    
    h_eff0_data->SetLineColor(kBlue);
    h_eff0_data->SetMarkerColor(kBlue);
    h_eff0_data->SetMarkerStyle(20);
    h_eff0_data->SetMarkerSize(0.5);
    h_eff0_data->GetXaxis()->SetRangeUser(pitch_xmin, pitch_xmax);
    h_eff0_data->GetYaxis()->SetRangeUser(eff_ymin, eff_ymax);
    h_eff0_data->GetXaxis()->SetTitleSize(0.05);
    h_eff0_data->GetYaxis()->SetTitleSize(0.05);
    h_eff0_data->Draw("P");
    
    h_eff1_data->SetLineColor(kRed);
    h_eff1_data->SetMarkerColor(kRed);
    h_eff1_data->SetMarkerStyle(21);
    h_eff1_data->SetMarkerSize(0.5);
    h_eff1_data->Draw("P SAME");
    
    h_eff2_data->SetLineColor(kGreen+2);
    h_eff2_data->SetMarkerColor(kGreen+2);
    h_eff2_data->SetMarkerStyle(22);
    h_eff2_data->SetMarkerSize(0.5);
    h_eff2_data->Draw("P SAME");
    
    h_eff0_mc->SetLineColor(kCyan);
    h_eff0_mc->SetMarkerColor(kCyan);
    h_eff0_mc->SetMarkerStyle(24);
    h_eff0_mc->SetMarkerSize(0.5);
    h_eff0_mc->Draw("P SAME");
    
    h_eff1_mc->SetLineColor(kMagenta);
    h_eff1_mc->SetMarkerColor(kMagenta);
    h_eff1_mc->SetMarkerStyle(25);
    h_eff1_mc->SetMarkerSize(0.5);
    h_eff1_mc->Draw("P SAME");
    
    h_eff2_mc->SetLineColor(kYellow+2);
    h_eff2_mc->SetMarkerColor(kYellow+2);
    h_eff2_mc->SetMarkerStyle(26);
    h_eff2_mc->SetMarkerSize(0.5);
    h_eff2_mc->Draw("P SAME");
    
    TLegend* leg1 = new TLegend(0.65, 0.2, 0.85, 0.5);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.04);
    leg1->AddEntry(h_eff0_data, "Plane 0 (Data)", "lp");
    leg1->AddEntry(h_eff1_data, "Plane 1 (Data)", "lp");
    leg1->AddEntry(h_eff2_data, "Plane 2 (Data)", "lp");
    leg1->AddEntry((TObject*)0, "", "");
    leg1->AddEntry(h_eff0_mc, "Plane 0 (MC)", "lp");
    leg1->AddEntry(h_eff1_mc, "Plane 1 (MC)", "lp");
    leg1->AddEntry(h_eff2_mc, "Plane 2 (MC)", "lp");
    leg1->Draw();
    
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.025);
    latex->DrawLatex(0.6, 0.77, "Data Statistics:");
    latex->DrawLatex(0.6, 0.74, Form("Total Tracks: %d", total_tracks_data));
    latex->DrawLatex(0.6, 0.71, Form("Pitches - Min: %.6f, Max: %.5f", min_pitch_data, max_pitch_data));
    latex->DrawLatex(0.6, 0.68, Form("Mean: %.6f", mean_pitch_data));
    latex->DrawLatex(0.6, 0.65, Form("Efficiency - Mean: %.5f", mean_eff_data));
    
    latex->DrawLatex(0.6, 0.59, "MC Statistics:");
    latex->DrawLatex(0.6, 0.56, Form("Total Tracks: %d", total_tracks_mc));
    latex->DrawLatex(0.6, 0.53, Form("Pitches - Min: %.6f, Max: %.5f", min_pitch_mc, max_pitch_mc));
    latex->DrawLatex(0.6, 0.50, Form("Mean: %.6f", mean_pitch_mc));
    latex->DrawLatex(0.6, 0.47, Form("Efficiency - Mean: %.5f", mean_eff_mc));
    
    c1->SaveAs(Form("%s/hit_efficiency_vs_pitch.png", output_dir));
    
    // ========== PLOT 2: Mean Hit Efficiency vs Pitch (Binned) ==========
    float bin_width = (binned_xmax - binned_xmin) / nbins;
    std::vector<float> x_vals_data, y_vals_data, x_errs_data, y_errs_data;
    std::vector<float> x_vals_mc, y_vals_mc, x_errs_mc, y_errs_mc;
    
    // Data binning
    for (int b = 0; b < nbins; ++b) {
        float x_low = binned_xmin + b * bin_width;
        float x_high = x_low + bin_width;
        float x_center = 0.5 * (x_low + x_high);
        
        std::vector<float> eff_bin;
        for (size_t i = 0; i < all_pitches_data.size(); ++i) {
            if (all_pitches_data[i] >= x_low && all_pitches_data[i] < x_high)
                eff_bin.push_back(all_effs_data[i]);
        }
        
        if (!eff_bin.empty()) {
            float sum = std::accumulate(eff_bin.begin(), eff_bin.end(), 0.0);
            float mean = sum / eff_bin.size();
            float err = 0.0;
            for (auto v : eff_bin)
                err += (v - mean) * (v - mean);
            err = std::sqrt(err / eff_bin.size()) / std::sqrt(eff_bin.size());
            
            x_vals_data.push_back(x_center);
            y_vals_data.push_back(mean);
            x_errs_data.push_back(bin_width / 2);
            y_errs_data.push_back(err);
        }
    }
    
    // MC binning
    for (int b = 0; b < nbins; ++b) {
        float x_low = binned_xmin + b * bin_width;
        float x_high = x_low + bin_width;
        float x_center = 0.5 * (x_low + x_high);
        
        std::vector<float> eff_bin;
        for (size_t i = 0; i < all_pitches_mc.size(); ++i) {
            if (all_pitches_mc[i] >= x_low && all_pitches_mc[i] < x_high)
                eff_bin.push_back(all_effs_mc[i]);
        }
        
        if (!eff_bin.empty()) {
            float sum = std::accumulate(eff_bin.begin(), eff_bin.end(), 0.0);
            float mean = sum / eff_bin.size();
            float err = 0.0;
            for (auto v : eff_bin)
                err += (v - mean) * (v - mean);
            err = std::sqrt(err / eff_bin.size()) / std::sqrt(eff_bin.size());
            
            x_vals_mc.push_back(x_center);
            y_vals_mc.push_back(mean);
            x_errs_mc.push_back(bin_width / 2);
            y_errs_mc.push_back(err);
        }
    }
    
    TGraphErrors* graph_data = new TGraphErrors(x_vals_data.size(), x_vals_data.data(), y_vals_data.data(), x_errs_data.data(), y_errs_data.data());
    graph_data->SetTitle("Mean Hit Efficiency vs Pitch (All planes);Average Pitch [cm];Efficiency");
    graph_data->SetLineColor(kBlack);
    graph_data->SetMarkerStyle(20);
    graph_data->SetMarkerSize(0.8);
    graph_data->SetLineWidth(2);
    
    TGraphErrors* graph_mc = new TGraphErrors(x_vals_mc.size(), x_vals_mc.data(), y_vals_mc.data(), x_errs_mc.data(), y_errs_mc.data());
    graph_mc->SetLineColor(kBlue);
    graph_mc->SetMarkerStyle(24);
    graph_mc->SetMarkerSize(0.8);
    graph_mc->SetLineWidth(2);
    
    TCanvas* c3 = new TCanvas("c3", "Line Chart of Hit Efficiency", 800, 600);
    c3->SetLeftMargin(0.12);
    c3->SetBottomMargin(0.12);
    c3->SetGrid(1, 1);
    
    graph_data->GetYaxis()->SetRangeUser(binned_ymin, binned_ymax);
    graph_data->Draw("AP");
    graph_mc->Draw("P SAME");
    
    TLegend* leg2 = new TLegend(0.15, 0.15, 0.35, 0.25);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    leg2->AddEntry(graph_data, "Data", "lp");
    leg2->AddEntry(graph_mc, "MC", "lp");
    leg2->Draw();
    
    c3->SaveAs(Form("%s/mean_hit_efficiency.png", output_dir));
    
    // Print summary
    std::cout << "\nPlots saved in " << output_dir << "/:" << std::endl;
    std::cout << "1. hit_efficiency_vs_pitch.png" << std::endl;
    std::cout << "2. mean_hit_efficiency.png" << std::endl;
}