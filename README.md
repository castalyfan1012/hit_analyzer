# SBND Hit Efficiency Study Tools

**Important Notice**  
This repository was **originally designed for internal SBND collaboration use only**.  
It is not intended as a general-purpose, user-friendly public analysis package.  
Many parts assume specific file paths, SBN software environment (https://github.com/SBNSoftware/sbndcode), and internal access to samples.

Use with caution if you are not already familiar with the SBND offline software stack.

## Quick Overview

This set of scripts is used to study **hit efficiency** vs average pitch angle  
(mainly comparing Data vs Monte Carlo) using **SBND Ntuples**.

Typical workflow:
1. Setup environment  
2. Identify dead channels/wires  
3. Run hit efficiency analyzer (data + MC)  
4. Make comparison plots  
5. (Optional) Detailed event/wire inspection

## Prerequisites

- SL7 (Scientific Linux 7) container/environment  
- Working SBND offline software setup (usually `sbndcode` + dependencies)  
- Access to appropriate Ntuple files (both data and MC)  
- ROOT6 (usually comes with sbndcode)

## Typical Directory Structure & Important Files

```bash
/exp/sbnd/app/users/${USER}/hiteff_XXXX/
├── setup_simple.sh               # ← basic environment setup
├── get_xrootd.sh                 # (optional) prepare xrootd file lists
├── analyze_ntuple.C              # quick inspection of ntuple structure
├── dead_wires.C                  # find dead channels/wires → hit_wires.root, dead_channels.csv
├── hit_analyzer.C                # main hit efficiency analyzer → hiteff_data.csv, hiteff_mc.csv
├── hit_split_regions_data.C      # hit efficiency – split by TPC regions (data)
├── hit_split_regions_mc.C        # hit efficiency – split by TPC regions (MC)
├── hit_plotter.C                 # main comparison plots (hit eff vs avg pitch)
├── plot_split_regions.C          # comparison plots for split TPC regions
└── event_info_viewer.C           # interactive event/wire/timestamp browser
```

## Step-by-Step Quick Start (most common path)

1. Go to your working directory  
   `cd /exp/sbnd/app/users/castalyf/hiteff_2512`    # ← change to your actual working directory!

2. Setup environment  
   `source setup_simple.sh`

3. (Recommended) Quick check of ntuple structure  
   `root -l analyze_ntuple.C`  
   → make sure the input file path(s) inside the macro are correct!

4. Identify dead channels/wires (needed for masking)  
   `root -l dead_wires.C`  
   → produces `hit_wires.root` and `dead_channels.csv`

5. Run main hit efficiency calculation  
   (this step can take a long time → recommended to use nohup / screen / tmux)  
   `nohup root -l -b -q 'hit_analyzer.C' > log_analyzer.txt 2>&1 &`  
   If you want to study different TPC regions separately:  
   `nohup root -l -b -q 'hit_split_regions_data.C' > log_data_split.txt 2>&1 &`  
   `nohup root -l -b -q 'hit_split_regions_mc.C'   > log_mc_split.txt   2>&1 &`

6. Make standard hit efficiency vs average pitch plots  
   `root -l hit_plotter.C`

7. Make plots for different TPC regions (after split-region analysis)  
   `root -l plot_split_regions.C`

8. (Optional) Interactive viewer for debugging specific events/wires/timestamps  
    ```
   root 
   .L event_info_viewer.C 
   show_event_info("path/to/your/file.root")     # ← remember to change the path!
    ```

## Tips & Common Practices

- Most users run **data** and **MC** analysis in **separate terminals**  
  → easier to manage long jobs and compare results
- Always keep log files (use `> log_*.txt 2>&1` redirection)
- Double-check input file paths inside each `.C` macro before running
- The analyzer scripts (`hit_analyzer.C`, `hit_split_regions_*.C`) can take **several hours** depending on statistics and sample size
- Output CSV files (`hiteff_data.csv`, `hiteff_mc.csv`, etc.) are automatically read by the plotting macros

## Final Reminder

This is **collaboration-internal style code** — very much work-in-progress.  
File paths, ntuple branch names, and assumptions can change without notice.

