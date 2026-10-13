// ============================================================================
//  Alpha-gamma correlation heatmap and selection of correlated events
//
//  Usage :  root -l -b -q 'correlation.C+'
//
//  Reads every  coincidence_params/coincidence_s*.csv  (all columns are kept),
//  converts pulse heights to energies  (E = k * V, from the calibration),
//  plots the alpha-gamma heatmap (linear and log colour scale) and, once you
//  set the selection box, writes the events inside it to nice_correlated.csv.
//
//  Output: heatmap_linear.pdf/png, heatmap_log.pdf/png, nice_correlated.csv
// ============================================================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <filesystem>
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TBox.h"
#include "TAxis.h"

namespace fs = std::filesystem;

// ----------------------------- USER SETTINGS --------------------------------
const std::string COINC_FOLDER = "coincidence_params";   // folder with the csv files
const std::string FILE_PREFIX  = "coincidence_s";        // coincidence_s1.csv, _s2.csv ...

// column names in the csv files holding the pulse heights (V)
const std::string COL_ALPHA = "alpha_height_V";
const std::string COL_GAMMA = "gamma_height_V";

// calibration  E = k * V
const double K_ALPHA = 7.626;  const std::string UNIT_ALPHA = "MeV";   // MeV per V
const double K_GAMMA = 20.06;  const std::string UNIT_GAMMA = "keV";   // keV per V

// Heatmap binning and plotted range (0,0 = automatic range from the data)
const int    NBINS_X = 100, NBINS_Y = 100;
const double PLOT_A_MIN = 0, PLOT_A_MAX = 0;     // alpha energy range
const double PLOT_G_MIN = 0, PLOT_G_MAX = 0;     // gamma energy range

// Selection box = "correlated events" (alpha and gamma energy ranges).
// Leave at 0,0 for the first run (only the heatmaps are made); then fill in.
const double SEL_A_MIN = 0, SEL_A_MAX = 0;       // alpha energy (UNIT_ALPHA)
const double SEL_G_MIN = 0, SEL_G_MAX = 0;       // gamma energy (UNIT_GAMMA)
// ----------------------------------------------------------------------------

struct Row {
  std::string set;
  std::map<std::string, std::string> v;     // column name -> text, as in the file
  double Ea, Eg;
};

std::vector<std::string> Split(const std::string& line) {
  std::vector<std::string> out; std::string cur; std::stringstream ss(line);
  while (std::getline(ss, cur, ',')) {
    while (!cur.empty() && (cur.back() == '\r' || cur.back() == ' ')) cur.pop_back();
    out.push_back(cur);
  }
  return out;
}

void StyleSetup() {
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetTextFont(42);
  for (const char* ax : {"X", "Y", "Z"}) {
    gStyle->SetLabelFont(42, ax); gStyle->SetTitleFont(42, ax);
    gStyle->SetLabelSize(0.045, ax); gStyle->SetTitleSize(0.05, ax);
  }
  gStyle->SetTitleOffset(1.0, "X"); gStyle->SetTitleOffset(1.2, "Y");
  gStyle->SetTitleOffset(1.2, "Z");
  gStyle->SetPadLeftMargin(0.14); gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.17); gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255);
}

void correlation() {
  // ---------------- read all coincidence files --------------------------------
  std::vector<fs::path> files;
  for (auto& e : fs::directory_iterator(COINC_FOLDER)) {
    std::string name = e.path().filename().string();
    if (e.path().extension() == ".csv" && name.rfind(FILE_PREFIX, 0) == 0) files.push_back(e.path());
  }
  std::sort(files.begin(), files.end());
  if (files.empty()) { std::cerr << "No " << FILE_PREFIX << "*.csv in " << COINC_FOLDER << "\n"; return; }

  std::vector<std::string> columns;          // union of all columns, in order of appearance
  std::vector<Row> rows;
  for (auto& fp : files) {
    std::ifstream f(fp);
    std::string line;
    if (!std::getline(f, line)) continue;
    std::vector<std::string> hdr = Split(line);
    if (std::find(hdr.begin(), hdr.end(), COL_ALPHA) == hdr.end() ||
        std::find(hdr.begin(), hdr.end(), COL_GAMMA) == hdr.end()) {
      std::cerr << "Skipping " << fp.filename() << ": no " << COL_ALPHA << " / " << COL_GAMMA << " column\n";
      continue;
    }
    for (auto& c : hdr) if (std::find(columns.begin(), columns.end(), c) == columns.end()) columns.push_back(c);
    size_t n0 = rows.size();
    while (std::getline(f, line)) {
      auto fld = Split(line);
      if (fld.size() != hdr.size()) continue;
      Row r; r.set = fp.stem().string();
      for (size_t i = 0; i < hdr.size(); ++i) r.v[hdr[i]] = fld[i];
      r.Ea = K_ALPHA * std::atof(r.v[COL_ALPHA].c_str());
      r.Eg = K_GAMMA * std::atof(r.v[COL_GAMMA].c_str());
      rows.push_back(r);
    }
    std::cout << fp.filename().string() << ": " << rows.size() - n0 << " events\n";
  }
  std::cout << "Total events: " << rows.size() << "\n";
  if (rows.empty()) return;

  // ---------------- heatmaps ---------------------------------------------------
  double amin = PLOT_A_MIN, amax = PLOT_A_MAX, gmin = PLOT_G_MIN, gmax = PLOT_G_MAX;
  if (amax <= amin || gmax <= gmin) {
    amin = gmin = 1e30; amax = gmax = -1e30;
    for (auto& r : rows) {
      amin = std::min(amin, r.Ea); amax = std::max(amax, r.Ea);
      gmin = std::min(gmin, r.Eg); gmax = std::max(gmax, r.Eg);
    }
    double pa = 0.02 * (amax - amin > 0 ? amax - amin : 1), pg = 0.02 * (gmax - gmin > 0 ? gmax - gmin : 1);
    amin -= pa; amax += pa; gmin -= pg; gmax += pg;
  }
  StyleSetup();
  TH2D h("h2", "", NBINS_X, amin, amax, NBINS_Y, gmin, gmax);
  for (auto& r : rows) h.Fill(r.Ea, r.Eg);
  h.GetXaxis()->SetTitle(Form("Alpha energy (%s)", UNIT_ALPHA.c_str()));
  h.GetYaxis()->SetTitle(Form("Gamma energy (%s)", UNIT_GAMMA.c_str()));
  h.GetZaxis()->SetTitle("Counts / bin");
  h.SetMinimum(0.5);                          // empty bins stay white

  const bool haveSel = (SEL_A_MAX > SEL_A_MIN && SEL_G_MAX > SEL_G_MIN);
  for (int logz = 0; logz < 2; ++logz) {
    TCanvas c(Form("c%d", logz), "", 800, 700);
    c.SetLogz(logz);
    h.Draw("COLZ");
    if (haveSel) {
      TBox b(SEL_A_MIN, SEL_G_MIN, SEL_A_MAX, SEL_G_MAX);
      b.SetFillStyle(0); b.SetLineColor(kRed); b.SetLineWidth(2); b.SetLineStyle(2);
      b.DrawClone();
    }
    std::string name = logz ? "heatmap_log" : "heatmap_linear";
    c.SaveAs((name + ".pdf").c_str());
    c.SaveAs((name + ".png").c_str());
  }

  // ---------------- correlated events -> csv -----------------------------------
  if (!haveSel) {
    std::cout << "\nSelection box not set (SEL_* = 0): only heatmaps made.\n"
                 "Look at heatmap_*.png, set SEL_A_MIN/MAX and SEL_G_MIN/MAX, and rerun.\n";
    return;
  }
  std::ofstream o("nice_correlated.csv");
  o << "set";
  for (auto& c : columns) o << "," << c;
  o << ",alpha_E_" << UNIT_ALPHA << ",gamma_E_" << UNIT_GAMMA << "\n";
  std::map<std::string, int> perSet; int nSel = 0;
  for (auto& r : rows) {
    if (r.Ea < SEL_A_MIN || r.Ea > SEL_A_MAX || r.Eg < SEL_G_MIN || r.Eg > SEL_G_MAX) continue;
    o << r.set;
    for (auto& c : columns) { auto it = r.v.find(c); o << "," << (it == r.v.end() ? "" : it->second); }
    o << "," << r.Ea << "," << r.Eg << "\n";
    ++perSet[r.set]; ++nSel;
  }
  std::cout << "\nSelection: alpha " << SEL_A_MIN << "-" << SEL_A_MAX << " " << UNIT_ALPHA
            << ",  gamma " << SEL_G_MIN << "-" << SEL_G_MAX << " " << UNIT_GAMMA << "\n";
  for (auto& [s, n] : perSet) std::cout << "  " << s << ": " << n << "\n";
  std::cout << "Wrote nice_correlated.csv (" << nSel << " of " << rows.size() << " events)\n";
}
