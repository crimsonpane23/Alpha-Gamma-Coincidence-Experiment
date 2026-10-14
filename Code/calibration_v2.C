// ============================================================================
//  One-point energy calibration (alpha or gamma) from LeCroy WaveSurfer CSVs
//
//  Usage :  root -l -b -q 'calibration_v2.C+'
//
//  Folder layout expected:  BASE_PATH/<any 5 folders>/Autosave/C3agc_00000.csv ...
//  Steps : 1) pulse height (peak - baseline) of every event in every folder
//          2) pulse-height spectrum
//          3) Gaussian fit in [FIT_MIN, FIT_MAX]  (whole range if not set)
//          4) one-point calibration through the origin:  E = k * V,  k = E0 / mu
//  Output: spectrum_<signal>.pdf/png          (pulse height + Gaussian fit, publication plot)
//          calibration_<signal>.txt           (fit result + equation)
//          (EXTRA_PLOTS = true also makes rms_<signal> and spectrum_<signal>_energy)
//          pulses_<signal>.csv                (cache, see USE_CACHE)
// ============================================================================
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TColor.h"
#include "TAxis.h"

namespace fs = std::filesystem;

// ----------------------------- USER SETTINGS --------------------------------
const std::string SIGNAL    = "alpha";        // "alpha" (CH3)  or  "gamma" (CH4)
const std::string BASE_PATH = "/home/aryan-shrivastava/Documents/Integrated_OpenLAB_Alpha_gamma_coincidence/Data_analysis_agc_calib/agc_final_calib_data/alpha_calib";   // contains the 5 folders

const std::string PREFIX_ALPHA = "C3", PREFIX_GAMMA = "C4";
const int POLARITY_ALPHA = +1, POLARITY_GAMMA = +1;          // +1 positive pulses, -1 negative

// Known energy of the peak you will fit (241Am: alpha 5.486 MeV, gamma 59.54 keV)
const double      KNOWN_ENERGY_ALPHA = 5.486;   const std::string UNIT_ALPHA = "MeV";
const double      KNOWN_ENERGY_GAMMA = 59.54;   const std::string UNIT_GAMMA = "keV";

// Time windows (s): baseline and peak search. Sample file: pulse starts at ~0,
// peaks at ~1 us, decays over several us; file spans -3 us ... +7 us.
const double BASE_T0 = -2.0e-6, BASE_T1 = -1.0e-6;
const double PEAK_T0 = -1.0e-6, PEAK_T1 = +4.0e-6;

// Ignore events whose baseline RMS (V) exceeds this (<= 0 : no cut)
const double RMS_CUT = 35e-3;

// Pulse-height threshold (V): events with height below it are ignored completely
// (not in the spectrum, not in the RMS plot, not counted).  <= 0 : no threshold
const double HEIGHT_MIN_ALPHA = 0.5;
const double HEIGHT_MIN_GAMMA = 2.3;

// Spectrum binning (V). HIST_MAX <= HIST_MIN : automatic range
const double HIST_MIN = 0.035, HIST_MAX = 1.81;
const int    NBINS    = 50;

// Gaussian fit range (V). FIT_MAX <= FIT_MIN : fit the whole histogram range
const double FIT_MIN = 0.6, FIT_MAX = 0.95;

// Reading ~5000 files is slow; the pulse heights are cached in pulses_<signal>.csv.
// Set true to reuse the cache when it exists (delete the file if you change the
// time windows or the signal files; the RMS cut and fit range can change freely).
const bool USE_CACHE = true;

// ---- plot style (does not affect any number) ----
const std::string COLOR_ALPHA = "#0379bd";   // histogram colour, alpha (hex)
const std::string COLOR_GAMMA = "#fb0633";   // histogram colour, gamma (hex)
const bool ERROR_BARS  = true;               // Poisson error bars on the histogram
const bool EXTRA_PLOTS = false;              // true: also RMS distribution + energy spectrum
// ----------------------------------------------------------------------------

struct Pulse { double height, rms; };

// Pulse height = peak - baseline, from one CSV file
bool ReadPulse(const std::string& path, int pol, Pulse& p) {
  std::ifstream f(path);
  if (!f) return false;
  std::string line;
  for (int i = 0; i < 5; ++i) std::getline(f, line);        // 5 header lines
  double s = 0, s2 = 0, vmax = -1e30; long nb = 0; bool havePeak = false;
  const double tEnd = std::max(BASE_T1, PEAK_T1);
  while (std::getline(f, line)) {
    double t, v;
    if (std::sscanf(line.c_str(), "%lf,%lf", &t, &v) != 2) continue;
    if (t > tEnd) break;
    v *= pol;                                               // pulse positive
    if (t >= BASE_T0 && t <= BASE_T1) { s += v; s2 += v * v; ++nb; }
    if (t >= PEAK_T0 && t <= PEAK_T1) { vmax = std::max(vmax, v); havePeak = true; }
  }
  if (nb < 2 || !havePeak) return false;
  double base = s / nb;
  p.rms = std::sqrt(std::max(0.0, (s2 - nb * base * base) / (nb - 1)));
  p.height = vmax - base;
  return true;
}

void StyleSetup() {
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetTextFont(42);
  for (const char* ax : {"X", "Y"}) {
    gStyle->SetLabelFont(42, ax); gStyle->SetTitleFont(42, ax);
    gStyle->SetLabelSize(0.045, ax); gStyle->SetTitleSize(0.05, ax);
  }
  gStyle->SetTitleOffset(1.0, "X"); gStyle->SetTitleOffset(1.2, "Y");
  gStyle->SetPadLeftMargin(0.14); gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.05); gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  gStyle->SetHistLineWidth(2); gStyle->SetEndErrorSize(0);
}

void calibration_v2() {
  const bool isAlpha = (SIGNAL == "alpha");
  if (!isAlpha && SIGNAL != "gamma") { std::cerr << "SIGNAL must be \"alpha\" or \"gamma\"\n"; return; }
  const std::string prefix = isAlpha ? PREFIX_ALPHA : PREFIX_GAMMA;
  const int pol            = isAlpha ? POLARITY_ALPHA : POLARITY_GAMMA;
  const double E0          = isAlpha ? KNOWN_ENERGY_ALPHA : KNOWN_ENERGY_GAMMA;
  const std::string unit   = isAlpha ? UNIT_ALPHA : UNIT_GAMMA;
  const double hMin        = isAlpha ? HEIGHT_MIN_ALPHA : HEIGHT_MIN_GAMMA;
  const std::string cache  = "pulses_" + SIGNAL + ".csv";

  // ---------------- 1) pulse heights (cache or read all folders) -------------
  std::vector<Pulse> all;
  if (USE_CACHE && fs::exists(cache)) {
    std::ifstream c(cache); std::string line; std::getline(c, line);
    Pulse p;
    while (std::getline(c, line))
      if (std::sscanf(line.c_str(), "%lf,%lf", &p.height, &p.rms) == 2) all.push_back(p);
    std::cout << "Loaded " << all.size() << " pulses from cache " << cache << "\n";
  } else {
    std::vector<fs::path> dirs;
    for (auto& d : fs::directory_iterator(BASE_PATH))
      if (d.is_directory() && fs::exists(d.path() / "Autosave")) dirs.push_back(d.path() / "Autosave");
    std::sort(dirs.begin(), dirs.end());
    if (dirs.empty()) { std::cerr << "No <folder>/Autosave found in " << BASE_PATH << "\n"; return; }
    for (auto& dir : dirs) {
      size_t n0 = all.size();
      for (auto& e : fs::directory_iterator(dir)) {
        std::string name = e.path().filename().string();
        if (e.path().extension() != ".csv" || name.rfind(prefix, 0) != 0) continue;
        Pulse p;
        if (ReadPulse(e.path().string(), pol, p)) all.push_back(p);
      }
      std::cout << dir.string() << ": " << all.size() - n0 << " events\n";
    }
    std::ofstream c(cache);
    c << "height_V,baseline_rms_V\n";
    for (auto& p : all) c << p.height << "," << p.rms << "\n";
  }

  // height threshold: these events are dropped everywhere below
  if (hMin > 0) {
    size_t n0 = all.size();
    all.erase(std::remove_if(all.begin(), all.end(),
              [&](const Pulse& p) { return p.height < hMin; }), all.end());
    std::cout << "Height threshold " << hMin << " V removed " << n0 - all.size()
              << " of " << n0 << " events\n";
  }

  std::vector<double> h;
  for (auto& p : all) if (RMS_CUT <= 0 || p.rms <= RMS_CUT) h.push_back(p.height);
  std::cout << "Events used: " << h.size() << " of " << all.size() << "\n";
  if (h.empty()) return;

  // ---------------- 2) spectrum ----------------------------------------------
  double lo = HIST_MIN, hi = HIST_MAX;
  if (hi <= lo) {
    lo = *std::min_element(h.begin(), h.end());
    hi = *std::max_element(h.begin(), h.end());
    double pad = 0.02 * (hi - lo > 0 ? hi - lo : 1);
    lo -= pad; hi += pad;
  }
  StyleSetup();

  // ---- baseline-RMS distribution (all events, before the RMS cut) ----
  if (EXTRA_PLOTS) {
    double rmax = 0;
    for (auto& p : all) rmax = std::max(rmax, p.rms * 1e3);
    TH1D hr("hr", "", 100, 0, 1.05 * (rmax > 0 ? rmax : 1));
    for (auto& p : all) hr.Fill(p.rms * 1e3);
    const int rc = TColor::GetColor((isAlpha ? COLOR_ALPHA : COLOR_GAMMA).c_str());
    hr.GetXaxis()->SetTitle(Form("%s baseline RMS (mV)", isAlpha ? "Alpha" : "Gamma"));
    hr.GetYaxis()->SetTitle("Counts");
    hr.SetLineColor(rc); hr.SetFillColorAlpha(rc, 0.25);
    TCanvas cr("cr", "", 800, 600);
    //cr.SetLogy();
    hr.Draw("HIST");
    if (RMS_CUT > 0) {
      TLine l; l.SetLineColor(kBlack); l.SetLineStyle(2); l.SetLineWidth(2);
      l.DrawLine(RMS_CUT * 1e3, hr.GetMinimum(), RMS_CUT * 1e3, 1.05 * hr.GetMaximum());
    }
    cr.SaveAs(("rms_" + SIGNAL + ".pdf").c_str());
    cr.SaveAs(("rms_" + SIGNAL + ".png").c_str());
  }

  TH1D hist("hist", "", NBINS, lo, hi);
  for (double v : h) hist.Fill(v);
  const double bw = (hi - lo) / NBINS;
  hist.GetXaxis()->SetTitle(Form("%s pulse height (V)", isAlpha ? "Alpha" : "Gamma"));
  hist.GetYaxis()->SetTitle(Form("Counts / %.3g V", bw));
  const int col = TColor::GetColor((isAlpha ? COLOR_ALPHA : COLOR_GAMMA).c_str());
  hist.SetLineColor(col); hist.SetFillColorAlpha(col, 0.30);
  hist.SetMarkerStyle(20); hist.SetMarkerSize(0.6); hist.SetMarkerColor(col);
  hist.SetMaximum(1.3 * hist.GetMaximum());

  // ---------------- 3) Gaussian fit -------------------------------------------
  double flo = FIT_MIN, fhi = FIT_MAX;
  if (fhi <= flo) { flo = lo; fhi = hi; }
  hist.GetXaxis()->SetRangeUser(flo, fhi);
  double mean0 = hist.GetMean(), rms0 = hist.GetRMS(), amp0 = hist.GetMaximum();
  hist.GetXaxis()->SetRange(0, 0);                          // restore full range
  TF1 f("gfit", "gaus", flo, fhi);
  f.SetParameters(amp0 / 1.3, mean0, rms0);
  int status = hist.Fit(&f, "RLQ0");                        // Poisson likelihood
  if (status != 0) std::cerr << "Warning: fit status " << status << " (check range)\n";
  double mu = f.GetParameter(1), dmu = f.GetParError(1);
  double sig = std::fabs(f.GetParameter(2)), dsig = f.GetParError(2);
  double res = 2.3548 * sig / mu * 100;                     // FWHM / mu in %
  double chi2 = f.GetChisquare(), ndof = f.GetNDF();

  // ---------------- 4) one-point calibration  E = k * V -----------------------
  double k = E0 / mu, dk = k * dmu / mu;

  std::ostringstream out;
  out << "Signal            : " << SIGNAL << "\n"
      << "Events in spectrum: " << h.size() << "\n"
      << "Fit range         : " << flo << " - " << fhi << " V\n"
      << "Peak mean mu      : " << mu << " +- " << dmu << " V\n"
      << "Peak sigma        : " << sig << " +- " << dsig << " V\n"
      << "Resolution FWHM/mu: " << res << " %\n"
      << "Known energy E0   : " << E0 << " " << unit << "\n"
      << "Calibration       : E [" << unit << "] = (" << k << " +- " << dk << ") * V [V]\n"
      << "Gain              : " << 1.0 / k << " V per " << unit << "\n";
  std::cout << "\n" << out.str();
  std::ofstream("calibration_" + SIGNAL + ".txt") << out.str();

  // ---------------- plot: spectrum + Gaussian fit -----------------------------
  const int grey = TColor::GetColor("#9a9a9a");
  hist.SetMinimum(0);
  const double ymin = 0, ymax = hist.GetMaximum();
  TCanvas c("c", "", 800, 600);
  hist.Draw("AXIS");
  TBox shade(flo, ymin, fhi, ymax);              // lightly shaded fitted region
  shade.SetFillColorAlpha(grey, 0.15); shade.SetLineColor(0);
  shade.DrawClone();
  TLine edge;                                     // light lines marking the fit range
  edge.SetLineColor(grey); edge.SetLineWidth(2);
  edge.DrawLine(flo, ymin, flo, ymax);
  edge.DrawLine(fhi, ymin, fhi, ymax);
  hist.Draw("HIST SAME");
  if (ERROR_BARS) hist.Draw("E SAME");
  f.SetLineColor(kBlack); f.SetLineWidth(3);
  f.DrawCopy("SAME");
  gPad->RedrawAxis();
  TLatex tx; tx.SetNDC(); tx.SetTextFont(42);
  tx.SetTextSize(0.05);
  tx.DrawLatex(0.18, 0.86, isAlpha ? "Alpha spectrum" : "Gamma spectrum");
  tx.SetTextSize(0.042);
  tx.DrawLatex(0.60, 0.86, Form("Fit results"));
  tx.DrawLatex(0.60, 0.79, Form("#mu = %.4g #pm %.2g V", mu, dmu));
  tx.DrawLatex(0.60, 0.72, Form("#sigma = %.3g #pm %.2g V", sig, dsig));
  tx.DrawLatex(0.60, 0.65, Form("#chi^{2}/ndof = %.2f / %.0f", chi2, ndof));
  tx.DrawLatex(0.60, 0.58, Form("E = %.3g #times V  (%s)", k, unit.c_str()));
  c.SaveAs(("spectrum_" + SIGNAL + ".pdf").c_str());
  c.SaveAs(("spectrum_" + SIGNAL + ".png").c_str());

  if (EXTRA_PLOTS) {
  TH1D hE("hE", "", NBINS, lo * k, hi * k);
  for (double v : h) hE.Fill(v * k);
  hE.GetXaxis()->SetTitle(Form("%s energy (%s)", isAlpha ? "Alpha" : "Gamma", unit.c_str()));
  hE.GetYaxis()->SetTitle(Form("Counts / %.3g %s", bw * k, unit.c_str()));
  hE.SetLineColor(col); hE.SetFillColorAlpha(col, 0.25);
  TCanvas c2("c2", "", 800, 600);
  hE.Draw("HIST");
  c2.SaveAs(("spectrum_" + SIGNAL + "_energy.pdf").c_str());
  c2.SaveAs(("spectrum_" + SIGNAL + "_energy.png").c_str());
  }
}
