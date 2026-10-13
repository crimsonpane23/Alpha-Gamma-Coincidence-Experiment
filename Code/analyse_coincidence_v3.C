// ============================================================================
//  Alpha-gamma coincidence analysis  (LeCroy WaveSurfer 3034z CSV files)
//  CH3 = alpha, CH4 = gamma
//
//  Usage :  root -l -b -q 'analyse_coincidence.C+'
//  Output:  alpha_CH3.csv, gamma_CH4.csv, coincidence.csv, outliers.csv,
//           PDF/PNG plots, event_checks/event_XXXXX.png (waveforms to inspect)
//
//  Per event: peak height, peak time, trigger-crossing time,
//             constant-fraction (CF) time, baseline RMS, saturation flag
//  Delta t is built three ways: peak, fixed-trigger crossing, CF crossing.
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
#ifndef NO_ROOT
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#endif

// ----------------------------- USER SETTINGS --------------------------------
const int         SETTING        = 3;                             // 1, 2 or 3
const std::string FOLDER         = "/home/aryan-shrivastava/Documents/Integrated_OpenLAB_Alpha_gamma_coincidence/Data_analysis_agc_final/Raw_Data/Auto7_40kS";
const std::string PREFIX_ALPHA   = "C3";                          // channel 3
const std::string PREFIX_GAMMA   = "C4";                          // channel 4

const double TRIGGER_ALPHA = 0.60;                                // V
const double TRIGGER_GAMMA = 2.30;                                // V
const int    POLARITY_ALPHA = +1;                                 // +1 positive pulses, -1 negative
const int    POLARITY_GAMMA = +1;

// Baseline-RMS cutoffs in V (events with RMS > cutoff are ignored; <= 0 : no cut)
const double RMS_CUT_ALPHA = 20e-3;  // 20 mV
const double RMS_CUT_GAMMA = 30e-3;  // 50 mV

// Constant-fraction timing: time at which the baseline-subtracted pulse reaches
// this fraction of its own peak height (leading edge, linear interpolation)
const double CF_FRACTION = 0.20;

// Saturation: raw peak (V) at/above which a pulse is flagged as clipped.
// Set from the flat top you see in the waveforms (alpha: 1e9 = never flagged)
const double SAT_V_ALPHA = 1e9;
const double SAT_V_GAMMA = 6.2;
const bool   REMOVE_SATURATED = true;                            // true : drop flagged events

// Delta-t histograms
const double DT_MIN = -1500, DT_MAX = 3000;                       // ns
const int    DT_BINS = 90;

// Outlier inspection: events with |dt_peak| > OUTLIER_NS are listed in
// outliers.csv and their waveforms drawn (at most MAX_EVENT_PLOTS of them).
// Add any event numbers you want to look at to EVENTS_TO_VIEW.
const double OUTLIER_NS     = 2000;
const int    MAX_EVENT_PLOTS = 20;
const std::vector<int> EVENTS_TO_VIEW = {};
// ----------------------------------------------------------------------------

const double ROI_MIN = -4.0e-6, ROI_MAX = +4.0e-6;   // s (same for all settings)

struct Setting { long nSamples; double rate; double window; };
Setting GetSetting(int s) {
  switch (s) {
    case 1: return {50000,  500e6, 100e-6};
    case 2: return {200000, 2e9,   100e-6};
    case 3: return {40000,  2e9,   20e-6};
  }
  std::cerr << "SETTING must be 1, 2 or 3\n"; std::exit(1);
}

struct Result {
  int event = -1;
  double height = 0, tPeak = 0, tTrig = -999, tCF = -999, baseRMS = 0, vmax = 0;
  bool ok = false, okCF = false, sat = false;
  std::string path;
};

// Read the samples inside the ROI (voltage flipped so the pulse is positive)
bool ReadROI(const std::string& path, int pol,
             std::vector<double>& t, std::vector<double>& v) {
  std::ifstream f(path);
  if (!f) return false;
  std::string line;
  for (int i = 0; i < 5; ++i) std::getline(f, line);   // 5 header lines
  while (std::getline(f, line)) {
    double tt, vv;
    if (std::sscanf(line.c_str(), "%lf,%lf", &tt, &vv) != 2) continue;
    if (tt < ROI_MIN) continue;
    if (tt > ROI_MAX) break;
    t.push_back(tt); v.push_back(pol * vv);
  }
  return t.size() >= 2;
}

bool AnalyseFile(const std::string& path, double trig, int pol, double satV,
                 double expectedDt, Result& r) {
  std::vector<double> t, v;
  if (!ReadROI(path, pol, t, v)) return false;
  r.path = path;

  double dt = (t.back() - t.front()) / (t.size() - 1);
  if (std::fabs(dt - expectedDt) > 0.02 * expectedDt)
    std::cerr << "Warning: " << path << " sample spacing " << dt
              << " s differs from setting (" << expectedDt << " s)\n";

  // baseline (average of -4us to -2us) and its RMS
  double base = 0; int nBase = 0; double baserms = 0;
  for (size_t i = 0; i < t.size(); ++i)
    if (t[i] >= -4e-6 && t[i] <= -2e-6) { base += v[i]; ++nBase; }
  if (nBase > 0) base /= nBase;
  if (nBase > 1) {
    for (size_t i = 0; i < t.size(); ++i)
      if (t[i] >= -4e-6 && t[i] <= -2e-6) baserms += (v[i] - base) * (v[i] - base);
    baserms = std::sqrt(baserms / (nBase - 1));
  }
  r.baseRMS = baserms;

  // peak
  size_t ip = std::max_element(v.begin(), v.end()) - v.begin();
  r.height = pol * (v[ip] - base);  r.tPeak = t[ip];
  r.vmax = v[ip];  r.sat = (v[ip] >= satV);

  // constant-fraction time: walk back from the peak to the leading-edge crossing
  // of base + CF_FRACTION*(peak-base), then interpolate linearly
  double amp = v[ip] - base;
  double lvl = base + CF_FRACTION * amp;
  r.tCF = -999; r.okCF = false;
  if (amp > 0) {
    size_t i = ip;
    while (i > 0 && v[i-1] >= lvl) --i;
    if (i > 0) {
      r.tCF = t[i-1] + (lvl - v[i-1]) / (v[i] - v[i-1]) * (t[i] - t[i-1]);
      r.okCF = true;
    }
  }

  // first positive-edge crossing of the fixed trigger level
  const double th = pol * trig;
  r.tTrig = -999; r.ok = false;
  for (size_t i = 1; i < v.size(); ++i)
    if (v[i-1] < th && v[i] >= th) {
      //r.tTrig = t[i-1] + (th - v[i-1]) / (v[i] - v[i-1]) * (t[i] - t[i-1]);
      r.tTrig = t[i];
      r.ok = true; break;
    }
  return true;
}

// event number = last group of digits in the file name (e.g. C3agc_00010 -> 10)
int EventNumber(const std::string& stem) {
  size_t e = stem.find_last_of("0123456789");
  if (e == std::string::npos) return -1;
  size_t b = e; while (b > 0 && std::isdigit(stem[b-1])) --b;
  return std::atoi(stem.substr(b, e - b + 1).c_str());
}

std::map<int, Result> ReadChannel(const std::string& prefix, double trig,
                                  int pol, double satV, double dt) {
  std::map<int, Result> out;
  for (auto& p : std::filesystem::directory_iterator(FOLDER)) {
    std::string name = p.path().filename().string();
    if (p.path().extension() != ".csv" || name.rfind(prefix, 0) != 0) continue;
    Result r; r.event = EventNumber(p.path().stem().string());
    if (AnalyseFile(p.path().string(), trig, pol, satV, dt, r)) out[r.event] = r;
  }
  return out;
}

// Quality cuts: baseline RMS, no trigger crossing found, (optionally) saturated
void ApplyCuts(std::map<int, Result>& m, double rmsCut, bool dropSat,
               const std::string& label) {
  size_t n0 = m.size(), nRMS = 0, nNoTrig = 0, nSat = 0;
  for (auto it = m.begin(); it != m.end();) {
    const Result& r = it->second;
    bool bad = true;
    if (rmsCut > 0 && r.baseRMS > rmsCut) ++nRMS;
    else if (!r.ok)                       ++nNoTrig;
    else if (dropSat && r.sat)            ++nSat;
    else bad = false;
    it = bad ? m.erase(it) : std::next(it);
  }
  std::cout << label << ": kept " << m.size() << " of " << n0
            << "  (removed: RMS " << nRMS << ", no trigger crossing " << nNoTrig
            << ", saturated " << nSat << ")\n";
}

void WriteCSV(const std::string& fn, const std::map<int, Result>& m) {
  std::ofstream o(fn);
  o << "event,peak_height_V,peak_time_ns,trigger_cross_time_ns,cf_time_ns,baseline_rms_V,saturated\n";
  for (auto& [ev, r] : m)
    o << ev << "," << r.height << "," << r.tPeak * 1e9 << ","
      << (r.ok ? r.tTrig * 1e9 : -999) << "," << (r.okCF ? r.tCF * 1e9 : -999) << ","
      << r.baseRMS << "," << r.sat << "\n";
  std::cout << "Wrote " << fn << " (" << m.size() << " events)\n";
}

#ifndef NO_ROOT
// ------------------------------ plotting ------------------------------------
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

void Plot(const std::vector<double>& x, const std::string& xt,
          const std::string& file, int color, int nbins = 100) {
  if (x.empty()) return;
  double lo = *std::min_element(x.begin(), x.end());
  double hi = *std::max_element(x.begin(), x.end());
  double pad = 0.05 * (hi - lo > 0 ? hi - lo : 1);
  TH1D h(file.c_str(), "", nbins, lo - pad, hi + pad);
  for (double v : x) h.Fill(v);
  h.GetXaxis()->SetTitle(xt.c_str());
  h.GetYaxis()->SetTitle("Counts");
  h.SetLineColor(color); h.SetFillColorAlpha(color, 0.25);
  TCanvas c("c", "", 800, 600);
  h.Draw("HIST");
  c.SaveAs((file + ".pdf").c_str());
  c.SaveAs((file + ".png").c_str());
}

// Baseline-RMS distribution (mV), dashed line at the cutoff if one is set
void PlotRMS(const std::map<int, Result>& m, const std::string& xt,
             const std::string& file, int color, double cutV) {
  if (m.empty()) return;
  double hi = 0;
  for (auto& [ev, r] : m) hi = std::max(hi, r.baseRMS * 1e3);
  TH1D h(file.c_str(), "", 50, 0, 1.05 * (hi > 0 ? hi : 1));
  for (auto& [ev, r] : m) h.Fill(r.baseRMS * 1e3);
  h.GetXaxis()->SetTitle(xt.c_str());
  h.GetYaxis()->SetTitle("Counts");
  h.SetLineColor(color); h.SetFillColorAlpha(color, 0.25);
  TCanvas c("c_rms", "", 800, 600);
  c.SetLogy();
  h.Draw("HIST");
  TLine line;
  if (cutV > 0) {
    line.SetLineColor(kBlack); line.SetLineStyle(2); line.SetLineWidth(2);
    line.DrawLine(cutV * 1e3, h.GetMinimum(), cutV * 1e3, 1.05 * h.GetMaximum());
  }
  c.SaveAs((file + ".pdf").c_str());
  c.SaveAs((file + ".png").c_str());
}

// Delta-t histogram, saved twice:  <file>  (data only)  and  <file>_fit
// (data + fit with an exponential convolved with a Gaussian, i.e. an
//  exponentially modified Gaussian:  tau = lifetime, sigma = timing resolution)
void PlotDeltaT(const std::vector<double>& x, const std::string& xt,
                const std::string& file, int color) {
  if (x.empty()) return;
  TH1D h(file.c_str(), "", DT_BINS, DT_MIN, DT_MAX);
  int nOut = 0;
  for (double v : x) { if (v < DT_MIN || v >= DT_MAX) ++nOut; h.Fill(v); }
  std::cout << file << ": " << nOut << " of " << x.size()
            << " events outside [" << DT_MIN << ", " << DT_MAX << "] ns (not shown)\n";
  double bw = (DT_MAX - DT_MIN) / DT_BINS;
  h.GetXaxis()->SetTitle(xt.c_str());
  h.GetYaxis()->SetTitle(Form("Counts / %.0f ns", bw));
  h.SetLineColor(color); h.SetFillColorAlpha(color, 0.25);
  h.SetMaximum(1.25 * h.GetMaximum());

  TCanvas c("c_dt", "", 800, 600);
  h.Draw("HIST");
  c.SaveAs((file + ".pdf").c_str());
  c.SaveAs((file + ".png").c_str());

  TF1 f("emg",
        "[0]/(2*[3])*exp([2]*[2]/(2*[3]*[3])-(x-[1])/[3])"
        "*TMath::Erfc(([2]*[2]/[3]-(x-[1]))/(sqrt(2)*[2]))", DT_MIN, DT_MAX);
  f.SetParNames("N", "mu", "sigma", "tau");
  f.SetParameters(h.Integral() * bw, h.GetBinCenter(h.GetMaximumBin()) - 100, 200, 150);
  f.SetParLimits(2, 5, 2000);
  f.SetParLimits(3, 5, 5000);
  h.Fit(&f, "RLQ0");                       // Poisson likelihood fit, no drawing
  std::cout << "  fit " << file << ":  mu = " << f.GetParameter(1)
            << " +- " << f.GetParError(1) << " ns,  sigma = " << f.GetParameter(2)
            << " +- " << f.GetParError(2) << " ns,  tau = " << f.GetParameter(3)
            << " +- " << f.GetParError(3) << " ns\n";

  h.Draw("HIST");
  f.SetLineColor(kBlack); f.SetLineWidth(2);
  f.DrawCopy("same");
  TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.042);
  tx.DrawLatex(0.60, 0.86, Form("#tau = %.0f #pm %.0f ns", f.GetParameter(3), f.GetParError(3)));
  tx.DrawLatex(0.60, 0.80, Form("#sigma = %.0f #pm %.0f ns", f.GetParameter(2), f.GetParError(2)));
  c.SaveAs((file + "_fit.pdf").c_str());
  c.SaveAs((file + "_fit.png").c_str());
}

// Waveforms of one event (alpha top, gamma bottom) for manual checking.
// dashed black = peak time, dotted red = CF time
void PlotEvent(const Result& a, const Result& g) {
  gSystem->mkdir("event_checks", true);
  TCanvas c(Form("cev%d", a.event), "", 900, 700);
  c.Divide(1, 2);
  for (int k = 0; k < 2; ++k) {
    const Result& r = (k == 0) ? a : g;
    std::vector<double> t, v;
    if (!ReadROI(r.path, k == 0 ? POLARITY_ALPHA : POLARITY_GAMMA, t, v)) continue;
    for (auto& x : t) x *= 1e9;
    c.cd(k + 1);
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.16);
    auto* gr = new TGraph(t.size(), t.data(), v.data());
    gr->SetLineColor(k == 0 ? kBlue + 1 : kRed + 1);
    gr->GetXaxis()->SetTitle("Time (ns)");
    gr->GetYaxis()->SetTitle(k == 0 ? "Alpha (V)" : "Gamma (V)");
    gr->GetXaxis()->SetTitleSize(0.07); gr->GetYaxis()->SetTitleSize(0.07);
    gr->GetXaxis()->SetLabelSize(0.06); gr->GetYaxis()->SetLabelSize(0.06);
    gr->GetYaxis()->SetTitleOffset(0.7);
    gr->Draw("AL");
    gPad->Update();
    double y0 = gPad->GetUymin(), y1 = gPad->GetUymax();
    TLine l;
    l.SetLineStyle(2); l.SetLineColor(kBlack); l.DrawLine(r.tPeak * 1e9, y0, r.tPeak * 1e9, y1);
    if (r.okCF) { l.SetLineStyle(3); l.SetLineColor(kRed); l.DrawLine(r.tCF * 1e9, y0, r.tCF * 1e9, y1); }
  }
  c.SaveAs(Form("event_checks/event_%05d.png", a.event));
}
#endif

// -------------------------------- main --------------------------------------
void analyse_coincidence_v3() {
  Setting s = GetSetting(SETTING);
  double dt = 1.0 / s.rate;
  std::cout << "Setting " << SETTING << ": " << s.nSamples << " samples, "
            << s.rate / 1e6 << " MS/s, window " << s.window * 1e6 << " us\n";

  auto alpha = ReadChannel(PREFIX_ALPHA, TRIGGER_ALPHA, POLARITY_ALPHA, SAT_V_ALPHA, dt);
  auto gamma = ReadChannel(PREFIX_GAMMA, TRIGGER_GAMMA, POLARITY_GAMMA, SAT_V_GAMMA, dt);

#ifndef NO_ROOT
  StyleSetup();
  // RMS distributions are made BEFORE the cuts
  PlotRMS(alpha, "Alpha baseline RMS (mV)", "alpha_baseline_rms", kBlue + 1, RMS_CUT_ALPHA);
  PlotRMS(gamma, "Gamma baseline RMS (mV)", "gamma_baseline_rms", kRed + 1,  RMS_CUT_GAMMA);
#endif

  ApplyCuts(alpha, RMS_CUT_ALPHA, REMOVE_SATURATED, "Alpha");
  ApplyCuts(gamma, RMS_CUT_GAMMA, REMOVE_SATURATED, "Gamma");
  WriteCSV("alpha_CH3.csv", alpha);
  WriteCSV("gamma_CH4.csv", gamma);

  // combined file for events present in both channels (after the cuts)
  std::ofstream o("coincidence.csv");
  o << "event,alpha_height_V,alpha_peak_ns,alpha_trig_ns,alpha_cf_ns,"
       "gamma_height_V,gamma_peak_ns,gamma_trig_ns,gamma_cf_ns,gamma_saturated\n";
  std::ofstream oo("outliers.csv");
  oo << "event,dt_peak_ns,alpha_peak_ns,gamma_peak_ns,alpha_height_V,gamma_height_V\n";

  std::vector<double> ha, hg, dPeak, dTrig, dCF;
  std::vector<int> outliers;
  for (auto& [ev, a] : alpha) {
    auto it = gamma.find(ev); if (it == gamma.end()) continue;
    const Result& g = it->second;
    bool cfOK = a.okCF && g.okCF;
    double dp = (g.tPeak - a.tPeak) * 1e9;
    o << ev << "," << a.height << "," << a.tPeak*1e9 << "," << a.tTrig*1e9 << ","
      << (a.okCF ? a.tCF*1e9 : -999) << ","
      << g.height << "," << g.tPeak*1e9 << "," << g.tTrig*1e9 << ","
      << (g.okCF ? g.tCF*1e9 : -999) << "," << g.sat << "\n";
    ha.push_back(a.height); hg.push_back(g.height);
    dPeak.push_back(dp);
    dTrig.push_back((g.tTrig - a.tTrig) * 1e9);
    if (cfOK) dCF.push_back((g.tCF - a.tCF) * 1e9);
    if (std::fabs(dp) > OUTLIER_NS) {
      outliers.push_back(ev);
      oo << ev << "," << dp << "," << a.tPeak*1e9 << "," << g.tPeak*1e9 << ","
         << a.height << "," << g.height << "\n";
    }
  }
  o.close(); oo.close();
  std::cout << "Wrote coincidence.csv (" << ha.size() << " matched events)\n";
  std::cout << "Wrote outliers.csv (" << outliers.size() << " events with |dt_peak| > "
            << OUTLIER_NS << " ns)\n";

#ifndef NO_ROOT
  Plot(ha, "Alpha pulse height (V)", "alpha_height", kBlue + 1);
  Plot(hg, "Gamma pulse height (V)", "gamma_height", kRed + 1);
  PlotDeltaT(dPeak, "t_{peak}(#gamma) #minus t_{peak}(#alpha) (ns)", "delta_t_peak",    kBlack);
  PlotDeltaT(dTrig, "t_{trig}(#gamma) #minus t_{trig}(#alpha) (ns)", "delta_t_trigger", kGreen + 2);
  PlotDeltaT(dCF,   "t_{CF}(#gamma) #minus t_{CF}(#alpha) (ns)",     "delta_t_cf",      kMagenta + 1);

  // waveforms for manual inspection: outliers (capped) + user-chosen events
  std::vector<int> toView(outliers.begin(),
        outliers.begin() + std::min<size_t>(outliers.size(), MAX_EVENT_PLOTS));
  toView.insert(toView.end(), EVENTS_TO_VIEW.begin(), EVENTS_TO_VIEW.end());
  for (int ev : toView) {
    auto ia = alpha.find(ev), ig = gamma.find(ev);
    if (ia != alpha.end() && ig != gamma.end()) PlotEvent(ia->second, ig->second);
    else std::cout << "Event " << ev << " not in both channels (removed by cuts?)\n";
  }
#endif
}

#ifdef NO_ROOT
int main() { analyse_coincidence(); }
#endif
