// ============================================================================
//  Final analysis: walk-corrected Delta t of the correlated events,
//  exponential fit of the decay tail -> publication plots
//
//  Usage :  root -l -b -q 'lifetime.C+'
//
//  Reads nice_correlated.csv. For each timing method (peak, fixed trigger
//  crossing, constant fraction):
//    1) dt = t(gamma) - t(alpha)
//    2) walk correction: dt = c + f(V_gamma) + h(V_alpha) (polynomials, linear
//       least squares); the amplitude-dependent part is subtracted, the mean
//       position of dt is kept.
//    3) plots (PDF = vector graphics for the report, PNG for quick viewing):
//         <out>walk_<method>.pdf/png          dt vs heights, before / after
//         <out>delta_t_<method>_fit.pdf/png   dt histogram + exponential fit
//    <out> = "final_walkcorr_" (correction on) or "final_" (correction off)
//  Fit results: <out>fit_results.txt
//  Method order in all per-method tables below:  peak, trigger, cf
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
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TBox.h"
#include "TLine.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TAxis.h"

// ----------------------------- USER SETTINGS --------------------------------
const std::string NICE_FILE = "nice_correlated.csv";

// pulse-height columns (V) in nice_correlated.csv
const std::string COL_AH = "alpha_height_V";
const std::string COL_GH = "gamma_height_V";

// which methods to plot (peak is only kept for comparison)
const bool PLOT_PEAK = true, PLOT_TRIGGER = true, PLOT_CF = true;

// ---- walk correction ----
const bool   WALK_CORRECTION = true;
const int    WALK_DEG_GAMMA  = 2;      // polynomial degree in gamma height (0 : none)
const int    WALK_DEG_ALPHA  = 2;      // polynomial degree in alpha height (0 : none)
// only events with dt inside this window (ns) are used to determine the walk
const double WALK_FIT_MIN = -1000, WALK_FIT_MAX = 2000;
const int    WALK_PROFILE_BINS = 25;   // bins of the red mean-dt markers

// ---- plot ranges (ns), per method:        peak     trigger     cf
const double WALK_PLOT_MIN[3] = {          -1000,      -1000,   -1000 };   // walk plots, dt axis
const double WALK_PLOT_MAX[3] = {           2000,       2000,    2000 };
const double DT_PLOT_MIN[3]   = {          -1000,       -500,    -500 };   // dt histograms
const double DT_PLOT_MAX[3]   = {           2000,       1500,    1500 };
const double DT_BIN_WIDTH     = 25;    // ns

// ---- exponential fit: start (ns) per method, common end ----
const double FIT_START[3] = {  500,      500,      500 };                  // peak, trigger, cf
const double FIT_END      = 0;         // <= start : fit up to the end of the plot range
const double TAU_GUESS_NS = 100;       // starting value for the fit

// ---- plot style ----
const bool LOGY       = true;          // log y-axis on the dt plots
const bool ERROR_BARS = true;          // Poisson error bars on the histogram
// ----------------------------------------------------------------------------

const std::string OUT = WALK_CORRECTION ? "final_walkcorr_" : "final_";

struct Method {
  int idx; std::string name, label, colA, colG, xtitle;
  int color; bool sentinel;                    // sentinel: -999 marks "no time found"
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
  gStyle->SetPadLeftMargin(0.14); gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.05); gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  gStyle->SetHistLineWidth(2); gStyle->SetEndErrorSize(0);
  gStyle->SetPalette(kViridis); gStyle->SetNumberContours(255);
}

// ------------------------------ walk correction -----------------------------
struct Walk {
  int dg = 0, da = 0;
  double mg = 0, sg = 1, ma = 0, sa = 1;      // standardisation of the heights
  std::vector<double> c;                       // c[0] constant, then gamma terms, then alpha terms
  double meanAmp = 0;                          // mean of the amplitude-dependent part
};

std::vector<double> Basis(const Walk& w, double vg, double va) {
  std::vector<double> b(1, 1.0);
  double x = (vg - w.mg) / w.sg, p = 1;
  for (int k = 1; k <= w.dg; ++k) { p *= x; b.push_back(p); }
  x = (va - w.ma) / w.sa; p = 1;
  for (int k = 1; k <= w.da; ++k) { p *= x; b.push_back(p); }
  return b;
}

double Amp(const Walk& w, double vg, double va) {      // amplitude-dependent part of dt
  auto b = Basis(w, vg, va);
  double s = 0;
  for (size_t i = 1; i < b.size(); ++i) s += w.c[i] * b[i];
  return s;
}

bool Solve(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double>& x) {
  int n = b.size();
  for (int i = 0; i < n; ++i) {
    int p = i;
    for (int r = i + 1; r < n; ++r) if (std::fabs(A[r][i]) > std::fabs(A[p][i])) p = r;
    if (std::fabs(A[p][i]) < 1e-12) return false;
    std::swap(A[i], A[p]); std::swap(b[i], b[p]);
    for (int r = i + 1; r < n; ++r) {
      double f = A[r][i] / A[i][i];
      for (int c = i; c < n; ++c) A[r][c] -= f * A[i][c];
      b[r] -= f * b[i];
    }
  }
  x.assign(n, 0);
  for (int i = n - 1; i >= 0; --i) {
    double s = b[i];
    for (int c = i + 1; c < n; ++c) s -= A[i][c] * x[c];
    x[i] = s / A[i][i];
  }
  return true;
}

bool FitWalk(const std::vector<double>& vg, const std::vector<double>& va,
             const std::vector<double>& dt, Walk& w) {
  w.dg = WALK_DEG_GAMMA; w.da = WALK_DEG_ALPHA;
  std::vector<size_t> sel;
  for (size_t i = 0; i < dt.size(); ++i) if (dt[i] >= WALK_FIT_MIN && dt[i] <= WALK_FIT_MAX) sel.push_back(i);
  int n = 1 + w.dg + w.da;
  if ((int)sel.size() < 10 * n) return false;
  double sg = 0, sg2 = 0, sa = 0, sa2 = 0;
  for (size_t i : sel) { sg += vg[i]; sg2 += vg[i] * vg[i]; sa += va[i]; sa2 += va[i] * va[i]; }
  double N = sel.size();
  w.mg = sg / N; w.ma = sa / N;
  w.sg = std::sqrt(std::max(0.0, sg2 / N - w.mg * w.mg)); if (w.sg <= 0) w.sg = 1;
  w.sa = std::sqrt(std::max(0.0, sa2 / N - w.ma * w.ma)); if (w.sa <= 0) w.sa = 1;
  std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
  std::vector<double> y(n, 0.0);
  for (size_t i : sel) {
    auto b = Basis(w, vg[i], va[i]);
    for (int r = 0; r < n; ++r) { y[r] += b[r] * dt[i]; for (int c = 0; c < n; ++c) A[r][c] += b[r] * b[c]; }
  }
  if (!Solve(A, y, w.c)) return false;
  w.meanAmp = 0;
  for (size_t i : sel) w.meanAmp += Amp(w, vg[i], va[i]);
  w.meanAmp /= N;
  return true;
}

double RmsInWindow(const std::vector<double>& x) {
  double s = 0, s2 = 0; int n = 0;
  for (double v : x) if (v >= WALK_FIT_MIN && v <= WALK_FIT_MAX) { s += v; s2 += v * v; ++n; }
  if (n < 2) return 0;
  double m = s / n; return std::sqrt(std::max(0.0, s2 / n - m * m));
}

// ------------------------------ walk plot -----------------------------------
// dt vs gamma / alpha height, before (left) and after (right) the correction
void WalkPlot(const Method& m, const std::vector<double>& vg, const std::vector<double>& va,
              const std::vector<double>& before, const std::vector<double>& after) {
  const double ylo = WALK_PLOT_MIN[m.idx], yhi = WALK_PLOT_MAX[m.idx];
  TCanvas c("c_walk", "", 1400, 1100);
  c.Divide(2, 2);
  int pad = 1;
  for (int row = 0; row < 2; ++row) {
    const std::vector<double>& v = (row == 0) ? vg : va;
    double lo = *std::min_element(v.begin(), v.end()), hi = *std::max_element(v.begin(), v.end());
    double p = 0.02 * (hi - lo > 0 ? hi - lo : 1); lo -= p; hi += p;
    for (int col = 0; col < 2; ++col, ++pad) {
      const std::vector<double>& y = (col == 0) ? before : after;
      c.cd(pad);
      gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.17);
      gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.09);
      auto* h2 = new TH2D(Form("w2_%s_%d", m.name.c_str(), pad), "", 60, lo, hi, 60, ylo, yhi);
      auto* pf = new TProfile(Form("wp_%s_%d", m.name.c_str(), pad), "", WALK_PROFILE_BINS, lo, hi, ylo, yhi);
      for (size_t i = 0; i < v.size(); ++i) { h2->Fill(v[i], y[i]); pf->Fill(v[i], y[i]); }
      h2->GetXaxis()->SetTitle(row == 0 ? "Gamma pulse height (V)" : "Alpha pulse height (V)");
      h2->GetYaxis()->SetTitle(Form("#Deltat_{%s} (ns)", m.label.c_str()));
      h2->GetZaxis()->SetTitle("Counts / bin");
      h2->SetMinimum(0.5);                      // empty bins stay white
      h2->Draw("COLZ");
      pf->SetMarkerStyle(20); pf->SetMarkerSize(0.9);
      pf->SetMarkerColor(kRed + 1); pf->SetLineColor(kRed + 1); pf->SetLineWidth(2);
      pf->Draw("SAME");
      TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.055);
      tx.DrawLatex(0.14, 0.935, col == 0 ? "Before walk correction" : "After walk correction");
    }
  }
  c.SaveAs((OUT + "walk_" + m.name + ".pdf").c_str());
  c.SaveAs((OUT + "walk_" + m.name + ".png").c_str());
}

// ------------------------- delta t histogram + exponential fit ---------------
void PlotMethod(const std::vector<double>& x, const Method& m, bool walkCorrected, std::ostream& res) {
  const std::string file = OUT + "delta_t_" + m.name + "_fit";
  const double lo = DT_PLOT_MIN[m.idx];
  const int nb = std::max(1, (int)std::lround((DT_PLOT_MAX[m.idx] - lo) / DT_BIN_WIDTH));
  const double hi = lo + nb * DT_BIN_WIDTH;
  TH1D h(file.c_str(), "", nb, lo, hi);
  int nOut = 0;
  for (double v : x) { if (v < lo || v >= hi) ++nOut; h.Fill(v); }
  std::cout << m.name << ": " << x.size() << " events, " << nOut
            << " outside [" << lo << ", " << hi << "] ns (not shown)\n";

  h.GetXaxis()->SetTitle(m.xtitle.c_str());
  h.GetYaxis()->SetTitle(Form("Counts / %.0f ns", DT_BIN_WIDTH));
  h.SetLineColor(m.color); h.SetFillColorAlpha(m.color, 0.30);
  h.SetMarkerStyle(20); h.SetMarkerSize(0.6); h.SetMarkerColor(m.color);
  const double hmax = h.GetMaximum();
  const double ymin = LOGY ? 0.5 : 0.0, ymax = LOGY ? 8 * hmax : 1.35 * hmax;
  h.SetMinimum(ymin); h.SetMaximum(ymax);

  // ---- exponential fit: N(dt) = A exp(-(dt - start)/tau) ----
  const double fs = FIT_START[m.idx];
  const double fe = (FIT_END > fs) ? std::min(FIT_END, hi) : hi;
  TF1 f("fexp", Form("[0]*exp(-(x-%.6f)/[1])", fs), fs, fe);
  f.SetParNames("A", "tau");
  double a0 = h.GetBinContent(h.FindBin(fs));
  if (a0 <= 0) a0 = std::max(1.0, 0.5 * hmax);
  f.SetParameters(a0, TAU_GUESS_NS);
  f.SetParLimits(0, 1e-3, 1e9);
  f.SetParLimits(1, 1, 1e5);
  int status = h.Fit(&f, "RLIQ0");            // Poisson likelihood, integral over each bin
  double tau = f.GetParameter(1), dtau = f.GetParError(1);
  const double ln2 = std::log(2.0);
  std::cout << "  fit [" << fs << ", " << fe << "] ns:  tau = " << tau << " +- " << dtau
            << " ns,  T1/2 = " << tau * ln2 << " +- " << dtau * ln2 << " ns   (status " << status << ")\n";
  res << m.name << ": fit range " << fs << " - " << fe << " ns, tau = " << tau << " +- " << dtau
      << " ns, T1/2 = " << tau * ln2 << " +- " << dtau * ln2 << " ns, events " << x.size()
      << ", fit status " << status << "\n";

  // ---- drawing ----
  const int grey = TColor::GetColor("#9a9a9a");
  TCanvas c("c_dt", "", 800, 600);
  c.SetLogy(LOGY);
  h.Draw("AXIS");
  TBox shade(fs, ymin, fe, ymax);              // lightly shaded fitted region
  shade.SetFillColorAlpha(grey, 0.15); shade.SetLineColor(0);
  shade.DrawClone();
  TLine start(fs, ymin, fs, ymax);             // light line marking the fit start
  start.SetLineColor(grey); start.SetLineWidth(2);
  start.DrawClone();
  h.Draw("HIST SAME");
  if (ERROR_BARS) h.Draw("E SAME");
  f.SetLineColor(kBlack); f.SetLineWidth(3);
  f.DrawCopy("SAME");
  gPad->RedrawAxis();

  TLatex tx; tx.SetNDC(); tx.SetTextFont(42);
  tx.SetTextSize(0.05);
  tx.DrawLatex(0.18, 0.86, m.label.c_str());
  tx.SetTextSize(0.040);
  if (walkCorrected) tx.DrawLatex(0.18, 0.80, "walk corrected");
  tx.SetTextSize(0.045);
  tx.DrawLatex(0.60, 0.86, Form("#tau = %.0f #pm %.0f ns", tau, dtau));
  tx.DrawLatex(0.60, 0.79, Form("T_{1/2} = %.0f #pm %.0f ns", tau * ln2, dtau * ln2));
  c.SaveAs((file + ".pdf").c_str());
  c.SaveAs((file + ".png").c_str());
}

// -------------------------------- main --------------------------------------
void lifetime() {
  std::ifstream f(NICE_FILE);
  if (!f) { std::cerr << "Cannot open " << NICE_FILE << "\n"; return; }
  std::string line;
  std::getline(f, line);
  std::vector<std::string> hdr = Split(line);
  std::map<std::string, int> idx;
  for (size_t i = 0; i < hdr.size(); ++i) idx[hdr[i]] = i;
  std::vector<std::vector<std::string>> rows;
  while (std::getline(f, line)) {
    auto fl = Split(line);
    if (fl.size() == hdr.size()) rows.push_back(fl);
  }
  std::cout << "Read " << rows.size() << " correlated events from " << NICE_FILE << "\n";
  if (rows.empty()) return;

  const bool haveHeights = idx.count(COL_AH) && idx.count(COL_GH);
  bool doWalk = WALK_CORRECTION;
  if (doWalk && !haveHeights) {
    std::cerr << "Columns " << COL_AH << " / " << COL_GH << " not found: walk correction off\n";
    doWalk = false;
  }

  StyleSetup();
  // colour-blind-safe palette (Okabe-Ito)
  std::vector<Method> methods = {
    {0, "peak",    "Peak",             "alpha_peak_ns", "gamma_peak_ns", "t_{peak}(#gamma) #minus t_{peak}(#alpha) (ns)", TColor::GetColor("#009E73"), false},
    {1, "trigger", "Fixed threshold",  "alpha_trig_ns", "gamma_trig_ns", "t_{trig}(#gamma) #minus t_{trig}(#alpha) (ns)", TColor::GetColor("#D55E00"), true},
    {2, "cf",      "Constant fraction","alpha_cf_ns",   "gamma_cf_ns",   "t_{CF}(#gamma) #minus t_{CF}(#alpha) (ns)",     TColor::GetColor("#0072B2"), true}};
  const bool plotIt[3] = {PLOT_PEAK, PLOT_TRIGGER, PLOT_CF};

  std::ofstream res(OUT + "fit_results.txt");
  for (auto& m : methods) {
    if (!plotIt[m.idx]) continue;
    if (!idx.count(m.colA) || !idx.count(m.colG)) {
      std::cerr << "Skipping " << m.name << ": columns " << m.colA << " / " << m.colG << " not found\n";
      continue;
    }
    int ia = idx[m.colA], ig = idx[m.colG];
    std::vector<double> dt, vg, va;
    for (auto& r : rows) {
      double a = std::atof(r[ia].c_str()), g = std::atof(r[ig].c_str());
      if (m.sentinel && (std::fabs(a + 999) < 1e-6 || std::fabs(g + 999) < 1e-6)) continue;
      dt.push_back(g - a);
      if (haveHeights) { vg.push_back(std::atof(r[idx[COL_GH]].c_str())); va.push_back(std::atof(r[idx[COL_AH]].c_str())); }
    }
    if (dt.empty()) { std::cerr << "No valid events for " << m.name << "\n"; continue; }

    std::vector<double> dtUse = dt;
    bool corrected = false;
    if (doWalk) {
      Walk w;
      if (FitWalk(vg, va, dt, w)) {
        for (size_t i = 0; i < dt.size(); ++i) dtUse[i] = dt[i] - (Amp(w, vg[i], va[i]) - w.meanAmp);
        corrected = true;
        std::cout << m.name << ": walk correction (gamma deg " << w.dg << ", alpha deg " << w.da
                  << ");  dt RMS in [" << WALK_FIT_MIN << ", " << WALK_FIT_MAX << "] ns: "
                  << RmsInWindow(dt) << " -> " << RmsInWindow(dtUse) << " ns\n";
        res << m.name << ": walk correction, dt RMS " << RmsInWindow(dt) << " -> " << RmsInWindow(dtUse) << " ns\n";
        WalkPlot(m, vg, va, dt, dtUse);
      } else {
        std::cerr << m.name << ": walk fit failed (too few events in window?), no correction\n";
      }
    }
    PlotMethod(dtUse, m, corrected, res);
  }
}
