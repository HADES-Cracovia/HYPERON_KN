
void format_epem_with_qcdcorr(TH1 * h) {
  h->SetLineColor(14);
  h->SetLineStyle(9);
  h->SetLineWidth(4);
  h->GetXaxis()->SetTitle("M_{#Lambda e^{+}e^{-}}[MeV c^{-2}]");
  h->GetXaxis()->SetRange(3,49);
  h->GetXaxis()->SetNdivisions(505);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.93);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitle("d#sigma/dM_{e^{+}e^{-}} [#mub/MeV c^{-2}]");
  h->GetYaxis()->SetRange(1,1);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.007);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.29);
  h->GetYaxis()->SetTitleFont(42);
  h->GetZaxis()->SetLabelFont(42);
  h->GetZaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetTitleSize(0.035);
  h->GetZaxis()->SetTitleFont(42);
}

void format_epem(TH1 * h) {
  h->SetTitle("Y #rightarrow #Lambda e^{+}e^{-}");
  h->SetLineColor(1);
  h->SetLineWidth(4);
  //    h->SetLineStyle(9);
  h->GetXaxis()->SetTitle("M_{#Lambdae^{+}e^{-}} [MeV c^{-2}]");
  h->GetXaxis()->SetRange(3,49);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetNdivisions(505);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitle("d#sigma/dM_{e^{+}e^{-}} [#mub/MeV c^{-2}]");
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleFont(42);
  h->GetZaxis()->SetLabelFont(42);
  h->GetZaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetTitleSize(0.035);
  h->GetZaxis()->SetTitleFont(42);
}

void format_cb(TH1 * h)
{
  h->SetLineWidth(3);
  h->SetLineColor(kRed);
  h->SetMarkerColor(kRed);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.5);
//h->GetXaxis()->SetRange(1,13);
  //h->GetXaxis()->SetLabelFont(42);
  //h->GetXaxis()->SetLabelSize(0.035);
  //h->GetXaxis()->SetTitleSize(0.035);
  //h->GetXaxis()->SetTitleFont(42);
  //h->GetYaxis()->SetLabelFont(42);
  //h->GetYaxis()->SetLabelSize(0.035);
  //h->GetYaxis()->SetTitleSize(0.035);
  //h->GetYaxis()->SetTitleFont(42);
  //h->GetZaxis()->SetLabelFont(42);
  //h->GetZaxis()->SetLabelSize(0.035);
  //h->GetZaxis()->SetTitleSize(0.035);
  //h->GetZaxis()->SetTitleFont(42);
}

void format_l1520(TH1 * h) {
  //h->SetTitle("#Lambda(1520) #rightarrow #Lambda e^{+}e^{-}");
  h->SetLineStyle(1);
  h->SetLineColor(kYellow-6);
  h->SetLineWidth(3);
  //h->GetXaxis()->SetRange(3,49);
  //h->GetXaxis()->SetLabelFont(42);
  //h->GetXaxis()->SetLabelSize(0.035);
  //h->GetXaxis()->SetTitleSize(0.035);
  //h->GetXaxis()->SetTitleFont(42);
  //h->GetYaxis()->SetLabelFont(42);
  //h->GetYaxis()->SetLabelSize(0.035);
  //h->GetYaxis()->SetTitleSize(0.035);
  //h->GetYaxis()->SetTitleFont(42);
  //h->GetZaxis()->SetLabelFont(42);
  //h->GetZaxis()->SetLabelSize(0.035);
  //h->GetZaxis()->SetTitleSize(0.035);
  //h->GetZaxis()->SetTitleFont(42);
}

void format_l1405(TH1 * h) {
  //h->SetTitle("#Lambda(1405) #rightarrow #Lambda e^{+}e^{-}");
  h->SetLineColor(30);
  h->SetLineStyle(1);
  h->SetLineWidth(3);
  //h->GetXaxis()->SetRange(3,49);
  //h->GetXaxis()->SetLabelFont(42);
  //h->GetXaxis()->SetLabelSize(0.035);
  //h->GetXaxis()->SetTitleSize(0.035);
  //h->GetXaxis()->SetTitleFont(42);
  //h->GetYaxis()->SetLabelFont(42);
  //h->GetYaxis()->SetLabelSize(0.035);
  //h->GetYaxis()->SetTitleSize(0.035);
  //h->GetYaxis()->SetTitleFont(42);
  //h->GetZaxis()->SetLabelFont(42);
  //h->GetZaxis()->SetLabelSize(0.035);
  //h->GetZaxis()->SetTitleSize(0.035);
  //h->GetZaxis()->SetTitleFont(42);
}

void format_s1385(TH1 * h) {
  //h->SetTitle("#Sigma(1385) #rightarrow #Lambda e^{+}e^{-}");
  h->SetLineColor(kGray+2);
  h->SetLineStyle(1);
  h->SetLineWidth(3);
  //h->GetXaxis()->SetRange(3,49);
  //h->GetXaxis()->SetLabelFont(42);
  //h->GetXaxis()->SetLabelSize(0.035);
  //h->GetXaxis()->SetTitleSize(0.035);
  //h->GetXaxis()->SetTitleFont(42);
  //h->GetYaxis()->SetLabelFont(42);
  //h->GetYaxis()->SetLabelSize(0.035);
  //h->GetYaxis()->SetTitleSize(0.035);
  //h->GetYaxis()->SetTitleFont(42);
  //h->GetZaxis()->SetLabelFont(42);
  //h->GetZaxis()->SetLabelSize(0.035);
  //h->GetZaxis()->SetTitleSize(0.035);
  //h->GetZaxis()->SetTitleFont(42);
}

void format_pi0(TH1 *h)
{
  h->SetLineColor(kCyan-6);
  h->SetLineStyle(1);
  h->SetLineWidth(3);
  h->SetMarkerColor(kCyan-6);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.5);

}

void format_delta(TH1 *h)
{
  h->SetLineColor(38);
  h->SetLineStyle(1);
  h->SetLineWidth(3);
  //h->SetMarkerColor(30);
  //h->SetMarkerStyle(8);
  //h->SetMarkerSize(1.4);

}


void format_sumSignals(TH1 *h)
{
  h->SetLineColor(kMagenta-3);
  h->SetLineWidth(3);
  h->SetFillStyle(3144);
}
