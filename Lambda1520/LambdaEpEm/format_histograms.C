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

void format_cb(TH1 * h) {
  h->SetLineColor(0);
  h->SetLineStyle(3);
  h->SetMarkerColor(kRed);
  h->SetMarkerStyle(8);
  h->SetMarkerSize(2.4);
  h->GetXaxis()->SetRange(1,13);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.035);
  h->GetXaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleFont(42);
  h->GetZaxis()->SetLabelFont(42);
  h->GetZaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetTitleSize(0.035);
  h->GetZaxis()->SetTitleFont(42);
}

void format_l1520(TH1 * h) {
  //h->SetTitle("#Lambda(1520) #rightarrow #Lambda e^{+}e^{-}");
  //    h->SetLineStyle(9);
  h->SetLineColor(46);
  h->SetLineWidth(2);
  //h->GetXaxis()->SetRange(3,49);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.035);
  h->GetXaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleFont(42);
  h->GetZaxis()->SetLabelFont(42);
  h->GetZaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetTitleSize(0.035);
  h->GetZaxis()->SetTitleFont(42);
}

void format_l1405(TH1 * h) {
  //h->SetTitle("#Lambda(1405) #rightarrow #Lambda e^{+}e^{-}");
  h->SetLineColor(30);
  //    h->SetLineStyle(2);
  h->SetLineWidth(2);
  //h->GetXaxis()->SetRange(3,49);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.035);
  h->GetXaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleFont(42);
  h->GetZaxis()->SetLabelFont(42);
  h->GetZaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetTitleSize(0.035);
  h->GetZaxis()->SetTitleFont(42);
}

void format_s1385(TH1 * h) {
  //h->SetTitle("#Sigma(1385) #rightarrow #Lambda e^{+}e^{-}");
  h->SetLineColor(38);
  //    h->SetLineStyle(10);
  h->SetLineWidth(2);
  //h->GetXaxis()->SetRange(3,49);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.035);
  h->GetXaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleFont(42);
  h->GetZaxis()->SetLabelFont(42);
  h->GetZaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetTitleSize(0.035);
  h->GetZaxis()->SetTitleFont(42);
}

/* For legend use code below */
/*
  TLegend *leg = new TLegend(0.17,0.6,0.6,0.95,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.027);
  TLegendEntry *entry=leg->AddEntry(h,"Y #rightarrow #Lambda e^{+}e^{-}","lpf");
  entry->SetFillStyle(1001);
  entry->SetTextFont(42);
  entry=leg->AddEntry(h,"CB","lpf");
  entry->SetTextFont(42);
  entry=leg->AddEntry(h,"#Lambda(1520) #rightarrow #Lambda e^{+}e^{-}","lpf");
  entry->SetTextFont(42);
  entry=leg->AddEntry(h,"#Lambda(1405) #rightarrow #Lambda e^{+}e^{-}","lpf");
  entry->SetFillStyle(1001);
  entry->SetTextFont(42);
  entry=leg->AddEntry(h,"#Sigma(1385) #rightarrow #Lambda e^{+}e^{-}","lpf");
  entry->SetTextFont(42);
  leg->Draw();
*/

/* Text on canvas for e+e-  */
/*
  tex = new TLatex();
  tex->SetTextSize(0.03);
  tex->SetNDC();
  tex->DrawLatex(0.16, 0.96,"HADES preliminary");

  tex->DrawLatex(0.75, 0.35,"#Lambda(1520) #rightarrow #Lambda e^{+}e^{-}");
  tex->DrawLatex(0.18, 0.29,"#Lambda(1405) #rightarrow #Lambda e^{+}e^{-}");
  tex->DrawLatex(0.18, 0.39,"#Sigma(1385) #rightarrow #Lambda e^{+}e^{-}");

  TArrow *ar1 = 0;
  ar1 = new TArrow(1284,17,1390,2.5, 0.01, "|>");
  ar1->Draw();
  ar1 = new TArrow(1283,26,1385,18.5, 0.01, "|>");
  ar1->Draw();
  ar1 = new TArrow(1565,22.5,1480,14.5, 0.01, "|>");
  ar1->Draw();
*/

/* Text on canvas for Ye+e- */
/*
  tex = new TLatex();
  tex->SetTextSize(0.03);
  tex->SetNDC();
  tex->DrawLatex(0.11, 0.96,"HADES preliminary");

  tex->DrawLatex(0.20, 0.88,"#pi^{0}");
  tex->DrawLatex(0.15, 0.50,"#Lambda(1520) #rightarrow #Lambda e^{+}e^{-}");
  tex->DrawLatex(0.40, 0.40,"#Sigma(1385) #rightarrow #Lambda e^{+}e^{-}");
  tex->DrawLatex(0.40, 0.26,"#Lambda(1405) #rightarrow #Lambda e^{+}e^{-}");

  TArrow *ar1 = 0;
  ar1 = new TArrow(190,0.35,220, 0.88, 0.01, "|>");
  ar1->Draw();
  ar1 = new TArrow(190,3.2,220, 6.9, 0.01, "|>");
  ar1->Draw();
  ar1 = new TArrow(60,13.9, 90, 36.0, 0.01, "|>");
  ar1->Draw();
*/
