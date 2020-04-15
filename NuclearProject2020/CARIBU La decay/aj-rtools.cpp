#include <stdlib.h>
#include <fcntl.h>
#include <stddef.h>
#include <sys/stat.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <TSpectrum.h>
#include <Riostream.h>
#define DIMX 2048
#define DIMY 2048
#define NCHX 2048
#define NCHY 2048
#define MAXGE 110

TFile *f1;
TF1   *f2x=new TF1("f2x","[0]+[1]*x+[2]*x^2",0,5);
TF1   *f2y=new TF1("f2y","[0]+[1]*x+[2]*x^2",0,5);
TCanvas *c1;
TSpectrum *s1;
Axis_t axmin=0,aymin=0;
Float_t numchx=DIMX,numchy=DIMY;
Float_t fwx=3.0,gwx=1.0,hwx=0.0;
Float_t fwy=3.0,gwy=1.0,hwy=0.0;
Char_t rootfilename[255];
TH1D *histx;
TH1D *histy;
struct{
  Float_t swpars[3];
  Int_t maxch;
  Float_t spec[16384];
}gf3gd;

/*//FUNCTIUON LIST
// load a root file e.g. dload("file.root")
void dload(const Char_t *fn, const Char_t *options="READ")
// draw 2D histogram e.g. d2d("EhiCn")
void d2d(const Char_t *matn, Int_t n=1)
// do x projection and get background spectrum e.g. pjx("CEhiCln")
Int_t pjx(const Char_t *matn, Int_t w=45)
// do y projection and get background spectrum e.g. pjy("EvsTape")
Int_t pjy(const Char_t *matn, Int_t w=45)
// do x and y projection and get background spectrum e.g. pj2d("EvsTape")
Int_t pj2d(const Char_t *matn, Float_t w=1.0)
// do an x-projection of a 2D histogram and subtracts a backgroud with scaling factor 'w' i.e. dx("EhiCln", 10).
// do a y-projection of a 2D histogram and subtracts a backgroud with scaling factor 'w' i.e. dy("EhiCln", 10).
//set a gate on a gamma-gamma 2D matrix i.e. GX("GG",181,3,1)
Float_t GX(const Char_t *matn,Float_t peak,Float_t width=0.0,Int_t cy=0)
//set a gate on a gamma-gamma 2D matrix i.e. GY("GG",181,3,1)
Float_t GY(const Char_t *matn,Float_t peak,Float_t width=0.0, Int_t cy=0)
// write out 2D matrix e.g. write("gg")
void writemat(const Char_t *filename)
// write out 1D spectrum e.g. writespe("px","spectrum")
void writespe(const Char_t *hisname, Char_t *spename)
// delete a histogram e.g. del hist"GamE_bgx")
void delhis(const Char_t *hn)
// fit a peak in a spectrum
void gaus1peakfit(const Char_t *s, Float_t x1, Float_t x2, Float_t x3, Float_t x4)
// fit two peaks in a spectrum
void gaus2peakfit(const Char_t *s, Float_t x1, Float_t x2, Float_t x3, Float_t x4)
*/

/*=========================================================*/
/*=========================================================*/
// load a root file e.g. dload("file.root")
void dload(const Char_t *fn, const Char_t *options="READ")
{
  FILE *in;
  Double_t parx[3],pary[3];

  sprintf(rootfilename,"%s",fn); 
  in=fopen("fwhm.out","r"); 
  if(in==NULL){
    printf("fwhm.out is not found\n");
  }
  else{
    printf("reading width parameters\n");
    fscanf(in,"ffx=%f,ggx=%f,hhx=%f\n",&fwx,&gwx,&hwx);
    fscanf(in,"ffy=%f,ggy=%f,hhy=%f\n",&fwy,&gwy,&hwy);
  }
  parx[0]=fwx*fwx;
  parx[1]=gwx*gwx;
  parx[2]=hwx*hwx;
  pary[0]=fwy*fwy;
  pary[1]=gwy*gwy;
  pary[2]=hwy*hwy;

  f2x->SetParameters(parx);
  f2y->SetParameters(pary);
  f1 = new TFile(fn,options);
  //setcanvas(1);
  return; //dload end
}
/*=========================================================*/
void setcanvas(Int_t ny=1, Int_t nx=1)
{
  c1=(TCanvas*)gROOT->FindObject("c1");
  if(c1==NULL)
    c1=new TCanvas("c1","c1");
  else
    c1->Clear();
  c1->SetFillColor(0);
  c1->SetCrosshair(1);
  if(!c1->GetShowEventStatus())
    c1->ToggleEventStatus();
  c1->Divide(nx,ny);
  return; //setcanvas end
}
/*=========================================================*/
// draw 2D histogram e.g. d2d("EhiCn")
void d2d(const Char_t *matn, Int_t n=1)
{
  TH2F *hist;
  Axis_t axmax=axmin+numchx;
  Axis_t aymax=aymin+numchy;

  c1=(TCanvas *)gROOT->FindObject("c1");
  if(c1!=NULL)c1->cd(n);
  else setcanvas(1);
  gStyle->SetPalette(1);
  gPad->SetLogz(1);
  hist=(TH2F*)f1->Get(matn);
  hist->SetAxisRange(axmin,axmax,"X");
  hist->SetAxisRange(aymin,aymax,"Y");
  hist->Draw("COLZ");
  return; //d2d end
}
/*=========================================================*/
// do x projection and get background spectrum e.g. pjx("CEhiCln")
Int_t pjx(const Char_t *matn, Int_t w=45){
  TH2 *hist2;
  TH1 *hist1[2];
  Char_t str[255];
  s1=new TSpectrum();

  hist2 = (TH2 *) f1->Get(matn);

  sprintf(str,"%s_px",matn);
  hist1[0]=(TH1*)gROOT->FindObject(str);
  if(hist1[0]==NULL) hist1[0] = (TH1*)hist2->ProjectionX();
  sprintf(str,"%s_bgx",matn);
  hist1[1]=(TH1*)gROOT->FindObject(str);
  if(hist1[1]==NULL){
    hist1[1]=s1->Background(hist1[0],w);
    hist1[1]->SetName(str);
    hist1[1]->SetTitle(str);
  }
  return 0;
}
/*=========================================================*/
// do y projection and get background spectrum e.g. pjy("EvsTape")
Int_t pjy(const Char_t *matn, Int_t w=45){
  TH2 *hist2;
  TH1 *hist1[2];
  Char_t str[255];
  s1=new TSpectrum();

  hist2 = (TH2 *) f1->Get(matn);

  sprintf(str,"%s_py",matn);
  hist1[0]=(TH1*)gROOT->FindObject(str);
  if(hist1[0]==NULL) hist1[0] = (TH1*)hist2->ProjectionY();
  sprintf(str,"%s_bgy",matn);
  hist1[1]=(TH1*)gROOT->FindObject(str);
  if(hist1[1]==NULL){
    hist1[1]=s1->Background(hist1[0],w);
    hist1[1]->SetName(str);
    hist1[1]->SetTitle(str);
  }
  return 0;
}
/*=========================================================*/
// set up find_bg2 function for 'autobkgnd'
int find_bg2 (int loch, int hich, float *x1, float *y1, float *x2, float *y2)
{
  float a, y, bg;
  int i, mid;

  mid = (loch + hich) / 2;
  bg = 0.f;

  /* find channel with least counts on each side of middle */
  *y1 =
    (gf3gd.spec[loch] + gf3gd.spec[loch + 1] + gf3gd.spec[loch + 2]) / 3.f;
  *x1 = (float) (loch + 1);
  for(i = loch + 2; i <= mid - 2; ++i){
    a = (gf3gd.spec[i - 1] + gf3gd.spec[i] + gf3gd.spec[i + 1]) / 3.f;
    if(*y1 > a){
      *y1 = a;
      *x1 = (float) i;
    }
  }
  *y2 =
    (gf3gd.spec[mid + 2] + gf3gd.spec[mid + 3] + gf3gd.spec[mid + 4]) / 3.f;
  *x2 = (float) (mid + 3);
  for(i = mid + 4; i <= hich - 1; ++i){
    a = (gf3gd.spec[i - 1] + gf3gd.spec[i] + gf3gd.spec[i + 1]) / 3.f;
    if(*y2 > a){
      *y2 = a;
      *x2 = (float) i;
    }
  }
  /* check that there are no channels between loch and hich
     that are below the background. */
  if(*y2 > *y1){
    for(i = (int)ceil(*x1) + 1; i <= mid - 2; ++i){
      y = *y1 - (*x1 - (float) i) * (*y1 - *y2) / (*x1 - *x2);
      a = (gf3gd.spec[i - 1] + gf3gd.spec[i] + gf3gd.spec[i + 1]) / 3.f;
      if(y > a){
        *y1 = a;
        *x1 = (float) i;
      }
    }
    for(i = (int)ceil(*x2) + 1; i <= hich - 1; ++i){
      y = *y1 - (*x1 - (float) i) * (*y1 - *y2) / (*x1 - *x2);
      a = (gf3gd.spec[i - 1] + gf3gd.spec[i] + gf3gd.spec[i + 1]) / 3.f;
      if(y > a){
        *y2 = a;
        *x2 = (float) i;
      }
    }
  }
  else{
    for(i = (int)ceil(*x1) - 1; i >= loch + 1; --i){
      y = *y1 - (*x1 - (float) i) * (*y1 - *y2) / (*x1 - *x2);
      a = (gf3gd.spec[i - 1] + gf3gd.spec[i] + gf3gd.spec[i + 1]) / 3.f;
      if(y > a){
        *y1 = a;
        *x1 = (float) i;
      }
    }
    for(i = (int)ceil(*x2) - 1; i >= mid + 3; --i){
      y = *y1 - (*x1 - (float) i) * (*y1 - *y2) / (*x1 - *x2);
      a = (gf3gd.spec[i - 1] + gf3gd.spec[i] + gf3gd.spec[i + 1]) / 3.f;
      if(y > a){
        *y2 = a;
        *x2 = (float) i;
      }
    }
  }
  if(*y1 < 1.f)
    *y1 = 1.f;
  if(*y2 < 1.f)
    *y2 = 1.f;
  return 0;
}                               /* find_bg2_ */
/*=========================================================*/
// set up autobackground finder
Int_t autobkgnd(Float_t w1)
{
  static float w2 = 6.f;
  float fwhm, x, x1, y1, x2, y2, bg, r1;
  int i, hi, lo, ich, jch, ihi = 0, ilo = 0, ihi2, ilo2;
  ilo = 0;
  ihi = gf3gd.maxch;
  ilo2 = ihi;
  ihi2 = ilo;
  for(i = 1; i <= 3; ++i)
    for(ich = ilo + 2; ich <= ihi - 2; ich += 5){
      x = (float) ich;
      fwhm=sqrt(gf3gd.swpars[0]+
                gf3gd.swpars[1]*x +
                gf3gd.swpars[2]*x*x);
      r1 = w1 * w2 * fwhm;
      lo = ich-(int)floor(r1);
      r1 = w1 * w2 * fwhm;
      hi = ich + (int)floor(r1);
      if(lo < ilo)
        lo = ilo;
      if(hi > ihi)
        hi = ihi;
      find_bg2 (lo, hi, &x1, &y1, &x2, &y2);
      if(ilo2 > (int)floor(x1))
        ilo2 = (int)floor(x1);
      if(ihi2 < (int)floor(x2))
        ihi2 = (int)floor(x2);
      for(jch = (int)floor(x1); jch <= (int)floor(x2); ++jch){
        x = (float) jch;
        bg = y1 - (x1 - x) * (y1 - y2) / (x1 - x2);
        gf3gd.spec[jch] = bg;
      }
    }
  for(ich = ilo2; ich <= ihi2; ++ich){
    gf3gd.spec[ich] += sqrt (gf3gd.spec[ich]) * 1.6f;
  }
  return 0;
}                               /* autobkgnd */
/*=========================================================*/
// do x and y projection and get background spectrum e.g. pj2d("EvsTape")
Int_t pj2d(const Char_t *matn, Float_t w=1.0)
{
  TH2F *hist2;
  TH1 *hist1[4];
  Char_t str[255];
  Int_t i;
  hist2 = (TH2F *) f1->Get(matn);

  hist1[0] = hist2->ProjectionX();
  sprintf(str,"%s_bgx",matn);
  hist1[1] = new TH1D(str,matn,NCHX,0,DIMX);
  gf3gd.swpars[0]=fwx*fwx;
  gf3gd.swpars[1]=gwx*gwx/1e3;
  gf3gd.swpars[2]=hwx*hwx/1e6;
  gf3gd.maxch=NCHX;
  for(i=1;i<=NCHX;i++)
    gf3gd.spec[i-1]=hist1[0]->GetBinContent(i);
  autobkgnd(w);
  for(i=1;i<=NCHX;i++)
    hist1[1]->SetBinContent(i,gf3gd.spec[i-1]);

  hist1[2] = hist2->ProjectionY();
  sprintf(str,"%s_bgy",matn);
  hist1[3] = new TH1D(str,matn,NCHY,0,DIMY);
  gf3gd.swpars[0]=fwy*fwy;
  gf3gd.swpars[1]=gwy*gwy/1e3;
  gf3gd.swpars[2]=hwy*hwy/1e6;
  gf3gd.maxch=NCHY;
  for(i=1;i<=NCHY;i++)
    gf3gd.spec[i-1]=hist1[2]->GetBinContent(i);
  autobkgnd(w);
  for(i=1;i<=NCHY;i++)
    hist1[3]->SetBinContent(i,gf3gd.spec[i-1]);

  return 0;
}
/*=========================================================*/
// do an x-projection of a 2D histogram and subtracts a backgroud with scaling factor 'w' i.e. dx("EhiCln", 10).
void dx(const Char_t *matn, Int_t w)
{
  TH1D *hist[3];
  TH2F *h2;
  Axis_t axmax=axmin+numchx;
  Char_t str[255];
  Int_t NX,NY;
  Double_t XMIN,XMAX;
  Double_t YMIN,YMAX;
  s1=new TSpectrum();

  h2 = (TH2F*)gROOT->FindObject(matn);
  NX = h2->GetNbinsX();
  XMIN = h2->GetXaxis()->GetXmin();
  XMAX = h2->GetXaxis()->GetXmax();
  NY = h2->GetNbinsY();
  YMIN = h2->GetYaxis()->GetXmin();
  YMAX = h2->GetYaxis()->GetXmax();
  h2->SetAxisRange(XMIN,XMAX-1,"X");
  h2->SetAxisRange(YMIN,YMAX-1,"Y");

  sprintf(str,"%sbx",matn);
  hist[2]=(TH1D*)gROOT->FindObject(str);
  if(hist[2]!=NULL)hist[2]->Delete();
  hist[2]=new TH1D(str,str,NX,XMIN,XMAX);

  sprintf(str,"%s_px",matn);
  hist[0]=(TH1D*)gROOT->FindObject(str);
  if(hist[0]==NULL) hist[0] = (TH1D*)h2->ProjectionX();
  //(Int_t)N = hist[0]->GetNbinsX();
  //printf("N=%i\n",N);
  sprintf(str,"%s_bgx",matn);
  hist[1]=(TH1D*)gROOT->FindObject(str);
  if(hist[1]!=NULL)hist[1]->Delete();
  hist[1]=(TH1D*)s1->Background(hist[0],w);
  hist[1]->SetName(str);
  hist[1]->SetTitle(str);

  setcanvas(2);
  c1->cd(1);
  hist[0]->Draw();
  hist[1]->SetLineColor(3);
  hist[1]->Draw("same");
  c1->cd(2);
  hist[2]->Add(hist[0],hist[1],1.0,-1.0);
  hist[2]->SetLineColor(2);
  hist[2]->SetAxisRange(axmin+1,axmax-1);
  hist[2]->Draw();
  return; //end DX
}
/*=========================================================*/
// do a y-projection of a 2D histogram and subtracts a backgroud with scaling factor 'w' i.e. dy("EhiCln", 10).
void dy(const Char_t *matn, Int_t w)
{
  TH1D *hist[3];
  TH2F *h2;
  Axis_t axmax=axmin+numchx;
  Char_t str[255];
  Int_t NX,NY;
  Double_t XMIN,XMAX;
  Double_t YMIN,YMAX;
  s1=new TSpectrum();

  h2 = (TH2F*)gROOT->FindObject(matn);
  NX = h2->GetNbinsX();
  XMIN = h2->GetXaxis()->GetXmin();
  XMAX = h2->GetXaxis()->GetXmax();
  NY = h2->GetNbinsY();
  YMIN = h2->GetYaxis()->GetXmin();
  YMAX = h2->GetYaxis()->GetXmax();
  h2->SetAxisRange(XMIN,XMAX-1,"X");
  h2->SetAxisRange(YMIN,YMAX-1,"Y");

  sprintf(str,"%sby",matn);
  hist[2]=(TH1D*)gROOT->FindObject(str);
  if(hist[2]!=NULL)hist[2]->Delete();
  hist[2]=new TH1D(str,str,NY,YMIN,YMAX);

  sprintf(str,"%s_py",matn);
  hist[0]=(TH1D*)gROOT->FindObject(str);
  if(hist[0]==NULL) hist[0] = (TH1D*)h2->ProjectionY();

  sprintf(str,"%s_bgy",matn);
  hist[1]=(TH1D*)gROOT->FindObject(str);
  if(hist[1]!=NULL)hist[1]->Delete();
  hist[1]=(TH1D*)s1->Background(hist[0],w);
  hist[1]->SetName(str);
  hist[1]->SetTitle(str);

  setcanvas(2);
  c1->cd(1);
  hist[0]->Draw();
  hist[1]->SetLineColor(3);
  hist[1]->Draw("same");
  c1->cd(2);
  hist[2]->Add(hist[0],hist[1],1.0,-1.0);
  hist[2]->SetLineColor(2);
  hist[2]->SetAxisRange(axmin+1,axmax-1);
  hist[2]->Draw();
  return; //end DY
}
/*=========================================================*/
//set a gate on a gamma-gamma 2D matrix i.e. GX("GG",181,3,1)
Float_t GX(const Char_t *matn,Float_t peak,Float_t width=0.0,Int_t cy=0)
{
  TH2 *hist2;
  TH1 *hist[5];
  Char_t str[255];
  Int_t nxmin,nxmax;
  Float_t ratio,a,b;
  Float_t xmin,xmax;
  Axis_t aymax=aymin+numchy;
  s1=new TSpectrum();
  Int_t NX,NY;
  Double_t XMIN,XMAX;
  Double_t YMIN,YMAX;

  hist2 = (TH2*)f1->Get(matn);
  if(hist2==NULL){
    printf("Couldn't find %s\n",matn);
    return -1.0;
  }
  NX = hist2->GetNbinsX();
  XMIN = hist2->GetXaxis()->GetXmin();
  XMAX = hist2->GetXaxis()->GetXmax();
  NY = hist2->GetNbinsY();
  YMIN = hist2->GetYaxis()->GetXmin();
  YMAX = hist2->GetYaxis()->GetXmax();
  hist2->SetAxisRange(XMIN,XMAX-1,"X");
  hist2->SetAxisRange(YMIN,YMAX-1,"Y");

  sprintf(str,"%s_px",matn);
  hist[0]=(TH1*)gROOT->FindObject(str);
  sprintf(str,"%s_py",matn);
  hist[1]=(TH1*)gROOT->FindObject(str);
  sprintf(str,"%s_bgx",matn);
  hist[2]=(TH1*)gROOT->FindObject(str);

  if(hist[0]==NULL || hist[1]==NULL || hist[2]==NULL){
    dx(matn,20);
    dy(matn,20);
    return 0;
  }
  hist[0]->SetAxisRange(XMIN,XMAX);
  hist[1]->SetAxisRange(YMIN,YMAX);
  hist[2]->SetAxisRange(XMIN,XMAX);

  if(width==0.0){
    width = (Float_t)f2x->Eval(peak/1000.);
    width=sqrtf(width);
  }
  printf("width=%f\n",width);
  xmin=peak-width;
  xmax=peak+width;
  nxmin=hist[0]->FindBin(xmin);
  nxmax=hist[0]->FindBin(xmax);
  a=hist[2]->Integral(nxmin,nxmax);
  b=hist[0]->Integral();

  sprintf(str,"%sg%4.4iy",matn,(Int_t)peak);
  hist[3]=(TH1*)gROOT->FindObject(str);
  if(hist[3]!=NULL)hist[3]->Delete();
  hist[3]=hist2->ProjectionY(str, nxmin, nxmax);
  ratio=a/b;
  hist[3]->Add(hist[1],-ratio);
  hist[4]=s1->Background(hist[3],45);
  hist[3]->Add(hist[4],-1);
  hist[3]->SetLineColor(2);
  hist[3]->SetAxisRange(YMIN,YMAX);
  //histy->Add(histy,hist[3],0.0,1.0);
  if(cy>=1){
    c1=(TCanvas *)gROOT->FindObject("c1");
    if(c1==NULL)setcanvas(cy);
    c1->cd(cy);
    c1->SetCrosshair(1);
    hist[3]->Draw();
  }
  return ratio; //end GX
}
/*=========================================================*/
//set a gate on a gamma-gamma 2D matrix i.e. GY("GG",181,3,1)
Float_t GY(const Char_t *matn,Float_t peak,Float_t width=0.0, Int_t cy=0)
{
  TH2 *hist2;
  TH1 *hist[5];
  Char_t str[255];
  Int_t nymin,nymax;
  Float_t ratio,a,b;
  Float_t ymin,ymax;
  Axis_t axmax=axmin+numchx;
  s1=new TSpectrum();
  Int_t NX,NY;
  Double_t XMIN,XMAX;
  Double_t YMIN,YMAX;

  hist2 = (TH2*)f1->Get(matn);
  if(hist2==NULL){
    printf("Couldn't find %s\n",matn);
    return -1.0;
  }
  NX = hist2->GetNbinsX();
  XMIN = hist2->GetXaxis()->GetXmin();
  XMAX = hist2->GetXaxis()->GetXmax();
  NY = hist2->GetNbinsY();
  YMIN = hist2->GetYaxis()->GetXmin();
  YMAX = hist2->GetYaxis()->GetXmax();
  hist2->SetAxisRange(XMIN,XMAX-1,"X");
  hist2->SetAxisRange(YMIN,YMAX-1,"Y");

  sprintf(str,"%s_py",matn);
  hist[0]=(TH1*)gROOT->FindObject(str);
  sprintf(str,"%s_px",matn);
  hist[1]=(TH1*)gROOT->FindObject(str);
  sprintf(str,"%s_bgy",matn);
  hist[2]=(TH1*)gROOT->FindObject(str);

  if(hist[0]==NULL || hist[1]==NULL || hist[2]==NULL){
    dx(matn,20);
    dy(matn,20);
    return 0;
  }
  hist[0]->SetAxisRange(YMIN,YMAX);
  hist[1]->SetAxisRange(XMIN,XMAX);
  hist[2]->SetAxisRange(YMIN,YMAX);

  if(width==0.0){
    width = f2y->Eval(peak/1000.);
    width=sqrt(width);
  }
  printf("width=%f\n",width);
  ymin=peak-width;
  ymax=peak+width;
  nymin=hist[0]->FindBin(ymin);
  nymax=hist[0]->FindBin(ymax);
  a=hist[2]->Integral(nymin,nymax);
  b=hist[0]->Integral();

  sprintf(str,"%sg%4.4ix",matn,(Int_t)peak);
  hist[3]=(TH1*)gROOT->FindObject(str);
  if(hist[3]!=NULL)hist[3]->Delete();
  hist[3]=hist2->ProjectionX(str, nymin, nymax);
  ratio=a/b;
  hist[3]->Add(hist[1],-ratio);
  hist[4]=s1->Background(hist[3],45);
  hist[3]->Add(hist[4],-1);
  hist[3]->SetLineColor(2);
  hist[3]->SetAxisRange(YMIN,YMAX);

  if(cy>=1){
    c1=(TCanvas *)gROOT->FindObject("c1");
    if(c1==NULL)setcanvas(cy);
    c1->cd(cy);
    c1->SetCrosshair(1);
    hist[3]->Draw();
  }
  return ratio;
}
/*=========================================================*/
// read in a 2D matrix e.g. readmat("xx")
void readmat(const Char_t *filename){
  Int_t ix,iy,nzi;
  TH2F *hist=new TH2F(filename,filename,NCHX,0,DIMX,NCHY,0,DIMY);
  Char_t str[16];
  Int_t zmat[NCHY];
  FILE *zin;

  zin = fopen(filename,"r+");
  if(zin!=NULL){
    nzi = 4*NCHY;
    for(ix=0;ix<NCHX;ix++){
      fseek(zin, ix*nzi, SEEK_SET);
      fread(zmat,nzi,1,zin);
      for(iy=0;iy<NCHY;iy++)
        hist->SetBinContent(ix+1,iy+1,zmat[iy]);
    }
    fclose(zin);
  }else
    printf("spectrum %s not found\n", filename);
  hist->Write(filename);
  return; //end readmat
}
/*=========================================================*/
// write out 2D matrix e.g. write("gg")
void writemat(const Char_t *filename)
{
  Int_t ix,iy,nx,ny,zi;
  TH2F *hist;
  Char_t str[16];
  Int_t zmat[4096][4096];
  Int_t NX,NY;

  hist=(TH2F*)f1->Get(filename);
  FILE *zout;
  NX = hist->GetNbinsX();
  NY = hist->GetNbinsY();
  if(NX>4096)NX=4096;
  if(NY>4096)NY=4096;

  if(hist!=NULL){
    for(iy=0; iy<4096; iy++)
      for(ix=0; ix<4096; ix++)
        zmat[ix][iy] = 0;

    for(iy=1; iy<=NY; iy++)
      for(ix=1; ix<=NX; ix++)
        zmat[ix-1][iy-1]=(Int_t)hist->GetBinContent(ix,iy);

    sprintf(str, "%s.m4b", filename);
    zout = fopen(str,"wb+");
    for(zi=0;zi<4096;zi++)
      fwrite(zmat[zi],4,4096,zout);
    fclose(zout);
  } else
    printf("spectrum %s not found\n", filename);
  return; //end writemat
}
/*=========================================================*/
// write out 1D spectrum e.g. writespe("px","spectrum")
void writespe(const Char_t *hisname, Char_t *spename)
{
  TH1     *hist;
  Int_t   i,j,NN,size;
  Int_t   i1,i2;
  Char_t  str[32];
  float *sp;
  FILE *out;
  hist=(TH1*)gROOT->FindObject(hisname);
  if (hist!=NULL){
    NN = hist->GetNbinsX();
    if (!(sp = (float*) malloc(NN*sizeof(float)))) {
      printf("\007  ERROR: Could not malloc data buffer.\n");
      exit(-1);
    }
    for(i=1;i<NN+1;i++)
      sp[i-1]=(Float_t)hist->GetBinContent(i);
    sprintf(str, "%s.spe", spename);
    out=fopen(str, "wb+");
    i=1;
    j=24;
    fwrite(&j,4,1,out);
    fwrite(str,8,1,out);
    fwrite(&NN,4,1,out);
    fwrite(&i,4,1,out);
    fwrite(&i,4,1,out);
    fwrite(&i,4,1,out);
    fwrite(&j,4,1,out);
    size=sizeof(float)*NN;
    fwrite(&size,4,1,out);
    fwrite(sp,4,NN,out);
    fwrite(&size,4,1,out);
    fclose(out);
    printf("wrote %i channels to %s\n", NN, str);
  } else
    printf("spectrum %s not found\n", hisname);
  free(sp);
  return; //end writespe
}

/*=========================================================*/
// delete a histogram e.g. del hist"GamE_bgx")
void delhis(const Char_t *hn)
{
  f1->Delete(hn);
}
/*=========================================================*/
// fit a peak in a spectrum
void gaus1peakfit(const Char_t *s, Float_t x1, Float_t x2, Float_t x3, Float_t x4)
{
  Double_t par[5],epar[5],x[4],y[4];
  TH1 *hist;
  hist = (TH1 *) gROOT->FindObject(s);
  setcanvas(1);
  TCanvas *c1=(TCanvas*) gROOT->FindObject("c1");
  if(c1==NULL)setcanvas(1);
  c1->Clear();
  hist->SetAxisRange(x1-30,x4+30);
  hist->Draw();

  //--**-- Linear background estimation --**--//
  x[0] = x1;
  x[1] = x2;
  x[2] = x3;
  x[3] = x4;
  Int_t bin1 = hist->FindBin(x1);
  y[0] = hist->GetBinContent(bin1);
  Int_t bin2 = hist->FindBin(x2);
  y[1] = hist->GetBinContent(bin2);
  Int_t bin3 = hist->FindBin(x3);
  y[2] = hist->GetBinContent(bin3);
  Int_t bin4 = hist->FindBin(x4);
  y[3] = hist->GetBinContent(bin4);
  TGraph *g = new TGraph(4,x,y);
  TF1 *fpol1 = new TF1("POL1","pol1",x1,x4);
  g->Fit(fpol1,"RQN");
  par[3]=fpol1->GetParameter(0);
  par[4]=fpol1->GetParameter(1);

  //--**-- Gaussian Peak estimation without background --**--//
  TF1 *fgaus = new TF1("GAUS","gaus",x2,x3);
  hist->Fit(fgaus,"RQN");
  fgaus->GetParameters(&par[0]);

  //--**-- Final Peak Fit with Background --**--//
  TF1 *func = new TF1("FGAUS","gaus(0)+pol1(3)",x1,x4);
  func->SetParameters(par);
  hist->Fit(func,"R+QN");
  func->GetParameters(par);
  epar[0]=func->GetParError(0);
  epar[1]=func->GetParError(1);
  epar[2]=func->GetParError(2);
  Double_t fwhm = par[2]*TMath::Sqrt(8*TMath::Log(2));
  Double_t efwhm = epar[2]*TMath::Sqrt(8*TMath::Log(2));
  Double_t N0 = par[0]*(TMath::Sqrt(TMath::TwoPi())*par[2]);
  Double_t r0 = epar[0]/par[0];
  Double_t r2 = epar[2]/par[2];
  Double_t eN0= N0*TMath::Sqrt(r0*r0+r2*r2);
  func->Print();
  printf("Peak = %f +- %f; FFHM = %f +- %f; Area = %f +- %f\n",
          par[1],epar[1],fwhm,efwhm,N0,eN0);
  func->SetLineWidth(1);
  func->SetLineStyle(1);
  func->SetLineColor(4);
  func->SetFillColor(4);
  func->Draw("same");
}
/*============================================================================*/
// fit two peaks in a spectrum
void gaus2peakfit(const Char_t *s, Float_t x1, Float_t x2, Float_t x3, Float_t x4)
{
  Double_t par[8],epar[8],x[2],y[2];
  TH1 *hist;
  hist = (TH1 *) gROOT->FindObject(s);
  TCanvas *c1=(TCanvas*) gROOT->FindObject("c1");
  if(c1==NULL)setcanvas(1);
  c1->Clear();
  hist->SetAxisRange(x1-30,x4+30);
  hist->Draw();

  //--**-- Linear background estimation --**--//
  x[0] = x1;
  x[1] = x4;
  Int_t bin1 = hist->FindBin(x1);
  y[0] = hist->GetBinContent(bin1);
  Int_t bin2 = hist->FindBin(x4);
  y[1] = hist->GetBinContent(bin2);
  TGraph *g = new TGraph(2,x,y);
  TF1 *fpol1 = new TF1("POL1","pol1",x1,x4);
  g->Fit(fpol1,"RQN");
  par[6]=fpol1->GetParameter(0);
  par[7]=fpol1->GetParameter(1);

  //--**-- Two Gaussian Peak estimation without background --**--//
  TF1 *fgaus1 = new TF1("m1","gaus",x1,x2);
  TF1 *fgaus2 = new TF1("m2","gaus",x3,x4);
  hist->Fit(fgaus1,"R+QN");
  hist->Fit(fgaus2,"R+QN");
  fgaus1->GetParameters(&par[0]);
  fgaus2->GetParameters(&par[3]);

  //--**-- Final Peak Fit with Background --**--//
  TF1 *func = new TF1("m","gaus(0)+gaus(3)+pol1(6)",x1,x4);
  func->SetParameters(par);

  hist->Fit(func,"R+QN");
  func->SetLineWidth(1);
  func->SetLineStyle(1);
  func->SetLineColor(4);
  func->SetFillColor(4);
  func->Draw("same");
  func->GetParameters(par);
  epar[0]=func->GetParError(0);
  epar[1]=func->GetParError(1);
  epar[2]=func->GetParError(2);
  epar[3]=func->GetParError(3);
  epar[4]=func->GetParError(4);
  epar[5]=func->GetParError(5);
 
  Double_t fwhm1 = par[2]*TMath::Sqrt(8*TMath::Log(2));
  Double_t efwhm1 = epar[2]*TMath::Sqrt(8*TMath::Log(2));
  Double_t N10 = par[0]*(TMath::Sqrt(TMath::TwoPi())*par[2]);
  Double_t r10 = epar[0]/par[0];
  Double_t r12 = epar[2]/par[2];
  Double_t eN10= N10*TMath::Sqrt(r10*r10+r12*r12);

  Double_t fwhm2 = par[5]*TMath::Sqrt(8*TMath::Log(2));
  Double_t efwhm2 = epar[5]*TMath::Sqrt(8*TMath::Log(2));
  Double_t N20 = par[3]*(TMath::Sqrt(TMath::TwoPi())*par[5]);
  Double_t r20 = epar[3]/par[3];
  Double_t r22 = epar[5]/par[5];
  Double_t eN20= N20*TMath::Sqrt(r20*r20+r22*r22);

  printf("%11.4f %11.4f %11.0f %11.0f\n",
          par[1],epar[1],N10,eN10);
  printf("%11.4f %11.4f %11.0f %11.0f\n",
          par[4],epar[4],N20,eN20);
}

/*===========================================================================*/
// ignore from here
/*============================================================================*/
/*void dgate(Float_t peak, Float_t width=0.0, const Char_t XY="X")
{
  Float_t xx[4], yy[4];
  if(width==0.0){
    if(XY=="Y" || XY=="y"){
      width = f2y->Eval(peak/1000.);
    }
    else{
      width = f2x->Eval(peak/1000.);
    }
    width=sqrt(width);
  }
  yy[1]=0;yy[0]=10000000;yy[2]=0;yy[3]=10000000;
  Float_t ymin=peak-width;
  Float_t ymax=peak+width;
  xx[0]=ymin;
  xx[1]=xx[0];
  xx[2]=ymax;
  xx[3]=xx[2];
  TGraph *g1=(TGraph*)gROOT->FindObject("win");
  if(g1!=NULL)g1->Delete();
  g1=new TGraph(4,xx,yy);
  g1->SetName("win");
  g1->SetLineColor(6);
  g1->Draw("samel");
}
*/
/*============================================================================*/
void c2dwin()
{
  TCutG          *mycutg;

  mycutg = (TCutG *) gROOT->GetListOfSpecials()->FindObject("CUTG");
  if (mycutg != NULL)
    mycutg->Delete();
  gROOT->SetEditorMode("CutG");
  c1->Update();
  return;
}
/*=========================================================*/
void s2dwin(const char *winam, const Char_t *cnam)
{
  TCutG *mycutg;
  mycutg=(TCutG *)gROOT->GetListOfSpecials()->FindObject("CUTG");

  if (mycutg != NULL)
    {
      mycutg->SetName(cnam);
      mycutg->Print();
      TFile *fcutg=new TFile(winam,"UPDATE");
      mycutg->Write(cnam,TObject::kOverwrite);
      fcutg->ls();
      fcutg->Close();
    }
  mycutg->Delete();
  return;
}
/*=========================================================*/
void d2dwin(const Char_t *fnam, const Char_t *wnam, Int_t n, const Char_t *c)
{
  c1=(TCanvas *)gROOT->FindObject(c);
  c1->cd(n);
  TCutG *mycutg;
  TFile *fcutg;
  fcutg = new TFile(fnam);
  mycutg=(TCutG *)fcutg->Get(wnam);
  fcutg->Close();
  if (mycutg != NULL)
    {
      mycutg->SetLineColor(2);
      mycutg->Draw();
    }
  f1 = new TFile(rootfilename,"READ");
  return;
}
/*=========================================================*/
Float_t AGX(const Char_t *matn,Float_t peak,Float_t width=0.0)
{
  Axis_t aymax=aymin+numchy;
  setcanvas(2);
  TH1D *hist;
  hist=(TH1D*)gROOT->FindObject("ytemp");
  if(hist!=NULL)hist->Delete();
  hist=new TH1D("ytemp","ytemp",NCHY,0,DIMY);
  hist->Add(histy,1.0);
  Float_t a=GX(matn,peak,width);
  histy->Add(hist,1.0);
  c1->cd(2);
  histy->SetAxisRange(aymin,aymax);
  histy->SetLineColor(3);
  histy->Draw();
  return a;
}
/*=========================================================*/
Float_t AGY(const Char_t *matn,Float_t peak,Float_t width=0.0)
{
  Axis_t axmax=axmin+numchx;
  setcanvas(2);
  TH1D *hist;
  hist=(TH1D*)gROOT->FindObject("xtemp");
  if(hist!=NULL)hist->Delete();
  hist=new TH1D("xtemp","xtemp",NCHX,0,DIMX);
  hist->Add(histx,1.0);
  Float_t a=GY(matn,peak,width);
  histx->Add(hist,1.0);
  c1->cd(2);
  histx->SetAxisRange(axmin,axmax);
  histx->SetLineColor(3);
  histx->Draw();
  return a;
}
/*=========================================================*/
void fwhmfit(const Char_t *iname)
{
  Double_t f[40],c[40];
  ifstream in;
  ofstream out;
  in.open(iname);
  Int_t n=0;
  while(n<40)
    {
      in >> f[n] >> c[n];
      f[n]=f[n]*f[n];
      if (!in.good()) break;
      n++;
    }
  in.close();
  TF1 *f2 = new TF1("f2","[0]+[1]*x+[2]*x^2",0,10);
  TGraph *g1=new TGraph(n,c,f);
  g1->SetMarkerStyle(21);
  g1->SetMarkerSize(0.8);
  g1->SetMarkerColor(2);
  g1->SetLineColor(2);
  g1->SetTitle("fwhm");
  g1->Draw("alp");
  g1->Fit("f2","qn");
  f2->SetLineWidth(1);
  f2->SetLineStyle(1);
  f2->SetLineColor(4);
  f2->SetFillColor(4);
  f2->Draw("same");
  Double_t ff=f2->GetParameter(0);
  Double_t gg=f2->GetParameter(1);
  Double_t hh=f2->GetParameter(2);
  ff = sqrt(ff);
  gg = sqrt(gg);
  hh = sqrt(hh);
  out.open("fwhm.out");
  out<<"ffx="<<ff<<",ggx="<<gg<<",hhx="<<hh<< endl;
  out.close();
  printf("ff=%f,gg=%f,hh=%f\n",ff,gg,hh);
  return;
}
/*===========================================================================*/
void
adjdiff(Int_t M1LOW, Int_t M1HIGH, Int_t M2LOW, Int_t M2HIGH, const Char_t *iname, const Char_t *oname)
{
  Char_t          str[126];
  Float_t         off[110], gain[110];
  Float_t       pz[110],basecal[110],baseThr[110],M1Run[110],M2Run[110];
  TH2F           *hist2;
  TH1D           *hist1;
  Float_t         dp, *xpeaks;
  Int_t           nfound;
  Int_t           i,j,n,k;
  string          OneLine1,OneLine2;
  TSpectrum      *s1=new TSpectrum();
  Float_t         b1,a,b;
  Int_t           bin;

  //ofstream BACKCALFILE(backup, ios::out);
  //if(!BACKCALFILE.is_open()) {
  //  cerr << "Error opening Calibration file: " << backup << endl;
  //  exit(1);
  //}
  ifstream DGSCALFILE(iname, ios::in);
  if(!DGSCALFILE.is_open()) {
    cerr << "Error opening Calibration file: " << iname << endl;
    exit(1);
  }
  getline(DGSCALFILE,OneLine1);
  cout << OneLine1 << endl;
  //BACKCALFILE << OneLine1 << endl;
  getline(DGSCALFILE,OneLine2);
  cout << OneLine2 << endl;
  //BACKCALFILE << OneLine2 << endl;
  for(i=0;i<110;i++){
    DGSCALFILE >> k >> off[i] >> gain[i] >> pz[i] >> basecal[i] >> baseThr[i] >> M1Run[i] >> M2Run[i];
    cout << setw(4) << k << setw(15) << off[i] << setw(15) << gain[i] << setw(15) << pz[i]
    << setw(15) << basecal[i] << setw(15) << baseThr[i] << setw(15) << M1Run[i] << setw(15) << M2Run[i] << endl;
    //BACKCALFILE.precision(8);
    //BACKCALFILE.setf(ios::fixed,ios::floatfield);
    //BACKCALFILE << setw(4) << k << setw(15) << off[i] << setw(15) << gain[i] << setw(15) << pz[i]
    //<< setw(15) << basecal[i] << setw(15) << M1Cal[i] << setw(15) << M2Cal[i] << setw(15) << baserun[i] << setw(15) << M1Run[i]
    //<< setw(15) << M2Run[i] << endl;
  }
  DGSCALFILE.close();
  //BACKCALFILE.close();

/*
  hist2=(TH2F*)f1->Get("BaseLine");
  if(hist2==NULL){
    printf("please load the root file\n");
    return;
  }
  for(i=0;i<110;i++){
    sprintf(str,"BaseLine%3.3i",i+1);
    hist1=hist2->ProjectionX(str,i+1,i+1);
    hist1->SetAxisRange(M1LOW,M1HIGH);
    b = (Float_t)hist1->GetMaximum();
    bin = hist1->GetMaximumBin();
      a = (Float_t)hist1->GetXaxis()->GetBinCenter(bin);
    if(b<10) a = 0;
    baserun[i] = a;
    printf("ge %i: %i; %f; %f\n",i+1,bin,a,M1Run[i]);
  }
*/

  hist2=(TH2F*)f1->Get("M1OFF");
  if(hist2==NULL){
    printf("please load the root file\n");
    return;
  }
  for(i=0;i<110;i++){
    sprintf(str,"M1OFF%3.3i",i+1);
    hist1=hist2->ProjectionX(str,i+1,i+1);
    hist1->SetAxisRange(M1LOW,M1HIGH);
    b = (Float_t)hist1->GetMaximum();
    bin = hist1->GetMaximumBin();
      a = (Float_t)hist1->GetXaxis()->GetBinCenter(bin);
    if(b<10) a = 0;
    printf("ge %i: %i; %f; %f\n",i+1,bin,a,M1Run[i]);
    M1Run[i] = a;
    printf("%f\n",M1Run[i]);
  }

  hist2=(TH2F*)f1->Get("M2OFF");
  if(hist2==NULL){
    printf("please load the root file\n");
    return;
  }
  for(i=0;i<110;i++){
    sprintf(str,"M2OFF%3.3i",i+1);
    hist1=hist2->ProjectionX(str,i+1,i+1);
    hist1->SetAxisRange(M2LOW,M2HIGH);
    b = (Float_t)hist1->GetMaximum();
    bin = hist1->GetMaximumBin();
      a = (Float_t)hist1->GetXaxis()->GetBinCenter(bin);
    if(b<10) {a = 0;}
    printf("ge %i: %i; %f; %f\n",i+1,bin,a,M2Run[i]);
    M2Run[i] = a;
    printf("%f\n",M2Run[i]);
  }

  ofstream OUTCALFILE(oname, ios::out);
  if(!OUTCALFILE.is_open()){
    cerr << "Error opening Calibration file: " << oname << endl;
    exit(1);
  }
  OUTCALFILE << OneLine1 << endl;
  OUTCALFILE << OneLine2 << endl;
  for(i=0;i<110;i++){
    OUTCALFILE.precision(8);
    OUTCALFILE.setf(ios::fixed,ios::floatfield);
    OUTCALFILE << setw(4) << i+1 << setw(18) << off[i] << setw(18) << gain[i] << setw(18) << pz[i]
    << setw(18) << basecal[i] << setw(18) << baseThr[i] << setw(18)
    << M1Run[i] << setw(18) << M2Run[i] << endl;
    cout << setw(4) << i+1 << setw(18) << off[i] << setw(18) << gain[i] << setw(18) << pz[i]
    << setw(15) << basecal[i] << setw(18) << baseThr[i] << setw(18)
    << M1Run[i] << setw(18) << M2Run[i] << endl;
  }
  OUTCALFILE.close();
}
/*===========================================================================*/
void
adjoff(Int_t LOW, Int_t HIGH, Float_t Peak, const Char_t *iname, const Char_t *oname)
{
  Char_t          str[126];
  Float_t         off[110], gain[110];
  Float_t       pz[110],base[110],baseThr[110],M1Run[110],M2Run[110],offcor[110];
  TH2F           *hist2;
  TH1D           *hist1;
  Float_t         dp, *xpeaks;
  Int_t           nfound;
  Int_t           i,j,n,k;
  string          OneLine1,OneLine2;
  TSpectrum      *s1=new TSpectrum();
  Float_t         b1,a,b;
  Int_t           bin;

  //ofstream BACKCALFILE(backup, ios::out);
  //if(!BACKCALFILE.is_open()) {
  //  cerr << "Error opening Calibration file: " << backup << endl;
  //  exit(1);
  //}
  ifstream DGSCALFILE(iname, ios::in);
  if(!DGSCALFILE.is_open()) {
    cerr << "Error opening Calibration file: " << iname << endl;
    exit(1);
  }
  getline(DGSCALFILE,OneLine1);
  cout << OneLine1 << endl;
  //BACKCALFILE << OneLine1 << endl;
  getline(DGSCALFILE,OneLine2);
  cout << OneLine2 << endl;
  //BACKCALFILE << OneLine2 << endl;
  for(i=0;i<110;i++){
    DGSCALFILE >> k >> off[i] >> gain[i] >> pz[i] >> base[i] >> baseThr[i] >> M1Run[i] >> M2Run[i];
    cout << setw(4) << k << setw(18) << off[i] << setw(18) << gain[i] << setw(18) << pz[i]
    << setw(18) << base[i] << setw(18) << baseThr[i] << M1Run[i] << setw(18)
    << M2Run[i] << endl;
    //BACKCALFILE.precision(8);
    //BACKCALFILE.setf(ios::fixed,ios::floatfield);
    //BACKCALFILE << setw(4) << k << setw(15) << off[i] << setw(15) << gain[i] << setw(15) << pz[i]
    //<< setw(15) << base[i] << setw(15) << M1Cal[i] << setw(15) << M2Cal[i] << setw(15) << M1Run[i] << setw(15)
    //<< M2Run[i] << endl;
  }
  DGSCALFILE.close();
  //BACKCALFILE.close();

  hist2=(TH2F*)f1->Get("eGe");
  if(hist2==NULL){
    printf("please load the root file\n");
    return;
  }
  for(i=0;i<110;i++){
    sprintf(str,"eGe%3.3i",i+1);
    hist1=hist2->ProjectionX(str,i+1,i+1);
    hist1->SetAxisRange(LOW,HIGH);
    b = (Float_t)hist1->GetMaximum();
    bin = hist1->GetMaximumBin();
    a = (Float_t)hist1->GetXaxis()->GetBinCenter(bin);
    if(b<10){
      off[i] = 0.0;;
    }
    else{
      hist1->SetAxisRange(a-10,a+10);
      off[i] += Peak/0.333333 - (Float_t)hist1->GetMean()/0.333333;
    }
  }
  ofstream OUTCALFILE(oname, ios::out);
  if(!OUTCALFILE.is_open()){
    cerr << "Error opening Calibration file: " << iname << endl;
    exit(1);
  }
  OUTCALFILE << OneLine1 << endl;
  OUTCALFILE << OneLine2 << endl;
  for(i=0;i<110;i++){
    OUTCALFILE.precision(8);
    OUTCALFILE.setf(ios::fixed,ios::floatfield);
    OUTCALFILE << setw(4) << i+1 << setw(18) << off[i] << setw(18) << gain[i] << setw(18) << pz[i]
    << setw(18) << base[i] << setw(18) << baseThr[i] << setw(18) << M1Run[i] << setw(18) << M2Run[i] << endl;
    cout << setw(4) << i+1 << setw(18) << off[i] << setw(18) << gain[i] << setw(18) << pz[i]
    << setw(18) << base[i] << setw(18) << baseThr[i] << setw(18) << M1Run[i] << setw(18) << M2Run[i] << endl;
  }
  OUTCALFILE.close();
}
/*===========================================================================*/
void
talign(const Char_t *histname, Int_t LOW, Int_t HIGH, const Char_t *oname)
{
  Char_t          str[126];
  Float_t         off[110], gain[110];
  TH2F           *hist2;
  TH1D           *hist1;
  Float_t         dp, *xpeaks;
  Int_t           nfound;
  Int_t           i,j,n,k;
  string          OneLine1,OneLine2;
  TSpectrum      *s1=new TSpectrum();
  Int_t           bin,binc;

  ofstream OUTCALFILE(oname, ios::out);
  if(!OUTCALFILE.is_open()){
    cerr << "Error opening Calibration file: " << oname << endl;
    exit(1);
  }
  hist2=(TH2F*)f1->Get(histname);
  if(hist2==NULL){
    printf("please load the root file\n");
    return;
  }
  n = hist2->GetNbinsY();
  for(i=0;i<n;i++){
    sprintf(str,"%s%3.3i",histname,i+1);
    hist1=hist2->ProjectionX(str,i+1,i+1);
    hist1->SetAxisRange(LOW,HIGH);
    bin = hist1->GetMaximumBin();
    binc = (Int_t)hist1->GetXaxis()->GetBinCenter(bin);
    //hist1->SetAxisRange(binc-100,binc+100);
    //bin = (Int_t)hist1->GetMean();
    cout << "toff[" << i+1 <<"] = " << binc <<";"<< endl;
    OUTCALFILE << binc << endl;
  }
  OUTCALFILE.close();
}
/*===========================================================================*/
/*void
gainmatch(Int_t LOW1, Int_t HIGH1, Float_t Peak1, Int_t LOW2, Int_t HIGH2, Float_t Peak2, const Char_t *oname)
{
  Char_t          str[126];
  Float_t         off[110], gain[110];
  Float_t       pz[110],base[110],baseThr[110],M1Run[110],M2Run[110],offcor[110];
  TH2F           *hist2;
  TH1D           *hist1;
  Float_t         dp, *xpeaks;
  Int_t           nfound;
  Int_t           i,j,n,k;
  string          OneLine1,OneLine2;
  TSpectrum      *s1=new TSpectrum();
  Float_t         b1,a,b;
  Int_t           bin;
  Float_t         p1,p2;

  hist2=(TH2F*)f1->Get("gsingle");
  if(hist2==NULL){
    printf("please load the root file\n");
    return;
  }
  for(i=0;i<110;i++){
    sprintf(str,"gsingle%3.3i",i+1);
    hist1=hist2->ProjectionX(str,i+1,i+1);
    hist1->SetAxisRange(LOW1,HIGH1);
    b = (Float_t)hist1->GetMaximum();
    bin = hist1->GetMaximumBin();
    a = (Float_t)hist1->GetXaxis()->GetBinCenter(bin);
    if(b<10){
      p1 = 0.0;;
    }
    else{
      hist1->SetAxisRange(a-13,a+13);
      p1 = (Float_t)hist1->GetMean();
    }
    hist1->SetAxisRange(LOW2,HIGH2);
    b = (Float_t)hist1->GetMaximum();
    bin = hist1->GetMaximumBin();
    a = (Float_t)hist1->GetXaxis()->GetBinCenter(bin);
    if(b<10){
      p2 = 0.0;;
    }
    else{
      hist1->SetAxisRange(a-13,a+13);
      p2 = (Float_t)hist1->GetMean();
    }
    if(p1!=0.0 && p2!=0.0){
      gain[i] = (Peak1-Peak2)/(p1-p2);
      off[i] = Peak1 - gain[i]*p1;
    }
    else{
      gain[i] = 1.0;
      off[i] = 0.0;
    }
  }
  ofstream OUTCALFILE(oname, ios::out);
  if(!OUTCALFILE.is_open()){
    cerr << "Error opening Calibration file: " << iname << endl;
    exit(1);
  }
  for(i=0;i<110;i++){
    OUTCALFILE.precision(8);
    OUTCALFILE.setf(ios::fixed,ios::floatfield);
    OUTCALFILE << "off[" << i+1 << "] = " << off[i] << "; gain[" << i+1 << "] = " << gain[i] << ";" << endl;
    cout << "off[" << i+1 << "] = " << off[i] << "; gain[" << i+1 << "] = " << gain[i] << ";" << endl;
  }
  OUTCALFILE.close();
}
 */
/*===========================================================================*/
void ds(Int_t ch)
{
  Char_t          str[126];
  TH2F           *hist2;
  TH1D           *hist1;
  hist2=(TH2F*)f1->Get("Egam");
  sprintf(str,"g%i",ch);
  hist1=hist2->ProjectionY(str,ch,ch);
  hist1->SetAxisRange(200,3000);
  hist1->Draw();
}
/*===========================================================================*/
void dsall(Int_t ch0)
{
  Int_t i;
  Char_t          str[126];
  TH2F           *hist2;
  TH1D           *hist1;
  setcanvas(3,3);
  hist2=(TH2F*)f1->Get("Egam");
  for(i=ch0;i<ch0+9;i++){
    c1->cd(i-ch0+1);
    sprintf(str,"g%i",i);
    hist1=(TH1D*)f1->Get(str);
    if(hist1!=NULL)hist1->Delete();
    hist1=hist2->ProjectionY(str,i,i);
    hist1->SetAxisRange(1500,2500);
    hist1->Draw();
  }
}
/*===========================================================================*/
void dscc1(Int_t Qnum, Int_t ch0)
{
  Int_t i;
  Char_t          str1[126];
  Char_t          str2[126];
  TH2F           *hist2;
  TH1D           *hist1;
  setcanvas(4,1);
  sprintf(str1,"Q%1.1iEgam",Qnum);
  hist2=(TH2F*)f1->Get(str1);
  for(i=1;i<=4;i++){
    c1->cd(i);
    sprintf(str2,"Q%1.1iCC%1.1i_%1.1i",Qnum,ch0,i);
    hist1=(TH1D*)f1->Get(str2);
    if(hist1!=NULL)hist1->Delete();
    hist1=hist2->ProjectionY(str2,i*10+(ch0-1)*40,i*10+(ch0-1)*40);
    hist1->SetAxisRange(100,7500);
    hist1->Draw();
  }
}
/*===========================================================================*/
void CheckFWHM(Float_t p1, Float_t p2, Float_t fwhm)
{
  Float_t a = (1332.5-1173.2)/(p2-p1);
  printf("fwhm = %f\n",fwhm*a);
}
/*===========================================================================*/
void his2txt(const Char_t *h1, const Char_t *txtname, Int_t binlo, Int_t binhi)
{
  TH1 *hist;
  ofstream out;
  Int_t i,j;
  hist=(TH1*)gROOT->FindObject(h1);
  if (hist!=NULL){
    out.open(txtname);
    for(i=binlo;i<=binhi;i++){
      j=(Int_t)hist->GetBinContent(i);
      out<<i<<"\t"<<j<<endl;
    }
    out.close();
  } else
    printf("histgram not found");
  return;
}
/*===========================================================================*/
void addhist(const Char_t *h1, const Char_t *h2, Int_t scale)
{
  TH2F *hist1,*hist2;
  Axis_t axmax=axmin+numchx;
  Axis_t aymax=aymin+numchy;

  c1=(TCanvas *)gROOT->FindObject("c1");
  setcanvas(1);
  gStyle->SetPalette(1);
  gPad->SetLogz(1);
  hist1=(TH2F*)f1->Get(h1);
  hist1->SetAxisRange(axmin,axmax,"X");
  hist1->SetAxisRange(aymin,aymax,"Y");
  hist1->Draw("COLZ");
  hist2=(TH2F*)f1->Get(h2);
  hist1->Add(hist2,scale);
  return;
}
/*===========================================================================*/
