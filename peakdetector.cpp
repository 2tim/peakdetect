//______________________________________________________________________________
Int_t TSpectrum::SearchHighRes(Double_t *source,Double_t *destVector, int ssize,
                                     Double_t sigma, Double_t threshold,
                                     bool backgroundRemove,int deconIterations,
                                     bool markov, int averWindow)
{
   /* Begin_Html
   <b>One-dimensional high-resolution peak search function</b>
   <p>
   This function searches for peaks in source spectrum. It is based on
   deconvolution method. First the background is removed (if desired), then
   Markov smoothed spectrum is calculated (if desired), then the response
   function is generated according to given sigma and deconvolution is
   carried out. The order of peaks is arranged according to their heights in
   the spectrum after background elimination. The highest peak is the first in
   the list. On success it returns number of found peaks.
   <p>
   <b>Function parameters:</b>
   <ul>
   <li> source: pointer to the vector of source spectrum.
   <li> destVector: pointer to the vector of resulting deconvolved spectrum.
   <li> ssize: length of source spectrum.
   <li> sigma: sigma of searched peaks, for details we refer to manual.
   <li> threshold: threshold value in % for selected peaks, peaks with
        amplitude less than threshold*highest_peak/100
        are ignored, see manual.
   <li> backgroundRemove: logical variable, set if the removal of
        background before deconvolution is desired.
   <li> deconIterations-number of iterations in deconvolution operation.
   <li> markov: logical variable, if it is true, first the source spectrum
        is replaced by new spectrum calculated using Markov
        chains method.
   <li> averWindow: averanging window of searched peaks, for details
        we refer to manual (applies only for Markov method).
   </ul>
   <p>
   <b>Peaks searching:</b>
   <p>
   The goal of this function is to identify automatically the peaks in spectrum
   with the presence of the continuous background and statistical
   fluctuations - noise.
   <p>
   The common problems connected with correct peak identification are:
   <ul>
   <li> non-sensitivity to noise, i.e., only statistically
     relevant peaks should be identified.
   <li> non-sensitivity of the algorithm to continuous
     background.
   <li> ability to identify peaks close to the edges of the
     spectrum region. Usually peak finders fail to detect them.
   <li> resolution, decomposition of Double_tts and multiplets.
     The algorithm should be able to recognize close positioned peaks.
   <li> ability to identify peaks with different sigma.
   </ul>
   <img width=600 height=375 src="gif/TSpectrum_Searching1.jpg">
   <p>
   Fig. 27 An example of one-dimensional synthetic spectrum with found peaks
   denoted by markers.
   <p>
   <b>References:</b>
   <ol>
   <li> M.A. Mariscotti: A method for identification of peaks in the presence of
   background and its application to spectrum analysis. NIM 50 (1967),
   309-320.
   <li> M. Morhá&#269;, J. Kliman, V.  Matoušek, M. Veselský,
   I. Turzo.:Identification of peaks in
   multidimensional coincidence gamma-ray spectra. NIM, A443 (2000) 108-125.
   <li> Z.K. Silagadze, A new algorithm for automatic photopeak searches. NIM
   A 376 (1996), 451.
   </ol>
   <p>
   <b>Examples of peak searching method:</b>
   <p>
   The SearchHighRes function provides users with the possibility to vary the
   input parameters and with the access to the output deconvolved data in the
   destination spectrum. Based on the output data one can tune the parameters.
   <p>
   Example 15 - script SearchHR1.c:
   <img width=600 height=321 src="gif/TSpectrum_Searching1.jpg">
   <p>
   Fig. 28 One-dimensional spectrum with found peaks denoted by markers, 3
   iterations steps in the deconvolution.
   <p>
   <img width=600 height=323 src="gif/TSpectrum_Searching2.jpg">
   Fig. 29 One-dimensional spectrum with found peaks denoted by markers, 8
   iterations steps in the deconvolution.
   <p>
   Script:
   <pre>
   // Example to illustrate high resolution peak searching function (class TSpectrum).
   // To execute this example, do
   // root > .x SearchHR1.C

   #include <TSpectrum>

   void SearchHR1() {
      Double_t fPositionX[100];
      Double_t fPositionY[100];
      Int_t fNPeaks = 0;
      Int_t i,nfound,bin;
      Double_t nbins = 1024,a;
      Double_t xmin  = 0;
      Double_t xmax  = nbins;
      Double_t * source = new Double_t[nbins];
      Double_t * dest = new Double_t[nbins];
      TH1F *h = new TH1F("h","High resolution peak searching, number of iterations = 3",nbins,xmin,xmax);
      TH1F *d = new TH1F("d","",nbins,xmin,xmax);
      TFile *f = new TFile("spectra\\TSpectrum.root");
      h=(TH1F*) f->Get("search2;1");
      for (i = 0; i < nbins; i++) source[i]=h->GetBinContent(i + 1);
      TCanvas *Search = gROOT->GetListOfCanvases()->FindObject("Search");
      if (!Search) Search = new TCanvas("Search","Search",10,10,1000,700);
      h->SetMaximum(4000);
      h->Draw("L");
      TSpectrum *s = new TSpectrum();
      nfound = s->SearchHighRes(source, dest, nbins, 8, 2, kTRUE, 3, kTRUE, 3);
      Double_t *xpeaks = s->GetPositionX();
      for (i = 0; i < nfound; i++) {
         a=xpeaks[i];
         bin = 1 + Int_t(a + 0.5);
         fPositionX[i] = h->GetBinCenter(bin);
         fPositionY[i] = h->GetBinContent(bin);
      }
      TPolyMarker * pm = (TPolyMarker*)h->GetListOfFunctions()->FindObject("TPolyMarker");
      if (pm) {
         h->GetListOfFunctions()->Remove(pm);
         delete pm;
      }
      pm = new TPolyMarker(nfound, fPositionX, fPositionY);
      h->GetListOfFunctions()->Add(pm);
      pm->SetMarkerStyle(23);
      pm->SetMarkerColor(kRed);
      pm->SetMarkerSize(1.3);
      for (i = 0; i < nbins; i++) d->SetBinContent(i + 1,dest[i]);
      d->SetLineColor(kRed);
      d->Draw("SAME");
      printf("Found %d candidate peaks\n",nfound);
      for(i=0;i<nfound;i++)
         printf("posx= %d, posy= %d\n",fPositionX[i], fPositionY[i]);
      }
   </pre>
   <p>
   Example 16 - script SearchHR3.c:
   <p>
   <table border=solid>
   <tr><td> Peak # </td><td> Position </td><td> Sigma </td></tr>
   <tr><td> 1      </td><td> 118      </td><td> 26    </td></tr>
   <tr><td> 2      </td><td> 162      </td><td> 41    </td></tr>
   <tr><td> 3      </td><td> 310      </td><td> 4     </td></tr>
   <tr><td> 4      </td><td> 330      </td><td> 8     </td></tr>
   <tr><td> 5      </td><td> 482      </td><td> 22    </td></tr>
   <tr><td> 6      </td><td> 491      </td><td> 26    </td></tr>
   <tr><td> 7      </td><td> 740      </td><td> 21    </td></tr>
   <tr><td> 8      </td><td> 852      </td><td> 15    </td></tr>
   <tr><td> 9      </td><td> 954      </td><td> 12    </td></tr>
   <tr><td> 10     </td><td> 989      </td><td> 13    </td></tr>
   </table>
   <p>
   Table 4 Positions and sigma of peaks in the following examples.
   <p>
   <img width=600 height=328 src="gif/TSpectrum_Searching3.jpg">
   <p>
   Fig. 30 Influence of number of iterations (3-red, 10-blue, 100- green,
   1000-magenta), sigma=8, smoothing width=3.
   <p>
   <img width=600 height=321 src="gif/TSpectrum_Searching4.jpg">
   <p>
   Fig. 31 Influence of sigma (3-red, 8-blue, 20- green, 43-magenta),
   num. iter.=10, sm. width=3.
   <p>
   <img width=600 height=323 src="gif/TSpectrum_Searching5.jpg"></p>
   <p>
   Fig. 32 Influence smoothing width (0-red, 3-blue, 7- green, 20-magenta), num.
   iter.=10, sigma=8.
   <p>
   Script:
   <pre>
   // Example to illustrate the influence of number of iterations in deconvolution in high resolution peak searching function (class TSpectrum).
   // To execute this example, do
   // root > .x SearchHR3.C

   #include <TSpectrum>

   void SearchHR3() {
      Double_t fPositionX[100];
      Double_t fPositionY[100];
      Int_t fNPeaks = 0;
      Int_t i,nfound,bin;
      Double_t nbins = 1024,a;
      Double_t xmin  = 0;
      Double_t xmax  = nbins;
      Double_t * source = new Double_t[nbins];
      Double_t * dest = new Double_t[nbins];
      TH1F *h = new TH1F("h","Influence of # of iterations in deconvolution in peak searching",nbins,xmin,xmax);
      TH1F *d1 = new TH1F("d1","",nbins,xmin,xmax);
      TH1F *d2 = new TH1F("d2","",nbins,xmin,xmax);
      TH1F *d3 = new TH1F("d3","",nbins,xmin,xmax);
      TH1F *d4 = new TH1F("d4","",nbins,xmin,xmax);
      TFile *f = new TFile("spectra\\TSpectrum.root");
      h=(TH1F*) f->Get("search3;1");
      for (i = 0; i < nbins; i++) source[i]=h->GetBinContent(i + 1);
      TCanvas *Search = gROOT->GetListOfCanvases()->FindObject("Search");
      if (!Search) Search = new TCanvas("Search","Search",10,10,1000,700);
      h->SetMaximum(1300);
      h->Draw("L");
      TSpectrum *s = new TSpectrum();
      nfound = s->SearchHighRes(source, dest, nbins, 8, 2, kTRUE, 3, kTRUE, 3);
      Double_t *xpeaks = s->GetPositionX();
      for (i = 0; i < nfound; i++) {
         a=xpeaks[i];
         bin = 1 + Int_t(a + 0.5);
         fPositionX[i] = h->GetBinCenter(bin);
         fPositionY[i] = h->GetBinContent(bin);
      }
      TPolyMarker * pm = (TPolyMarker*)h->GetListOfFunctions()->FindObject("TPolyMarker");
      if (pm) {
         h->GetListOfFunctions()->Remove(pm);
         delete pm;
      }
      pm = new TPolyMarker(nfound, fPositionX, fPositionY);
      h->GetListOfFunctions()->Add(pm);
      pm->SetMarkerStyle(23);
      pm->SetMarkerColor(kRed);
      pm->SetMarkerSize(1.3);
      for (i = 0; i < nbins; i++) d1->SetBinContent(i + 1,dest[i]);
      h->Draw("");
      d1->SetLineColor(kRed);
      d1->Draw("SAME");
      for (i = 0; i < nbins; i++) source[i]=h->GetBinContent(i + 1);
      s->SearchHighRes(source, dest, nbins, 8, 2, kTRUE, 10, kTRUE, 3);
      for (i = 0; i < nbins; i++) d2->SetBinContent(i + 1,dest[i]);
      d2->SetLineColor(kBlue);
      d2->Draw("SAME");
      for (i = 0; i < nbins; i++) source[i]=h->GetBinContent(i + 1);
      s->SearchHighRes(source, dest, nbins, 8, 2, kTRUE, 100, kTRUE, 3);
      for (i = 0; i < nbins; i++) d3->SetBinContent(i + 1,dest[i]);
      d3->SetLineColor(kGreen);
      d3->Draw("SAME");
      for (i = 0; i < nbins; i++) source[i]=h->GetBinContent(i + 1);
      s->SearchHighRes(source, dest, nbins, 8, 2, kTRUE, 1000, kTRUE, 3);
      for (i = 0; i < nbins; i++) d4->SetBinContent(i + 1,dest[i]);
      d4->SetLineColor(kMagenta);
      d4->Draw("SAME");
      printf("Found %d candidate peaks\n",nfound);
   }
   </pre>
   End_Html */

   int i, j, numberIterations = (Int_t)(7 * sigma + 0.5);
   Double_t a, b, c;
   int k, lindex, posit, imin, imax, jmin, jmax, lh_gold, priz;
   Double_t lda, ldb, ldc, area, maximum, maximum_decon;
   int xmin, xmax, l, peak_index = 0, size_ext = ssize + 2 * numberIterations, shift = numberIterations, bw = 2, w;
   Double_t maxch;
   Double_t nom, nip, nim, sp, sm, plocha = 0;
   Double_t m0low=0,m1low=0,m2low=0,l0low=0,l1low=0,detlow,av,men;
   if (sigma < 1) {
      Error("SearchHighRes", "Invalid sigma, must be greater than or equal to 1");
      return 0;
   }

   if(threshold<=0 || threshold>=100){
      Error("SearchHighRes", "Invalid threshold, must be positive and less than 100");
      return 0;
   }

   j = (Int_t) (5.0 * sigma + 0.5);
   if (j >= PEAK_WINDOW / 2) {
      Error("SearchHighRes", "Too large sigma");
      return 0;
   }

   if (markov == true) {
      if (averWindow <= 0) {
         Error("SearchHighRes", "Averanging window must be positive");
         return 0;
      }
   }

   if(backgroundRemove == true){
      if(ssize < 2 * numberIterations + 1){
         Error("SearchHighRes", "Too large clipping window");
         return 0;
      }
   }

   k = int(2 * sigma+0.5);
   if(k >= 2){
      for(i = 0;i < k;i++){
         a = i,b = source[i];
         m0low += 1,m1low += a,m2low += a * a,l0low += b,l1low += a * b;
      }
      detlow = m0low * m2low - m1low * m1low;
      if(detlow != 0)
         l1low = (-l0low * m1low + l1low * m0low) / detlow;

      else
         l1low = 0;
      if(l1low > 0)
         l1low=0;
   }

   else{
      l1low = 0;
   }

   i = (Int_t)(7 * sigma + 0.5);
   i = 2 * i;
   Double_t *working_space = new Double_t [7 * (ssize + i)];
   for (j=0;j<7 * (ssize + i);j++) working_space[j] = 0;
   for(i = 0; i < size_ext; i++){
      if(i < shift){
         a = i - shift;
         working_space[i + size_ext] = source[0] + l1low * a;
         if(working_space[i + size_ext] < 0)
            working_space[i + size_ext]=0;
      }

      else if(i >= ssize + shift){
         a = i - (ssize - 1 + shift);
         working_space[i + size_ext] = source[ssize - 1];
         if(working_space[i + size_ext] < 0)
            working_space[i + size_ext]=0;
      }

      else
         working_space[i + size_ext] = source[i - shift];
   }

   if(backgroundRemove == true){
      for(i = 1; i <= numberIterations; i++){
         for(j = i; j < size_ext - i; j++){
            if(markov == false){
               a = working_space[size_ext + j];
               b = (working_space[size_ext + j - i] + working_space[size_ext + j + i]) / 2.0;
               if(b < a)
                  a = b;

               working_space[j]=a;
            }

            else{
               a = working_space[size_ext + j];
               av = 0;
               men = 0;
               for (w = j - bw; w <= j + bw; w++){
                  if ( w >= 0 && w < size_ext){
                     av += working_space[size_ext + w];
                     men +=1;
                  }
               }
               av = av / men;
               b = 0;
               men = 0;
               for (w = j - i - bw; w <= j - i + bw; w++){
                  if ( w >= 0 && w < size_ext){
                     b += working_space[size_ext + w];
                     men +=1;
                  }
               }
               b = b / men;
               c = 0;
               men = 0;
               for (w = j + i - bw; w <= j + i + bw; w++){
                  if ( w >= 0 && w < size_ext){
                     c += working_space[size_ext + w];
                     men +=1;
                  }
               }
               c = c / men;
               b = (b + c) / 2;
               if (b < a)
                  av = b;
               working_space[j]=av;
            }
         }
         for(j = i; j < size_ext - i; j++)
            working_space[size_ext + j] = working_space[j];
      }
      for(j = 0;j < size_ext; j++){
         if(j < shift){
                  a = j - shift;
                  b = source[0] + l1low * a;
                  if (b < 0) b = 0;
            working_space[size_ext + j] = b - working_space[size_ext + j];
         }

         else if(j >= ssize + shift){
                  a = j - (ssize - 1 + shift);
                  b = source[ssize - 1];
                  if (b < 0) b = 0;
            working_space[size_ext + j] = b - working_space[size_ext + j];
         }

         else{
            working_space[size_ext + j] = source[j - shift] - working_space[size_ext + j];
         }
      }
      for(j = 0;j < size_ext; j++){
         if(working_space[size_ext + j] < 0) working_space[size_ext + j] = 0;
      }
   }

   for(i = 0; i < size_ext; i++){
      working_space[i + 6*size_ext] = working_space[i + size_ext];
   }

   if(markov == true){
      for(j = 0; j < size_ext; j++)
         working_space[2 * size_ext + j] = working_space[size_ext + j];
      xmin = 0,xmax = size_ext - 1;
      for(i = 0, maxch = 0; i < size_ext; i++){
         working_space[i] = 0;
         if(maxch < working_space[2 * size_ext + i])
            maxch = working_space[2 * size_ext + i];
         plocha += working_space[2 * size_ext + i];
      }
      if(maxch == 0) {
         delete [] working_space;
         return 0;
      }

      nom = 1;
      working_space[xmin] = 1;
      for(i = xmin; i < xmax; i++){
         nip = working_space[2 * size_ext + i] / maxch;
         nim = working_space[2 * size_ext + i + 1] / maxch;
         sp = 0,sm = 0;
         for(l = 1; l <= averWindow; l++){
            if((i + l) > xmax)
               a = working_space[2 * size_ext + xmax] / maxch;

            else
               a = working_space[2 * size_ext + i + l] / maxch;

            b = a - nip;
            if(a + nip <= 0)
               a=1;

            else
               a = TMath::Sqrt(a + nip);

            b = b / a;
            b = TMath::Exp(b);
            sp = sp + b;
            if((i - l + 1) < xmin)
               a = working_space[2 * size_ext + xmin] / maxch;

            else
               a = working_space[2 * size_ext + i - l + 1] / maxch;

            b = a - nim;
            if(a + nim <= 0)
               a = 1;

            else
               a = TMath::Sqrt(a + nim);

            b = b / a;
            b = TMath::Exp(b);
            sm = sm + b;
         }
         a = sp / sm;
         a = working_space[i + 1] = working_space[i] * a;
         nom = nom + a;
      }
      for(i = xmin; i <= xmax; i++){
         working_space[i] = working_space[i] / nom;
      }
      for(j = 0; j < size_ext; j++)
         working_space[size_ext + j] = working_space[j] * plocha;
      for(j = 0; j < size_ext; j++){
         working_space[2 * size_ext + j] = working_space[size_ext + j];
      }
      if(backgroundRemove == true){
         for(i = 1; i <= numberIterations; i++){
            for(j = i; j < size_ext - i; j++){
               a = working_space[size_ext + j];
               b = (working_space[size_ext + j - i] + working_space[size_ext + j + i]) / 2.0;
               if(b < a)
                  a = b;
               working_space[j] = a;
            }
            for(j = i; j < size_ext - i; j++)
               working_space[size_ext + j] = working_space[j];
         }
         for(j = 0; j < size_ext; j++){
            working_space[size_ext + j] = working_space[2 * size_ext + j] - working_space[size_ext + j];
         }
      }
   }
//deconvolution starts
   area = 0;
   lh_gold = -1;
   posit = 0;
   maximum = 0;
//generate response vector
   for(i = 0; i < size_ext; i++){
      lda = (Double_t)i - 3 * sigma;
      lda = lda * lda / (2 * sigma * sigma);
      j = (Int_t)(1000 * TMath::Exp(-lda));
      lda = j;
      if(lda != 0)
         lh_gold = i + 1;

      working_space[i] = lda;
      area = area + lda;
      if(lda > maximum){
         maximum = lda;
         posit = i;
      }
   }
//read source vector
   for(i = 0; i < size_ext; i++)
      working_space[2 * size_ext + i] = TMath::Abs(working_space[size_ext + i]);
//create matrix at*a(vector b)
   i = lh_gold - 1;
   if(i > size_ext)
      i = size_ext;

   imin = -i,imax = i;
   for(i = imin; i <= imax; i++){
      lda = 0;
      jmin = 0;
      if(i < 0)
         jmin = -i;
      jmax = lh_gold - 1 - i;
      if(jmax > (lh_gold - 1))
         jmax = lh_gold - 1;

      for(j = jmin;j <= jmax; j++){
         ldb = working_space[j];
         ldc = working_space[i + j];
         lda = lda + ldb * ldc;
      }
      working_space[size_ext + i - imin] = lda;
   }
//create vector p
   i = lh_gold - 1;
   imin = -i,imax = size_ext + i - 1;
   for(i = imin; i <= imax; i++){
      lda = 0;
      for(j = 0; j <= (lh_gold - 1); j++){
         ldb = working_space[j];
         k = i + j;
         if(k >= 0 && k < size_ext){
            ldc = working_space[2 * size_ext + k];
            lda = lda + ldb * ldc;
         }

      }
      working_space[4 * size_ext + i - imin] = lda;
   }
//move vector p
   for(i = imin; i <= imax; i++)
      working_space[2 * size_ext + i - imin] = working_space[4 * size_ext + i - imin];
//initialization of resulting vector
   for(i = 0; i < size_ext; i++)
      working_space[i] = 1;
//START OF ITERATIONS
   for(lindex = 0; lindex < deconIterations; lindex++){
      for(i = 0; i < size_ext; i++){
         if(TMath::Abs(working_space[2 * size_ext + i]) > 0.00001 && TMath::Abs(working_space[i]) > 0.00001){
            lda=0;
            jmin = lh_gold - 1;
            if(jmin > i)
               jmin = i;

            jmin = -jmin;
            jmax = lh_gold - 1;
            if(jmax > (size_ext - 1 - i))
               jmax=size_ext-1-i;

            for(j = jmin; j <= jmax; j++){
               ldb = working_space[j + lh_gold - 1 + size_ext];
               ldc = working_space[i + j];
               lda = lda + ldb * ldc;
            }
            ldb = working_space[2 * size_ext + i];
            if(lda != 0)
               lda = ldb / lda;

            else
               lda = 0;

            ldb = working_space[i];
            lda = lda * ldb;
            working_space[3 * size_ext + i] = lda;
         }
      }
      for(i = 0; i < size_ext; i++){
         working_space[i] = working_space[3 * size_ext + i];
      }
   }
//shift resulting spectrum
   for(i=0;i<size_ext;i++){
      lda = working_space[i];
      j = i + posit;
      j = j % size_ext;
      working_space[size_ext + j] = lda;
   }
//write back resulting spectrum
   maximum = 0, maximum_decon = 0;
   j = lh_gold - 1;
   for(i = 0; i < size_ext - j; i++){
      if(i >= shift && i < ssize + shift){
         working_space[i] = area * working_space[size_ext + i + j];
         if(maximum_decon < working_space[i])
            maximum_decon = working_space[i];
         if(maximum < working_space[6 * size_ext + i])
            maximum = working_space[6 * size_ext + i];
      }

      else
         working_space[i] = 0;
   }
   lda=1;
   if(lda>threshold)
      lda=threshold;
   lda=lda/100;

//searching for peaks in deconvolved spectrum
   for(i = 1; i < size_ext - 1; i++){
      if(working_space[i] > working_space[i - 1] && working_space[i] > working_space[i + 1]){
         if(i >= shift && i < ssize + shift){
            if(working_space[i] > lda*maximum_decon && working_space[6 * size_ext + i] > threshold * maximum / 100.0){
               for(j = i - 1, a = 0, b = 0; j <= i + 1; j++){
                  a += (Double_t)(j - shift) * working_space[j];
                  b += working_space[j];
               }
               a = a / b;
               if(a < 0)
                  a = 0;

               if(a >= ssize)
                  a = ssize - 1;
               if(peak_index == 0){
                  fPositionX[0] = a;
                  peak_index = 1;
               }

               else{
                  for(j = 0, priz = 0; j < peak_index && priz == 0; j++){
                     if(working_space[6 * size_ext + shift + (Int_t)a] > working_space[6 * size_ext + shift + (Int_t)fPositionX[j]])
                        priz = 1;
                  }
                  if(priz == 0){
                     if(j < fMaxPeaks){
                        fPositionX[j] = a;
                     }
                  }

                  else{
                     for(k = peak_index; k >= j; k--){
                        if(k < fMaxPeaks){
                           fPositionX[k] = fPositionX[k - 1];
                        }
                     }
                     fPositionX[j - 1] = a;
                  }
                  if(peak_index < fMaxPeaks)
                     peak_index += 1;
               }
            }
         }
      }
   }

   for(i = 0; i < ssize; i++) destVector[i] = working_space[i + shift];
   delete [] working_space;
   fNPeaks = peak_index;
   if(peak_index == fMaxPeaks)
      Warning("SearchHighRes", "Peak buffer full");
   return fNPeaks;
}