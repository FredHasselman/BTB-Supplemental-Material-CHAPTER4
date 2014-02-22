%% Appendix 4.1: EXTRACTING MEASURES FROM .WAV FILES
%
%%% Introduction
% This is a demonstration script accompanying the fourth chapter of my dissertation (Beyond the Boundary). Its purpose is to
% provide an example of how to use various freely available MATLAB sources on the web to extract variables from speech
% stimuli. This script should thus be used as an example to build your own scripts only. It is not a function or toolbox, not
% optimized for speed or functionality! Evaluate the code per Cell (using MATLAB's Cell Mode) and inspect the workspace to
% see what is going on.
%
% Extracts spectral and pitch data, formant sweeps and amplitude envelopes and various measures of temporal evolution of the
% signal.
%
% Results should be comparable to those obtained by <http://www.praat.org Praat> (Boersma & Weeninck).

%% Toolboxes / Scripts used
%
% * Signal processing toolbox (<http://www.mathworks.com MathWorks>)
% * HILBERT2 and DERIVATIVE by Scott McKinney (<http://www.mathworks.com/matlabcentral/fileexchange/authors/110216
% FileExchange>)
% * The Sound Processing Toolbox by Naotoshi Seo (<http://note.sonots.com/SciSoftware/Pitch.html SciSoftware>)
% * MIR toolbox created by departement of Music at the University of Jyv?skyl? Finland
% (<https://www.jyu.fi/hum/laitokset/musiikki/en/research/coe/materials/mirtoolbox/mirtoolbox MIR>)

%% Author / Version
%
% Created by: <http://www.fredhasselman.com/site/Creations/Creations.html *Fred Hasselman*> / January 2011
%
% Affiliations: <http://www.ru.nl/bsi _Behavioural Science Institute_> - <http://www.ru.nl _Radboud University Nijmegen_>

%% PREP: SETTINGS
% Reads .WAV data and performs spectral, envelope and pitch analyses.

% Change these variables to reflect your datafile names and locations. 
% Here there are four different manipulations (manis) of a 
% 10-step continuum (stims), with files named accordingly.

manis     = {'BAKDAK', 'BAKDAKV', 'BAKDAKA', 'BAKDAKB'};
stims     = 10;

%Change ... to the path on your machine where you stored the files
wavPath   = '...';
wavPath   = '/Users/Fred/Documents/Matlab/PAPER - LOGIT/';
txtPath   = '/Users/Fred/Documents/Matlab/PAPER - LOGIT/FORMANT TABLES/';
ext       = '.WAV';

DELIMITER = '\t';
SKIP      = '--undefined--';

%Spectrogram settings put into structure SPEC
SPEC.fs        = 44100;     % Or read from file, see for loop below
SPEC.ts        = .002;      % Praat settings: Time step         (s)
SPEC.wl        = .005;      % Praat settings: Window length     (s)
SPEC.mf        = 5500;      % Praat settings: Maximum frequency (Hz) 
SPEC.sf        = 20;        % Praat settings: Frequency step    (Hz)


%%% Praat default behaviour
%
% Praat changes Frequencey step (|sf|) and Time step (|ts|) to the most 
% sensible minimal values given current Window length (wl) in the case they
% were set too low by the user.

if SPEC.ts < SPEC.wl/(8*sqrt(pi))
SPEC.ts = SPEC.wl/(8*sqrt(pi));
end
if SPEC.sf < (sqrt(pi)/8/SPEC.wl)
SPEC.sf = (sqrt(pi)/8/SPEC.wl);
end

SPEC.f       = 1:SPEC.sf:SPEC.mf; % Spectrogram will be est. at freqs in SPEC.f
SPEC.wintype = 'gausswin';        % Praat settings: Gaussian window is used

SPEC.nfft    = round(SPEC.wl*SPEC.fs); % Convert s to samples
SPEC.noverlap= round(SPEC.ts*SPEC.fs); % Convert s to samples
SPEC.window  = eval(sprintf('%s(SPEC.nfft)', SPEC.wintype)); % Create window

%Formant tracking with LPC, female speaker of 5 formants
SPEC.fl = 25;                     % Frame length  (ms)
SPEC.fo = 5;                      % Frame Overlap (ms)
SPEC.po = 10;                     % Prediction order (LPC poles / coefs)
SPEC.d  = gcd(2*SPEC.mf,SPEC.fs); % Get smallest integer conversion ratio
SPEC.p  = (2*SPEC.mf)/SPEC.d;     % Resample to   factor for LPC
SPEC.q  = SPEC.fs/SPEC.d;         % Resample from factor for LPC


% CHECK SETTINGS
%
% * *Check frequency vector:* [f(1) f(end)]
% * *Check the window properties:* wvtool(SPEC.window)
%

%% Load the .WAV files into MATLAB

cnt=0;
cnt2=0;
for mani = 1:length(manis)
 for stim = 1:stims
  
  wavFile = [wavPath char(manis(mani)) num2str(stim) ext];
  txtFile = [txtPath char(manis(mani)) num2str(stim) '.txt'];
  
  if exist(wavFile,'file')
   cnt=cnt+1;
   
   % Get wavefile and put into structure [stimuli(cnt)._]
   [stimuli(cnt).y stimuli(cnt).fs stimuli(cnt).nbits,...
    stimuli(cnt).opts]=wavread(wavFile);
   stimuli(cnt).name = [char(manis(mani)) num2str(stim)];
   
  else
   disp(wavFile)
  end
   
   % Read Praat Formant track files
   if exist(txtFile,'file')
    cnt2=cnt2+1;
    fid=fopen(txtFile,'r');
    Formants(cnt2).tracks = textscan(fid,'%*s %f %f %f %f %f %f %f %f','TreatAsEmpty',SKIP,'HeaderLines',1);
    Formants(cnt2).name = txtFile;
    fclose(fid);
   else
    disp(txtFile);
   end
   
 end
end

clear cnt cnt2 wavFile wavPath ext txtPath txtFile DELIMITER SKIP fid

save('QDAmeasures.mat');

%% Get Envelope and Formant measures

cnt=0;
for mani = 1:length(manis)
 for stim = 1:stims
  
  cnt=cnt+1;
  
  % Perform spectrogram analysis (SP toolbox) and put results into structure
  [STIM(cnt).S STIM(cnt).F STIM(cnt).T, STIM(cnt).P] = spectrogram(stimuli(cnt).y,SPEC.window,SPEC.noverlap,SPEC.f,SPEC.fs);
  STIM(cnt).name = [char(manis(mani)) num2str(stim)];
  
  % Get immediate amplitude envelope and frequency by Hilbert transform
  % (HILBERT2 by Scott McKinney) and put results into structure STIM(cnt).IA .IF
  [STIM(cnt).IA, STIM(cnt).IF] = hilbert2(stimuli(cnt).y,SPEC.fs);
  
  %Smooth the Envelope
  STIM(cnt).IAsm = smooth(STIM(cnt).IA,.1,'loess');
  STIM(cnt).IAT  = [1:length(STIM(cnt).IAsm)]./stimuli(cnt).fs;
  STIM(cnt).IAdf = diff(STIM(cnt).IAT);
  
  %Get MAX envelope amplitude at formant
  [mxIA mxIn] = max(STIM(cnt).IAsm);
  STIM(cnt).IAmx = mxIA; STIM(cnt).IATmx = STIM(cnt).IAT(mxIn);
  
  %Slope till MAX amplitude envelope (from stimulus onset)
  STIM(cnt).IASmxO = (STIM(cnt).IAmx-STIM(cnt).IAsm(1))/(STIM(cnt).IATmx-STIM(cnt).IAT(1));
  
  % Get formant tracks using LPC (Sound Processing Toolbox by Naotoshi Seo).
  % Resample so Nyquist = max frequency before LPC
  % (this means sample rate should be 2*SPEC.mf)
  y_rs = resample(stimuli(cnt).y,SPEC.p,SPEC.q);
  [STIM(cnt).FT, STIM(cnt).FTt] = spFormantsTrackLpc(y_rs,2*SPEC.mf,SPEC.po,...
   SPEC.fl,SPEC.fo,SPEC.wintype,0);
  
  % Get pitch track by using autocorrelation method
  % (Sound Processing Toolbox by Naotoshi Seo)
  [STIM(cnt).F0, STIM(cnt).F0t, STIM(cnt).F0r] = spPitchTrackCorr(stimuli(cnt).y,...
   SPEC.fs,SPEC.fl,SPEC.fo,[],0);
  STIM(cnt).name = [char(manis(mani)) num2str(stim)];
  
  clear mxIA mxIn y_rs
  
 end % for stims
end % for mani

clear cnt mani manis stim stims


save('QDAmeasures.mat');

%% Create TS with equal samples, normalized to [-1 1] for HNR and RQA

for i = 1:40
 
 fs = stimuli(i).fs;
 x  = stimuli(i).y;
 
 % Find transition part
 indb = find(x.^2>=.2,1,'first');
 xlr  = flipud(x);
 inde = find(xlr.^2>=.2,1,'first');
 xc   = x(indb:(length(x)-inde)); ts = indb/fs;
 t    = ([1:length(xc)]'./fs)+ts;
 
 % Create TS of size 4096
 step = floor(length(xc)/4096);
 xx   = decimate(xc,step);
 tx   = decimate(t,step);
 
 % Shrink or pad TS if resampling wasn't exactly 4096
 if length(xx)<=4095
  padsize(i)=(4096-length(xx));
  if padsize(i) == 1
   b = 1; en = 0;
  elseif mod(padsize(i),2) ~= 0
   b = floor(padsize(i)/2); en = b+1;
  else
   b = padsize(i)/2; en = b;
  end
  xx=padarray(xx,b,'pre');xx=padarray(xx,en,'post');
  tx=padarray(tx,b,'pre');tx=padarray(tx,en,'post');
 else
  shrink(i) = length(xx)-4095;
  if shrink(i) == 1
   b = 1; en = 0;
  elseif mod(shrink(i),2) ~= 0
   b = floor(shrink(i)/2); en = b+1;
  else
   b = shrink(i)/2; en = b;
  end
  xx = xx(b:length(xx)-en);
  tx = tx(b:length(tx)-en);
 end
 
 %Normalize waveform to [-1 1]
 rng=(max(xx)-min(xx));
 mid=(max(xx)+min(xx))/2;
 xu = (xx-(mid))/(rng/2);
 
 %Store into structure rpTS
 rpTS(i).ts = [tx xu];
 
 clear t ts tx fs x xc xx xn xlr xu en b step indb inde rng mid j
 
end

keep SPEC rpTS STIM stimuli

save('QDAmeasures.mat');

%% Calculate HNR based on resampled waveforms 

for i= 1:40
 
  %Create miraudio (MIR toolbox University of Jyv?skyl?)
 [mirAFile] = miraudio(rpTS(i).ts(:,2));
 [mirSFile] = mirspectrum(mirAFile,'Frame',SPEC.ts,'Min',...
  SPEC.f(1),'Max',max(SPEC.f),'Window',SPEC.wintype,'Res',SPEC.sf);
 
 % Get Pitch track by using MIR toolbox (also ac). Pitch floor
 % and ceiling frequencies are Praat defaults.
 [mirF0] = mirpitch(mirAFile,'Mono','Min',75,'Max',600);
 HNR(i).F0mir = mirgetdata(mirF0);
 
 % Get inharmonicity based on pitch estimates
 % mirSTATS(i).F0acm = median(mirSTATS(i).F0);
 [hnr] = mirinharmonicity(mirAFile,'f0',mirF0);
 HNR(i).HNR = mirgetdata(hnr);
 
 clear mirF0 mirAFile mirSFile hnr;
 
end

save('QDAmeasures.mat');


%% Calculate RECURRENCE MEASURES
% WARNING - this takes a long time!

rpSTATS = zeros(40,14);

for i=1:40
 
 xu = rpTS(i).ts(:,2);
 
 %RQA settings (based on Mutual Information and Nearest Neigbour analysis
 tau=6;m=3; e=.01;W=[];WS=[];LMIN=2;VMIN=2;TW=0;thr= 'rr';
 
 % Get RP for plotting and threshold of Fixed RR using an adaptation of the
 % crp.m code in Marwan's toolbox. The edited version is included in the zip
 % file available at http:/fredhasselman.com/site/BTB.html
 % Uncomment if you want the RP matrix and threshold
 % [rpMTRX(i).rp rpSTATS(i,1)] = crp(xu,m,tau,e,thr,'nonormalize','silent');
    
 %RQA measures
 rpSTATS(i,2:14) = crqa(xu,m,tau,e,W,WS,LMIN,VMIN,TW,thr,'nonormalize','nogui');
 
 clear xu
 
end

keep SPEC stimuli rpTS HNR rpSTATS STIM

save('QDAmeasures.mat');
