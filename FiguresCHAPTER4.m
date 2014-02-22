%% Appendix 4.3: CODE TO GENERATE FIGURES USED IN CHAPTER 4

%%% Introduction
% There is a lot of redundant code here.
% To create Figure 4.8 see Appendix 4.4

%% Data / Toolboxes / Scripts used
%
% * Signal processing toolbox (<http://www.mathworks.com MathWorks>)
% * HILBERT2 and DERIVATIVE by Scott McKinney (<http://www.mathworks.com/matlabcentral/fileexchange/authors/110216
% FileExchange>)
% * The Sound Processing Toolbox by Naotoshi Seo (<http://note.sonots.com/SciSoftware/Pitch.html SciSoftware>)
% * MIR toolbox created by departement of Music at the University of Jyv?skyl? Finland
% (<https://www.jyu.fi/hum/laitokset/musiikki/en/research/coe/materials/mirtoolbox/mirtoolbox MIR>)
% * Several helper scripts and data files created and authored by Fred Hasselman included in the .zip file

%% Author / Version
%
% Created by: <http://www.fredhasselman.com/site/Creations/Creations.html *Fred Hasselman*> / January 2011
%
% Affiliations: <http://www.ru.nl/bsi _Behavioural Science Institute_> - <http://www.ru.nl _Radboud University Nijmegen_>

%% PREP
%Change ... to the path on your machine where you stored the files
path='...';

%% Figure 4.1: F2 Slope
%
% Redundant code: Calculate MIR Spectrum again

% Uncomment next line to clear and close everything... detergent grade!
% omo

cd(path)
load('QDAmeasures.mat')
SPEC = getspec;

h0=figure;
maximize(h0);

cm=flipud(gray);
v = [-200:10:0];
icm = 1:floor(length(cm)/length(v)):length(cm-20);
cmp = cm(icm,:);
colormap(cmp);

cnt=0;
for man=1:4
 for stim = 1:10
  cnt=cnt+1;
  
  y = stimuli(cnt).y;
  t = (1:length(y))./stimuli(cnt).fs;
  [mirAFile] = miraudio(y);
  [mirSFile] = mirspectrum(mirAFile,'Frame',SPEC.ts,'Min',SPEC.f(1),'Max',max(SPEC.f),'Window',SPEC.wintype,'Res',SPEC.sf);
  SPECgram(cnt).gram = mirgetdata(mirSFile);
  [f tl] = size(SPECgram(cnt).gram);
  
  step = floor(length(t)/tl);
  tx = decimate(t,step);
  
  if length(tx)<=(tl-1)
   padsize(cnt)=(tl-length(tx));
   if padsize(cnt) == 1
    b = 1; en = 0;
   elseif mod(padsize(cnt),2) ~= 0
    b = floor(padsize(cnt)/2); en = b+1;
   else
    b = padsize(cnt)/2; en = b;
   end
   tx=padarray(tx,b,'pre');tx=padarray(tx,en,'post');
  else
   shrink(cnt) = length(tx)-(tl-1);
   if shrink(cnt) == 1
    b = 1; en = 0;
   elseif mod(shrink(cnt),2) ~= 0
    b = floor(shrink(cnt)/2); en = b+1;
   else
    b = shrink(cnt)/2; en = b;
   end
   tx = tx(b:length(tx)-en);
  end
  
  [t,m,n] = unique(Formants(1,cnt).tracks{1,1});
  
  IN = Formants(1,cnt).tracks{1,2};
  F1 = Formants(1,cnt).tracks{1,4};
  F2 = Formants(1,cnt).tracks{1,5};
  F3 = Formants(1,cnt).tracks{1,6};
  
  IN = IN(m); F1 = F1(m); F2 = F2(m); F3 = F3(m);
  
  lvl = .1;
  
  Sweep(cnt).TI =  t(IN>=lvl);
  Sweep(cnt).F1 = F1(IN>=lvl);
  Sweep(cnt).F2 = F2(IN>=lvl);
  Sweep(cnt).F3 = F3(IN>=lvl);
  
  dsF1 = smooth(Sweep(cnt).F1,.6,'rloess');
  dsF2 = smooth(Sweep(cnt).F2,.6,'rloess');
  dsF3 = smooth(Sweep(cnt).F3,.6,'rloess');
  
  subplot(4,10,cnt);
  
  [c h_spec] = contourf(tx,SPEC.f(2:end),20*log10(abs(SPECgram(cnt).gram)+eps),v,'LineColor','none');hold on;
  ax0 = gca;
  
  h_F1 = plot(ax0,Sweep(cnt).TI,dsF1,'-','Color',[.9 .9 .9],'LineWidth',3); hold on;
  h_F2 = plot(ax0,Sweep(cnt).TI,dsF2,'-','Color',[.9 .9 .9],'LineWidth',3); hold on;
  h_F3 = plot(ax0,Sweep(cnt).TI,dsF3,'-','Color',[.9 .9 .9],'LineWidth',3); hold on;
  
  set(ax0,'Ytick',[0:1000:3000],'YtickLabel',{'','1','2',''});
  
  ylim([1 3000]);
  
  grid on;
  
  [mxF2 In] = max(dsF2);
  Sweep(cnt).F2mx = mxF2; Sweep(cnt).tF2mx = Sweep(cnt).TI(In);
  [mnF2 mIn] = min(dsF2);
  Sweep(cnt).F2mn = mnF2; Sweep(cnt).F2mn = Sweep(cnt).TI(mIn); Sweep(cnt).swpF2 = (mxF2-mnF2)/(Sweep(cnt).TI(In)-Sweep(cnt).TI(mIn));
  
  if ismember(cnt,[1:10])
   title(num2str(cnt));
   xlim([0.01 .6]);
  end
  if ismember(cnt,[11:20])
   title(num2str(cnt-10));
   xlim([0.01 .9]);
  end
  if ismember(cnt,[21:30])
   title(num2str(cnt-20));
   xlim([0.01 .6]);
  end
  if ismember(cnt,[31:40])
   title(num2str(cnt-30));
   xlim([0.01 .9]);
  end
  
  Tpos = [-.1 -0.06 0 0];
  
  if cnt==1
   Opos=get(ax0,'Position');
   h_1= annotation('textbox',[Opos+Tpos],'String','None','EdgeColor','none','FontSize',16);
   set(h_1,'FitBoxToText','on');
   ylabel('Frquency (kHz)');
   xlabel('Time (s)');
  end
  if cnt==11
   Opos=get(ax0,'Position');
   h_2= annotation('textbox',[Opos+Tpos],'String','Slowed Down','EdgeColor','none','FontSize',16);
   set(h_2,'FitBoxToText','on');
   ylabel('Frquency (kHz)');
   xlabel('Time (s)');
  end
  if cnt==21
   Opos=get(ax0,'Position');
   h_3= annotation('textbox',[Opos+Tpos],'String','Amplified','EdgeColor','none','FontSize',16);
   set(h_3,'FitBoxToText','on');
   ylabel('Frquency (kHz)');
   xlabel('Time (s)');
  end
  if cnt==31
   Opos=get(ax0,'Position');
   h_4=annotation('textbox',[Opos+Tpos],'String','Both','EdgeColor','none','FontSize',16);
   set(h_4,'FitBoxToText','on');
   h_5=annotation('textbox',[.146 .065 0 0],'String','/bAk/','EdgeColor','none','FontSize',16);
   set(h_5,'FitBoxToText','on');
   ylabel('Frquency (kHz)');
   xlabel('Time (s)');
  end
  if cnt==40
   h_6=annotation('textbox',[.862 .065 0 0],'String','/dAk/','EdgeColor','none','FontSize',16);
   set(h_6,'FitBoxToText','on');
  end
  
  h_7 = annotation('arrow',[.2 .84],[.05 .05],'HeadStyle','cback2','LineWidth',1.5);
  h_t = text(.2,150,['\DeltaF2 = ',deblank(num2str(Sweep(cnt).swpF2/1000,'%#1.2f'))],'BackgroundColor',[.9 .9 .9],'Margin',0.01); %[.8 .8 .8] ,'EdgeColor','k'
  
  clear IN F1 F2 F3 t m n c h_spec h_s
  
 end
end

grab('sweepFORM1',0)

%% Figure 4.2: maxENVELOPE Slope
%
% Slope till max formant amplitude from stimulus onset and formant onset

% Uncomment next line to clear and close everything... detergent grade!
% omo

cd(path)
load('justwav.mat');

h0=figure;
maximize(h0);

cnt=0;
for man=1:4
 for stim = 1:10
  cnt=cnt+1;
  
  dsY   = stimuli(cnt).y;
  
  %Create miraudio for further use
  [mirAFile] = miraudio(dsY);
  
  %Get MIRenvelope
  [mirENV]     = mirenvelope(mirAFile);
  envSTATS(cnt).ENV = mirgetdata(mirENV);
  dsENV = envSTATS(cnt).ENV;
  lENV  = length(dsENV);
  t = (1:length(dsY))./stimuli(cnt).fs;
  
  step = floor(length(dsY)/lENV);
  xx = decimate(dsY,step);
  tx = decimate(t,step);
  
  if length(xx)<=(lENV-1)
   padsize(cnt)=(lENV-length(xx));
   if padsize(cnt) == 1
    b = 1; en = 0;
   elseif mod(padsize(cnt),2) ~= 0
    b = floor(padsize(cnt)/2); en = b+1;
   else
    b = padsize(cnt)/2; en = b;
   end
   xx=padarray(xx,b,'pre');xx=padarray(xx,en,'post');
   tx=padarray(tx,b,'pre');tx=padarray(tx,en,'post');
  else
   shrink(cnt) = length(xx)-(lENV-1);
   if shrink(cnt) == 1
    b = 1; en = 0;
   elseif mod(shrink(cnt),2) ~= 0
    b = floor(shrink(cnt)/2); en = b+1;
   else
    b = shrink(cnt)/2; en = b;
   end
   xx = xx(b:length(xx)-en);
   tx = tx(b:length(tx)-en);
  end
  
  xx = xx./5;
  
  subplot(4,10,cnt);
  ax1 = gca;
  
  plot(tx,[xx],'LineWidth',.1,'Color',[.7 .7 .7]);hold on;
  axis tight;
  plot(tx,[dsENV],'LineWidth',2,'Color',[.3 .3 .3]);hold on;
  axis tight;
  plot(tx,[-dsENV],'LineWidth',2,'Color',[.3 .3 .3]);hold on;
  axis tight;
  
  
  %Plot min to max AMP line
  [mxENV mxI] = max(dsENV);
  envSTATS(cnt).mxENV = (mxENV-dsENV(1))/tx(mxI);
  IASmxO = plot([tx(1) tx(mxI)],[dsENV(1) mxENV],'Color','k');hold on
  plot([tx(1) tx(mxI)],[dsENV(1) mxENV],'o','MarkerSize',4,'MarkerEdgeColor',[.3 .3 .3],'MarkerFaceColor',[.8 .8 .8]);
  
  axis tight;
  set(ax1,'Ytick',[-.3:.1:.3],'YtickLabel',{'','','','0','','',''});
  
  ylim([-.3 .3]);
  xlim([0 .9]);
  grid on;
  
  %Print slope in figure
  text(.5,-.25,['\Delta = ',num2str(envSTATS(cnt).mxENV,'%1.2f')]);
  
  %Garnish
  if ismember(cnt,[1:10])
   title(num2str(cnt));
  end
  if ismember(cnt,[11:20])
   title(num2str(cnt-10));
  end
  if ismember(cnt,[21:30])
   title(num2str(cnt-20));
  end
  if ismember(cnt,[31:40])
   title(num2str(cnt-30));
  end
  
  Tpos = [-.1 -0.06 0 0];
  
  if cnt==1
   Opos=get(ax1,'Position');
   h_1= annotation('textbox',[Opos+Tpos],'String','None','EdgeColor','none','FontSize',16);
   set(h_1,'FitBoxToText','on');
   ylabel('Amplitude (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==11
   Opos=get(ax1,'Position');
   h_2= annotation('textbox',[Opos+Tpos],'String','Slowed Down','EdgeColor','none','FontSize',16);
   set(h_2,'FitBoxToText','on');
   ylabel('Amplitude (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==21
   Opos=get(ax1,'Position');
   h_3= annotation('textbox',[Opos+Tpos],'String','Amplified','EdgeColor','none','FontSize',16);
   set(h_3,'FitBoxToText','on');
   ylabel('Amplitude (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==31
   Opos=get(ax1,'Position');
   h_4=annotation('textbox',[Opos+Tpos],'String','Both','EdgeColor','none','FontSize',16);
   set(h_4,'FitBoxToText','on');
   h_5=annotation('textbox',[.146 .065 0 0],'String','/bAk/','EdgeColor','none','FontSize',16);
   set(h_5,'FitBoxToText','on');
   ylabel('Amplitude (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==40
   h_6=annotation('textbox',[.862 .065 0 0],'String','/dAk/','EdgeColor','none','FontSize',16);
   set(h_6,'FitBoxToText','on');
  end
  
  h_7 = annotation('arrow',[.2 .84],[.05 .05],'HeadStyle','cback2','LineWidth',1.5);
  
 end
end

grab('sweepENV',0);

%% Figure 4.3: Rise and Fall Time Entropy
%
% Redundant code: Calculate Envelopes again

% Uncomment next line to clear and close everything... detergent grade!
% omo

cd(path);
load('bakdakdata.mat');

h0=figure;
maximize(h0);

cnt=0;
for man=1:4
 for stim = 1:10
  cnt=cnt+1;
  
  dsY   = stimuli(cnt).y;
  
  % Create miraudio for further use
  [mirAFile] = miraudio(dsY);
  
  % Get MIRenvelope
  [mirENV]     = mirenvelope(mirAFile);
  envSTATS(cnt).ENV = mirgetdata(mirENV);
  sENV = envSTATS(cnt).ENV;
  t = (1:length(dsY))./stimuli(cnt).fs;
  
  dsENV = derivative(sENV);
  lENV  = length(dsENV);
  
  % Resample
  step = floor(length(dsY)/lENV);
  xx = decimate(dsY,step);
  tx = decimate(t,step);
  
  if length(xx)<=(lENV-1)
   padsize(cnt)=(lENV-length(xx));
   if padsize(cnt) == 1
    b = 1; en = 0;
   elseif mod(padsize(cnt),2) ~= 0
    b = floor(padsize(cnt)/2); en = b+1;
   else
    b = padsize(cnt)/2; en = b;
   end
   xx=padarray(xx,b,'pre');xx=padarray(xx,en,'post');
   tx=padarray(tx,b,'pre');tx=padarray(tx,en,'post');
  else
   shrink(cnt) = length(xx)-(lENV-1);
   if shrink(cnt) == 1
    b = 1; en = 0;
   elseif mod(shrink(cnt),2) ~= 0
    b = floor(shrink(cnt)/2); en = b+1;
   else
    b = shrink(cnt)/2; en = b;
   end
   xx = xx(b:length(xx)-en);
   tx = tx(b:length(tx)-en);
  end
  
  % Use to create ghost image
  xx = xx./5;
  
  scaleENV = 1.5*10^2;
  dsENV =  dsENV.*scaleENV;
  
  subplot(4,10,cnt);
  ax1 = gca;
  
  plot(tx,xx,'LineWidth',.1,'Color',[.7 .7 .7]);hold on;
  axis tight;
  plot(tx,dsENV,'LineWidth',2,'Color',[.3 .3 .3]);hold on;
  axis tight;
  
  [zeroX,t0,s0] = crossing(dsENV,tx);
  
  t0 = [0 t0]; s0 = [0 s0];
  RTent(cnt) = entropy(nonzeros(sort(round(diff(t0)*1000),'ascend')));
  
  % Print entropy in figure
  text(.1,-.25,['RFTe = ',num2str(RTent(cnt),'%1.2f')]);
  
  plot(t0,s0,'.','MarkerSize',6,'MarkerEdgeColor','w','MarkerFaceColor',[.8 .8 .8]);
  
  set(ax1,'Ytick',[-.3:.1:.3],'YtickLabel',{'','','','0','','',''});
  
  ylim([-.3 .3]);
  xlim([0 .9]);
  grid on;
  
  % Garnish
  if ismember(cnt,[1:10])
   title(num2str(cnt));
  end
  if ismember(cnt,[11:20])
   title(num2str(cnt-10));
  end
  if ismember(cnt,[21:30])
   title(num2str(cnt-20));
  end
  if ismember(cnt,[31:40])
   title(num2str(cnt-30));
  end
  
  Tpos = [-.1 -0.06 0 0];
  
  if cnt==1
   Opos=get(ax1,'Position');
   h_1= annotation('textbox',[Opos+Tpos],'String','None','EdgeColor','none','FontSize',16);
   set(h_1,'FitBoxToText','on');
   ylabel('Amplitude change (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==11
   Opos=get(ax1,'Position');
   h_2= annotation('textbox',[Opos+Tpos],'String','Slowed Down','EdgeColor','none','FontSize',16);
   set(h_2,'FitBoxToText','on');
   ylabel('Amplitude change (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==21
   Opos=get(ax1,'Position');
   h_3= annotation('textbox',[Opos+Tpos],'String','Amplified','EdgeColor','none','FontSize',16);
   set(h_3,'FitBoxToText','on');
   ylabel('Amplitude change (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==31
   Opos=get(ax1,'Position');
   h_4=annotation('textbox',[Opos+Tpos],'String','Both','EdgeColor','none','FontSize',16);
   set(h_4,'FitBoxToText','on');
   h_5=annotation('textbox',[.146 .065 0 0],'String','/bAk/','EdgeColor','none','FontSize',16);
   set(h_5,'FitBoxToText','on');
   ylabel('Amplitude change (a.u.)');
   xlabel('Time (s)');
  end
  if cnt==40
   h_6=annotation('textbox',[.862 .065 0 0],'String','/dAk/','EdgeColor','none','FontSize',16);
   set(h_6,'FitBoxToText','on');
  end
  
  clear pks trgh pk tr dsENV
  
 end
end

h_7 = annotation('arrow',[.2 .84],[.05 .05],'HeadStyle','cback2','LineWidth',1.5);

keep RTent

save('RFTe.mat');
grab('sweepENV_d',0);


%% Figure 4.4: Phase Space Reconstruction Example

% Uncomment next line to clear and close everything... detergent grade!
% omo

cd(path);

load('eq_ssize.mat');
load('rpde_dfa3.mat','rpSTATS')

cnt=10; ds=1;
x = rpTS(cnt).ts(1:1024,2);
t = rpTS(cnt).ts(1:1024,1);
tau=6; m=3;

x=minplus1(x);

x1=x(1+(m-3)*tau:end-((m-1)*tau)+1); ts1=t(1+(m-3)*tau:end-((m-1)*tau)+1); %  t1=minplus1(t1);
x2=x(0+(m-2)*tau:end-((m-2)*tau));   ts2=t(0+(m-2)*tau:end-((m-2)*tau));   %  t2=minplus1(t2);
x3=x(0+(m-1)*tau:end-((m-3)*tau));   ts3=t(0+(m-1)*tau:end-((m-3)*tau));   %  t3=minplus1(t3);

t = -1:2/(length(x1)-1):1;

t1=t;
t2=t;
t3=t;

z1=-1+zeros(length(t1),1);
z2= 1+zeros(length(t2),1);
z3= 1+zeros(length(t3),1);

h0=figure;
maximize(h0);

ax_PS =axes('NextPlot','add');

h_ts1=plot3(x1(:),x2(:),z1,'-k');axis square;
xlim([-1 1]),ylim([-1 1]);zlim([-1 1]); view(3);
h_ts2=plot3(z2,x2(:),x3(:),'-k');axis square;
h_ts3=plot3(x1(:),z3,x3(:),'-k');axis square;

set(ax_PS,'XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0 1]);
set([h_ts1 h_ts2 h_ts3],'Color',[.7 .7 .7]);

h_ps=plot3(x1(:), x2(:), x3(:),'-ko'); grid on; axis square;
set(h_ps,'MarkerFaceColor',[.5 .5 .5]);
xlabel('X_m_=_1'); ylabel('X_m_=_2'); zlabel('X_m_=_3');
title({'Reconstructed Phase Space of Stimulus 10 (first 1024 samples)', ['Delay Embedding with m = 3, \tau = 6, \epsilon = ',num2str(rpSTATS(cnt,1),1),' of maximum norm distance']});

s=.18;
point=525;
xc=[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]*s;xc=xc+x1(point)-(s/2);
yc=[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]*s;yc=yc+x2(point)-(s/2);
zc=[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]*s;zc=zc+x3(point)-(s/2);

axes(ax_PS);
for i=1:6
 h=patch(xc(:,i),yc(:,i),zc(:,i),[.8 .8 .8]);
 set(h,'edgecolor','k','FaceAlpha',.5,'LineWidth',1)
end

xc1 = xc(:,1); yc1 = yc(:,1); zc1 = zc(:,1);
xc2 = xc(:,2); yc2 = yc(:,2); zc2 = zc(:,2);
xc3 = xc(:,3); yc3 = yc(:,3); zc3 = zc(:,3);
xc4 = xc(:,4); yc4 = yc(:,4); zc4 = zc(:,4);
xc5 = xc(:,5); yc5 = yc(:,5); zc5 = zc(:,5);
xc6 = xc(:,6); yc6 = yc(:,6); zc6 = zc(:,6);

axes(ax_PS);
h1=patch(xc5(:),yc5(:),[-1 -1 -1 -1]',[.8 .8 .8]);
h2=patch(xc6(:),[1 1 1 1]',zc1(:),[.8 .8 .8]);
h3=patch([1 1 1 1]',yc2(:),zc1(:),[.8 .8 .8]);

set([h1 h2 h3],'edgecolor','k','FaceAlpha',.5,'LineWidth',1)

text(xc(1,1)-(s/3),yc(1,1)+(s/3),zc(1,1)+.15,'\epsilon','FontSize',20);

axes(ax_PS);
hc1 = plot3(x1(point),x2(point),-1,'ok');
hc2 = plot3(1,x2(point),x3(point),'ok');
hc3 = plot3(x1(point),1,x3(point),'ok');
set([hc1 hc2 hc3],'MarkerFaceColor',[.2 .2 .2],'MarkerSize',8);

point3=795;
hc7 = plot3(x1(point3),x2(point3),-1,'sk');
hc8 = plot3(1,x2(point3),x3(point3),'sk');
hc9 = plot3(x1(point3),1,x3(point3),'sk');
set([hc7 hc8 hc9],'MarkerFaceColor',[.7 .7 .7],'MarkerSize',8);

keep point rpTS s stimuli
grab('psEXAMPLE',0)

%% Figure 4.5 (LEFT): RP example

% Uncomment next line to clear and close everything... detergent grade!
% omo

cd(path)
load('rpde_dfa3.mat')
cnt=10; ds=1;m=3;tau=6;

MTRX = rpMTRX(cnt).rp;
[r c] = size(MTRX);

scrsz = get(0,'ScreenSize');
h0 = figure('Position',[scrsz(3)/4 0 scrsz(3)/2 scrsz(3)/2],'NextPlot','add');

RPdims  =[.35 .35 .6 .6];

TSHdims1=[.35 .06 .6 .06];
TSHdims2=[.35 .15 .6 .06];
TSHdims3=[.35 .24 .6 .06];

TSVdims1=[.06 .35 .06 .6];
TSVdims2=[.15 .35 .06 .6];
TSVdims3=[.24 .35 .06 .6];

%Recurrence Plot
ax_RP =axes('Position',RPdims);
spy(MTRX,'.k',1);
axis square; axis xy;
xlabel('Recurrent values in m-dimensional phase space');
title({'(Auto) Recurrence Plot of Stimulus 10:',['\tau = 6, m = 3, \epsilon = ',num2str(rpSTATS(cnt,1),1)]});
loi=line([0 r],[0 c],'Color','k','LineWidth',3);

x1=rpTS(cnt).ts(1+(m-3)*tau:end-((m-1)*tau),2); y1=rpTS(cnt).ts(1+(m-3)*tau:end-((m-1)*tau),1);
x2=rpTS(cnt).ts(0+(m-2)*tau:end-((m-2)*tau),2); y2=rpTS(cnt).ts(0+(m-2)*tau:end-((m-2)*tau),1);
x3=rpTS(cnt).ts(0+(m-1)*tau:end-((m-3)*tau),2); y3=rpTS(cnt).ts(0+(m-1)*tau:end-((m-3)*tau),1);

%Horizontal TS
ax_TSH1=axes('Position',TSHdims1);%drawnow
h_TSH1=line(y1,x1); axis tight
xlabel('Surrogate Dimensions: m Time Series Delayed by m * \tau');
ylabel('X_m_=_1');

ax_TSH2=axes('Position',TSHdims2);%drawnow
h_TSH2=line(y2,x2); axis tight
ylabel('X_m_=_2');

ax_TSH3=axes('Position',TSHdims3);%drawnow
h_TSH3=line(y3,x3); axis tight
ylabel('X_m_=_3');

%Vertical TS
ax_TSV1=axes('Position',TSVdims1);%drawnow
h_TSV1=line(x1,y1); axis tight
ylabel('Surrogate Dimensions: m Time Series Delayed by m * \tau');
xlabel('X_m_=_1');

ax_TSV2=axes('Position',TSVdims2);%drawnow
h_TSV2=line(x2,y2); axis tight
xlabel('X_m_=_2');

ax_TSV3=axes('Position',TSVdims3);%drawnow
h_TSV3=line(x3,y3); axis tight
xlabel('X_m_=_3');

set([h_TSH1,h_TSH2,h_TSH3,h_TSV1,h_TSV2,h_TSV3],'Color',[.5 .5 .5]);

set([ax_RP],...
 'XTick',[1 round(r/2) r],...
 'YTick',[1 round(c/2) c],...
 'XTickLabel',{'1' num2str(round(r/2)) num2str(r)},...
 'YTickLabel',{'1' num2str(round(r/2)) num2str(r)});

set([ax_TSH1],...
 'XTick',[y1(1) y1(round(length(y1)/2)) y1(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y1(1),4) num2str(y1(round(length(y1)/2)),4) num2str(y1(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSH2],...
 'XTick',[y2(1) y2(round(length(y2)/2)) y2(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y2(1),4) num2str(y2(round(length(y2)/2)),4) num2str(y2(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSH3],...
 'XTick',[y3(1) y3(round(length(y3)/2)) y3(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y3(1),4) num2str(y3(round(length(y1)/2)),4) num2str(y3(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSV1],...
 'YTick',[y1(1) y1(round(length(y1)/2)) y1(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' '0' '1'},'Box','on');

set([ax_TSV2],...
 'YTick',[y2(1) y2(round(length(y2)/2)) y2(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' 0'' '1'},'Box','on');

set([ax_TSV3],...
 'YTick',[y3(1) y3(round(length(y3)/2)) y3(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' '0' '1'},'Box','on');

rpST = {['RQA measures:'],...
 [' '],...
 ['REC = ',num2str(rpSTATS(cnt,2),2)],...
 ['DET = ',num2str(rpSTATS(cnt,3),2)],...
 ['Lmn = ',num2str(rpSTATS(cnt,4),3)],...
 ['ENT = ',num2str(rpSTATS(cnt,6),3)],...
 ['LAM = ',num2str(rpSTATS(cnt,7),2)],...
 ['Vmn = ',num2str(rpSTATS(cnt,8),3)]};

h_s = annotation('textbox',[.11 .25 0 0],'String',rpST,'EdgeColor','none','FontName','Courier','FontSize',12);
set(h_s,'FitBoxToText','on');

grab('rpEXAMPLE')

%% Figure 4.5 (RIGHT): RP RANDOM example

% Uncomment next line to clear and close everything... detergent grade!
% omo

cd(path);
load('randBD10.mat')
load('rpde_dfa3.mat','rpTS')

cnt=10;

%RQA settings (based on Mutual Information and Nearest Neigbour analysis
tau=6; m=3; e=.01;
[r c] = size(MTRX);

scrsz = get(0,'ScreenSize');
h0 = figure('Position',[scrsz(3)/4 0 scrsz(3)/2 scrsz(3)/2],'NextPlot','add');

RPdims  =[.35 .35 .6 .6];

TSHdims1=[.35 .06 .6 .06];
TSHdims2=[.35 .15 .6 .06];
TSHdims3=[.35 .24 .6 .06];

TSVdims1=[.06 .35 .06 .6];
TSVdims2=[.15 .35 .06 .6];
TSVdims3=[.24 .35 .06 .6];

%Recurrence Plot
ax_RP =axes('Position',RPdims);
spy(MTRX,'.k',1);
axis square; axis xy;
xlabel('Recurrent values in m-dimensional phase space');
title({'(Auto) Recurrence Plot of Stimulus 10 (RANDOMISED):',['\tau = 6, m = 3, \epsilon = ',num2str(STATS(1),1)]});
loi=line([0 r],[0 c],'Color','k','LineWidth',3);

x1=xus(1+(m-3)*tau:end-((m-1)*tau)); y1=rpTS(cnt).ts(1+(m-3)*tau:end-((m-1)*tau),1);
x2=xus(0+(m-2)*tau:end-((m-2)*tau)); y2=rpTS(cnt).ts(0+(m-2)*tau:end-((m-2)*tau),1);
x3=xus(0+(m-1)*tau:end-((m-3)*tau)); y3=rpTS(cnt).ts(0+(m-1)*tau:end-((m-3)*tau),1);

%Horizontal TS
ax_TSH1=axes('Position',TSHdims1);%drawnow
h_TSH1=line(y1,x1); axis tight
xlabel('Surrogate Dimensions: m Time Series Delayed by m * \tau');
ylabel('X_m_=_1');

ax_TSH2=axes('Position',TSHdims2);%drawnow
h_TSH2=line(y2,x2); axis tight
ylabel('X_m_=_2');

ax_TSH3=axes('Position',TSHdims3);%drawnow
h_TSH3=line(y3,x3); axis tight
ylabel('X_m_=_3');

%Vertical TS
ax_TSV1=axes('Position',TSVdims1);%drawnow
h_TSV1=line(x1,y1); axis tight
ylabel('Surrogate Dimensions: m Time Series Delayed by m * \tau');
xlabel('X_m_=_1');

ax_TSV2=axes('Position',TSVdims2);%drawnow
h_TSV2=line(x2,y2); axis tight
xlabel('X_m_=_2');

ax_TSV3=axes('Position',TSVdims3);%drawnow
h_TSV3=line(x3,y3); axis tight
xlabel('X_m_=_3');

set([h_TSH1,h_TSH2,h_TSH3,h_TSV1,h_TSV2,h_TSV3],'Color',[.5 .5 .5]);

set([ax_RP],...
 'XTick',[1 round(r/2) r],...
 'YTick',[1 round(c/2) c],...
 'XTickLabel',{'1' num2str(round(r/2)) num2str(r)},...
 'YTickLabel',{'1' num2str(round(r/2)) num2str(r)});

set([ax_TSH1],...
 'XTick',[y1(1) y1(round(length(y1)/2)) y1(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y1(1),4) num2str(y1(round(length(y1)/2)),4) num2str(y1(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSH2],...
 'XTick',[y2(1) y2(round(length(y2)/2)) y2(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y2(1),4) num2str(y2(round(length(y2)/2)),4) num2str(y2(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSH3],...
 'XTick',[y3(1) y3(round(length(y3)/2)) y3(end)],...
 'YTick',[-1 0 1],...
 'XTickLabel',{num2str(y3(1),4) num2str(y3(round(length(y1)/2)),4) num2str(y3(end),4)},...
 'YTickLabel',{'' '' ''},'Box','on');

set([ax_TSV1],...
 'YTick',[y1(1) y1(round(length(y1)/2)) y1(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' '0' '1'},'Box','on');

set([ax_TSV2],...
 'YTick',[y2(1) y2(round(length(y2)/2)) y2(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' 0'' '1'},'Box','on');

set([ax_TSV3],...
 'YTick',[y3(1) y3(round(length(y3)/2)) y3(end)],...
 'XTick',[-1 0 1],...
 'YTickLabel',{'' '' ''},...
 'XTickLabel',{'-1' '0' '1'},'Box','on');

rpST = {['RQA measures:'],...
 [' '],...
 ['REC = ',num2str(STATS(2),2)],...
 ['DET = ',num2str(STATS(3),1)],...
 ['Lmn = ',num2str(STATS(4),3)],...
 ['ENT = ',num2str(STATS(6),1)],...
 ['LAM = ',num2str(STATS(7),1)],...
 ['Vmn = ',num2str(STATS(8),3)]};

h_s = annotation('textbox',[.11 .25 0 0],'String',rpST,'EdgeColor','none','FontName','Courier','FontSize',12);
set(h_s,'FitBoxToText','on');

grab('randRP',0);

keep MTRX STATS xus
save('randBD10.mat')

%% Figure 4.6: RP plots

% Uncomment next line to clear and close everything... detergent grade!
% omo

cd(path)
load('rpde_dfa3.mat')

ds = 2;
h0=figure;
maximize(h0);

for cnt=1:40
 
 subplot(4,10,cnt)
 tssz=length(rpTS(cnt).ts(:,1));
 
 rr = downsample(rpMTRX(cnt).rp,ds); cc = downsample(transpose(rr),ds);
 [r c] = size(cc);
 spy(cc,'.k',1);
 ax0 = gca;
 axis square; axis xy;
 xlabel(''); ylabel('');title(num2str(cnt));
 xlim([0 r(end)]);ylim([0 c(end)]);
 
 set(ax0,'XTick',[0 (round(r(end)/2)) r(end)],...
  'YTick',[0 (round(c(end)/2)) c(end)],...
  'XTickLabel',{'','',''},...
  'YTickLabel',{'','',''});
 
 if ismember(cnt,[1:10])
  title([num2str(cnt),'- \epsilon =',num2str(rpSTATS(cnt,1),1)]);
 end
 if ismember(cnt,[11:20])
  title([num2str(cnt-10),'- \epsilon =',num2str(rpSTATS(cnt,1),1)]);
 end
 if ismember(cnt,[21:30])
  title([num2str(cnt-20),'- \epsilon =',num2str(rpSTATS(cnt,1),1)]);
 end
 if ismember(cnt,[31:40])
  title([num2str(cnt-30),'- \epsilon =',num2str(rpSTATS(cnt,1),1)]);
 end
 
 Tpos = [-.1 -0.06 0 0];
 Opos = get(ax0,'Position');
 
 ax_TS = axes('Position',[Opos(1),Opos(2)-.03,Opos(3),Opos(4)/4]);
 h_TSH = line(rpTS(cnt).ts(1:end,1),rpTS(cnt).ts(1:end,2),'Color',[.5 .5 .5]); axis tight
 set(ax_TS,'Visible','off');
 
 if cnt==1
  h_1= annotation('textbox',[Opos+Tpos],'String','None','EdgeColor','none','FontSize',16);
  set(h_1,'FitBoxToText','on');
  
 end
 if cnt==11
  h_2= annotation('textbox',[Opos+Tpos],'String','Slowed Down','EdgeColor','none','FontSize',16);
  set(h_2,'FitBoxToText','on');
  
 end
 if cnt==21
  h_3= annotation('textbox',[Opos+Tpos],'String','Amplified','EdgeColor','none','FontSize',16);
  set(h_3,'FitBoxToText','on');
  
 end
 if cnt==31
  h_4=annotation('textbox',[Opos+Tpos],'String','Both','EdgeColor','none','FontSize',16);
  set(h_4,'FitBoxToText','on');
  h_5=annotation('textbox',[.146 .065 0 0],'String','/bAk/','EdgeColor','none','FontSize',16);
  set(h_5,'FitBoxToText','on');
 end
 if cnt==40
  h_6=annotation('textbox',[.862 .065 0 0],'String','/dAk/','EdgeColor','none','FontSize',16);
  set(h_6,'FitBoxToText','on');
 end
 
 h_7 = annotation('arrow',[.2 .84],[.05 .05],'HeadStyle','cback2','LineWidth',1.5);
 h_8 = annotation('textbox',[.39 .99 0 0],'String','(Auto) Recurrence Plots for all Stimuli','EdgeColor','none','FontSize',16);
 set(h_8,'FitBoxToText','on');
 
end

grab('rpMTRX',0)

%% Figure 4.7: LOGIT predictions

% Uncomment next line to clear and close everything... detergent grade!
% omo

cd(path);
load('predicted odds CI.mat');
%%
f0=figure;
maximize(f0);
subplot(2,2,1)
ax0=gca; pos=get(ax0,'Position');

none  = errorbar([1:10]-0.1,data(:,1),[data(:,1)-data(:,2)],[data(:,3)-data(:,1)],'-');
hold on;

noned = errorbar([1:10]+0.1,data2(:,1),[data2(:,1)-data2(:,2)],[data2(:,3)-data2(:,1)],'-.');
hold on;
xlim([0.7 10.3]);
ylim([-0.05 1.05]);
threshold = line(xlim,[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle','--');

set([none, noned], ...
 'Marker'          , 'o'        , ...
 'MarkerSize'      , 5          , ...
 'MarkerEdgeColor' , 'k'     , ...
 'MarkerFaceColor' , [.5 .5 .5] , ...
 'Color'       , 'k');

set(none, ...
 'MarkerFaceColor' , [0 0 0]);

set(gca, ...
 'TickDir'     , 'out'     , ...
 'TickLength'  , [.01 .01] , ...
 'YTick'       , 0:.1:1, ...
 'XTickLabel'  , {'/bAk/','2','3','4','5','6','7','8','9','/dAk/'},...
 'Box','off', ...
 'LineWidth'   , 1         );

legend('Average reader','Dyslexic reader','Location','SouthEast');
title('A. None');
xlabel('Stimulus');
ylabel('\fontsize{12}\it\pi\fontsize{10}\rm with \itCI\fontsize{8}_{.95}\fontsize{10}\rm  for Perceiving /dAk/');


subplot(2,2,2)
none=errorbar([1:10]-0.1,data(:,4),[data(:,4)-data(:,5)],[data(:,6)-data(:,4)],'-');
hold on;
noned=errorbar([1:10]+0.1,data2(:,4),[data2(:,4)-data2(:,5)],[data2(:,6)-data2(:,4)],'-.');
hold on;
xlim([0.7 10.3]);
ylim([-0.05 1.05]);
threshold = line(xlim,[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle','--');

set([none, noned], ...
 'Marker'          , 'v'        , ...
 'MarkerSize'      , 5          , ...
 'MarkerEdgeColor' , 'k'     , ...
 'MarkerFaceColor' , [.5 .5 .5] , ...
 'Color'       , 'k');

set(none, ...
 'MarkerFaceColor' , [0 0 0]);

set(gca, ...
 'TickDir'     , 'out'     , ...
 'TickLength'  , [.01 .01] , ...
 'YTick'       , 0:.1:1, ...
 'XTickLabel'  , {'/bAk/','2','3','4','5','6','7','8','9','/dAk/'},...
 'Box','off', ...
 'LineWidth'   , 1         );


legend('Average reader','Dyslexic reader','Location','SouthEast');
title('B. Slowed down');
xlabel('Stimulus');
ylabel('\fontsize{12}\it\pi\fontsize{10}\rm with \itCI\fontsize{8}_{.95}\fontsize{10}\rm  for Perceiving /dAk/');


subplot(2,2,3)
none=errorbar([1:10]-0.1,data(:,7),[data(:,7)-data(:,8)],[data(:,9)-data(:,7)],'-');
hold on;
noned=errorbar([1:10]+0.1,data2(:,7),[data2(:,7)-data2(:,8)],[data2(:,9)-data2(:,7)],'-.');
hold on;
xlim([0.7 10.3]);
ylim([-0.05 1.05]);
threshold = line(xlim,[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle','--');

set([none, noned], ...
 'Marker'          , '^'        , ...
 'MarkerSize'      , 5          , ...
 'MarkerEdgeColor' , 'k'     , ...
 'MarkerFaceColor' , [.5 .5 .5] , ...
 'Color'       , 'k');

set(none, ...
 'MarkerFaceColor' , [0 0 0]);


set(gca, ...
 'TickDir'     , 'out'     , ...
 'TickLength'  , [.01 .01] , ...
 'YTick'       , 0:.1:1, ...
 'XTickLabel'  , {'/bAk/','2','3','4','5','6','7','8','9','/dAk/'},...
 'Box','off', ...
 'LineWidth'   , 1         );

legend('Average reader','Dyslexic reader','Location','SouthEast');
title('C. Amplified');
xlabel('Stimulus');
ylabel('\fontsize{12}\it\pi\fontsize{10}\rm with \itCI\fontsize{8}_{.95}\fontsize{10}\rm  for Perceiving /dAk/');


subplot(2,2,4)
none=errorbar([1:10]-0.1,data(:,10),[data(:,10)-data(:,11)],[data(:,12)-data(:,10)],'-');
hold on;
noned=errorbar([1:10]+0.1,data2(:,10),[data2(:,10)-data2(:,11)],[data2(:,12)-data2(:,10)],'-.');
hold on;
xlim([0.7 10.3]);
ylim([-0.05 1.05]);
threshold = line(xlim,[0.5 0.5],'Color',[0.5 0.5 0.5],'LineStyle','--');

set([none, noned], ...
 'Marker'          , 'd'        , ...
 'MarkerSize'      , 5          , ...
 'MarkerEdgeColor' , 'k'     , ...
 'MarkerFaceColor' , [.5 .5 .5] , ...
 'Color'       , 'k');

set(none, ...
 'MarkerFaceColor' , [0 0 0]);

set(gca, ...
 'TickDir'     , 'out'     , ...
 'TickLength'  , [.01 .01] , ...
 'YTick'       , 0:.1:1, ...
 'XTickLabel'  , {'/bAk/','2','3','4','5','6','7','8','9','/dAk/'},...
 'Box','off', ...
 'LineWidth'   , 1         );

legend('Average reader','Dyslexic reader','Location','SouthEast');
title('D. Both');
xlabel('Stimulus');
ylabel('\fontsize{12}\it\pi\fontsize{10}\rm with \itCI\fontsize{8}_{.95}\fontsize{10}\rm  for Perceiving /dAk/');

% Print to EPS
grab('ident',0)

%% Figure 4.x: Equal sample TS -- NOT USED

% Uncomment next line to clear and close everything... detergent grade!
% omo

cd(path)
load('eq_ssize.mat');

for cnt=1:40
 
 subplot(4,10,cnt)
 plot([1:4096],rpTS(cnt).ts(:,2),'-k');
 ax0 = gca; set(ax0,'XTick',[1 2048 4096]);
 xlim([1 4096]); ylim([-1 1]);
 xlabel('Sample #');
 axis square;
 
 %Garnish
 if ismember(cnt,[1:10])
  title(num2str(cnt));
 end
 if ismember(cnt,[11:20])
  title(num2str(cnt-10));
 end
 if ismember(cnt,[21:30])
  title(num2str(cnt-20));
 end
 if ismember(cnt,[31:40])
  title(num2str(cnt-30));
 end
 
 Tpos = [-.1 -0.06 0 0];
 
 if cnt==1
  Opos=get(ax0,'Position');
  h_1= annotation('textbox',[Opos+Tpos],'String','None','EdgeColor','none','FontSize',16);
  set(h_1,'FitBoxToText','on');
  ylabel('Amplitude (a.u.)');
 end
 if cnt==11
  Opos=get(ax0,'Position');
  h_2= annotation('textbox',[Opos+Tpos],'String','Slowed Down','EdgeColor','none','FontSize',16);
  set(h_2,'FitBoxToText','on');
  ylabel('Amplitude (a.u.)');
 end
 if cnt==21
  Opos=get(ax0,'Position');
  h_3= annotation('textbox',[Opos+Tpos],'String','Amplified','EdgeColor','none','FontSize',16);
  set(h_3,'FitBoxToText','on');
  ylabel('Amplitude (a.u.)');
 end
 if cnt==31
  Opos=get(ax0,'Position');
  h_4=annotation('textbox',[Opos+Tpos],'String','Both','EdgeColor','none','FontSize',16);
  set(h_4,'FitBoxToText','on');
  h_5=annotation('textbox',[.146 .065 0 0],'String','/bAk/','EdgeColor','none','FontSize',16);
  set(h_5,'FitBoxToText','on');
  ylabel('Amplitude (a.u.)');
 end
 if cnt==40
  h_6=annotation('textbox',[.862 .065 0 0],'String','/dAk/','EdgeColor','none','FontSize',16);
  set(h_6,'FitBoxToText','on');
 end
 
 h_7 = annotation('arrow',[.2 .84],[.05 .05],'HeadStyle','cback2','LineWidth',1.5);
 
end

clear Opos Tpos ax0 ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 cnt h_1 h_2 h_3 h_4 h_5 h_6 h_7

grab('equalsamples',0)
