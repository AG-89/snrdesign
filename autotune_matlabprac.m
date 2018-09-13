%snr design autotune practice
%Team: Ashlin G, Alex G, Brandon Y
clc, clear variables, close all, clear sound, clear sounds;

global count; count = 420; %temp for quadratic interp graphs

%CONTROL VARIABLES
    %constants
    samplerate = 44100; %sample rate
    global A4freq; A4freq = 440; %Freq (Hz) for note A4
    C4 = A4freq * 2^(-9/12); %middle c (C4) from A4

    %control
    downsample_rate = 8; %down/up sample by these factors, new pitch = f * d/u
    upsample_rate = 1;
    lags = 110; %upper limit to range of lags used, (1 to 'lags')
    readfile = false; %read file or generate wave
        filename = 'File.wav';
        %if readfile=false: 
        frequency = C4; %f to generate (Hz)
        fixedlength = false; %generate fixed seconds
            secondslength = 1; %number of seconds length to generate
            pcount = 6; %number of periods to generate otherwise
    playsounds = false;
    volume = 10; %0-100 please

    %graphing stuff
    xTimeUnits = false; %x units as samples or time
    xLimTo_pcount = false; %only show pcount periods in graph (for many samples)
    graphAnything = true; %graph... anything
        graphY = true; %graph y
        graphRy = true; %graph R(y)
        %graphRyM = true; %graph R(y)M with specified lag
        %graphCF = false; %graph curve fitting
        %graphResample = false; %understand that LP filter is needed
        graphFFT = false; %graph FFT(y)
        doHE = true; %compute the H & E part
        graphHE = true; %graph H & E
%END CVARS

resample_rate = upsample_rate / downsample_rate;
%wave input section
if(readfile == false)
    f = frequency; %known frequency
    p = f^-1; %period = 1/f
    phase_samples = 0; %how many samples to shift the wave forward
    a = samplerate * secondslength; %number of samples needed for X seconds
    c = samplerate/f * pcount; %number of samples to take for x
    samplesperperiod = samplerate/f;
    if(fixedlength)
        if(c < a)
            c = a; %if c<a, make it enough samples for X seconds
        end
    end
    fprintf("Generated wave properties:\nfreq = %.2f Hz\nperiod = %.5g s\nsamples/period = %.3f samples\n",f,p,samplesperperiod);
    fprintf('Resampling Rate = %d up / %d down, RS freq = %.2f Hz, RS samples/period = %.3f samples\n',upsample_rate,downsample_rate,f/resample_rate,samplesperperiod*resample_rate); %print detected pitch
    fetchnote_print(f);
    x = 0:c-1; %x = c # of samples
    %x = (x ./ c)  .* (p * pcount);  %(( x = 0:1) * period length) * # of periods
    %changed to x = samples only
    if(fixedlength)
        x = x(1:a); 
    end
    y = sin(2*pi*f.*x/samplerate + phase_samples);%,1/2); %make a wave y, ex: y = sin(w*x + O), w = 2*pi*f
    %base period is samplerate
else
    y = audioread(filename); %have y be a sound file of X samples, may not know the frequency?
    x = 1:length(y); %x = X samples to match y sample count
    f = 0;
    samplesperperiod = length(x);
    fprintf("number of samples in file = %d\n",length(x));
end
%prep stuff
%y_normalised = normalize(y); %comment this out if licensing complaint
if(~isrow(y))
    y = y';
end
resample_rate_s = num2str(upsample_rate) + "/" + num2str(downsample_rate);
%[x_rs, y_rs] = resampleY(y,upsample_rate,downsample_rate,false); %resample
[x_rs, y_rs] = resampleY(y,upsample_rate,downsample_rate,true); %resample+LPF
lag_s = num2str(lags);
spp = samplesperperiod;
if(volume < 1) %don't play sound if volume is low
    playsounds = false;
end
%graphing stuff
if(xTimeUnits)
    xTimeUnits_modifier = 1/samplerate;
else
    xTimeUnits_modifier = 1;
end
graphY = graphAnything && graphY;
graphRy = graphAnything && graphRy;
%graphRyM = graphAnything && graphRyM;
% graphCF = graphAnything && graphCF;
graphHE = graphAnything && graphHE;
graphFFT = graphAnything && graphFFT;
SubplotSize = 2;
screensize = get(groot, 'ScreenSize'); screensize = [screensize(3) screensize(4)]; % get screensize
screencenter = screensize./2;
ylimits = [-1.25 1.25];
graphsize = [560 420]; %size of figure windows
graphcolors = [[0 0.2 0.6]; [0.4 0.4 0]; [0 0.6 0.2]];
MarkerSizeSmall = 4;
MarkerSizeNormal = 6;
MarkerSizeBig = 10;
MarkerSizeResampled = MarkerSizeBig;
if(playsounds)
    soundedStr = "(sounded) ";
else
    soundedStr = "";
end
fprintf('\n');

%graph wave section
xlim_x = (length(y)-1) * (xTimeUnits_modifier); %x axis limits
xlim_x_rs = xlim_x * resample_rate; %for resampled
if(xLimTo_pcount) %limit to pcount
    xlim_x = (samplesperperiod*pcount) * (xTimeUnits_modifier);
    if(xlim_x > length(y))
        xlim_x = length(y);
    end
end
xlim_x_rs = xlim_x * resample_rate;
if(graphY)
    yplots = figure('Name','Input waves');
    gs = graphsize; %alias, below sets position
    set(yplots,'Position', [screencenter(1)-gs(1)/2-gs(1) screencenter(2)-gs(2)/2+100 gs(1) gs(2)]);
    subplot(SubplotSize,1,1) %normal
        plot(x.*(xTimeUnits_modifier),y,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
        xlim([0 xlim_x])
        title("y"  + soundedStr)
        ylim(ylimits)
    subplot(SubplotSize,1,SubplotSize) %w/ filter
        plot(x_rs.*(xTimeUnits_modifier),y_rs,'.-.','MarkerSize',MarkerSizeResampled,'Color',graphcolors(3,:))
        xlim([0 xlim_x_rs])
        if(playsounds)
            soundedStr = "(sounded 2nd) ";
        end
        title("y rs " + soundedStr + resample_rate_s)
        ylim(ylimits)
end
y_fft = abs(fft(y)/length(y)); %fft(y)
y_fft_lefthalf = y_fft(1:round(length(y_fft)/2)+1);
if(graphFFT)
    FFTplots = figure('Name','Fourier Transform');
    gs = graphsize; %alias, below sets position
    set(FFTplots,'Position', [screencenter(1)-gs(1)/2 screencenter(2)-gs(2)/2-300 gs(1) gs(2)]);
    subplot(1,1,1) %FFT(y)
     plot(samplerate * (0:length(y_fft_lefthalf)-1),y_fft_lefthalf,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
     title("fft(y)")
     
end

%no lag autocorr part
%xlim_R = xlim_x;%int32(1000); % how far to scale the x axis
fprintf("R(y):\n");
    [R_x, R, R_f, R_fQI] = handle_R(y,length(y),readfile,f,samplerate,resample_rate,0,0);
    R_p = R_f^-1; %periods
    R_pQI = R_fQI^-1;
    R_spp = R_p * samplerate; %samples per period
    R_sppQI = R_pQI * samplerate;
%fprintf("R(y rs):\n");
    %[R_rs_x, R_rs, R_rs_f, R_rs_fQI] = handle_R(y_rs,length(y_rs)-1,readfile,f,samplerate,resample_rate,1,0);
xlim_x_R = (length(R)-1) * (xTimeUnits_modifier);
if(xLimTo_pcount) %limit to pcount
    xlim_x_R = (samplesperperiod*pcount*2) * (xTimeUnits_modifier);
    if(xlim_x_R > length(R))
        xlim_x_R = length(R);
    end
end
xlim_x_R_rs = xlim_x_R*resample_rate; %for resampled
if(graphRy)
  Rplots = figure('Name','Full Autocorrelation');
  gs = graphsize; %alias, below sets position
  set(Rplots,'Position', [screencenter(1)-gs(1)/2 screencenter(2)-gs(2)/2+100 gs(1) gs(2)]);
  subplot(1,1,1) %R(y)
     plot(R_x.*(xTimeUnits_modifier),normalize(R),'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
     title("R(y)")
     xlim([0 xlim_x_R])
     ylim(ylimits)
%   subplot(SubplotSize,1,SubplotSize) %R(y_rs)
% 	 plot(R_rs_x.*(xTimeUnits_modifier),normalize(R_rs),'.-.','MarkerSize',MarkerSizeResampled,'Color',graphcolors(3,:));
% 	 title("R(y rs) " + resample_rate_s);
%      xlim([0 xlim_x_R_rs])
%      ylim(ylimits)
end

%_M_anually specified lag autocorr part, no longer used
fprintf("Lags used: [0-%d]\n",lags-1);
% fprintf("R(y)M:\n");
%     [RM_x, RM, RM_f, RM_fQI] = handle_R(y,lags,readfile,f,samplerate,resample_rate,0,0);
%     xlim_x_RM = (length(RM)-1) * (xTimeUnits_modifier);
% fprintf("R(y rs)M:\n");
%     [RM_rs_x, RM_rs, RM_rs_f, RM_rs_fQI] = handle_R(y_rs,lags,readfile,f,samplerate,resample_rate,1,0);
%     xlim_x_RM_rs = (length(RM_rs)-1) * (xTimeUnits_modifier);
% if(graphRyM)
%   RMplots = figure('Name','Manual Autocorrelations');
%   gs = graphsize; %alias, below sets position
%   set(RMplots,'Position', [screencenter(1)-gs(1)/2+gs(1) screencenter(2)-gs(2)/2+100 gs(1) gs(2)]);
% %   subplot(SubplotSize,1,1) %RM(y)
% %       plot(RM_x.*(xTimeUnits_modifier),RM,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
% %       title("R(y)M" + " L=" + num2str(lags))
% %       axis tight
% %       xlim([0 xlim_x_RM])
% %       %ylim(ylimits)
%   subplot(1,1,1) %RM(y_rs)
%       plot(RM_rs_x.*(xTimeUnits_modifier),RM_rs,'.-.','MarkerSize',MarkerSizeResampled,'Color',graphcolors(3,:));
%       title("R(y rs)M " + resample_rate_s + " L=" + num2str(lags));
%       axis tight
%       xlim([0 xlim_x_RM_rs])
%       %ylim(ylimits)
% end
    
if(doHE) %H,E equations part
    %E: energy equation
    EM_rs = calcE(y_rs,lags);
    %H is similar to RM except it only Acorr's a range of points 0-L, not every point 
    HM_rs = calcH(y_rs,lags);
    %handle_R(HM_rs,lags,readfile,f,samplerate,resample_rate,1,1);
    HEMcompare_rs = EM_rs(1:length(HM_rs)) - 2.*HM_rs;
    %these printing false again for low numbers
    %idk why hopefully fixable. process still worked regardless...
    if(all(HEMcompare_rs >= 0)) %E >= 2H
        fprintf("EM,HM comparison true.\n");
    else
        fprintf("! EM,HM comparison false, failed. fix this somehow?\n");
        disp(find(HEMcompare_rs < 0));
    end
    epsilon = 0.4; %arbitrary small value
    fprintf("epsilon = %.2f\n",epsilon);
    %E - 2H > epsilon * E
    HEMcompare_rs_epsilon = HEMcompare_rs <= (epsilon * EM_rs(1:length(HEMcompare_rs)));
    %Lmin
    [~, Lmin1_x_rs] = findpeaks(-HEMcompare_rs(3:end)); %only check 2:lags
    if(isempty(Lmin1_x_rs)) %error check
        fprintf("! No peak found for Lmin1_x_rs. Set as right bound (more data needed?)...\n")
        Lmin1_x_rs = length(HEMcompare_rs);
    else
        Lmin1_x_rs = Lmin1_x_rs(1)-1+2; %first valley, add 2 offset from above
    end
    [~, Lmin2_x_rs] = findpeaks(-HEMcompare_rs(Lmin1_x_rs+1:end));
    if(isempty(Lmin2_x_rs))
        fprintf("! No peak found for Lmin2_x_rs. Set as Lmin1_x_rs\n")
        Lmin2_x_rs = Lmin1_x_rs;
    else
        Lmin2_x_rs = Lmin2_x_rs(1)-1 + Lmin1_x_rs(1); %2nd valley
    end
    Lmin1_x = (Lmin1_x_rs)/resample_rate;
    Lmin2_x = (Lmin2_x_rs)/resample_rate;
    %add code for no peaks found later (error check)
    %recalcuate E&H for original function (not rs)
    E = calcE(y,Lmin2_x+4+4); %prob E too
    H = calcH(y,Lmin2_x+4+4); %H coming out wrong cuz of lag stuff, fix
    HEcompare = E(1:length(E)) - 2.*H;
    if(all(HEcompare >= 0)) %E >= 2H
        fprintf("E,H comparison true.\n");
    else
        fprintf("! E,H comparison false, failed. fix this somehow?\n");
        disp(find(HEcompare < 0));
    end
    HEcompare_epsilon = HEcompare <= (epsilon * E(1:length(HEcompare)));
    Lmin = Lmin1_x; %select final Lmin to use in H
    %can give wrong harmonic? try comparing local minimums instead
    %fix E, H comparison stuff first
    if(Lmin1_x_rs(1) < lags/2 && HEcompare(Lmin2_x) < HEcompare(Lmin1_x))
        Lmin = Lmin2_x;
    end
    %16 point range is temporary, will use 8 (+- 4) later
    Lmin_points = (Lmin+1)-3:(Lmin+1)+4; %get range of 8 points
    Lmin_points_ext = (Lmin+1)-7:(Lmin+1)+8; %range of 16 points
    Lmin_points_ext(Lmin_points_ext > length(HEcompare)) = length(HEcompare);
    %determine HE_f:
    [~, peak] = findpeaks(-HEcompare(Lmin_points_ext));
    if(isempty(peak)) %error handle
        fprintf("! No peak found for HEcompare Lmin set evaluation. Set as right bound...\n")
        peak = Lmin_points_ext(end);
    else
        peak = Lmin_points_ext(peak(1));
    end
    peaknfriends = [peak-1,peak,peak+1]; %array of 3 points
    if(~isrow(HEcompare)) %to row vector
        HEcompare = HEcompare';
    end
    %shift left with -1 to get samples per period (from QI x):
    HE_spp = QInterp_peak(peaknfriends,-HEcompare(peaknfriends))-1; %QI
    HE_p = HE_spp/samplerate;
    HE_f = samplerate/(HE_spp); %HE detected f
    fprintf('\nHE p detected pitch: HE_f = ') %print out HE_f results
    report_ptof(HE_f,readfile,f)
    
    if(graphHE)
    E_pmult = 2; %2*x range for E
    EM_lag = lags*E_pmult; %2L
    EH_plot_x = Lmin_points_ext-1;
    HEplots = figure('Name','E and H functions');
    gs = graphsize; %alias, below sets position
    set(HEplots,'Position', [screencenter(1)-gs(1)/2+gs(1) screencenter(2)-gs(2)/2+100 gs(1) gs(2)]);
    subplot(2,4,1) %E(y)
      plot((0:2:length(E)*E_pmult-2).*(xTimeUnits_modifier),E,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
      title("E(y)M" + " 2L=" + num2str(EM_lag))
      axis tight
      hold on
      plot(EH_plot_x.*E_pmult,E(Lmin_points_ext),'.-','MarkerSize',MarkerSizeNormal+1,'Color',[0.9,0,0]);
      plot((peak-1)*E_pmult,E(peak),'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
    subplot(2,4,2) %H(y)
      plot((0:length(H)-1).*(xTimeUnits_modifier),H,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
      title("H(y)M")
      axis tight
      hold on
      plot(EH_plot_x,H(Lmin_points_ext),'.-','MarkerSize',MarkerSizeNormal+1,'Color',[0.9,0,0]);
      plot(peak-1,H(peak),'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
    subplot(2,4,3) %HEdiffs(y)
      plot((0:length(H)-1).*(xTimeUnits_modifier),HEcompare,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
      title("E-2H(y)M")
      axis tight
      hold on
      plot(EH_plot_x,HEcompare(Lmin_points_ext),'.-','MarkerSize',MarkerSizeNormal+1,'Color',[0.9,0,0]);
      plot(peak-1,HEcompare(peak),'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
    subplot(2,4,4) %HEdiffs_epsilon(y)
      plot((0:length(H)-1).*(xTimeUnits_modifier),HEcompare_epsilon,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
      title("E-2H<=eE(y)M")
      axis tight
      hold on
      plot(EH_plot_x,HEcompare_epsilon(Lmin_points_ext),'.-','MarkerSize',MarkerSizeNormal+1,'Color',[0.9,0,0]);
      plot(peak-1,HEcompare_epsilon(peak),'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
    subplot(2,4,5) %E(y_rs)
      plot((0:E_pmult:EM_lag-E_pmult).*(xTimeUnits_modifier),EM_rs,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
      title("E(y rs)M" + " 2L=" + num2str(EM_lag))
      axis tight
      hold on
      plot([Lmin1_x_rs*E_pmult Lmin2_x_rs*E_pmult],[EM_rs(Lmin1_x_rs+1) EM_rs(Lmin2_x_rs+1)],'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
    subplot(2,4,6) %H(y_rs)
      plot((0:lags-1).*(xTimeUnits_modifier),HM_rs,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
      title("H(y rs)M")
      axis tight
      hold on
      plot([Lmin1_x_rs Lmin2_x_rs],[HM_rs(Lmin1_x_rs+1) HM_rs(Lmin2_x_rs+1)],'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
    subplot(2,4,7) %HEdiffs(y_rs)
      plot((0:lags-1).*(xTimeUnits_modifier),HEMcompare_rs,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
      title("E-2H(y rs)M")
      axis tight
      hold on
      plot([Lmin1_x_rs Lmin2_x_rs],[HEMcompare_rs(Lmin1_x_rs+1) HEMcompare_rs(Lmin2_x_rs+1)],'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
    subplot(2,4,8) %HEdiffs_epsilon(y_rs)
      plot((0:lags-1).*(xTimeUnits_modifier),HEMcompare_rs_epsilon,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
      title("E-2H<=eE(y rs)M")
      axis tight
      hold on
      plot([Lmin1_x_rs Lmin2_x_rs],[HEMcompare_rs_epsilon(Lmin1_x_rs+1) HEMcompare_rs_epsilon(Lmin2_x_rs+1)],'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
    end
end

%Curve Fitting part, not viable cuz the peak doesn't change
% if(graphCF)
%   fprintf("-\nCurve Fitting\n-\n");
%   fittype = 'cubicinterp'; %what type of curve to fit
%   Gfitplots = figure('Name','RM fits');
%   movegui(Gfitplots,'southwest')
%   subplot(SubplotSize,1,1)
%       Gfit = fit(RM_x',RM',fittype);
%       plot(Gfit,RM_x,RM);
%       title("R(y)M fit")
%       ylim(ylimits)
%   subplot(SubplotSize,1,SubplotSize)
%       Gfit_rs = fit(RM_rs_x',RM_rs',fittype);
%       plot(Gfit_rs,RM_rs_x,RM_rs);
%       title("R(y rs)M fit")
%       ylim(ylimits)
%   GfitplotsD = figure('Name','RM fits discrete');
%   movegui(GfitplotsD,'south')
%   subplot(SubplotSize,1,1)
%     GfitY = feval(Gfit,0:110)';
%     fprintf("R(y)CF:\n");
%     handle_R(GfitY,lags,readfile,f,samplerate,resample_rate,0,1);
%     plot(GfitY)
%     ylim(ylimits)
%   subplot(SubplotSize,1,SubplotSize)
%     Gfit_rsY = feval(Gfit_rs,0:110)';
%     fprintf("R(y rs)CF:\n");
%     handle_R(Gfit_rsY,lags,readfile,f,samplerate,resample_rate,1,1);
%     plot(Gfit_rsY,'.-.','MarkerSize',MarkerSizeResampled,'Color',graphcolors(3,:))
%     ylim(ylimits)
% end

%playing sounds through resampling
[~, ~, target_f] = fetchnote(R_fQI); %find needed f
samplerate_tuned_ratio = target_f/R_fQI; %ratio between needed/found f
samplerate_tuned = round(samplerate * samplerate_tuned_ratio); %new S.R.
%higher SR = higher pitch, lower SR = lower pitch
y_tuned = y; %with diff samplerates, wave will be the same

%uncomment these if you are using your own pitch-shifted wave
%y_tuned = YOUR_WAVE_TO_PLAY; %put wave here
%samplerate_tuned = samplerate; %restore samplerate

if(playsounds) %plays the sounds
    volpeakpoint = 100/volume;
    y_toplay = [volpeakpoint y];
    fprintf("Playing y...\n");
    soundsc(y_toplay,(samplerate))
    pause
    clear sound, clear sounds;
    y_toplay = [volpeakpoint y_tuned];
    fprintf("Playing y (pitch corrected w/ samplerate adjustment from R_fQI)...\n\n");
    soundsc(y_toplay,(samplerate_tuned)) 
    %pause
end
%clear sound, clear sounds;

%clean up, currently does nothing tho
if(graphAnything && (graphY || graphRy || graphFFT || graphHE))
    %pause
    %close all;
end
%clear sound, clear soundsc

%FUNCTIONS
function [x, y] = resampleY(y,upsample_rate,downsample_rate,useLPfilter)
%plot and resample/(decimate/interpolate).
    if(upsample_rate > 1) %upsample/interp
        if(useLPfilter)
            y = interp(y, upsample_rate);
        else
            y = upsample(y, upsample_rate);
        end
    end
    if(downsample_rate > 1) %downsample/decimate
        if(useLPfilter)
            y = decimate(y, downsample_rate);
        else
            y = downsample(y, downsample_rate);
        end
    end
    x = 0:length(y)-1;
end

function [x, R, fR, fR_QI] = handle_R(y,lags,readfile,f,samplerate,resample_rate,isResampled,noAcorr)
%do all the autocorrelation stuff
%outputs x, R, guessed frequency fR, and ^2 interpolated fR
global count;
%defaults:
    fR = 0;
    fR_QI = 0;
    if(lags > length(y))
        lags = length(y);
    end
    if(~noAcorr)%for passing an R directly
        %R = autocorr(y,lags-1); %divides out variance
        R = autoc(y,lags); %doesn't divide out variance
    else
        R = y;
    end
    if(~isResampled)
        resample_rate = 1;
    end
    x = 0:length(R)-1; %create x values
    [~,locs] = findpeaks(R); %find the peaks' x locations
    if(isempty(locs)) %check if no 2nd peak (invalid R reading)
        locs(1) = length(R);
        fprintf("! No valid R peak detected. Set as right bound (need more samples?)\n");%,locs(1));
    else %if 2nd peak found
        if(~isrow(locs)) %to row vector
            locs = locs';
        end
        locs = locs-1;%shift left cuz matlab
        meanp = mean(diff([0 locs])); %find the mean period, fixed calculation
        peaks_str = strtrim(sprintf('%d ', locs));
        if(length(peaks_str) > 30) %truncate if too long
            peaks_str = peaks_str(1:30) + " ...";
        end
        fprintf("Peaks at x: [%s], AVG Peak p = %.2f\n",peaks_str,meanp); %show peak x's
        fR = samplerate/(meanp)*resample_rate; %R(avg peak) detected f
        %quadratic interpolation stuff
        firstpeak = locs(1); %uses only the first peak
        peaknfriends = [firstpeak-1,firstpeak,firstpeak+1]; %array of 3 points
        if(~isrow(R)) %to row vector
            R = R';
        end
        %locs was shifted left, shift right with +1 to reverse:
        xweneed = QInterp_peak(peaknfriends,R(peaknfriends+1)); %get float x
        fR_QI = samplerate/(xweneed)*resample_rate; %QI detected f
        fprintf('AVG Peak p detected pitch: R_f = ') %print out f results
        report_ptof(fR,readfile,f)
        fprintf('Q. interpolated p detected pitch: R_fQI = ')
        report_ptof(fR_QI,readfile,f)
    end
    fprintf("\n");
end

function report_ptof(f,readfile,f_generated)
%reports frequency results. generatedf used if not reading a file
    %fprintf('detected pitch: f = %.2f Hz, p = %f s, samples/period = %.5g samples\n',fR,1/fR,samplerate/fR*resample_rate);
    fprintf('%.2f Hz, ',f)
    if(~readfile) %error from generated pitch
      %fprintf('Error from generated f = %+.2f Hz, cents = %+.1f\n', fR_QInter - f, log(double(fR_QInter/f))/log(2)*1200);
      errorc = log(double(f/f_generated))/log(2)*1200;
      fprintf('Error (cents) = %+.1f', errorc);
      if(errorc > 40)
          fprintf(', ouch');
      end
      fprintf('\n');
    end
    fetchnote_print(f); %find closest note
end

function Qx = QInterp_peak(x_3, y_3)
%finds peak of quadratic interpolation on 3 points
    p = polyfit(-1:1,y_3,2); %fit ^2, a = y, x = -1 to 1 (3 points effectively)
    %pval = polyval(p,-1:1); %show the 3 points
    pvalext = polyval(p,linspace(-1,1,500)); %full evaluation, spacing should be large number
    [~,plocs] = findpeaks(pvalext); %find the peak's x value
    if(isempty(plocs)) %error check
        plocs(1) = round(max(pvalext));
        fprintf("! No peak found in QI function. Set as max, bounds included\n");
    end
    Qx = plocs(1);
    %found x value of the peak between 0-500, turn back into 0-2:
    Qx = (Qx - 1) * (x_3(3) - x_3(1))/499 + x_3(1); %needed x value (non-integer)
end

function a = autoc(y,lags)
%autocorrelation with lag range 0 to (lags-1), doesn't divide out variance
    if(length(y) < lags*2+1)
        y = [y zeros(1,lags*2+1-length(y))];
    end
    y = y(1:lags*2+1); %use only 2L points
    y_len = length(y);
    n = 2^nextpow2(2*y_len-1); %fft runs fastest on 2^n sized data
    a = ifft(fft(y,n) .* conj(fft(y,n))); %acorr = iF(F * F*) = iF(S)
    a = a(1:lags); %only need up to L points
end

function E = calcE(y,lags)
%energy (no lag), with input range (0 to (lags-1)) * 2
    E_pmult = 2;% * L (= 2L used)
    E_lag = lags*E_pmult; %2L
    endy = E_lag*2;
    if(endy > length(y))
        endy = length(y);
    end
    y = y(1:endy);
    y = [y zeros(1,E_lag)+2]; %prep to stop array index errors
    E = zeros(1,lags); %prep, only doing limited lags
    y_squared = y .* y; %y^2 = y * y
%     for k = 2:lags %calculate E for all sets of 2L's, L=0 (k=1) is always 0
%         E(k) = sum(yE(1:E_pmult*k-E_pmult).^2); %sum up squares in range [0,2L]
%     end
    for c = 1:length(E)
        E(c) = sum(y_squared(1:c*E_pmult)); %only sum c points, not all points
    end
end

function H = calcH(y,lags)
%autocorrelation with lag range 0 to (lags-1), doesn't divide out variance
%only does up to L points for each lag
%check this function again sometime
    endy = lags*2;
    if(endy > length(y))
        endy = length(y);
    end
    y = y(1:endy); %use only y(0 to lags*2-1)
    y = [y zeros(1,lags+2)]; %append extra zeros to avoid index errors
    H = zeros(1,lags); %set up with zeros
    for c = 1:length(H) %check this later
        %m = y_i * y_i-L = y(t) * y(t+c)
        m = y .* [y(c:end) zeros(1,c-1)]; %y * y shifted left by L
        H(c) = sum(m(1:c)); %only sum c points, not all points
    end

    if(false) %attempt at FFT optimization, do later
    %y = [y zeros(1,lags+2)];
    x = [y zeros(1,lags+2)];
    y = y(1:lags);
    y_len = length(y);
    H = zeros(1,lags);
    for c = 1:lags %todo: fix me
        mm = y(c:c);
        nn = y(1:c-(c-1));
        m = mm .* nn;
        
        H(c) = sum(m);
        d = autoc(y(c:end),c-1);
        length(d);
        
        m = x .* [x(1+c-1:end) zeros(1,c-1)];
        b(c) = sum(m(1:c));
    end
    
    figure(69)
    subplot(3,1,1)
    plot(0:lags-1,H)
    subplot(3,1,2)
    plot(0:lags-1,b)
    subplot(3,1,3)
    plot(0:lags-1,H - b)
    H = b;
    end
end

function f = constrain(y,bl,bu)
%return bounded value clipped between bl and bu
	f = min(max(y,bl),bu);
end

function f = normalize(y)
%normalize a function to range [-1, 1]
    f = 2*(y - min(y))./(max(y) - min(y)) - 1;
end

function [note_name, error_cents, target_f] = fetchnote(x)
%return closest 12EDO note n of frequency x (Hz) and error e (cents)
%target_f is the frequency of the correct note
    global A4freq;
    notes = ["C" "C#/Db" "D" "D#/Eb" "E" "F" "F#/Gb" "G" "G#/Ab" "A" "A#/Bb" "B"];
    Czero = 2^(-(9+12*4)/12)*A4freq; %C0
    h = round(12*log(x/Czero)/log(2)); %convert to cents/100
    x_octave = idivide(int32(h),int32(12)); %get octave number
    note_index = mod(h,12)+1; %get array index
    
    x_c = log(x/Czero)/log(2) * 1200; %Hz to cents (logarithmic unit)
    x_cm = mod(x_c,1200); %octave reduce
    error = mod(x_cm,100); %calculate error in cents
    if(error >= 50)
        error = error - 100;
    end

    error_cents = error; %error from closest pitch in cents
    note_name = strcat(notes(note_index),num2str(x_octave)); %select from 'notes' + octave number
    target_f = 2^((note_index-1+12*double(x_octave))/12)*Czero;
end

function fetchnote_print(x)
%print out fetchnote information
    [n, e, f] = fetchnote(x);
    fprintf("Closest note = %s (%.2f Hz) Error (cents) = %+.1f\n",n,f,e);
end

%for MATLAB; ignore some warnings
    %#ok<*DEFNU>
    %#ok<*UNRCH>
    %#ok<*NASGU>
    %#ok<*NUSED>