%snr design autotune practice
%Team: Ashlin G, Alex G, Brandon Y
clc, clear variables, clear sound, clear sounds, close all;

%CONTROL VARIABLES
    %constants
    samplerate = 44100; %sample rate
    global A4freq; A4freq = 440; %Freq (Hz) for note A4
    C4 = A4freq * 2^(-9/12); %middle c (C4) from A4
    C0 = 2^(-(9+12*4)/12)*A4freq; %C0

    %control
    downsample_rate = 8; %down/up sample by these factors, new pitch = f * d/u
    upsample_rate = 1;
    lags = 110; %upper limit to range of lags used, (1 to 'lags')
    epsilon = 0.4; %arbitrary small value
    readfile = false; %read file or generate wave
        filename = 'File.wav'; %must be MONO, no stereo
        frequency = C4; %f to generate (Hz)
        fixedlength = true; %generate fixed seconds
            secondslength = 0.1; %number of seconds length to generate
            pcount = 30; %number of periods to generate otherwise
    playsounds = false;
    volume = 10; %0-100 please

    %graphing stuff
    xTimeUnits = false; %x units as samples or time
    xLimTo_pcount = false; %only show pcount periods in graph (for many samples)
    graphAnything = true; %graph... anything
        graphY = true; %graph y
        graphRy = true; %graph R(y)
        graphFFT = false; %graph FFT(y)
        doHEcont = true; %compute H & E continuous
        graphHEcont = true; %graph H & E continuous
%END CVARS

resample_rate = upsample_rate / downsample_rate;
resample_period = downsample_rate / upsample_rate;
RSR = resample_rate;
RSP = resample_period;
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
    fprintf('Resampling Rate = %d up / %d down, RS freq = %.2f Hz, RS samples/period = %.3f samples\n',upsample_rate,downsample_rate,f*resample_period,samplesperperiod*resample_rate); %print detected pitch
    fetchnote_print(f);
    x = 0:c-1; %x = c # of samples
    %x = (x ./ c)  .* (p * pcount);  %(( x = 0:1) * period length) * # of periods
    %changed to x = samples only
    if(fixedlength)
        x = x(1:a); 
    end
    %frequency modifier (modulator?)
    fmod = [ones(1,length(x)/2) linspace(1,4,length(x)/2)]; %linear
    %fmod_f = 20; %fmod oscillation f (Hz)
    %fmod = 1 + 0.5*sin(2*pi*fmod_f.*x/samplerate); %oscillate sine (can implement vibrato)
    [~, fmod_rs] = resampleY(fmod,upsample_rate,downsample_rate,false);
    %input y
    y = sin(2.*pi.*(f.*fmod).*x/samplerate + 2*pi*phase_samples/samplesperperiod); %make a wave y, ex: y = sin(w*x + O), w = 2*pi*f
    %base period is samplerate
    fprintf("Generated sample count = %d (%.3f seconds)\n",length(x),length(x)/samplerate);
else
    [y,SR_file] = audioread(filename); %import y as file, SR_file samplerate
    x = 1:length(y); %x = X samples to match y sample count
    f = 0;
    samplesperperiod = length(x);
    fprintf("File sample count = %d (%.3f seconds), samplerate = %d\n",length(x),length(x)/samplerate,SR_file);
    samplerate = SR_file; %experiment, override set samplerate with file's
end
%prep stuff
%y_normalised = normalize(y); %comment this out if licensing complaint
if(~isrow(y)) %make y a row vector cuz thats what i like
    y = y';
end
resample_rate_s = num2str(upsample_rate) + "/" + num2str(downsample_rate);
sample_rate = samplerate; %alias
resamplerate = resample_rate; %alias
%[x_rs, y_rs] = resampleY(y,upsample_rate,downsample_rate,false); %resample
[x_rs, y_rs] = resampleY(y,upsample_rate,downsample_rate,true); %resample+LPF
lag_s = num2str(lags);
spp = samplesperperiod; %alias
if(volume < 1) %don't play sound if volume is low
    playsounds = false;
end
%graphing stuff
if(xTimeUnits) %seconds instead of samples
    xTimeUnits_modifier = 1/samplerate;
else
    xTimeUnits_modifier = 1;
end
graphY = graphAnything && graphY;
graphRy = graphAnything && graphRy;
graphHEcont = graphAnything && graphHEcont;
graphFFT = graphAnything && graphFFT;
SubplotSize = 2;
screensize = get(groot, 'ScreenSize'); screensize = [screensize(3) screensize(4)]; % get screensize
screencenter = screensize./2;
ylimits = [-1.2 1.2];
graphsize = [560 420]; %size of figure windows
graphcolors = [[0 0.2 0.6]; [0.4 0.4 0]; [0 0.6 0.2]]; %some colors
MarkerSizeTiny = 2;
MarkerSizeSmall = 4;
MarkerSizeNormal = 6;
MarkerSizeBig = 10;
MarkerSizeResampled = 8;
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
        plot(x.*(xTimeUnits_modifier),y,'.-','MarkerSize',MarkerSizeTiny,'Color',graphcolors(1,:));
        xlim([0 xlim_x])
        title("y"  + soundedStr)
        ylim(ylimits)
    subplot(SubplotSize,1,SubplotSize) %w/ filter
        plot(x_rs.*(xTimeUnits_modifier),y_rs,'.-.','MarkerSize',MarkerSizeSmall,'Color',graphcolors(3,:))
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
     plot((0:length(y_fft_lefthalf)-1),y_fft_lefthalf,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
     title("fft(y)")
     
end

%full autocorr
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
if(graphRy && ~readfile)
  Rplots = figure('Name','Full Autocorrelation');
  gs = graphsize; %alias, below sets position
  set(Rplots,'Position', [screencenter(1)-gs(1)/2-gs(1) screencenter(2)-gs(2)/2+100-gs(2)/2 gs(1) gs(2)]);
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
figure(yplots) %bring to front

fprintf("Lags used: [0-%d]\n",lags-1);

debugHEgraphs = false; %debug graphs
if(debugHEgraphs)
    HEdebugplots = figure('Name','EH debug plots');
end
usetimesavers = false;%may change detection behavior, default off
if(doHEcont) %H,E continuous
    %prepare arrays. fill with zeros to avoid checks in the loop
    fC_tracker = zeros(1,length(y_rs)+0); %tracks detected frequency
    note_tracker = zeros(1,length(y_rs)+0); %tracks closest note frequency
    error_tracker = zeros(1,5); %track errors without spamming console
    %to eliminate computation time
    E_pmult = 2;% * L (= 2L used)
    E_lag = lags*E_pmult; %2L
    H_lag = lags; %L
    y_Lbuffer = zeros(1,E_lag*resample_period); %holds some y data
    y_Lbuffer_rs = zeros(1,E_lag);
    ynew_rsbuffer_diff = zeros(1,length(y_Lbuffer_rs)); %holds pitch shifting diffs
    ynew_rsbuffer_ratio = zeros(1,length(y_Lbuffer_rs)); %holds pitch shifting ratios
    y_smallbuffer = zeros(1,resample_period);
    EC_rs = zeros(1,H_lag); %H&E of rs and non-rs
    HC_rs = zeros(1,H_lag);
    EC = zeros(1,H_lag*resample_period);
    HC = zeros(1,H_lag*resample_period);
    %execution time saving checks, don't do things if result is the same
    LminC_last = -1;
    Lmin1C_rs_last = -1;
    Lmin2C_rs_last = -1;
    HEC_spp = 0;
    HEC_p = 0;
    HEC_f = 0;
    HEC_f_diff = 0; %difference in f for pitch shift buffer
    %HEC_f_diff_last = 0;
    HEC_f_ratio = 1; %ratio of needed/detected
    %keep indexes as matlab 1:end. only use 0 for period calculation
    s = "";
    i_rs = 1;
    print_samplei = true; %show sample i=#
    tic; %time this process
    for i = 1:length(y) %big ol' for loop, realtime will use a while loop
        y_Lbuffer = [y_Lbuffer(2:end) y(i)];
%         y_smallbuffer = [y_smallbuffer(2:end) y(i)]; %faster
        %iterative E & H
%         for a = 1:H_lag*RSP
%             EHC_newterms = [y_Lbuffer(E_lag*RSP) y_Lbuffer(E_lag*RSP-a+1) y_Lbuffer(E_lag*RSP-2*a+1)];
%             EC(a) = EC(a) + EHC_newterms(1).^2 - EHC_newterms(3).^2;
%             HC(a) = HC(a) + EHC_newterms(1).*EHC_newterms(2) - EHC_newterms(2).*EHC_newterms(3);
%         end
        a = 1:H_lag*RSP;
            EHC_newterm1 = ones(1,H_lag*RSP) .* y_Lbuffer(E_lag*RSP);
            EHC_newterm2 = y_Lbuffer(E_lag*RSP-a+1);
            EHC_newterm3 = y_Lbuffer(E_lag*RSP-2*a+1);
            EC = EC + EHC_newterm1.^2 - EHC_newterm3.^2;
            HC = HC + EHC_newterm1.*EHC_newterm2 - EHC_newterm2.*EHC_newterm3;
        if(mod(i,1*resample_period) == 0) %only check pitch every resampled point for now
            if(print_samplei) %only print in here too
                if(~debugHEgraphs)
                    fprintf(repmat(sprintf('\b'),1,strlength(s)));
                end
                s = "i = " + sprintf('%d',i);
                fprintf("%s",s);
                if(debugHEgraphs)
                    fprintf("\n");
                end
            end
%             y_Lbuffer = [y_Lbuffer(1+resample_period:end) y_smallbuffer]; %do this in here
            y_Lbuffer_rs = [y_Lbuffer_rs(2:end) y_rs(i_rs)];
            %iterative equations from patent
%             for a = 1:H_lag
%                 EHC_rs_newterms = [y_Lbuffer_rs(E_lag) y_Lbuffer_rs(E_lag-a+1) y_Lbuffer_rs(E_lag-2*a+1)];
%                 EC_rs(a) = EC_rs(a) + EHC_rs_newterms(1).^2 - EHC_rs_newterms(3).^2;
%                 HC_rs(a) = HC_rs(a) + EHC_rs_newterms(1).*EHC_rs_newterms(2) - EHC_rs_newterms(2).*EHC_rs_newterms(3);
%             end
            a = 1:H_lag;
                EHC_rs_newterm1 = ones(1,H_lag) .* y_Lbuffer_rs(E_lag);
                EHC_rs_newterm2 = y_Lbuffer_rs(E_lag-a+1);
                EHC_rs_newterm3 = y_Lbuffer_rs(E_lag-2*a+1);
                EC_rs(a) = EC_rs(a) + EHC_rs_newterm1.^2 - EHC_rs_newterm3.^2;
                HC_rs(a) = HC_rs(a) + EHC_rs_newterm1.*EHC_rs_newterm2 - EHC_rs_newterm2.*EHC_rs_newterm3;
            
            HECcompare_rs = EC_rs - 2.*HC_rs; %E - 2H
            if(~all(HECcompare_rs >= 0)) %E_rs >= 2H_rs?
                error_tracker(1) = error_tracker(1) + 1; %E < 2H detected
                %doesn't seem to matter
            end
            HECcompare_rs_epsilon = HECcompare_rs <= (epsilon * EC_rs);
            %make unqualifying points have high amplitude, excluding them from Lmin
            HECcompare_rs_epsilon = HECcompare_rs.*HECcompare_rs_epsilon + ((max(HECcompare_rs)+1) .* ~HECcompare_rs_epsilon);
            %Lmin
            [~, Lmin1C_rs_xarray] = findpeaks_fast(-HECcompare_rs_epsilon(3:end));
            Lmin1C_rs = Lmin1C_rs_xarray;
            if(isempty(Lmin1C_rs_xarray)) %No valley found for Lmin1C_rs. (bad)
                error_tracker(2) = error_tracker(2) + 1;
                Lmin1C_rs = length(HECcompare_rs_epsilon); %use right bound
            else
                Lmin1C_rs = Lmin1C_rs_xarray(1)+2; %first valley, add 2 offset from above
            end
            if(isempty(HECcompare_rs_epsilon(Lmin1C_rs+1:end)))
                Lmin2C_rs = Lmin1C_rs;
            else
                %[~, Lmin2C_rs_xarray] = findpeaks(-HECcompare_rs_epsilon(Lmin1C_rs+1:end));
                %Lmin2C_rs_xarray = Lmin1C_rs_xarray(2:end)-Lmin1C_rs+2;
                Lmin2C_rs_xarray = Lmin1C_rs_xarray(2:end);
                if(isempty(Lmin2C_rs_xarray)) %No valley found for Lmin2C_rs.
                    error_tracker(3) = error_tracker(3) + 1;
                    Lmin2C_rs = Lmin1C_rs; %No Lmin2
                else
                    Lmin2C_rs = Lmin2C_rs_xarray(1) + 2; %2nd valley
                end
            end 
            detectNewLmin = true;
            %don't do anything further if Lmin won't change
            if(usetimesavers && Lmin1C_rs == Lmin1C_rs_last && Lmin2C_rs == Lmin2C_rs_last)
                detectNewLmin = false;
            end
            Lmin1C_rs_last = Lmin1C_rs;
            Lmin2C_rs_last = Lmin2C_rs;
            if(detectNewLmin) %transfer to non-resampled H & E
                Lmin1C = (Lmin1C_rs)*resample_period;
                Lmin2C = (Lmin2C_rs)*RSP;
%                 EC = [calcE(y_Lbuffer,Lmin1C+RSP,Lmin1C-RSP) zeros(1,Lmin2C+RSP-Lmin1C-RSP)] + calcE(y_Lbuffer,Lmin2C+RSP,Lmin2C-RSP); %prob E too
%                 HC = [calcH(y_Lbuffer,Lmin1C+RSP,Lmin1C-RSP) zeros(1,Lmin2C+RSP-Lmin1C-RSP)] + calcH(y_Lbuffer,Lmin2C+RSP,Lmin2C-RSP); %H coming out wrong cuz of lag stuff, fix
                HECcompare = EC - 2.*HC;
                if(~all(HECcompare >= 0)) %E >= 2H?
                    error_tracker(4) = error_tracker(4) + 1;
                    %also doesn't seem to matter
                end
                LminC = Lmin1C; %select final Lmin to use in H
                %can give wrong harmonic? try comparing local minimums instead
                %check inequality signs
                if(Lmin1C_rs < lags/2 && HECcompare(Lmin2C) < HECcompare(Lmin1C))
                    LminC = Lmin2C;
                end
                detectNewPitch = true;
                if(usetimesavers && LminC == LminC_last) %don't run Qinterp if the pitch didn't change
                    detectNewPitch = false;
                end
                LminC_last = LminC;
                LminC_rs = LminC * resample_rate;
                if(detectNewPitch) %detect pitch from valley near Lmin
                    %use large range of points around Lmin to find valley
                    %later write function to find the valley using n points
                    LminC_points = (LminC+1)-resample_period/2+1:(LminC+1)+resample_period/2; %get range of 8 points
                    LminC_points_ext = (LminC+1)-resample_period*2+1:(LminC+1)+resample_period*2; %range of 32 points
                    %try and stop any errors
                    LminC_points_ext(LminC_points_ext > length(HECcompare)-1) = length(HECcompare)-1;
                    %determine HE_f:
                    [~, peakC] = findpeaks_fast(-HECcompare(LminC_points_ext));
                    if(isempty(peakC)) %no valley found for LminC
                        error_tracker(5) = error_tracker(5) + 1;
                        peakC = LminC_points_ext(end)-1; %-1 to stop indexing errors
                    else
                        peakC = LminC_points_ext(peakC(1));
                    end
                    peaknfriendsC = [peakC-1,peakC,peakC+1]; %array of 3 points
                    if(~isrow(HECcompare)) %to row vector
                        HECcompare = HECcompare';
                    end
                    %shift left with -1 to get samples per period (from QI x):
                    HEC_spp = QInterp_peak(peaknfriendsC,-HECcompare(peaknfriendsC))-1; %QI
                    HEC_p = HEC_spp/samplerate;
                    HEC_f = samplerate/(HEC_spp); %detected f
                    f_closest = fetchnote_fastf(HEC_f,A4freq,C0);
                    %if not updated, these will retain previous value
                    %update pitch diff = detected f - closest note f
                    HEC_f_diff = HEC_f - f_closest;
                    %ratio = closest note f / detected f
                    HEC_f_ratio = f_closest / HEC_f;
                    if(HEC_f_ratio < 0.001) %fix zeros
                        HEC_f_ratio = 1;
                    end
                end
            end
%             if(i_rs > length(fC_tracker)) %allocate more space if needed
%                 fC_tracker = [fC_tracker zeros(1,length(fC_tracker))];
%                 note_tracker = [note_tracker zeros(1,length(note_tracker))];
%             end
            ynew_rsbuffer_diff(i_rs) = HEC_f_diff; %add new pitch shift adjustment to buffer
            ynew_rsbuffer_ratio(i_rs) = HEC_f_ratio; %same but ratio instead of difference
            fC_tracker(i_rs) = HEC_f;
            f_for_note_tracker = f_closest; %closest note f
            note_tracker(i_rs) = f_for_note_tracker;
            i_rs = i_rs + 1;
            %debug graphs
            if(debugHEgraphs)
                fprintf("detected f = %.2f Hz",HEC_f);
                if(~readfile)
                    fprintf(", real f = %.2f Hz",f.*fmod_rs(i_rs));
                end
                fprintf("\n")
                %EHC_plot_x =(LminC_points_ext-1).*(xTimeUnits_modifier);
                EHC_plot_x = 0:length(EC)-1;
                figure(HEdebugplots)
                clf
                subplot(2,4,1) %E(y)
                  plot(EHC_plot_x,EC,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
                  title("E")
                  axis tight
                  hold on
                  plot(peakC-1,EC(peakC),'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
                subplot(2,4,2) %H(y)
                  plot(EHC_plot_x,HC,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
                  title("H")
                  axis tight
                  hold on
                  plot(peakC-1,HC(peakC),'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
                subplot(2,4,3) %HEdiffs(y)
                  plot(EHC_plot_x,HECcompare,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
                  title("E-2H")
                  axis tight
                  hold on
                  plot(peakC-1,HECcompare(peakC),'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
                subplot(2,4,5) %E(y_rs)
                  plot((1:lags),EC_rs,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
                  title("E rs")
                  axis tight
                  hold on
                  plot([Lmin1C_rs Lmin2C_rs],[EC_rs(Lmin1C_rs) EC_rs(Lmin2C_rs)],'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
                subplot(2,4,6) %H(y_rs)
                  plot((1:lags),HC_rs,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
                  title("H rs")
                  axis tight
                  hold on
                  plot([Lmin1C_rs Lmin2C_rs],[HC_rs(Lmin1C_rs) HC_rs(Lmin2C_rs)],'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
                subplot(2,4,7) %HEdiffs(y_rs)
                  plot((1:lags),HECcompare_rs,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
                  title("E-2H rs")
                  axis tight
                  hold on
                  plot([Lmin1C_rs Lmin2C_rs],[HECcompare_rs(Lmin1C_rs) HECcompare_rs(Lmin2C_rs)],'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
                subplot(2,4,8) %HEdiffs_epsilon(y_rs)
                  plot((1:lags),HECcompare_rs_epsilon,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
                  title("E-2H<=eE rs")
                  axis tight
                  hold on
                  plot([Lmin1C_rs Lmin2C_rs],[HECcompare_rs_epsilon(Lmin1C_rs) HECcompare_rs_epsilon(Lmin2C_rs)],'*','MarkerSize',MarkerSizeBig,'Color',[0.9,0,0]);
            end
        end
    end
    if(print_samplei)
        fprintf("\n");
    end
    HEcont_elapsedtime = toc; %end timing
    fprintf("Continuous detection elapsed time: %.3f seconds\nExec time < sample time?: %s\n",HEcont_elapsedtime,mat2str(HEcont_elapsedtime < length(y)/samplerate));
    if(graphHEcont)
        FCplot = figure('Name','EH continuous f');
            set(FCplot,'Position', [screencenter(1)-gs(1)/2 screencenter(2)-gs(2)/2+100 gs(1)*2 gs(2)]);
            if(~readfile)
                plot((0:length(fC_tracker)-1).*resample_period.*(xTimeUnits_modifier),ones(1,length(fC_tracker)).*f.*fmod_rs,'LineWidth',2,'Color',[0.5,0,0]);
            end
            hold on
            plot((0:length(fC_tracker)-1).*resample_period.*(xTimeUnits_modifier),fC_tracker,'.-','MarkerSize',MarkerSizeSmall,'Color',graphcolors(1,:));
            title("f detected")
            axis tight
        Noteplot = figure('Name','EH continuous notes');
            set(Noteplot,'Position', [screencenter(1)-gs(1)/2 screencenter(2)-gs(2)/2+100-gs(2)/2 gs(1)*2 gs(2)]);
            notes = ["C" "C#/Db" "D" "D#/Eb" "E" "F" "F#/Gb" "G" "G#/Ab" "A" "A#/Bb" "B"];
            note_tracker(note_tracker < 1) = 1; %remove zero f
            note_tracker_cents = round(12*log(note_tracker/C0)/log(2)); %convert to cents/100
            note_tracker_range = 19:max(note_tracker_cents)+1; %G1 : highest, +1 to avoid errors
            note_tracker_octave = idivide(int32(note_tracker_range),int32(12)); %get octave number
            note_tracker_index = mod(note_tracker_range,12)+1; %get array index
            note_tracker_name = strings(1,length(note_tracker_index)); %allocate
            for a = 1:length(note_tracker_name)
                note_tracker_name(a) = strcat(notes(note_tracker_index(a)),num2str(note_tracker_octave(a)));
            end
            noteplot = plot((0:length(note_tracker)-1).*resample_period.*(xTimeUnits_modifier),note_tracker_cents,'.-','MarkerSize',MarkerSizeSmall,'Color',graphcolors(1,:));
            set(gca,'Ytick',note_tracker_range,'YTickLabel',note_tracker_name);
            title("closest note detected")
            ylim([min(note_tracker_range) max(note_tracker_range)])
            set(gca,'YGrid','on')
         figure(FCplot)
    end
    %pitch shifting
    if(true) %graph the two pitch shifting buffers
    PSplot = figure('Name','Pitch shifting EH');
        PSplot_x = (0:length(ynew_rsbuffer_ratio)-1).*resample_period.*(xTimeUnits_modifier);
        subplot(2,1,1)
        plot(PSplot_x,zeros(1,length(ynew_rsbuffer_diff)),'LineWidth',1,'Color',[0.5,0.5,0.5]);
        hold on
        plot(PSplot_x,log(ynew_rsbuffer_ratio)/log(2)*1200,'.-','MarkerSize',MarkerSizeSmall,'Color',graphcolors(1,:));
        title("difference in cents")
        subplot(2,1,2)
        plot(PSplot_x,ones(1,length(ynew_rsbuffer_diff)),'LineWidth',1,'Color',[0.5,0.5,0.5]);
        hold on
        plot(PSplot_x,ynew_rsbuffer_ratio,'.-','MarkerSize',MarkerSizeSmall,'Color',graphcolors(1,:));
        title("ratios of f")
    end
    %samplerates for each section = ratios * original samplerate
    ynew_samplerates_modified = ones(1,length(ynew_rsbuffer_ratio)) .* samplerate .* ynew_rsbuffer_ratio;
end

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
%fully autocorrelates an input for specified lags, and prints info
%outputs x, R, guessed frequency fR, and ^2 interpolated fR
    %defaults:
    fR = 0;
    fR_QI = 0;
    if(lags > length(y))
        lags = length(y);
    end
    if(~noAcorr)%for passing an R directly
        %R = autocorr(y,lags-1); %divides out variance
        R = autoc(y,lags,0); %doesn't divide out variance
    else
        R = y;
    end
    if(~isResampled)
        resample_rate = 1;
    end
    x = 0:length(R)-1; %create x values
    [~,locs] = findpeaks_fast(R); %find the peaks' x locations
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
    [~,plocs] = findpeaks_fast(pvalext); %find the peak's x value
    if(isempty(plocs)) %error check
        plocs(1) = round(max(pvalext));
        %fprintf("! No peak found in QI function. Set as max, bounds included\n");
    end
    Qx = plocs(1);
    %found x value of the peak between 0-500, turn back into 0-2:
    Qx = (Qx - 1) * (x_3(3) - x_3(1))/499 + x_3(1); %needed x value (non-integer)
end

function a = autoc(y,lags,minlag)
%autocorrelation with lag range minlag to (lags-1), doesn't divide out variance
    minlag = minlag + 1; %ind0 -> ind1
    if(length(y) < lags*2+1)
        y = [y zeros(1,lags*2+1-length(y))];
    end
    y = y(1:lags*2+1); %use only 2L points
    y_len = length(y);
    n = 2^nextpow2(2*y_len-1); %fft runs fastest on 2^n sized data
    a = ifft(fft(y,n) .* conj(fft(y,n))); %acorr = iF(F * F*) = iF(S)
    a = [zeros(1,minlag-1) a(minlag:lags)]; %only need up to L points
end

function E = calcE(y,lags,minlag)
%energy (no lag), with input range (minlag to (lags-1)) * 2
    minlag = minlag + 1; %ind0 -> ind1
    E_pmult = 2;% * L (= 2L used)
    E_lag = lags*E_pmult; %2L
    endy = E_lag*2;
    if(endy > length(y))
        endy = length(y);
    end
    y = y(1:endy);
    y = [y zeros(1,E_lag)+2]; %zero padding to stop array index errors
    E = zeros(1,lags); %prep, only doing limited lags
    y_squared = y .* y; %y^2 = y * y
    E(1) = sum(y_squared(1:2)); %Error check this recursive optimization
    for c = (minlag+1):length(E)
        %only sum c points, not all points
        %E(c) = sum(y_squared(1:c*E_pmult));
        E(c) = E(c-1) + y_squared(c*E_pmult-1) + y_squared(c*E_pmult);
    end
end

function H = calcH(y,lags,minlag)
%autocorrelation with lag range minlag to (lags-1), doesn't divide out variance
%only does up to L points for each lag
%check this function again sometime
    minlag_set = 1;
    if(nargin > 2)
        minlag_set = minlag+1;
    end
    endy = lags*2;
    if(endy > length(y))
        endy = length(y);
    end
    y = y(1:endy); %use only y(0 to lags*2-1)
    y = [y zeros(1,lags+2)]; %zero padding to stop array index errors
    H = zeros(1,lags); %set up with zeros
    %shifts left instead of right, but may be ok?
%     for c = minlag_set:length(H) %check this later
%         %m = y_i * y_i-L = y(t) * y(t+c)
%         m = y .* [y(c:end) zeros(1,c-1)]; %y * y shifted left by L
%         H(c) = sum(m(1:c)); %only sum c points, not all points
%     end
    
    for c = minlag_set:length(H)
        %m = y_i * y_i-L = y(t) * y(t+c)
        m = y .* [zeros(1,c-1) y(1:end-c+1)]; %y * y shifted right by L
        H(c) = sum(m(1:c)); %only sum c points, not all points
    end
end

function [pks,locs] = findpeaks_fast(y)
%custom findpeaks function much faster than stock
%pks = y values, locs = x values (indices)
    pks = zeros(1,0);
    locs = zeros(1,0);
    if(length(y) < 3) %require at least 3 points
        return;
    end
    num = 0;
    %matches stock output well, may have extra peak at end
    for i = 2:length(y)-1 %exclude bounds like stock
        if((y(i) - y(i-1) > 0) && (y(i+1) - y(i) <= 0))
            num = num + 1;
            pks(num) = y(i);
            locs(num) = i;
        end
    end
%     pks = pks(1:num);
%     locs = locs(1:num);
end

function f = constrain(y,bl,bu)
%return bounded value clipped between lower and upper bound
	f = min(max(y,bl),bu);
end

function f = normalize(y)
%normalize a function to range [-1, 1]
    f = 2*(y - min(y))./(max(y) - min(y)) - 1;
end

function [note_name, error_cents, target_f] = fetchnote(x)
%return closest 12EDO note n of frequency x (Hz) and error e (cents)
%target_f is the frequency of the correct note
    if(x <= 0) %handle bad input
        fprintf("! fetchnote() only to be used with positive f input\n");
        note_name = "?";
        error_cents = 0;
        target_f = 0;
        return;
    end
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

function target_f = fetchnote_fastf(x,A4freq,Czero)
%return closest 12EDO note frequency (Hz), target_f
    if(x <= 0) %handle bad input
        %fprintf("! fetchnote_fastf() only to be used with positive f input\n");
        target_f = 0;
        return;
    end
    h = round(12*log(x/Czero)/log(2)); %convert to cents/100
    x_octave = idivide(int32(h),int32(12)); %get octave number
    note_index = mod(h,12)+1; %get array index

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