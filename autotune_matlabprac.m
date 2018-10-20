%snr design autotune practice%Team: Ashlin G, Alex G, Brandon Y
clc, clear variables, clear sound, clear sounds, close all;

%% CONTROL VARIABLES
    %constants
    samplerate = 44100; %sample rate
    global A4freq; A4freq = 440; %Freq (Hz) for note A4
    C4 = A4freq * 2^(-9/12); %middle c (C4) from A4
    C0 = 2^(-(9+12*4)/12)*A4freq; %C0
    %scales = 12edo; tbd: just5. Listing before graph y section
    scale_name = '12edo';

    %control
    downsample_rate = 8; %down/up sample by these factors, new pitch = f * d/u
    upsample_rate = 1;
    lags = 110; %upper limit to range of lags used, (2 to 'lags')
    epsilon = 0.4; %arbitrary small value
    readfile = false; %read file or generate wave
        filename = 'test.wav'; %stereo files converted to mono
        frequency = C4*1; %f to generate (Hz)
        fixedlength = true; %generate fixed seconds
            secondslength = 1; %number of seconds length to generate
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

%% Generate wave
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
    fetchnote_print(f,A4freq,scale_name,scalemap);
    x = 0:c-1; %x = c # of samples
    %x = (x ./ c)  .* (p * pcount);  %(( x = 0:1) * period length) * # of periods
    %changed to x = samples only
    if(fixedlength)
        x = x(1:a); 
    end
    %frequency modifier (modulator?):
    %constant
        fmod_constant = ones(1,length(x));
    %constant f then ramp
        fmod_constthenramp_constratio = 1/4;
        fmod_constthenramp = [ones(1,floor(length(x).*fmod_constthenramp_constratio)) linspace(1,2,floor(length(x).*(1-fmod_constthenramp_constratio)))];
    %stair steps
        fmod_stairstep_size = floor(length(x)/12);
        a = linspace(1,2,floor(length(x)/(fmod_stairstep_size)));
        b = repmat(a,fmod_stairstep_size,1);
        fmod_stairsteps = b(:)';
    %oscillator
        fmod_f = 20; %Hz
        fmod_osc = 1 + 0.5*sin(2*pi*fmod_f.*x/samplerate); %oscillate sine (can vibrato)
    %select which fmod
    fmod = fmod_constthenramp;
    fmod = [fmod ones(1,length(x)-length(fmod)).*fmod(end)]; %try to ensure same length
    [~, fmod_rs] = resampleY(fmod,upsample_rate,downsample_rate,false); %get rs fmod for graph
    %input y
    y = sin(2.*pi.*(f.*fmod).*x/samplerate + 2*pi*phase_samples/samplesperperiod); %make a wave y, ex: y = sin(w*x + O), w = 2*pi*f
    %base period is samplerate
    fprintf("Generated sample count = %d (%.3f seconds)\n",length(x),length(x)/samplerate);
else
    [y,SR_file] = audioread(filename); %import y as file, SR_file samplerate
    if(size(y,2) > 1)%multichannel file file
        y = sum(y,2);
    end
    x = 1:length(y); %x = X samples to match y sample count
    f = 0;
    p = 0;
    samplesperperiod = 0;
    fprintf("File sample count = %d (%.3f seconds), samplerate = %d\n",length(x),length(x)/samplerate,SR_file);
    samplerate = SR_file; %experiment, override set samplerate with file's
end
%% preparatory stuff
%y_normalised = normalize(y); %comment this out if licensing complaint
if(~isrow(y)) %make y a row vector cuz thats what i like
    y = y';
end
resample_rate_s = num2str(upsample_rate) + "/" + num2str(downsample_rate);
sample_rate = samplerate; %alias
resamplerate = resample_rate; %alias
spp = samplesperperiod; %alias
A4 = A4freq;
Czero = C0;
%[x_rs, y_rs] = resampleY(y,upsample_rate,downsample_rate,false); %resample
[x_rs, y_rs] = resampleY(y,upsample_rate,downsample_rate,true); %resample+LPF
lag_s = num2str(lags);
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
graphcolors = [[0 0.2 0.6]; [0.4 0.4 0]; [0 0.6 0.2]; [0.7 0.2 0.0]]; %some colors
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
%scales:
scalemap = containers.Map; %matlab hashmap
%For now: must repeat at the octave, and be 12 notes C -> note before repeat.
  %12edo (default)
    scalemap('12edo') = [2^(1/12) 2^(2/12) 2^(3/12) 2^(4/12) 2^(5/12) 2^(6/12) 2^(7/12) 2^(8/12) 2^(9/12) 2^(10/12) 2^(11/12) 2];
  %5-limit just intonation
    scalemap('just5') = [16/15 9/8 6/5 5/4 4/3 45/32 3/2 8/5 5/3 16/9 15/8 2];
  %C major scale JI5
    scalemap('cmajor') = [1 9/8 5/4 5/4 4/3 4/3 3/2 5/3 5/3 15/8 15/8 2];
scale = scalemap(scale_name); %load the scale array
scalename = scale_name; %alias
fprintf("Using scale: '%s'\n",scale_name);
%
fprintf('\n');

%% graph y
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
        title("y" + soundedStr)
        ylim(ylimits)
    subplot(SubplotSize,1,SubplotSize) %w/ filter
        plot(x_rs.*(xTimeUnits_modifier),y_rs,'.-.','MarkerSize',MarkerSizeSmall,'Color',graphcolors(3,:))
        xlim([0 xlim_x_rs])
        title("y rs " + resample_rate_s)
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

%% full fixed autocorrelation
fprintf("R(y):\n");
    [R_x, R, R_f, R_fQI] = handle_R(y,length(y),readfile,f,samplerate,resample_rate,0,0,A4freq,scale_name,scalemap);
    R_p = R_f^-1; %periods
    R_pQI = R_fQI^-1;
    R_spp = R_p * samplerate; %samples per period
    R_sppQI = R_pQI * samplerate;
%fprintf("R(y rs):\n");
    %[R_rs_x, R_rs, R_rs_f, R_rs_fQI] = handle_R(y_rs,length(y_rs)-1,readfile,f,samplerate,resample_rate,1,0,A4freq,scale_name,scalemap);
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

%% Continuous Pitch Detection w/ E & H
if(doHEcont) %H,E continuous
    %prepare arrays. fill with zeros to avoid checks in the loop
    fC_tracker = zeros(1,length(y_rs)+0); %tracks detected frequency
    note_tracker = zeros(1,length(y_rs)+0); %tracks closest note frequency
    error_tracker = zeros(1,5); %track errors without spamming console
    %to eliminate computation time
    E_pmult = 2;% * L (= 2L used)
    E_lag = lags*E_pmult; %2L
    H_lag = lags; %L
    y_Lbuffer = zeros(1,E_lag*resample_period+RSP); %holds some y data
    y_Lbuffer_rs = zeros(1,E_lag+1+1);
    ynew_rsbuffer_diff = zeros(1,length(y_Lbuffer_rs)); %holds pitch shifting diffs
    ynew_rsbuffer_ratio = zeros(1,length(y_Lbuffer_rs)); %holds pitch shifting ratios
    y_smallbuffer = zeros(1,resample_period);
    EC_rs = zeros(1,H_lag); %H&E of rs and non-rs
    HC_rs = zeros(1,H_lag);
    EC = zeros(1,H_lag*resample_period);
    HC = zeros(1,H_lag*resample_period);
    %values from previous detection cycle
    LminC_last = lags*RSP; %meaningless init values
    peakC_last = LminC_last-2;
    peaknfriendsC_last = [peakC_last-1,peakC_last,peakC_last+1];
    Lmin1C_rs_last = lags;
    Lmin2C_rs_last = lags;
    HEC_spp = 0;
    HEC_p = 0;
    HEC_f = 0;
    f_closest = 0;
    HEC_f_diff = 0; %difference in f for pitch shift buffer
    %HEC_f_diff_last = 0;
    HEC_f_ratio = 1; %ratio of needed/detected
    %keep indexes as matlab 1:end. only use 0 for period calculation
    s = "";
    i_rs = 1;
    print_samplei = true; %show sample i=#
    usetimesavers = false; %messes up detection, leave off
    debugHEgraphs = false; %HE debug mode
        debug_starti = 0; %sample to start debug graphing at
    
    if(debugHEgraphs)
        HEdebugplots = figure('Name','EH debug plots');
    end
    tic; %time this process
    for i = 1:length(y) %big ol' for loop, realtime will use a while loop
        y_Lbuffer = [y_Lbuffer(2:end) y(i)];
%         y_smallbuffer = [y_smallbuffer(2:end) y(i)]; %faster
            %iterative E & H
            a = 1:H_lag*RSP;
                EHC_newterm1 = ones(1,H_lag*RSP) .* y_Lbuffer(end);
                EHC_newterm2 = y_Lbuffer(end-a-0);
                EHC_newterm3 = y_Lbuffer(end-2*a);
                EHC_newterm4 = y_Lbuffer(end-a-0);
                EC = EC + EHC_newterm1.^2 - EHC_newterm3.^2;
                HC = HC + EHC_newterm1.*EHC_newterm2 - EHC_newterm4.*EHC_newterm3;

        if(mod(i,1*resample_period) == 0) %only check pitch every resampled point for now
            
            %can put breakpoint on this for debug
            '';
            
            if(print_samplei) %only print in here too
                if(~debugHEgraphs)
                    fprintf(repmat(sprintf('\b'),1,strlength(s)));
                end
                s = sprintf("i = %d (%.2f ms)",i,i/samplerate*1000);
                fprintf("%s",s);
                if(debugHEgraphs)
                    fprintf("\n");
                end
            end
            %y_Lbuffer = [y_Lbuffer(1+resample_period:end) y_smallbuffer]; %do this in here
            y_Lbuffer_rs = [y_Lbuffer_rs(2:end) y_rs(i_rs)];
            %iterative E & H resampled
              a = 1:H_lag;
                EHC_rs_newterm1 = ones(1,H_lag) .* y_Lbuffer_rs(end);
                EHC_rs_newterm2 = y_Lbuffer_rs(end-a-0);
                EHC_rs_newterm3 = y_Lbuffer_rs(end-2*a);
                EHC_rs_newterm4 = y_Lbuffer_rs(end-a-0);
                EC_rs = EC_rs + EHC_rs_newterm1.^2 - EHC_rs_newterm3.^2;
                HC_rs = HC_rs + EHC_rs_newterm1.*EHC_rs_newterm2 - EHC_rs_newterm4.*EHC_rs_newterm3;
            HECcompare_rs = EC_rs - E_pmult.*HC_rs; %E - 2H
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
                %Lmin1C_rs = length(HECcompare_rs_epsilon); %use right bound
                Lmin1C_rs = Lmin1C_rs_last; %use last value
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
                if(false) %non-iterative E & H calculation
                    %len: = E_lag*resample_period+RSP
                    EHC_calc_Lbuffer = y_Lbuffer(end-(E_lag*RSP-1):end);
                    EC = ones(1,length(EC));
                    HC = ones(1,length(EC));
                    %get ranges to calculate
                    EHC_rL1 = Lmin1C-RSP:Lmin1C+RSP;
                    EHC_rL1 = EHC_rL1(EHC_rL1 <= length(EC));
                    EHC_rL2 = Lmin2C-RSP:Lmin2C+RSP;
                    EHC_rL2 = EHC_rL2(EHC_rL2 <= length(EC));
                    %calc
                    EC(EHC_rL1) = calcE(EHC_calc_Lbuffer,EHC_rL1(end),EHC_rL1(1)-1);
                    EC(EHC_rL2) = calcE(EHC_calc_Lbuffer,EHC_rL2(end),EHC_rL2(1)-1);
                    HC(EHC_rL1) = calcH(EHC_calc_Lbuffer,EHC_rL1(end),EHC_rL1(1)-1);
                    HC(EHC_rL2) = calcH(EHC_calc_Lbuffer,EHC_rL2(end),EHC_rL2(1)-1);
                    %set uncalculated values (ones) to edges' values
                    EC(1:length(EC)-1 < EHC_rL1(1)) = EC(EHC_rL1(1));
                    EC(1:length(EC)-1 > EHC_rL2(end)) = EC(EHC_rL2(end));
                    EC(1:length(EC)-1 > EHC_rL1(end) & 1:length(EC)-1 < EHC_rL2(1)) = EC(EHC_rL2(1));
                    HC(1:length(EC)-1 < EHC_rL1(1)) = HC(EHC_rL1(1));
                    HC(1:length(EC)-1 > EHC_rL2(end)) = HC(EHC_rL2(end));
                    HC(1:length(EC)-1 > EHC_rL1(end) & 1:length(HC)-1 < EHC_rL2(1)) = HC(EHC_rL2(1));
                end
                %E-2H
                HECcompare = EC - E_pmult.*HC;
                if(~isrow(HECcompare)) %to row vector
                    HECcompare = HECcompare';
                end
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
                if(detectNewPitch) %detect pitch from valley near Lmin
                    %use large range of points around Lmin to find valley
                    %later write function to find the valley using n points
                    for k = 0:1 %do for Lmin 1 & 2 and then compare
                        if(k == 0)%1 & 2 part
                            LminC = Lmin2C;
                        else
                            LminC = Lmin1C;
                        end
                        LminC_points = (LminC)-RSP/2+1:(LminC)+RSP/2; %get range of 8 points
                        LminC_points_ext = (LminC)-RSP*2+1:(LminC)+RSP*2; %range of 32 points
                        %try and stop any errors
                        LminC_points_ext(LminC_points_ext > length(HECcompare)-1) = length(HECcompare)-1;
                        %determine HE_f:
                        [~, peakC] = findpeaks_fast(-HECcompare(LminC_points_ext));
                        if(isempty(peakC)) %no valley found for LminC
                            error_tracker(5) = error_tracker(5) + 1;
                            %peakC = LminC_points_ext(end)-1; %-1 to stop indexing errors
                            peakC = peakC_last;
                        else
                            peakC = LminC_points_ext(peakC(1)); %no shift -1 needed
                        end
                        peaknfriendsC = [peakC-1,peakC,peakC+1]; %array of 3 points
                        if(k == 0)%1 & 2 part
                            peakC2 = peakC;
                            peaknfriendsC2 = peaknfriendsC;
                            LminC_points_ext2 = LminC_points_ext;
                        else
                            peakC1 = peakC;
                            peaknfriendsC1 = peaknfriendsC;
                            LminC_points_ext1 = LminC_points_ext;
                        end
                    end
                    %now set the real Lmin and use its valley
                    if(Lmin1C_rs < lags/2 && HECcompare(peakC2) <= HECcompare(peakC1))
                        %use Lmin2 stuff if its valley is lower
                         LminC = Lmin2C;
                         peakC = peakC2;
                         peaknfriendsC = peaknfriendsC2;
                    end
                    %use first valley (higher f) if both are close
%                     if(abs(HECcompare(peakC2) - HECcompare(peakC1)) < 1)
%                         LminC = Lmin1C;
%                          peakC = peakC1;
%                          peaknfriendsC = peaknfriendsC1;
%                     end
                    LminC_last = LminC;
                    LminC_rs = LminC * resample_rate;
                    %no shift left -1 needed due to above peakC indexing
                    HEC_spp = QInterp_peak(peaknfriendsC,-HECcompare(peaknfriendsC))-0; %QI
                    HEC_p = HEC_spp/samplerate;
                    HEC_f = samplerate/(HEC_spp); %detected f
                    f_closest = fetchnote_fastf(HEC_f,C0,scale);
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
            if(debugHEgraphs && i >= debug_starti)
                fprintf("detected f = %.2f Hz",HEC_f);
                if(~readfile)
                    fprintf(", real f = %.2f Hz",f.*fmod_rs(i_rs));
                end
                fprintf("\n")
                %EHC_plot_x =(LminC_points_ext-1).*(xTimeUnits_modifier);
                EHC_plot_x = 0:length(EC)-1;
                %figure(HEdebugplots)
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
    HEcont_pd_elapsedtime = toc; %end timing
    fprintf("Continuous detection elapsed time: %.3f seconds\n",HEcont_pd_elapsedtime);
    if(graphHEcont)
        HEcont_plot_x = (0:length(fC_tracker)-1).*resample_period.*(xTimeUnits_modifier);
        FCplot = figure('Name','EH continuous f');
            set(FCplot,'Position', [screencenter(1)-gs(1)/2 screencenter(2)-gs(2)/2+100 gs(1)*2 gs(2)]);
            if(~readfile)
                plot(HEcont_plot_x,ones(1,length(fC_tracker)).*f.*fmod_rs,'LineWidth',2,'Color',[0.5,0,0]);
            end
            hold on
            plot(HEcont_plot_x,fC_tracker,'.-','MarkerSize',MarkerSizeSmall,'Color',graphcolors(1,:));
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
            noteplot = plot(HEcont_plot_x,note_tracker_cents,'.-','MarkerSize',MarkerSizeSmall,'Color',graphcolors(1,:));
            set(gca,'Ytick',note_tracker_range,'YTickLabel',note_tracker_name);
            title("closest note detected")
            ylim([min(note_tracker_range) max(note_tracker_range)])
            set(gca,'YGrid','on')
         figure(FCplot)
         %error plot if not reading file
         if(~readfile)
             fC_tracker_error = fC_tracker - ones(1,length(fC_tracker)).*f.*fmod_rs;
             fC_tracker_error_cents = Hz2Cents(fC_tracker,C0) - Hz2Cents(ones(1,length(fC_tracker)).*f.*fmod_rs,C0);
             fC_tracker_error_cents_mod = mod(fC_tracker_error_cents,1200);
             fC_tracker_error_cents_mod =  fC_tracker_error_cents_mod - (fC_tracker_error_cents_mod > 600).*1200;
             FCplot_error = figure('Name','EH continuous f error');
             subplot(3,1,1)
                plot(HEcont_plot_x,zeros(1,length(HEcont_plot_x)),'LineWidth',1,'Color',[0.5,0.5,0.5]);
                hold on
                plot(HEcont_plot_x,fC_tracker_error,'.-','MarkerSize',MarkerSizeSmall,'Color',graphcolors(1,:));
                title("f error in Hz")
                axis tight
             subplot(3,1,2)
                 plot(HEcont_plot_x,zeros(1,length(HEcont_plot_x)),'LineWidth',1,'Color',[0.5,0.5,0.5]);
                hold on
                plot(HEcont_plot_x,fC_tracker_error_cents,'.-','MarkerSize',MarkerSizeSmall,'Color',graphcolors(1,:));
                title("f error in cents")
                axis tight
             subplot(3,1,3)
                plot(HEcont_plot_x,zeros(1,length(HEcont_plot_x)),'LineWidth',1,'Color',[0.5,0.5,0.5]);
                hold on
                plot(HEcont_plot_x,fC_tracker_error_cents_mod,'.-','MarkerSize',MarkerSizeSmall,'Color',graphcolors(1,:));
                title("f error in cents octave reduced")
                axis tight
         end
    end
%% pitch shifting
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
    tic; %time
    %samplerates for each section = ratios * original samplerate
    ynew_samplerates_modified = ones(1,length(ynew_rsbuffer_ratio)) .* samplerate .* ynew_rsbuffer_ratio;
    %y to be retuned
    y_retunedHE = y;
    y_retunedHE = zeros(1,length(y));
    %make sure size of y is evenly divisible by resampling period
    if(length(y_retunedHE)/RSP < length(ynew_rsbuffer_ratio))
        %y_retunedHE = [y_retunedHE zeros(1,RSP-mod(length(y_retunedHE),RSP))]
        y_retunedHE = [y_retunedHE zeros(1,RSP)]; %just pad some zeros
    end
    c = 1; %counter
    length_changed_check = true;
    %set below to be: > than SPP of lowest f you want autotuned
    %128 * 8 (default spp = ds rate) = 1024 samples will be > spp of 50Hz (882).
    retune_sliceRatio = 128; %resize set of points (RSP*wR) (whole number)
    retune_ratio_transpose = 1; %transpose by additional ratio
    %for loop at the end for fixed-time
    for a = 1:floor(length(ynew_rsbuffer_ratio)/retune_sliceRatio)
        retune_RSP = RSP * retune_sliceRatio;
        retune_RSP_range = 1+retune_RSP*(a-1):retune_RSP+retune_RSP*(a-1); %RSP# point ranges
        %use pitch shift function for each set of points
        retune_ratio = ynew_rsbuffer_ratio(a*retune_sliceRatio) .* retune_ratio_transpose;
        if(retune_sliceRatio > 1) %try out averaging the PS ratio over the slice
            retune_ratio_sRrange = 1+retune_sliceRatio*(a-1):retune_sliceRatio+retune_sliceRatio*(a-1);
            if(retune_ratio_sRrange(end) > length(ynew_rsbuffer_ratio)) %stop errors
                retune_ratio_sRrange = retune_ratio_sRrange(1):length(ynew_rsbuffer_ratio);
            end
            retune_ratio = mean(ynew_rsbuffer_ratio(retune_ratio_sRrange) .* retune_ratio_transpose);
        end
        %perform the actual pitch shifting of the slice
        retune_part = pitch_shift_Npt(y(retune_RSP_range),retune_ratio);
        %notify if different # of points as input, once
        if(length(retune_part) ~= length(y(retune_RSP_range)) && length_changed_check)
            fprintf("! Using PS method that resizes slices\n");
            length_changed_check = false;
        end
        %y_retunedHE(retune_RSP_range) = retune_part;
        %supports differing lengths:
        y_retunedHE(c:c+length(retune_part)-1) = retune_part;
        c = c + length(retune_part);
    end
    HEcont_ps_elapsedtime = toc; %end timing
    fprintf("Continuous retuning elapsed time: %.3f seconds\nTotal time: %.3f seconds\nTotal Exec time < sample time?: %s\n",HEcont_ps_elapsedtime,HEcont_pd_elapsedtime+HEcont_ps_elapsedtime,mat2str(HEcont_pd_elapsedtime+HEcont_ps_elapsedtime < length(y)/samplerate));
    if(length(y_retunedHE) ~= length(y))
        fprintf("! Final PS wave different size than input. Sample difference = %d\n",length(y)-length(y_retunedHE));
    end
    %graph
    if(true)
      Retunedplot = figure('Name','Retuned waveform (HE)');
        if(playsounds)
            soundedStr = "(sounded 2nd) ";
        end
      subplot(1,1,1)
        set(Retunedplot,'Position', [screencenter(1)-gs(1)/2-gs(1) screencenter(2)-gs(2)/2+100-gs(2)/2 gs(1) gs(2)]);
        plot((0:length(y_retunedHE)-1).*(xTimeUnits_modifier),y_retunedHE,'.-','MarkerSize',MarkerSizeTiny,'Color',graphcolors(4,:));
        %xlim([0 xlim_x])
        title("y retuned" + soundedStr)
        ylim(ylimits)
    end
end

%% playback pitch shifted wave
%playing sounds through resampling
% [~, ~, target_f] = fetchnote(R_fQI,A4freq,scale_name,scalemap); %find needed f
% samplerate_tuned_ratio = target_f/R_fQI; %ratio between needed/found f
% samplerate_tuned = round(samplerate * samplerate_tuned_ratio); %new S.R.
% %higher SR = higher pitch, lower SR = lower pitch
% y_retuned = y; %with diff samplerates, wave will be the same

%uncomment these if you are using your own pitch-shifted wave
y_retuned = y_retunedHE; %put wave here
samplerate_tuned = samplerate; %restore samplerate

%save Y and Y_retuned to .wav files
[~,~] = mkdir('output/');
audiowrite('output/Y.wav',y,samplerate);
audiowrite('output/Y_retuned.wav',y_retuned,samplerate);
fprintf("\n")

if(playsounds) %plays the sounds
    volpeakpoint = 100/volume;
    y_toplay = [volpeakpoint y];
    fprintf("Playing y...\n");
    soundsc(y_toplay,(samplerate))
    pause
    clear sound, clear sounds;
    y_toplay = [volpeakpoint y_retuned];
    fprintf("Playing y (pitch corrected thru continuous HE method)...\n\n");
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

%% FUNCTIONS
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

function [x, R, fR, fR_QI] = handle_R(y,lags,readfile,f,samplerate,resample_rate,isResampled,noAcorr,A4freq,scale_name,scalemap)
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
        %locs = locs;%shift left cuz matlab
        meanp = mean(diff([1 locs])); %find the mean period, fixed calculation
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
        xweneed = QInterp_peak(peaknfriends,R(peaknfriends)); %QI
        fR_QI = samplerate/(xweneed-1)*resample_rate; %QI detected f
        fprintf('AVG Peak p detected pitch: R_f = ') %print out f results
        report_ptof(fR,readfile,f,A4freq,scale_name,scalemap)
        fprintf('Q. interpolated p detected pitch: R_fQI = ')
        report_ptof(fR_QI,readfile,f,A4freq,scale_name,scalemap)
    end
    fprintf("\n");
end

function report_ptof(f,readfile,f_generated,A4freq,scale_name,scalemap)
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
    fetchnote_print(f,A4freq,scale_name,scalemap); %find closest note
end

function Qx = QInterp_peak(x_3, y_3)
%finds peak of quadratic interpolation on 3 points
    p = polyfit_fast(-1:1,y_3,2); %fit ^2, a = y, x = -1 to 1 (3 points effectively)
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
%energy (no lag), with input range 2*[minlag, (lags-1)]
%only does up to L points for each lag
    minlag = minlag + 1; %calculate only from this lag onward (0 -> 1)
    E_pmult = 2;% * L (= 2L used)
    y_squared = y .* y; %y^2 = y * y
    E = zeros(1,lags); %set up with zeros
    %get E(minlag-1):
    prevEc_minlagpart = y_squared(end-(minlag-1)*E_pmult+1:end);
    prevEc = (minlag>1)*sum(prevEc_minlagpart);
    for c = minlag:length(E)
        %recursive sum + next two y^2 points
        E(c) = prevEc + y_squared(end-c*E_pmult+1) + y_squared(end-(c-1)*E_pmult);
        prevEc = E(c);
    end
    E = E(minlag:end); %only return range
end

function H = calcH(y,lags,minlag)
%autocorrelation with lag range [minlag, lags-1], doesn't divide out variance
%only does up to L points for each lag
    minlag = minlag + 1; %calculate only from this lag onward (0 -> 1)
    y = [zeros(1,lags+2) y]; %zero padding on left to stop array index errors
    H = zeros(1,lags); %set up with zeros
    for c = minlag:length(H)
        %sum of y_i * y_i-L
        H(c) = sum(y(end-c+1:end) .* y(end-c+1-c+1:end-c+1));
    end
    H = H(minlag:end); %only return range
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

function [p,S,mu] = polyfit_fast(x,y,n)
%same as polyfit.m but removed some warnings that take forever
if ~isequal(size(x),size(y))
    error(message('MATLAB:polyfit:XYSizeMismatch'))
end

x = x(:);
y = y(:);

if nargout > 2
   mu = [mean(x); std(x)];
   x = (x - mu(1))/mu(2);
end

% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end

% Solve least squares problem.
[Q,R] = qr(V,0);
p = R\(Q'*y);    % Same as p = V\y;
if size(R,2) > size(R,1)
   warning(message('MATLAB:polyfit:PolyNotUnique'))
end

if nargout > 1
    r = y - V*p;
    % S is a structure containing three elements: the triangular factor from a
    % QR decomposition of the Vandermonde matrix, the degrees of freedom and
    % the norm of the residuals.
    S.R = R;
    S.df = max(0,length(y) - (n+1));
    S.normr = norm(r);
end

p = p.';          % Polynomial coefficients are row vectors by convention.
end

function f = constrain(y,bl,bu)
%return bounded value clipped between lower and upper bound
	f = min(max(y,bl),bu);
end

function f = normalize(y)
%normalize a function to range [-1, 1]
    f = 2*(y - min(y))./(max(y) - min(y)) - 1;
end

%these need to be addressed for different scales now
% function cents = Hz2Cents(f,A4freq)
% %change Hz to cents based off C0 (from A4)
%     if(f <= 0) %handle bad input
%         cents = 0; %f > 0 Hz required
%         return;
%     end
%     Czero = 2^(-(9+12*4)/12)*A4freq; %C0
%     cents = 1200.*log(f./Czero)./log(2);
% end
% 
% function f = Cents2Hz(cents,A4freq)
% %change cents to Hz based off C0 (from A4)
%     Czero = 2^(-(9+12*4)/12)*A4freq; %C0
%     f = 2.^(cents./1200).*Czero;
% end

function f_OR = octave_reduce(f)
%octave reduce f (Hz) to range [1,2) Hz
    f_OR = exp(mod(log(f),log(2)));
end

function [note_name, error_cents, target_f] = fetchnote(x,A4freq,scale_name,scalemap)
%return closest 12EDO note n of frequency x (Hz) and error e (cents)
%target_f is the frequency of the correct note
    if(x <= 0) %handle bad input
        fprintf("! fetchnote() only to be used with positive f input\n");
        note_name = "?";
        error_cents = 0;
        target_f = 0;
        return;
    end
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

function target_f = fetchnote_fastf(x,A4freq,scale_ARRAY)
%return closest frequency (Hz) from scale ARRAY (12 notes), target_f.
%octave equivalency assumed
    if(x <= 0) %handle bad input
        %fprintf("! fetchnote_fastf() only to be used with positive f input\n");
        target_f = 0;
        return;
    end
    scale_ARRAY = [1 scale_ARRAY]; %append 1 at beginning
    C4 = A4freq/scale_ARRAY(9+1); %C4 from scale degree 9
    C0 = C4/16;
    SA_cents = Hz2Cents(scale_ARRAY,1);
    x_cents = Hz2Cents(x,C0);
    [~,note_index] = min(abs(SA_cents - mod(x_cents,1200)));

    target_f = Cents2Hz(x_cents + SA_cents(note_index) - mod(x_cents,1200),C0);
    
    function cents = Hz2Cents(f,Czero)
        cents = 1200.*log(f./Czero)./log(2);
    end
    function f = Cents2Hz(cents,Czero)
        f = 2.^(cents./1200).*Czero;
    end
end

function fetchnote_print(x,A4freq,scale_name,scalemap)
%print out fetchnote information
    [n, e, f] = fetchnote(x,A4freq,scale_name,scalemap);
    fprintf("Closest note = %s (%.2f Hz) Error (cents) = %+.1f\n",n,f,e);
end

function y_new = pitch_shift_Npt(y,ratio)
%the pitch shifting used on the groups of N points in HE
    [num, denom] = rat(ratio); %matlab's rational approximator
    
    y_new = pitch_shift_Npt_linterp(y,ratio,false); %change length
    %uncomment below to use
    %y_new = pitch_shift_ratio(y,num,denom);
    
end

function y_new = pitch_shift_ratio(y,num,denom)
%pitch shifts a data set by a integer ratio n/d
%n & d will usually be around 3-4 digits (base 10)
    y_new = y; %init
    %do the things in here
    
    
    %make sure y_new is set as output
    
end

function y_new = pitch_shift_Npt_linterp(y,ratio,retain_length)
%pitch shifts N points by a ratio through stretch or compress
%will create discontinuities/jumps esp with retain_length, iffy behavior
%try and clean this up later?
    if(ratio == 1 || ratio <= 0) %handle bad cases
        y_new = y;
        return;
    end
    y_to_interp = repmat(y,1,ceil(1/ratio)); %have enough periods avail
    if(ratio < 1)
        y_to_interp = y_to_interp(1:round(length(y)/ratio));
    end
    x_new = (0:length(y_to_interp)-1).*ratio; %new dx to interp
    y_new = interp1(0:length(x_new)-1,y_to_interp,x_new);
    y_new(isnan(y_new)) = [];
    if(length(y_new) < length(y))
        y_new = repmat(y_new,1,ceil(length(y)/length(y_new)));
        y_new = y_new(1:round(length(y)/ratio));
    end
    if(retain_length)
        y_new(isnan(y_new)) = [];
        y_new = repmat(y_new,1,ceil(length(y)/length(y_new)));
        y_new = y_new(1:length(y)); %only use first N points
    end
    y_new(isnan(y_new)) = []; %remove null points
end

%for MATLAB; ignore some warnings
    %#ok<*DEFNU>
    %#ok<*UNRCH>
    %#ok<*NASGU>
    %#ok<*NUSED>