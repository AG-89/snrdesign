%snr design autotune practice
%Team: Ashlin G, Alex G, Brandon Y
clc, clear variables, close all

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
            pcount = 20; %number of periods to generate otherwise
    playsounds = false;
    volume = 10; %0-100 please

    %graphing stuff
    xTimeUnits = false; %x units as samples or time
    xLimTo_pcount = false; %only show pcount periods in graph (for many samples)
    graphAnything = true; %graph... anything
        graphY = true; %graph y
        graphRy = true; %graph R(y)
        graphRyM = true; %graph R(y)M with specified lag
        graphCF = false; %graph curve fitting
        %graphResample = false; %understand that LP filter is needed
        doHE = true;
        graphHE = true; %graph E(L) & H(L)
%END CVARS

resample_rate = upsample_rate / downsample_rate;
%wave input section
if(readfile == false)
    f = frequency; %known frequency
    p = f^-1; %period = 1/f
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
    %x = (x ./ c)  .* (p * pcount);  %((x = 0:1) * period length) * # of periods
    %changed to x = samples only
    if(fixedlength)
        x = x(1:a); 
    end
    y = sin(2*pi*f.*x/samplerate);%,1/2); %makse y a wave, ex y = sin(w * x), w = 2 * pi * f (omega, radian frequency)
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
graphRyM = graphAnything && graphRyM;
graphCF = graphAnything && graphCF;
graphHE = graphAnything && graphHE;
SubplotSize = 2;
screensize = get(groot, 'ScreenSize'); screensize = [screensize(3) screensize(4)]; % get screensize
screencenter = screensize./2;
ylimits = [-1.25 1.25];
graphsize = [560 420]; %size of figure windows
graphcolors = [[0 0.2 0.6]; [0.4 0.4 0]; [0 0.6 0.2]];
MarkerSizeNormal = 6;
MarkerSizeResampled = 10;
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
    if(playsounds)%play sounds
        volpeakpoint = 100/volume;
        y_toplay = y;
        y_toplay(1) = volpeakpoint;
        fprintf("Playing y...\n");
        soundsc(y_toplay,(samplerate))
        pause
        y_toplay = y_rs;
        y_toplay(1) = volpeakpoint;
        %soundsc(y_vol,(samplerate / downsample_rate))
        fprintf("Playing y rs...\n\n");
        soundsc(y_toplay,(samplerate))
        pause
    end
end
clear sound, clear sounds;

%no lag autocorr part
%xlim_R = xlim_x;%int32(1000); % how far to scale the x axis
fprintf("R(y):\n");
    [R_x, R, R_f, R_fQI] = handle_R(y,length(y)-1,readfile,f,samplerate,resample_rate,0,0);
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
  Rplots = figure('Name','Autocorrelations');
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

%_M_anually specified lag autocorr part
fprintf("Lags used: [0-%d]\n",lags-1);
% fprintf("R(y)M:\n");
%     [RM_x, RM, RM_f, RM_fQI] = handle_R(y,lags,readfile,f,samplerate,resample_rate,0,0);
%     xlim_x_RM = (length(RM)-1) * (xTimeUnits_modifier);
fprintf("R(y rs)M:\n");
    [RM_rs_x, RM_rs, RM_rs_f, RM_rs_fQI] = handle_R(y_rs,lags,readfile,f,samplerate,resample_rate,1,0);
    xlim_x_RM_rs = (length(RM_rs)-1) * (xTimeUnits_modifier);
if(graphRyM)
  RMplots = figure('Name','Manual Autocorrelations');
  gs = graphsize; %alias, below sets position
  set(RMplots,'Position', [screencenter(1)-gs(1)/2+gs(1) screencenter(2)-gs(2)/2+100 gs(1) gs(2)]);
%   subplot(SubplotSize,1,1) %RM(y)
%       plot(RM_x.*(xTimeUnits_modifier),RM,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
%       title("R(y)M" + " L=" + num2str(lags))
%       axis tight
%       xlim([0 xlim_x_RM])
%       %ylim(ylimits)
  subplot(1,1,1) %RM(y_rs)
      plot(RM_rs_x.*(xTimeUnits_modifier),RM_rs,'.-.','MarkerSize',MarkerSizeResampled,'Color',graphcolors(3,:));
      title("R(y rs)M " + resample_rate_s + " L=" + num2str(lags));
      axis tight
      xlim([0 xlim_x_RM_rs])
      %ylim(ylimits)
end
    
if(doHE) %H,E equations part
    %E: energy equation
    E_pmult = 2;% * L (= 2L used)
    EM_lag = lags*E_pmult; %2L
    EM = zeros(1,lags-E_pmult); %prep, only doing limited lags
    EM_rs = EM;
    yEM = [y zeros(1,EM_lag)]; %prep to stop array index errors
    yEM_rs = [y_rs zeros(1,EM_lag)];
    for k = 2:lags %calculate E for all sets of 2L's, L=0 (k=1) is always 0
        EM(k) = sum(yEM(1:E_pmult*k-E_pmult).^2); %sum up squares in range [0,2L]
        EM_rs(k) = sum(yEM_rs(1:E_pmult*k-E_pmult).^2);
    end
    %H is similar to RM except it only Acorr's a range of points 0-L, 
    %not every point
    HM = calcH(y,lags);
    HM_rs = calcH(y_rs,lags);
    %handle_R(HM_rs,lags,readfile,f,samplerate,resample_rate,1,1);
    HEMcompare = EM(1:length(HM)) - 2.*HM;
    HEMcompare_rs = EM_rs(1:length(HM_rs)) - 2.*HM_rs;
    if(all(HEMcompare >= 0) && all(HEMcompare_rs >= 0)) %E >= 2H
        fprintf("EM,HM comparison true.\n");
    else
        fprintf("EM,HM comparison false, failed.\n");
    end
    epsilon = 0.4; %arbitrary small value
    fprintf("epsilon = %.2f\n",epsilon)
    %E - 2H > epsilon * E
    HEMcompare_epsilon = HEMcompare <= (epsilon * EM(1:length(HEMcompare)));
    HEMcompare_rs_epsilon = HEMcompare <= (epsilon * EM_rs(1:length(HEMcompare_rs)));
    %Lmin
    [Lmin_y, Lmin_x] = findpeaks(-HEMcompare,'SortStr','descend');
    [Lmin_y_rs, Lmin_x_rs] = findpeaks(-HEMcompare_rs,'SortStr','descend');
    %add code for no peaks found later (error check)
    Lmin_x = Lmin_x(1); %Lmin
    Lmin_x_rs = Lmin_x_rs(1);
    if(Lmin_x(1) < lags/2)
        [Lmin2_y, Lmin2_x] = findpeaks(-HEMcompare(lags/2+1:end),'SortStr','descend');
        [Lmin2_y(1), Lmin2_x(1)] %Lmin2_x(1) is Lmin2
        Lmin2_x = Lmin2_x(1) + lags/2;
    end
    if(Lmin_x_rs(1) < lags/2)
        [Lmin2_y_rs, Lmin2_x_rs] = findpeaks(-HEMcompare_rs(lags/2+1:end),'SortStr','descend');
        [Lmin2_y_rs(1), Lmin2_x_rs(1)];
        Lmin2_x_rs = Lmin2_x_rs(1) + lags/2;
    end
    %
    if(graphHE)
    figure('Name','E and H functions')
    subplot(2,4,1) %E(y)
      plot((0:E_pmult:EM_lag-E_pmult).*(xTimeUnits_modifier),EM,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
      title("E(y)M" + " 2L=" + num2str(EM_lag))
      axis tight
    subplot(2,4,5) %E(y_rs)
      plot((0:E_pmult:EM_lag-E_pmult).*(xTimeUnits_modifier),EM_rs,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
      title("E(y rs)M" + " 2L=" + num2str(EM_lag))
      axis tight
    subplot(2,4,2) %H(y)
      plot((0:lags-1).*(xTimeUnits_modifier),HM,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
      title("H(y)M")
      axis tight
    subplot(2,4,6) %H(y_rs)
      plot((0:lags-1).*(xTimeUnits_modifier),HM_rs,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
      title("H(y rs)M")
      axis tight
    subplot(2,4,3) %HEdiffs(y)
      plot((0:lags-1).*(xTimeUnits_modifier),HEMcompare,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
      title("E-2H(y)M")
      axis tight
    subplot(2,4,7) %HEdiffs(y_rs)
      plot((0:lags-1).*(xTimeUnits_modifier),HEMcompare_rs,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
      title("E-2H(y rs)M")
      axis tight
    subplot(2,4,4) %HEdiffs_epsilon(y)
      plot((0:lags-1).*(xTimeUnits_modifier),HEMcompare_epsilon,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(1,:));
      title("E-2H<=eE(y)M")
      axis tight
    subplot(2,4,8) %HEdiffs_epsilon(y_rs)
      plot((0:lags-1).*(xTimeUnits_modifier),HEMcompare_rs_epsilon,'.-','MarkerSize',MarkerSizeNormal,'Color',graphcolors(3,:));
      title("E-2H<=eE(y rs)M")
      axis tight
    end
end

%Curve Fitting part, not viable cuz the peak doesn't change
if(graphCF)
  fprintf("-\nCurve Fitting\n-\n");
  fittype = 'cubicinterp'; %what type of curve to fit
  Gfitplots = figure('Name','RM fits');
  movegui(Gfitplots,'southwest')
  subplot(SubplotSize,1,1)
      Gfit = fit(RM_x',RM',fittype);
      plot(Gfit,RM_x,RM);
      title("R(y)M fit")
      ylim(ylimits)
  subplot(SubplotSize,1,SubplotSize)
      Gfit_rs = fit(RM_rs_x',RM_rs',fittype);
      plot(Gfit_rs,RM_rs_x,RM_rs);
      title("R(y rs)M fit")
      ylim(ylimits)
  GfitplotsD = figure('Name','RM fits discrete');
  movegui(GfitplotsD,'south')
  subplot(SubplotSize,1,1)
    GfitY = feval(Gfit,0:110)';
    fprintf("R(y)CF:\n");
    handle_R(GfitY,lags,readfile,f,samplerate,resample_rate,0,1);
    plot(GfitY)
    ylim(ylimits)
  subplot(SubplotSize,1,SubplotSize)
    Gfit_rsY = feval(Gfit_rs,0:110)';
    fprintf("R(y rs)CF:\n");
    handle_R(Gfit_rsY,lags,readfile,f,samplerate,resample_rate,1,1);
    plot(Gfit_rsY,'.-.','MarkerSize',MarkerSizeResampled,'Color',graphcolors(3,:))
    ylim(ylimits)
end

%clean up
if(graphAnything && (graphY || graphRy || graphRyM || doHE))
    %pause
    %close all;
end
clear sound, clear soundsc

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
        fprintf("No valid 2nd peak detected\n");%,locs(1));
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
        %firstpeak_mean_diff = meanp - firstpeak; %unused
        peaknfriends = [firstpeak-1,firstpeak,firstpeak+1]; %array of 3 points
        peaknfriends_shift = peaknfriends + 1; %shift right to reverse SL
        if(~isrow(R)) %to row vector
            R = R';
        end
        a = R(peaknfriends_shift); %get R values
        xweneed = QInterp_peak(peaknfriends,a); %moved to function
        fR_QI = samplerate/(xweneed)*resample_rate; %QI detected f
        
        fprintf('AVG Peak p ') %print out f results
        report_ptof(fR,readfile,f)
        fprintf('Q. interpolated p ')
        report_ptof(fR_QI,readfile,f)
    end
    fprintf("\n");
end

function report_ptof(f,readfile,f_generated)
%reports frequency results. generatedf used if not reading a file
    %fprintf('Detected Pitch: f = %.2f Hz, p = %f s, samples/period = %.5g samples\n',fR,1/fR,samplerate/fR*resample_rate);
    fprintf('Detected Pitch: f = %.2f Hz, ',f)
    if(~readfile) %error from generated pitch
      %fprintf('Error from generated f = %+.2f Hz, cents = %+.1f\n', fR_QInter - f, log(double(fR_QInter/f))/log(2)*1200);
      fprintf('Error (cents) = %+.1f\n', log(double(f/f_generated))/log(2)*1200);
    end
    fetchnote_print(f); %find closest note
end

function Qx = QInterp_peak(x_3, y_3)
%finds peak of quadratic interpolation on 3 points
    p = polyfit(-1:1,y_3,2); %fit ^2, a = y, x = -1 to 1 (3 points effectively)
    %pval = polyval(p,-1:1); %show the 3 points
    pvalext = polyval(p,linspace(-1,1,500)); %full evaluation, spacing should be large number
    [~,plocs] = findpeaks(pvalext); %find the peak's x value
    Qx = plocs(1);
    %found x value of the peak between 0-500, turn back into 0-2:
    Qx = (Qx - 1) * (x_3(3) - x_3(1))/499 + x_3(1); %needed x value (non-integer)
end

function a = autoc(y,lags)
%autocorrelation with lag range 0 to (lags-1), doesn't divide out variance
    if(length(y) < lags+1)
        y = [y zeros(1,lags+1-length(y))];
    end
    y = y(1:lags+1); %use only L points
    y_len = length(y);
    n = 2^nextpow2(2*y_len-1); %fft runs fastest on 2^n sized data
    a = ifft(fft(y,n) .* conj(fft(y,n))); %acorr = iF(F * F*) = iF(S)
    a = a(1:lags); %only need up to L points
end

function H = calcH(y,lags)
%autocorrelation with lag range 0 to (lags-1), doesn't divide out variance
%only does up to L points for each lag
%check this function again sometime
    endy = lags;
    if(endy > length(y))
        endy = length(y);
    end
    y = y(1:endy); %use only y(0 to lags-1)
    y = [y zeros(1,lags+2)]; %append extra zeros to avoid index errors
    y_len = length(y);
    a = zeros(1,lags); %set up with zeros
    for c = 1:length(a)-1 %check this later
        m = y .* [y(c:end) zeros(1,c-1)]; %y(t) * y(t-c)
        a(c+1) = sum(m(1:c)); %only sum c points, not all points
    end
    H = a;

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