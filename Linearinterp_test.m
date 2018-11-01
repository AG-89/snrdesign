clc, close all, clear variables, clear sound
samplerate = 44100;
f = 500;
p = 1/f;
spp = samplerate*p;
N = ceil(spp); %# of points
x = 0:N-1;
y = square(2.*pi.*f.*x/samplerate);
y = repmat(y,1,1); 

ratio = 1.25;

%graph
figure()
plot(0:length(y)-1,y,'.-','MarkerSize',20)
% hold on
figure()
y_new = pitch_shift_Npt_linterp(y,ratio,false);
plot(0:length(y_new)-1,y_new,'.-','MarkerSize',20)

fprintf("ylen = %d, y_newlen = %d\n\n",length(y),length(y_new))


%play
volume = 2;
volpeakpoint = 100/volume;
y_toplay = [volpeakpoint repmat(y,1,ceil(samplerate*1.5/length(y)))];
fprintf("Playing original 1.5s...\n");
soundsc(y_toplay,(samplerate))
pause
clear sound, clear sounds;
y_toplay = [volpeakpoint repmat(y_new,1,ceil(samplerate*1.5/length(y_new)))];
fprintf("Playing pitch shifted 1.5s...\n\n");
soundsc(y_toplay,(samplerate))

function y_new = pitch_shift_Npt_linterp(y,ratio,retain_length)
%pitch shifts N points by a ratio through stretch or compress
%will create discontinuities/jumps with retain_length, iffy behavior
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