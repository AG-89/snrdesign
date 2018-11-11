clc, clear variables, close all
f = 261.6;
T = f^-1;
samplerate = 44100;
x = 1:2000;
y = sin(2*pi*f.*x/samplerate);
volume = 25;

pitch_scale = 1.25;
[n, d] = rat(pitch_scale);
N_D = [n, d];

figure()
plot(y)

y_new = doPSOLA(y, samplerate, N_D);

if(false) %plays the sounds
    volpeakpoint = 100/volume;
    y_toplay = [volpeakpoint y];
    fprintf("Playing y...\n");
    soundsc(y_toplay,(samplerate))
    pause
    clear sound, clear sounds;
    y_toplay = [volpeakpoint y_new];
    fprintf("Playing y pitch shifted...\n\n");
    soundsc(y_toplay,(samplerate_tuned))
    %pause
end

function y_out = doPSOLA(y, samplerate, N_D)
    N = N_D(1); %num and denom
    D = N_D(2);
    T = 261.6^-1; % period
    T_new = T / N * D;
    ylen = length(y);
    %later: find first peak offest and rest of them will be T apart
    [peaks, locs] = findpeaks(y); %update later, find pitch peaks
    
    T_s = round(samplerate * T); %period in samples (rounded)
    T_new_s = round(samplerate * T_new);
    
    first = true;
    windows = [];
    
    %overlap windows by 1/2
    %check rounding behavior very carefully later
    beginning_boundary = round(locs(1) + T_s/2);
    b = beginning_boundary;
    for i = 1:floor(length(y)/(T_s/2))
        lb = b;
        ub = b+T_s-1;
        zerolen = 0;
        %check start/end
        if(ub > ylen)
            ub = ylen;
        end
        if(lb < 1)
            lb = 1;
        end
        
        [lb,ub, ub-lb+1];
        w = y(lb:ub);
        
        lb = b;
        ub = b+T_s;
        if(ub > ylen)
            w = [w, zeros(1,T_s - length(w))];
        end
        if(lb < 1)
            w = [zeros(1,T_s - length(w)), w];
        end
        windows = [windows; w];
        b = b + round(T_s/2) - 1; %overlap end->start point
    end
    windows = [[zeros(1,T_s-beginning_boundary), y(1:beginning_boundary)]; windows]
    
    
    y_out = y;
end


%this doesnt overlap windows
%     beginning_boundary = round(locs(1) + T_s/2);
%     b = beginning_boundary;
%     for i = 1:length(locs)
%         lb = b;
%         ub = b+T_s-1;
%         zerolen = 0;
%         %check start/end
%         if(ub > ylen)
%             ub = ylen;
%         end
%         if(lb < 1)
%             lb = 1;
%         end
%         
%         w = y(lb:ub);
%         
%         lb = b;
%         ub = b+T_s;
%         if(ub > ylen)
%             w = [w, zeros(1,T_s - length(w))];
%         end
%         if(lb < 1)
%             w = [zeros(1,T_s - length(w)), w];
%         end
%         windows = [windows; w];
%         b = b + T_s - 1; %overlap end->start point
%     end
%     windows = [[zeros(1,T_s-beginning_boundary), y(1:beginning_boundary)]; windows]