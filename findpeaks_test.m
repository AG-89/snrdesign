clc,clear variables
y = [0 0 0 3 3 3 0 0 1 1 10 1 1 0 0 20 2 2 2];
%y = [0 1 2 3 4 5 4 3 7 3 2 3 4 5 4 3 2 0 3 1 3 4 2 5 3 3 4 2 9 2 3 3 5 3 2 3 4 5 4 3 2 0 3 1 3 4 6 2 5 3 3 4 2 2 3 3 5 4 3 2];

for a = 1:20000
% [Y,X] = findpeaks(y);
[Y,X] = findpeaks_fast(y);
[Yo,Xo] = findpeaks_fasto(y);
end
[Y,X] = findpeaks_fast(y)
[Yo,Xo] = findpeaks_fasto(y)

function [pks,locs] = findpeaks_fast(y)
%custom findpeaks function much faster than stock
%pks = y values, locs = x values (indices)
    pks = zeros(1,length(y)); %allocate (faster than growing)
    locs = pks;
    if(length(y) < 3) %require at least 3 points
        return;
    end
    num = 0;
    %I couldn't find a way to make this faster
    %matches stock output well, may have extra peaks at mesas/end
    for i = 2:length(y)-1 %exclude bounds like stock
        if((y(i) - y(i-1) > 0) && (y(i+1) - y(i) <= 0))
            num = num + 1;
            pks(num) = y(i);
            locs(num) = i;
        end
    end
    pks = pks(1:num); %output right sized array
    locs = locs(1:num);
end

%all matrix way. dxprev allocation takes massive time unfortunately
% function [pks,locs] = findpeaks_fasto(y)
% %custom findpeaks function much faster than stock
% %pks = y values, locs = x values (indices)
%     if(length(y) < 3) %require at least 3 points
%         return;
%     end
%     %matches stock output well, may have extra peak at end
%     dx = diff(y); %dx
%     dxprev = [0 dx(1:end-1)]; %dx previous
%     ispeak = (dx > 0) & (dxprev <= 0); %peak locations
%     pks = ispeak .* (y(2:end)); %zero out non peaks
%     locs = ispeak .* (2:length(y)); %zero out non locs
%     
%     pks = pks(pks > 0);
%     locs = locs(locs > 0);
% end
