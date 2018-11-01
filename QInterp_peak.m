function Qx = QInterp_peak(x_3, y_3)
%finds peak of quadratic interpolation on 3 points
    %V = [1 -1 1; 0 0 1; 1 1 1]; %V for below for x_3 = -1 to 1
    numpoints = 48; %use a decent amount of points
    %fit and evaluate ^2, a = y, x = -1 to 1 (3 points effectively)
    %pvalext = polyfitNval_fastQ(y_3,linspace(-1,1,numpoints),numpoints);
    x = linspace(-1,1,numpoints);
    
    % Solve least squares problem.
    p = [y_3(1)/2 - y_3(2) + y_3(3)/2,-(y_3(1)-y_3(3))/2,y_3(2)];
    %Horner's method:
    y_3 = x .* (x .* ones(1,numpoints) .* p(1) + p(2)) + p(3);
    
    plocs = findpeaks_fast(y_3); %find the peak's x value
    if(isempty(plocs)) %error check
        plocs(1) = round(max(y_3));
        %fprintf("! No peak found in QI function. Set as max, bounds included\n");
    end
    Qx = plocs(1);
    %found x value of the peak between 0-N, turn back into 0-2:
    Qx = (Qx - 1) * (x_3(3) - x_3(1))/(numpoints-1) + x_3(1); %needed x value (non-integer)
    
    function locs = findpeaks_fast(y)
    %custom findpeaks function much faster than stock
    %pks = y values, locs = x values (indices)
        %coder.inline('always') %force inline
        loc = zeros(1,length(y)); %allocate (faster than growing)
        if(length(y) < 3) %require at least 3 points
            return;
        end
        num = 0;
        %I couldn't find a way to make this faster
        %matches stock output well, may have extra peaks at mesas/end
        for i = 2:length(y)-1 %exclude bounds like stock
            if((y(i) - y(i-1) > 0) && (y(i+1) - y(i) <= 0))
                num = num + 1;
                loc(num) = i;
            end
        end
        locs = zeros(1,num);
        locs = loc(1:num); %output right sized array
    end
end