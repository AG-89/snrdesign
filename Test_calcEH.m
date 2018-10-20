clc, clear variables

y = -100000:100000;
calcH(y,10000,0)
calcE(y,10000,0)

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

    %E: BR triangle matrix, too slow to use
%     a = ones(length(y_squared)) .* y_squared; %repeat rows for square
%     t = flipud(triu(a)); %bottom right triangle matrix
%     t = t(E_pmult*minlag:E_pmult:end,:); %for L=2, use only even rows, starting at appropriate minlag
%     E = [zeros(1,minlag-1) sum(t,2)']; %pad any missing lags w/ zeros

%old functions that prob don't return correct results
% function E = calcE(y,lags,minlag)
% %energy (no lag), with input range (minlag to (lags-1)) * 2
%     minlag = minlag + 1; %ind0 -> ind1
%     E_pmult = 2;% * L (= 2L used)
%     E_lag = lags*E_pmult; %2L
%     endy = E_lag*2;
%     if(endy > length(y))
%         endy = length(y);
%     end
%     y = y(1:endy); %use only y(1 to lags*2)
%     y = [y zeros(1,E_lag)+2]; %zero padding to stop array index errors
%     E = zeros(1,lags); %prep, only doing limited lags
%     y_squared = y .* y; %y^2 = y * y
%     E(1) = sum(y_squared(1:2)); %Error check this recursive optimization
%     for c = (minlag+1):length(E)
%         %only sum c points, not all points
%         %E(c) = sum(y_squared(1:c*E_pmult));
%         E(c) = E(c-1) + y_squared(c*E_pmult-1) + y_squared(c*E_pmult);
%     end
% end
% 
% function H = calcH(y,lags,minlag)
% %autocorrelation with lag range minlag to (lags-1), doesn't divide out variance
% %only does up to L points for each lag
% %check this function again sometime
%     minlag = minlag + 1;
%     endy = lags*2;
%     if(endy > length(y))
%         endy = length(y);
%     end
%     y = [zeros(1,lags+2) y]; %zero padding on right to stop array index errors
%     H = zeros(1,lags); %set up with zeros
%     for c = minlag:length(H)
%         %m = y_i * y_i-L = y(t) * y(t+c)
%         m = y(end-c+1:end) .* y(end-c+1-c:end-c); %y * y shifted right by L
%         H(c) = sum(m); %only sum c points, not all points
%     end
% end