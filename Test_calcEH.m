clc, clear variables

lags = 880;
minlag = 0;
y = fliplr(1:1760);
for a = 1:1
    %autoc(y,lags,minlag);
    %calcHnew(y,lags,minlag);
    %calcH(y,lags,minlag);
     calcEnew(y,lags,minlag);
     calcE(y,lags,minlag);
end
all(zeros(1,lags-minlag) == calcEnew(y,lags,minlag) - calcE(y,lags,minlag))

function H = calcHnew(y,lags,minlag)
%autocorrelation with lag range [minlag, lags-1], doesn't divide out variance
%only does up to L points for each lag
    minlag = minlag + 1; %calculate only from this lag onward (0 -> 1)
    %y = [zeros(1,lags+2) y]; %zero padding on left to stop array index errors
    ylen = length(y);
    H = zeros(1,lags); %set up with zeros
    c = minlag:lags;
        %sum of y_i * y_i-L
        %todo: look at vectorization possibility
        ylen-c+1
        ylen-c+1:ylen
        a = y(ylen-c+1:ylen) .* y(ylen-c+1-c+1:ylen-c+1);
        H = sum(a);
    H = H(minlag:end); %only return range
end

function H = calcH(y,lags,minlag)
%autocorrelation with lag range [minlag, lags-1], doesn't divide out variance
%only does up to L points for each lag
    minlag = minlag + 1; %calculate only from this lag onward (0 -> 1)
    %y = [zeros(1,lags+2) y]; %zero padding on left to stop array index errors
    ylen = length(y);
    H = zeros(1,lags); %set up with zeros
    for c = minlag:lags
        %sum of y_i * y_i-L
        %todo: look at vectorization possibility
        H(c) = sum(y(ylen-c+1:ylen) .* y(ylen-c+1-c+1:ylen-c+1));
    end
    H = H(minlag:end); %only return range
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
    ylen = length(y);
    y_squared = y.^2; %y^2 = y * y
    E = zeros(1,lags); %set up with zeros
    %get E(minlag-1):
    prevEc_minlagpart = y_squared(ylen-(minlag-1)*E_pmult+1:ylen);
    prevEc = (minlag>1)*sum(prevEc_minlagpart);
    for c = minlag:lags
        %recursive sum + next two y^2 points
        E(c) = prevEc + y_squared(ylen-(c)*E_pmult+1) + y_squared(ylen-(c-1)*E_pmult);
        prevEc = E(c);
    end
    E = E(minlag:end); %only return range
end

%again slower than for loop. maybe try and eliminate fliplr
function E = calcEnew(y,lags,minlag)
%energy (no lag), with input range 2*[minlag, (lags-1)]
%only does up to L points for each lag
    minlag = minlag + 1; %calculate only from this lag onward (0 -> 1)
    E_pmult = 2;% * L (= 2L used)
    ylen = length(y);
    y_squared = y.^2; %y^2 = y * y
    E = zeros(1,lags); %set up with zeros
    %get E(minlag-1):
    prevEc_minlagpart = y_squared(ylen-(minlag-1)*E_pmult+1:ylen);
    prevEc = (minlag>1)*sum(prevEc_minlagpart);
    c = minlag:lags;
    a = y_squared(ylen-(lags)*E_pmult+1:ylen-(minlag-1)*E_pmult);
    cs = cumsum(a,'reverse');
    E = fliplr(cs(1:E_pmult:end)) + prevEc; %add last sum before minlag
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