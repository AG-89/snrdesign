clc, close all, clear variables
x = 0:7;
y = sin(2*pi*1/8*x);

figure()
plot(x,y,'.-','MarkerSize',20)
y_new = linear_interp(y,2);
hold on
plot(x,y_new,'.-','MarkerSize',20)

%not finished yet
function y_out = linear_interp(y,ratio)
    x = 1:length(y);
    x_adjusted = x.*ratio
    y_out = y;
    lin_eqs = ones(1,length(y)-1); %pt slope eqs
    k = 1:length(y)-1
        dy = y(k+1) - y(k);
        dx = x_adjusted(k+1) - x_adjusted(k);
        slope = (dy)./(dx)
        lin_eqs(k) = slope.*(dx) + y(k) %y values
    if(ratio < 1) %repeat
        return
    else %cutoff
        return
    end
end
