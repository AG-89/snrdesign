
close all
%Decimate by a fractional value
fs = 200;
t = 0:1/(20*fs):20/fs;
x = 0.75*sin(2*pi*fs*t);
p=1001;
q=1997;
y = resample(x,p,q);
tnew = 0:1/(20*fs):(20*p)/(fs*q);
plot(t,x,tnew,y);
figure
%Interpolate by a fractional value
y2 = resample(x,q,p);
tnew2 = 0:1/(20*fs):((20*q)/(fs*p)+1/(20*fs));
plot(t,x,tnew2,y2);
