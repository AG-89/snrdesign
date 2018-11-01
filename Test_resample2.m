clc, close all, clear variables
%Decimate by a fractional value
samplerate = 44100;
fs = 261.1;
t = 0:1/(20*fs):20/fs;
x = 0.75*sin(2*pi*fs*t);
p=100;
q=125;
y = resample(x,p,q);
lenx = length(x)
leny = length(y)
plot(t(1:leny),x(1:leny),t(1:leny),y);
figure
%Interpolate by a fractional value
p2=q;
q2=p;
y2 = resample(x,p2,q2);

plot(t,x,t,y2(1:lenx));

leny2 = length(y2)

[~,~] = mkdir('test/');
audiowrite('test/X.wav',y,samplerate);
audiowrite('test/Y.wav',y,samplerate);
audiowrite('test/Y2.wav',y2,samplerate);