close all
x = linspace(0,100,10000);
y = -0.045.*x.^2 + 4.45.*x - 0.35;
figure()
plot(x(1:length(x)/2), y(1:length(y)/2));
ylim ([0 120])
xlabel('cents')
ylabel('dissonance')
title('-0.045.*x.^2 + 4.45.*x - 0.35;')