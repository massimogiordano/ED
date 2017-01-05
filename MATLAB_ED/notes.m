tt = 1:500;
par = 0 + (tt/500 - 0.7).^4;
gan = -0.7./(tt+6);
tot = [par gan];
plot(tot)
Y = diff(tot,2)
hold on
plot(1*Y)