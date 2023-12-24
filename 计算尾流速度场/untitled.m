clc;clear;
t = 0:2*pi/1000:2*pi;
f = 10;
y = sin(2*pi*f*t);
Y = fft(y);
yy = ifft(Y);
abs_Y = abs(Y);