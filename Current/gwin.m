%function [ramp] = gwin(deg,winsize,duty,delay)
clc;
clear all;
close all;

winsize=4000;
deg=0.05;
duty=0.7;

slope=deg*winsize;

up=linspace(0,1,slope);
down=linspace(1,0,slope);
mid=ones(1,(round(duty*winsize)-(round(2*slope,0))));
z=zeros(1,winsize);
step=[up,mid,down];


ramp= z+[step,zeros(1,winsize-length(step))];
hold on
plot(ramp,'b');axis([-100 4100 -0.1 1.1])
hold off
%end

