%% rotation.m
% function which returns the planar (x-y plane) rotation given an angle t
function r = rotation(t)
%ROTATION Summary of this function goes here
%   Detailed explanation goes here
    r = [cos(t) -sin(t) 0;
         sin(t)  cos(t) 0;
         0        0     1];
end

