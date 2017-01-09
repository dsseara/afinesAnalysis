% Script for creating an interpolated velocity field from
% position of actin beads from AFiNeS simulation

% Daniel Seara 1/7/2017

clear;
close all;

% Returns simdata.mat if doesn't exist already
%run([pwd, 'read_data2.m']);
fname = 'simdata.mat';
vars = {'adata', 'timestep'};
load(fname, vars{:});

% Change this if you want to find velocities for steps longer than dt_frame
numTimeSteps = 1;

xyDisp = diff(adata(:,1:2, 1:numTimeSteps:end),1,3);
[numBeads, dim, numFrames] = size(xyDisp);
dt = uniquetol(diff(timestep));
dt = repmat(dt, size(xyDisp));

adata_vel = xyDisp ./ dt;