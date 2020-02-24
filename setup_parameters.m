%% Version
% Last revision: February 2020 (Matlab R2019b)
% Author: Anastasios Dimas
%
%% Purpose
% The purpose of this code is to support the published paper: 
%
% Anastasios Dimas, Dionysios S. Kalogerias, and Athina P. Petropulu,
% "Cooperative Beamforming With Predictive Relay Selection for Urban mmWave Communications", 
% IEEE Access, 7, 2019.
%
% Any part of this code used in your work should cite the above publication.
%
% This code is provided "as is" to support the ideals of reproducible research. Any issues with this
% code should be reported by email to tasos.dimas@rutgers.edu. However, no guarantees are being made
% that the reported issues will be eventually fixed.
%
% The code is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
% available at https://creativecommons.org/licenses/by-nc-sa/4.0/
%
%%
%%
%
%   This file contains the parameters used in the simulations
%
%%

%% Topology
d_full=100; %length of segment (meters)
x=(100:d_full:500); %x-cordinate grid dimension
y=(100:d_full:300); %y-cordinate grid dimension

%% Segment
%segment discritization
delta=50; %how many parts to divide a segment into

%discrete coordinates in segment
% discrete_points=((setup_parameters.d_full/setup_parameters.delta)/2:setup_parameters.d_full/setup_parameters.delta:setup_parameters.d_full);
discrete_points= linspace(1, d_full-1, delta);
%spacing=d_full/delta; %distance between antennas (in meters)

%% Source coodinate(s)
S_cord=[157,300];
% S_cord=[150,300];
S_orient=1; %Source orientation
I_S=6;     %Intersection index after source
%block source is located in
Spos= find(discrete_points==(mod(S_cord(:,S_orient),d_full)));

% Destination coordinate(s)
D_cord=[457,100];
%  D_cord=[450,100];
D_orient=1; %Destination orientation
I_D=10;     %Intersection index before destination
%block destination is located in
Dpos= find(discrete_points==(mod(D_cord(:,D_orient),d_full)));

% displacement from current position in neighborhood (meters)
% relay_pos=[-8*spacing:spacing: 8*spacing]';

%% Channel Parameters
% Path-loss exponents
% a_LoS=2.1;
% a_NLoS=2.1;

%         P_S=db2pow(100); %Source power (mag2db)
P_S=10^((80-30)/10); %dbm to power 30
P_C=10^((100-30)/10); %Total relay power
Delta=10;%10; %corner loss (in dB)

channel_parameters = struct('a_LoS',2.1,'a_NLoS',2.1,'etaSQ',40,'beta',10,'gamma',15,'sigmaSQ',1,'sigmaSQ_D',1,'sigmaSQ_xi',20);

% sigmaSQ=1; %relay noise
% sigmaSQ_D=1; %destination noise
% sigmaSQ_xi=20; %multipath variance

%shadowing parameters
% etaSQ=40; %shadowing power
% beta=10;%20; %correlation distance
% gamma=15; %correlation time


%% Simulation
T=50; %time slots to run the program
SAA_trials=500;
total_trials=5000;
%c=physconst('LightSpeed'); %speed of light

