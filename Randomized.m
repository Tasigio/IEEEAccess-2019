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
function [R_cord,round_pos] = Randomized(discrete_points,delta,R_cord,N_C,C_orient,C_bounds)
% This function implements the "randomized" relay selection approach. 
% Every relay is chosen randomly from the cluster withoug taking into
% considertion the observed CSI.

% Input parameters
%   - discrete_points: Coordinates of relays in cluster
%   - delta: Number of disrcete relay positions in cluster
%   - R_cord:
%   - N_C: Number of clusters
%   - C_orient: Orientation of cluster
%   - C_bounds: % Coordinate bounds of segment that cluster is placed

% Output parameters
%   - R_cord: the new relay coordinates for all clusters
%   - round_pos: the new relay positions in the segment for all clusters

%%

round_pos=zeros(N_C,1); %holds respective relay positions of all clusters
for r=1:N_C % for each cluster
    
    round_pos(r,1)=randi(delta); %choose a random relay position from the segment
    R_cord(r,C_orient(r,1))= C_bounds(r,1)-discrete_points(1,1)+discrete_points(1,round_pos(r,1)); %coordinate of new position
    %round_pos(r,1)=find(discrete_points==(mod( R_cord(r,R_orient(r,1)),d_full)));
    
end


end

