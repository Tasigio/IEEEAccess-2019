%% Version
% Last revision: February 2019 (Matlab R2019b)
% Author: Anastasios Dimas
%
%% Purpose
% The purpose of this code is to support the published paper:
% "Cooperative Beamforming with Predictive Relay Selection for Urban mmWave Communications"
%  by Anastasios Dimas, Dionysios S. Kalogerias, and Athina P. Petropulu.
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
function [SAA_samples] = Generate_SAAsamples(SAA_trials,t,Sigma,data_bld,sample_size)
% This function generates the Sample Average Approximation "shadowing+multipath" 
% samples for segments WITHOUT cluster, the segment the source is located,
% and the segment the destination is located on.


% Input parameters
%   - SAA_trials:
%   - t: Time slot
%   - Sigma
%   - data_bld
%   - sample_size

%
% Output parameters
%   - SAA_samples:

%%
    muData_bld=(Sigma(t+1,1:t)/Sigma(1:t,1:t))*data_bld;
    Sigma_bld=(Sigma(t+1,t+1))-((Sigma(t+1,1:t)/Sigma(1:t,1:t))*Sigma(t+1,1:t)');  %  [~,p] = chol(Sigma_bld);
 
    SAA_samples=mvnrnd(muData_bld, sample_size*Sigma_bld, SAA_trials); %shadowing samples (SAA trials x #segments)
    %Sigma_DL= C.etaSQ.*exp(-(t:-1:1)./C.gamma);


end

