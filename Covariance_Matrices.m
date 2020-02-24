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
function [Sigma,Sigma_R] = Covariance_Matrices(T,discrete_points,d_full,delta,channel_parameters)
%COVARIANCE_MATRICES Summary of this function goes here
%   Detailed explanation goes here
%Covariance matrices

etaSQ=channel_parameters.etaSQ;
sigmaSQ_xi=channel_parameters.sigmaSQ_xi;
gamma=channel_parameters.gamma;
beta=channel_parameters.beta;

TM=zeros(T,T);
for k=1:T
    for l=1:T
        TM(k,l)=exp(-abs(k-l)/gamma); %time lag
    end
end
%covariance matrix for segment WITHOUT relay
Sigma=etaSQ*TM+(sigmaSQ_xi*eye(T));

D=zeros(delta,delta);
for i=1:delta
    for j=1:delta
        D(i,j)=norm(discrete_points(i)-discrete_points(j));
    end
end

A1=etaSQ.*exp( -(D./beta) );

B1=zeros(delta,delta);
for i=1:delta
    for j=1:delta
        if (discrete_points(i)+discrete_points(delta-j+1))>d_full
            B1(i,j)=etaSQ*exp( (abs(discrete_points(i)-discrete_points(j))-98)/beta  );
        elseif (discrete_points(i)+discrete_points(delta-j+1))==d_full
            B1(i,j)=etaSQ*exp( (-98)/beta );
        elseif (discrete_points(i)+discrete_points(delta-j+1))<d_full
            B1(i,j)=etaSQ*exp( (-98-abs(discrete_points(i)-discrete_points(j)))/beta );
        end
    end
end
% issymmetric(B1)
% [~,p] = chol(B1) %if p==0 then positive definite

%Kernel
K=[A1,B1;B1',A1];
% issymmetric(K)
% [~,p] = chol(K)

%covariance matrix for segment WITH relay
Sigma_R=kron(TM,K)+(sigmaSQ_xi*eye(2*delta*T)); %this is positive SEMI definite, symmetric

% Sigma_R=Sigma_R+eye(5000)*(10^-12);
% [~,p] = chol(Sigma_R);
% issymmetric(Sigma_R)
% eig(Sigma_R);

%covariance matrix for segment WITHOUT relay
% SigmaSAA=Sigma+(sigmaSQ_xi*eye(T));

%Sigma_d=[etaSQ.*exp(-D./beta),(etaSQ.*exp((D-98)./beta));(etaSQ.*exp((D-98)./beta)),etaSQ.*exp(-D./beta)];
% [~,p] = chol(Sigma_d);
% A=mvnrnd(zeros(1,100),Sigma_d);
% plot(A,'b-*')

% Sigma_R=zeros(2*delta*T,2*delta*T);
% for k=1:T
%     for l=1:T
%         Sigma_R(1+2*delta*(k-1):(delta*2)+2*delta*(k-1),1+2*delta*(l-1):(delta*2)+2*delta*(l-1))=Sigma_d.*exp(-abs(k-l)/gamma);
%     end
% end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = zeros(length(long_grid_X), length(long_grid_X));
% 
% C = exp(-squareform(pdist([long_grid_X long_grid_Y]))/beta);
% C = eta_sq*C;
% % Add the Source and The Destination to .. (Correlation)
% C_SD = [C *exp(-norm(source-destination)/delta) ; C*exp(-norm(source-destination)/delta) C];
% 
% phi = exp(-1/gamma);
% C_SD_chol = chol((1-phi^2) * C_SD);
% %This is the initial shadowing term at time t=0;
% X_t = (randn(1, length(C_SD(1,:))) * chol(C_SD)).'; %(X_0 initial)
% 
% for i = 1 : T
%     % Innovation...
%     %W_t =  mvnrnd(zeros(1, length(C_SD(1,:))), (1-phi^2) * C_SD).';
%     W_t = (randn(1, length(C_SD(1,:))) *  C_SD_chol).'; %Generate the driving noise of the AR process
%     
%     X_t = phi * X_t + W_t; %Upldate the shadowing term for time t
%     % Add multipath fading...
%     temp  = X_t + sqrt(sigma_xi_sq) .* randn(length(C_SD(1,:)), 1);% Add the multipath fading term
%     % ...
%     temp_source = temp(1:length(temp)/2); %% "shadowing + multipath" term (in dB) of source to "points in grid" channels
%     temp_destination = temp(length(temp)/2+1:end); %% shadowing + multipath term (in dB) of "points in grid" to destination channels
%     sh_maps_source{i} = reshape(temp_source, grid_analysis, grid_analysis) + path_loss_source; % "F" channels (in dB) from source to "points in grid" (after adding the respective path loss)
%     sh_maps_destination{i} = reshape(temp_destination, grid_analysis, grid_analysis) + path_loss_destination;%"G" channels (in dB) from "points in grid" to destination (after adding the respective path loss)
% end
% 
% 
