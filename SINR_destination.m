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
function [Vopt] = SINR_destination(t,R_pos,R_cord,z,z_r,phi,phi_r,L_f,L_g,N_R,I_cord,full_seg,segidx_f,segidx_g,R_segments,channel_parameters,S_cord,D_cord,P_S,P_C,Delta,delta,I_S,d_full)
%This function implements the evaluation of the SINR.

% Input parameters
%   - C: the setup parameters
%   - 
%   -

%
% Output parameters
%   - Vopt:

a_LoS=channel_parameters.a_LoS;
a_NLoS=channel_parameters.a_NLoS;
sigmaSQ=channel_parameters.sigmaSQ;
sigmaSQ_D=channel_parameters.sigmaSQ_D;

f=zeros(max(L_f),N_R);
g=zeros(max(L_g),N_R);
for r=1:N_R     %for every relay
    
    
    %% Channel f: from p_S to p_r
    %dRpos=find(discrete_points==(mod(R_cord(r,R_orient(r,1)),d_full)));  %find discrete position of relay in segment
    
    dLoS_f=norm(S_cord-I_cord(I_S,:)); %LoS distance
    dNLoS_f=norm(I_cord(I_S,:)-R_cord(r,:),1); %Manhattan distance from first intersection after p_S to p_r (d=sum_1^n(|x_i-y_i|)
    aF= -(a_LoS*10*log10(dLoS_f)) - (a_NLoS*full_seg(r,1)*10*log10(d_full)) - (a_NLoS*10*log10(mod(dNLoS_f,d_full))) - (Delta*(full_seg(r,1)+1));
    
    for i=1:L_f(r,1)
        testidf=segidx_f(i,1:full_seg(r,1),r);
        Phi_F=phi(end-1,t)+sum(phi(testidf,t))+phi_r(r,R_pos(r,1)+(2*delta*(t-1)));%+phi(R_indseg(r,1),t); %segment with source + other segments + relay segment
        zF=z(end-1,t) + sum(z(testidf,t))+z_r(r,R_pos(r,1)+2*delta*(t-1)); %(2*delta-Spos)+2*delta*(t-1)
    
        F=aF+zF; %channel magnitude in dB scale
        f(i,r)=(10^(F/20))*exp(2*pi*1i*Phi_F);%convert dB to linear scale, =exp((chi/2)*F)*exp(2*pi*1i*Phi_F); 
      %  f_check(i,r)=10^(F/20)*exp(2*pi*1i*Phi_F);
    end
    
    %% Channel g: from p_r to p_D
    dLoS_g=norm(R_cord(r,:)-I_cord(R_segments(r,2),:)); % LoS distance
    dNLoS_g=norm(I_cord(R_segments(r,2),:)-D_cord,1); % Manhattan distance from first intersetion after p_r to p_D
    aG= -(a_LoS*10*log10(dLoS_g)) - (a_NLoS*full_seg(r,2)*10*log10(d_full)) - (a_NLoS*10*log10(mod(dNLoS_g,d_full))) - (Delta*(full_seg(r,2)+1));
    
    %for each segment in path that does not contain relay
    for i=1:L_g(r,1)
  
        testidg=segidx_g(i,1:full_seg(r,2),r);
        Phi_G=phi_r(r,(R_pos(r,1)+delta)+2*delta*(t-1))+sum(phi(testidg,t))+phi(end,t); %phi(R_indseg(r,1),t)
        zG=z_r(r,(R_pos(r,1)+delta)+2*delta*(t-1))+sum(z(testidg,t))+z(end,t);%+beta(end,Dpos+2*delta*(t-1));
        %cG=xi_r(r,(dRpos+delta)+2*delta*(t-1))+sum(xi(segidx_g(i,1:full_seg(r,2),r),t))+xi(end,t);%xi(R_indseg(r,1),t)
        
        G=aG+zG; %channel magnitude in dB scale
        g(i,r)=(10^(G/20))*exp(2*pi*1i*Phi_G);%convert magnitude to linear scale,g(i,r)=exp((chi/2)*G)*exp(2*pi*1i*Phi_G)
    end

    
end

%Defined vectors and matrices for the optimization
h_vector=zeros(N_R,1);
for r=1:N_R
    dbstop if error
    h_vector(r,1)=ones(L_g(r,1),1)'*g(1:L_g(r,1),r)*ones(L_f(r,1),1)'*f(1:L_f(r,1),r);
%         
%     VI_numerator=P_S*(ones(L_g(r,1),1)'*conj(g(1:L_g(r,1),r)))*(ones(L_f(r,1),1)'*conj(f(1:L_f(r,1),r)));
%     VI_denominator=(P_S*sigmaSQ_D*(abs(ones(L_f(r,1),1)'*f(1:L_f(r,1),r))^2)) + (P_C*sigmaSQ*(abs(ones(L_g(r,1),1)'*g(1:L_g(r,1),r))^2)) + (sigmaSQ*sigmaSQ_D) ;
%     All_VI(r) = VI_numerator/VI_denominator;    
end

% BFdenom=norm(All_VI);



D_matrix=zeros(N_R,N_R);
Q_matrix=zeros(N_R,N_R);

for r=1:N_R
    D_matrix(r,r)=P_S*abs(ones(L_f(r,1),1)'*f(1:L_f(r,1),r))^2 +sigmaSQ;
    Q_matrix(r,r)=sigmaSQ*abs(ones(L_g(r,1),1)'*g(1:L_g(r,1),r))^2;
    %check(r,1)=abs(ones(L_f(r,1),1)'*f(1:L_f(r,1),r))^2;
end

R_matrix=P_S*(h_vector*h_vector');
e=eig( ((sigmaSQ_D*eye(N_R)+P_C*D_matrix^(-0.5)*Q_matrix*D_matrix^(-0.5))^(-1))*D_matrix^(-0.5)*R_matrix*D_matrix^(-0.5) );

%Optimal value of SINR
Vopt=P_C*max(real(e)); %P_C*max(e);


end

