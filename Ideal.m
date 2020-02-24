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
function [R_cord,round_pos] = Ideal(t,R_cord,N_R,R_orient,R_bounds,z,z_r,phi,phi_r,L_f,L_g,I_cord,full_seg,segidx_f,segidx_g,R_segments,channel_parameters,S_cord,D_cord,P_S,P_C,Delta,delta,discrete_points,I_S,d_full)
%This function implements the "ideal" relay selection scheme.

% Input parameters
%   - C: the setup parameters
%   -
%   -

%
% Output parameters
%   - R_cord: the new relay coordinates for all clusters
%   - round_pos: the new relay positions in the segment for all clusters
%
%%

a_LoS=channel_parameters.a_LoS;
a_NLoS=channel_parameters.a_NLoS;
sigmaSQ=channel_parameters.sigmaSQ;
sigmaSQ_D=channel_parameters.sigmaSQ_D;

round_pos=zeros(N_R,1);  %holds antenna positions of all relays in segment
for r=1:N_R  %for every relay
    
    VI=zeros(delta,1); %holds respective values of V_I for every considered relay position
    R_cord_temp=R_cord(r,:); %coordinates of examined relay
    
    for relayPos=1:delta %every possible discrete position in segment
        
        R_cord_temp(1,R_orient(r,1))=R_bounds(r,1)-discrete_points(1,1)+ discrete_points(1,relayPos); %temporary examined relay coordiante
        
        %% Channel f: from p_S to p_r
        dLoS_f=norm(S_cord-I_cord(I_S,:)); %LoS distance
        dNLoS_f=norm(I_cord(I_S,:)-R_cord_temp(1,:),1); %dNLoS: Manhattan distance from first intersection after p_S to p_r
        aF=-(a_LoS*10*log10(dLoS_f)) - (a_NLoS*full_seg(r,1)*10*log10(d_full)) - (a_NLoS*10*log10(mod(dNLoS_f,d_full))) - (Delta*(full_seg(r,1)+1));
%         if pathlossF(r,relayPos)-aF
%             check=4;
%         end
            
        f=zeros(L_f(r,1),1); %channel f vector
        for i=1:L_f(r,1)
            testidf=segidx_f(i,1:full_seg(r,1),r);
            Phi_F=phi(end-1,t) + sum( phi(testidf,t)) + phi_r(r,relayPos+(2*delta*(t-1))); %phase of segment with p_S + phase of other segments in path + phase of segment with p_r
            zF= z(end-1,t) + sum(z(testidf,t)) + z_r(r,relayPos+2*delta*(t-1)); %beta(end-1,(2*delta-Spos)+2*delta*(t-1))
            %             cF=xi(end-1,t)+sum(xi(segidx_f(i,1:full_seg(r,1),r),t))+xi_r(r,dRpos+(2*delta*(t-1)));
            
            F=aF+zF; %channel magnitude in dB scale
            f(i,1)=(10^(F/20))*exp(2*pi*1i*Phi_F); %convert magnitude to linear scale, f(i,1)=exp((chi/2)*F)*exp(2*pi*1i*Phi_F)
        end
        
        %% Channel g: from p_r to p_D
        dLoS_g=norm(R_cord_temp(1,:)-I_cord(R_segments(r,2),:)); %LoS distance
        dNLoS_g=norm(I_cord(R_segments(r,2),:)-D_cord,1);  %dNLoS: Manhattan distance from first intersetion after p_r to p_D
        aG=-(a_LoS*10*log10(dLoS_g)) - (a_NLoS*full_seg(r,2)*10*log10(d_full)) -(a_NLoS*10*log10(mod(dNLoS_g,d_full))) - (Delta*(full_seg(r,2)+1)) ;
%         if pathlossG(r,relayPos)-aG
%             check=4;
%         end
        g=zeros(L_g(r,1),1); %channel g vector
        for i=1:L_g(r,1)
            testidg=segidx_g(i,1:full_seg(r,2),r);
            Phi_G=phi_r(r,(relayPos+delta)+2*delta*(t-1))+sum(phi(testidg,t))+phi(end,t);
            zG=z_r(r,(relayPos+delta)+2*delta*(t-1))+sum(z(testidg,t))+z(end,t);%+beta(end,Dpos+2*delta*(t-1));
            %             cG=xi_r(r,(dRpos+delta)+2*delta*(t-1))+sum(xi(segidx_g(i,1:full_seg(r,2),r),t))+xi(end,t);
            
            G=aG+zG;
            g(i,1)=(10^(G/20))*exp(2*pi*1i*Phi_G); %g(i,1)=exp((chi/2)*G)*exp(2*pi*1i*Phi_G);
        end
        
        VI_numerator=P_C*P_S*abs(ones(L_f(r,1),1)'*f(1:L_f(r,1),1))^2*abs(ones(L_g(r,1),1)'*g(1:L_g(r,1),1))^2;
        VI_denominator=(P_S*sigmaSQ_D*abs(ones(L_f(r,1),1)'*f(1:L_f(r,1),1))^2 ) + (P_C*sigmaSQ*abs(ones(L_g(r,1),1)'*g(1:L_g(r,1),1))^2) + (sigmaSQ*sigmaSQ_D);
        VI(relayPos,1)=VI_numerator/VI_denominator;
        
    end %end every possible relay positon
    
    [~, I]=sort(VI,'descend'); %sorted values and respective index
    %     if I(1)==0
    %         fprintf('Max VI in Oracle is zero!')
    %     end
    R_cord(r,R_orient(r,1))=  R_bounds(r,1)-discrete_points(1,1)+discrete_points(1,I(1));
    round_pos(r,1)=I(1);   %=find(discrete_points==(mod( R_cord(r,R_orient(r,1)),d_full))); %position of antenna in segment
end %end for every relay

end