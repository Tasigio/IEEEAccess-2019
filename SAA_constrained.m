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
function [R_cord,round_pos] = SAA_constrained(t,SAA_trials,R_cord,N_R,R_orient,R_bounds,L_f,L_g,I_cord,full_seg,segid_f,segid_g,R_segments,SAAz,SAAphi,z_r,Sigma_R,chi,SAApos_ALL,channel_parameters,S_cord,D_cord,P_S,P_C,Delta,delta,discrete_points,I_S,d_full,considered_neigh)
%UNTITLED2 Summary of this function goes here

a_LoS=channel_parameters.a_LoS;
a_NLoS=channel_parameters.a_NLoS;
sigmaSQ=channel_parameters.sigmaSQ;
sigmaSQ_D=channel_parameters.sigmaSQ_D;

round_pos=zeros(N_R,1);

% considered_neigh=[1:4,delta-4:delta]; %considered relay neighborhood
for r=1:N_R

    testidgA=1:full_seg(r,1);
    testidgB=1:full_seg(r,2);

    All_VI=zeros(SAA_trials,length(considered_neigh)); %keeps all VI values for all samples
    R_cord_temp=R_cord(r,:); %coordinates of examined relay
    
    for relaypos=1:length(considered_neigh)   % for every possible relay position in constrained neighborhood
        dRpos=considered_neigh(1,relaypos); %temporary examined discrete position in the segment
        R_cord_temp(1,R_orient(r,1))=  R_bounds(r,1)-discrete_points(1,1)+discrete_points(1,dRpos); %and its relay coordinate
        
        
        %% Channel f: from p_S to p_r
        dLoS_f=norm(S_cord-I_cord(I_S,:)); %LoS distance
        dNLoS_f=norm(I_cord(I_S,:)-R_cord_temp(1,:),1); %Manhattan distance from first intersetion after SRC to R_i
        aF=-(a_LoS*10*log10(dLoS_f)) - (a_NLoS*full_seg(r,1)*10*log10(d_full)) - (a_NLoS*10*log10(mod(dNLoS_f,d_full))) - (Delta*(full_seg(r,1)+1));
        
        %% Channel g: from p_r to p_D
        dLoS_g=norm(R_cord_temp(1,:)-I_cord(R_segments(r,2),:)); %LoS distance
        dNLoS_g=norm(I_cord(R_segments(r,2),:)-D_cord,1); %Manhattan distance from first intersetion after R_i to D
        aG=-(a_LoS*10*log10(dLoS_g)) - (a_NLoS*full_seg(r,2)*10*log10(d_full)) - (a_NLoS*10*log10(mod(dNLoS_g,d_full))) - (Delta*(full_seg(r,2)+1));
        
        %SAA samples for relay segment
        mFG_bld=zeros(2*t,1);
        cnt=1;
        for i=1:t
            mFG_bld(cnt:cnt+1,1)= [(z_r(r,SAApos_ALL(r,i)+2*delta*(i-1)))',(z_r(r,(SAApos_ALL(r,i)+delta)+2*delta*(i-1)))'];
            cnt=cnt+2;
        end
        
        positions=[SAApos_ALL(r,1:t),dRpos];
        Q={};
        for i=1:t+1
            Stemp=[Sigma_R(positions(i),positions(i)+2*delta*(i-1)), Sigma_R(positions(i),(positions(i)+delta)+2*delta*(i-1)) ;
                Sigma_R((delta+positions(i)),positions(i)+2*delta*(i-1)), Sigma_R((delta+positions(i)),(delta+positions(i))+2*delta*(i-1))];
            Q{i}=Stemp;
        end
        SigmaSAA = cell2mat(Q(toeplitz(1:t+1)));
        
        muFG_bld=(SigmaSAA(2*t+1:end,1:2*t)/SigmaSAA(1:2*t,1:2*t))*mFG_bld;
        SigmaFG_bld= SigmaSAA(2*t+1:end,2*t+1:end)-((SigmaSAA(2*t+1:end,1:2*t)/SigmaSAA(1:2*t,1:2*t))*SigmaSAA(2*t+1:end,1:2*t)');
        
        SAAz_r=mvnrnd(muFG_bld, SigmaFG_bld, SAA_trials); %SAA "shadowing+multipath" samples for relay segment
        
        for samples=1:SAA_trials
            
            total_F=zeros(L_f(r,1),1);
            phiF=zeros(L_f(r,1),full_seg(r,1));
            
            total_F(:,1)=SAAz(samples,end-1)+SAAz_r(samples,1);
            
            oneF_magsq=0;     
            for i=1:L_f(r,1)
                testidf=segid_f(i,testidgA,r);
                total_F(i,1)=total_F(i,1)+ sum(SAAz(samples,testidf));
                phiF(i,testidgA)=SAAphi(samples,testidf);
            end
            for i=1:L_f(r,1)
                for k=1:L_f(r,1)
                    total_Phi=sum(phiF(i,:)-phiF(k,:)); %all paths start from same segment and end to same segment
                    oneF_magsq=oneF_magsq+ ( exp((chi/2)*(aF+total_F(i,1)))*exp((chi/2)*(aF+total_F(k,1)))*cos(total_Phi) );
                end
            end
            
            %Channel g
            total_G=zeros(L_g(r,1),1);
            phiG=zeros(L_g(r,1),full_seg(r,2));
            
            total_G(:,1)=SAAz_r(samples,2)+SAAz(samples,end);
            
            oneG_magsq=0;         
            for i=1:L_g(r,1)
                testidg=segid_g(i,testidgB,r);
                total_G(i,1)=total_G(i,1)+sum(SAAz(samples,testidg));%+SAA_xi(samples,segid_g(i,1:full_seg(r,2),r)));
                phiG(i,testidgB)=SAAphi(samples,testidg);
            end
            for i=1:L_g(r,1)
                for k=1:L_g(r,1)
                    total_PhiG=sum(phiG(i,:)-phiG(k,:));
                    oneG_magsq=oneG_magsq+ ( exp((chi/2)*(aG+total_G(i,1)))*exp((chi/2)*(aG+total_G(k,1)))*cos(total_PhiG) );
                end
            end
            
            VI_numerator=P_S*P_C*oneF_magsq*oneG_magsq;
            VI_denominator=(P_S*sigmaSQ_D*oneF_magsq) + (P_C*sigmaSQ*oneG_magsq) + (sigmaSQ*sigmaSQ_D) ;
            All_VI(samples,relaypos) = VI_numerator/VI_denominator;
            
        end%end samples
        
    end %end relaypos
    
    Check=mean(All_VI);
    [~, I]=sort(Check,'descend'); %sorted values and respective index
    
    round_pos(r,1)= considered_neigh(1,I(1)); %=find(discrete_points==(mod( R_cord(r,R_orient(r,1)),d_full)));
    R_cord(r,R_orient(r,1))= R_bounds(r,1)-discrete_points(1,1)+discrete_points(1,round_pos(r,1));
    
    
end %end for every relay


end

