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

close all;
clear;
clc;

%%
setup_parameters;
a_LoS=channel_parameters.a_LoS;
a_NLoS=channel_parameters.a_NLoS;
sigmaSQ=channel_parameters.sigmaSQ;
sigmaSQ_D=channel_parameters.sigmaSQ_D;

rng(7);
%Initial antenna coordinate(s) for all relays
%Rcord_initial=[400,357;200,257;357,200;357,500;257,100;300,357];
%R_segments=[22,21;9,8;14,20;17,23;7,13;16,15]; %the intersections between which the relays are located
% Rcord_initial=[200,257; 100,257; 257,400];
% R_segments=[16,15;10,9;17,23];
% R_orient=[2;2;1];
% Rcord_initial=[350,300; 350,100; 250,200; 400,250; ];

% %6 relay topology
%  Rcord_initial=[357,300; 357,100; 257,200; 400,257; 200,257;257,100];
% C_segments=[9,12;7,10;5,8;12,11;6,5;4,7];
% C_orient=[1;1;1;2;2;1];

%Other topology
%  Rcord_initial=[257,300; 200,257; 357,100; 400,157];
%  R_segments=[6,9;6,5;7,10;11,10];
%  R_orient=[1;2;1;2];

%Paper topology
Rcord_initial=[357,300; 357,100; 257,200; 400,257; ];
C_segments=[9,12;7,10;5,8;12,11]; %Cluster intersection end points
C_orient=[1;1;1;2];
middlepos=29;

%make sure relays are on the discrete positions
% R_orient=[1;2;2;1];
% R_orient=[2;2;1;1;1;2]; %indicate coordinate index relay can move upon
N_C=size(Rcord_initial,1); % Number of clusters

show_topology=1; % Show topology figure
%Draw grid, return intersections and all paths from S to D
[all_paths,all_segments,I_cord]=Draw_Topology(x,y,S_cord,D_cord,Rcord_initial,show_topology,C_orient);

C_bounds=zeros(N_C,2); % Coordinate bounds of segment that cluster is placed
for r=1:N_C
    if C_orient(r,1)==2
        C_bounds(r,1) = I_cord(C_segments(r,2),C_orient(r,1))+discrete_points(1,1);
        C_bounds(r,2) = I_cord(C_segments(r,1),C_orient(r,1))-discrete_points(1,1);
    else
        C_bounds(r,1) = I_cord(C_segments(r,1),C_orient(r,1))+discrete_points(1,1);
        C_bounds(r,2) = I_cord(C_segments(r,2),C_orient(r,1))-discrete_points(1,1);
    end
end

%%
full_seg=zeros(N_C,2); % Stores number of full segments [from p_S to p_r, from p_r to p_D] (regardless of path)

L_f=zeros(N_C,1); % Stores the number of paths from p_S to every relay p_r
L_g=zeros(N_C,1); % Stores the number of paths from every relay p_r to p_D

relay_pathsF=zeros(6,5,N_C); % Stores every path from p_S to p_r
relay_pathsG=zeros(6,5,N_C); % Stores every path from p_r to p_D (10,7,N_c)

for r=1:N_C %for each relay p_r
    
    steps_right=floor(abs(Rcord_initial(r,1)-I_cord(I_S,1))/d_full);
    steps_down=floor(abs(Rcord_initial(r,2)-I_cord(I_S,2))/d_full);
    
    full_seg(r,1)=steps_right+steps_down;
    L_f(r,1)=nchoosek(steps_right+steps_down,steps_down); %number of possible paths from p_S to p_r
    
    A=unique( all_paths(:,[1:full_seg(r,1)+1]), 'rows');
    relay_pathsF(1:L_f(r,1),1:full_seg(r,1)+1,r)=A(A(:,full_seg(r,1)+1)==C_segments(r,1),:); %traversed nodes till before relay
    
    
    for i=1:L_f(r,1)
        if full_seg(r,1)==0 %if there is no full segment
            segid_f(i,1,r)=0;
        else
            for seg=1:full_seg(r,1)
                [ind_seg, ~]=find(ismember(all_segments,[relay_pathsF(i,seg,r),relay_pathsF(i,seg+1,r)],'rows'));
                segid_f(i,seg,r)=ind_seg;
            end
            
        end
    end
    
    steps_right=floor(abs(Rcord_initial(r,1)-I_cord(I_D,1))/d_full);
    steps_down=floor(abs(Rcord_initial(r,2)-I_cord(I_D,2))/d_full);
    
    full_seg(r,2)=steps_right+steps_down;
    L_g(r,1)=nchoosek(steps_right+steps_down,steps_down);
    
    A=unique( all_paths(:,[full_seg(r,1)+2:end]), 'rows');
    relay_pathsG(1:L_g(r,1),1:full_seg(r,2)+1,r)=A(A(:,1)==C_segments(r,2),:);
    
    for i=1:L_g(r,1)
        if full_seg(r,2)==0
            segid_g(i,1,r)=0;
        else
            for seg=1:full_seg(r,2)
                [ind_seg, ~]=find(ismember(all_segments,[relay_pathsG(i,seg,r),relay_pathsG(i,seg+1,r)],'rows'));
                segid_g(i,seg,r)=ind_seg;
            end
        end
    end
end

%%%
% pathlossG=zeros(N_C,delta)
% dLoS_f=norm(S_cord-I_cord(I_S,:)); %LoS distance
% for r=1:N_C
%     R_cord_temp=[C_bounds(r,1)-1, Rcord_initial(r,2)];
%     for relayPos=1:delta %every possible discrete position in segment
%         R_cord_temp(1,C_orient(r,1))=C_bounds(r,1)-discrete_points(1,1)+ discrete_points(1,relayPos); %temporary examined antenna coordiante
%         dNLoS_f=norm(I_cord(I_S,:)-R_cord_temp(1,:),1); %dNLoS: Manhattan distance from first intersection after p_S to p_r
%         pathlossF(r, relayPos)=-(a_LoS*10*log10(dLoS_f)) - (a_NLoS*full_seg(r,1)*10*log10(d_full)) - (a_NLoS*10*log10(mod(dNLoS_f,d_full))) - (Delta*(full_seg(r,1)+1));
%      
%     end
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Sigma,Sigma_R] = Covariance_Matrices(T,discrete_points,d_full,delta,channel_parameters);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chi=log(10)/(10); %constant used in calculations

% Stores optimal SINR value of all trials at every time slots
Vopt_IDEAL=zeros(total_trials,T);
Vopt_SAA=zeros(total_trials,T);
Vopt_RANDOMIZED=zeros(total_trials,T);

Vopt_SAAconst=zeros(total_trials,T);
Vopt_RANDOMIZEDconst=zeros(total_trials,T);

% Stores relay positions for all clusters
AllTrials_IDEALpos=zeros(total_trials,T,N_C);
AllTrials_SAApos=zeros(total_trials,T,N_C);
AllTrials_RANDOMIZEDpos=zeros(total_trials,T,N_C);

AllTrials_SAAconstpos=zeros(total_trials,T,N_C);
AllTrials_RANDOMIZEDconstpos=zeros(total_trials,T,N_C);

considered_neigh=[1:4,delta-4:delta]; % considered constrained neighborhood

tic
% WaitMessage = parfor_wait(C.total_trials, 'Waitbar', true);
parfor trial=1:total_trials % for each trial
    
    trial/total_trials
    %    WaitMessage.Send;
    z=mvnrnd(zeros(T,1),Sigma,size(all_segments,1)+2); %shadowing+multipath terms for every non-cluster segment, source, destination
    z_r=mvnrnd(zeros(T*2*delta,1), Sigma_R, N_C); %shadowing+multipath terms for the segments with cluster
    
    %     multipath=sqrt(C.sigmaSQ_xi).*randn(size(all_segments,1)+2,C.T); %multipath for every segment, source, destination
    %     multipath_r=sqrt(C.sigmaSQ_xi).*randn(N_c, C.T*2*C.delta);
    
    phase=rand(size(all_segments,1)+2, T); % phase terms for all segments, source, destination
    phase_r=rand(N_C, T*2*delta); % phase terms for segment with cluster
    
    %SAA samples
    SAAphi=rand(SAA_trials,size(all_segments,1));  %phase for SAA
    
    %Current relay coordinates for all clusters, depending on policy followed
    Rcord_ideal=Rcord_initial;
    Rcord_SAA=Rcord_initial;
    Rcord_randomized=Rcord_initial;
    Rcord_SAAconst=Rcord_initial;  %constrained
    Rcord_randomizedconst=Rcord_initial;
    
    %Relay positions of current trial for all clusters, depending on policy followed
    CurrTrial_IDEALpos=zeros(N_C,T);
    CurrTrial_SAApos=zeros(N_C,T);
    CurrTrial_RANDOMIZEDpos=zeros(N_C,T);
    CurrTrial_SAAconstpos=zeros(N_C,T); %constrained
    CurrTrial_RANDOMIZEDconstpos=zeros(N_C,T);
    
    %Initial relay positions for all clusters
    CurrTrial_IDEALpos(:,1)=middlepos*ones(N_C,1);
    CurrTrial_SAApos(:,1)=middlepos*ones(N_C,1);
    CurrTrial_RANDOMIZEDpos(:,1)=middlepos*ones(N_C,1);
    CurrTrial_SAAconstpos(:,1)=middlepos*ones(N_C,1);   %constrained
    CurrTrial_RANDOMIZEDconstpos(:,1)=middlepos*ones(N_C,1);
    
    % Vopt for all time slots of current trial (needed for parfor)
    Vopt_IDEAL_temp=zeros(1,T);
    Vopt_SAA_temp=zeros(1,T);
    Vopt_RANDOMIZED_temp=zeros(1,T);
    Vopt_SAAconst_temp=zeros(1,T);      %constrained
    Vopt_RANDOMIZEDconst_temp=zeros(1,T);
    
    for t=1:T-1 %for every time slot
        
        % Generate SAA "shadowing+multipath" samples for segments WITHOUT cluster, source, destination
        m_bld= (z(1:end,1:t))'; %measurements of segments without cluster
        [SAAz] = Generate_SAAsamples(SAA_trials,t,Sigma,m_bld,ones(1,size(all_segments,1)+2));
        
        
        %% IDEAL
        % Beamforming and SINR evaluation at destination (2nd stage problem of time slot t)
        [VT_ideal]=SINR_destination(t,CurrTrial_IDEALpos(:,t),Rcord_ideal,z,z_r,phase,phase_r,L_f,L_g,N_C,I_cord,full_seg,segid_f,segid_g,C_segments,channel_parameters,S_cord,D_cord,P_S,P_C,Delta,delta,I_S,d_full);
        Vopt_IDEAL_temp(1,t)=VT_ideal;
        
        % Relay selection (1st stage problem of time slot t+1)
        [Rcord_ideal,round_posIDEAL]= Ideal(t+1,Rcord_ideal,N_C,C_orient,C_bounds,z,z_r,phase,phase_r,L_f,L_g,I_cord,full_seg,segid_f,segid_g,C_segments,channel_parameters,S_cord,D_cord,P_S,P_C,Delta,delta,discrete_points,I_S,d_full); 
        
        %% SAA Unconstrained neighborhood
        % Beamforming and SINR evaluation at destination (2nd stage problem of time slot t)
        [VT_SAA]=SINR_destination(t,CurrTrial_SAApos(:,t),Rcord_SAA,z,z_r,phase,phase_r,L_f,L_g,N_C,I_cord,full_seg,segid_f,segid_g,C_segments,channel_parameters,S_cord,D_cord,P_S,P_C,Delta,delta,I_S,d_full);
        Vopt_SAA_temp(1,t)=VT_SAA;
        
        % Relay selection (1st stage problem of time slot t+1)
        [Rcord_SAA,round_posSAA]= Sample_Avg_Approx(t,SAA_trials,Rcord_SAA,N_C,C_orient,C_bounds,L_f,L_g,I_cord,full_seg,segid_f,segid_g,C_segments,SAAz,SAAphi,z_r,Sigma_R,chi,CurrTrial_SAApos,channel_parameters,S_cord,D_cord,P_S,P_C,Delta,delta,discrete_points,I_S,d_full); 
        
        %% RANDOM Unconstrained neighborhood
        % Beamforming and SINR evaluation at destination (2nd stage problem of time slot t)
        [VT_randomized]=SINR_destination(t,CurrTrial_RANDOMIZEDpos(:,t),Rcord_randomized,z,z_r,phase,phase_r,L_f,L_g,N_C,I_cord,full_seg,segid_f,segid_g,C_segments,channel_parameters,S_cord,D_cord,P_S,P_C,Delta,delta,I_S,d_full);
        Vopt_RANDOMIZED_temp(1,t)=VT_randomized;
        
        % Relay selection (1st stage problem of time slot t+1)
        [Rcord_randomized,round_posRANDOMIZED]= Randomized(discrete_points,delta,Rcord_randomized,N_C,C_orient,C_bounds); % Randomized
        
        %% SAA Constrained neighborhood
        % Beamforming and SINR evaluation at destination (2nd stage problem of time slot t)
        [VT_SAAconst]=SINR_destination(t,CurrTrial_SAAconstpos(:,t),Rcord_SAAconst,z,z_r,phase,phase_r,L_f,L_g,N_C,I_cord,full_seg,segid_f,segid_g,C_segments,channel_parameters,S_cord,D_cord,P_S,P_C,Delta,delta,I_S,d_full);
        Vopt_SAAconst_temp(1,t)=VT_SAAconst;
        
        % Relay selection (1st stage problem of time slot t+1)
        [Rcord_SAAconst,round_posSAAconst] = SAA_constrained(t,SAA_trials,Rcord_SAAconst,N_C,C_orient,C_bounds,L_f,L_g,I_cord,full_seg,segid_f,segid_g,C_segments,SAAz,SAAphi,z_r,Sigma_R,chi,CurrTrial_SAAconstpos,channel_parameters,S_cord,D_cord,P_S,P_C,Delta,delta,discrete_points,I_S,d_full,considered_neigh);
        
        %% RANDOM Constrained neighborhood
        % Beamforming and SINR evaluation at destination (2nd stage problem of time slot t)
        [VT_randomizedconst]=SINR_destination(t,CurrTrial_RANDOMIZEDconstpos(:,t),Rcord_randomizedconst,z,z_r,phase,phase_r,L_f,L_g,N_C,I_cord,full_seg,segid_f,segid_g,C_segments,channel_parameters,S_cord,D_cord,P_S,P_C,Delta,delta,I_S,d_full);
        Vopt_RANDOMIZEDconst_temp(1,t)=VT_randomizedconst;
        
        % Relay selection (1st stage problem of time slot t+1)
        [Rcord_randomizedconst,round_posRANDOMIZEDconst] = Randomized_constrained(discrete_points,Rcord_randomizedconst,N_C,C_orient,C_bounds,considered_neigh);
        
        %%
        
        %Current antenna positions in segment for all relays
        CurrTrial_IDEALpos(:,t+1)=round_posIDEAL;
        CurrTrial_SAApos(:,t+1)=round_posSAA;
        CurrTrial_RANDOMIZEDpos(:,t+1)=round_posRANDOMIZED;
        
        CurrTrial_SAAconstpos(:,t+1)=round_posSAAconst;
        CurrTrial_RANDOMIZEDconstpos(:,t+1)=round_posRANDOMIZEDconst;
        
    end %end for every time slot
    
    Vopt_IDEAL(trial,:)= Vopt_IDEAL_temp; %optimal SINR values of current trial for all time slots
    Vopt_SAA(trial,:)= Vopt_SAA_temp;
    Vopt_RANDOMIZED(trial,:)= Vopt_RANDOMIZED_temp;
    
    Vopt_SAAconst(trial,:)= Vopt_SAAconst_temp;
    Vopt_RANDOMIZEDconst(trial,:)= Vopt_RANDOMIZEDconst_temp;
    
    
    %Keep antenna positions
    for r=1:N_C
        AllTrials_IDEALpos(trial,:,r)=CurrTrial_IDEALpos(r,:);  %antenna positions of current trial for all time slots
        AllTrials_SAApos(trial,:,r)=CurrTrial_SAApos(r,:);
        AllTrials_RANDOMIZEDpos(trial,:,r)=CurrTrial_RANDOMIZEDpos(r,:);
        
        AllTrials_SAAconstpos(trial,:,r)=CurrTrial_SAAconstpos(r,:);
        AllTrials_RANDOMIZEDconstpos(trial,:,r)=CurrTrial_RANDOMIZEDconstpos(r,:);
    end
    
end %end trials

% WaitMessage.Destroy
endtime=toc;

minutes = endtime/60;
hours = floor(minutes/60);
minutes=mod(minutes, 60);

fprintf('The simulation took %3f hours and %f minutes\n', hours, minutes)


save test.mat

Draw_Plots(T,delta,total_trials,N_C,Vopt_IDEAL,Vopt_SAA,Vopt_RANDOMIZED,Vopt_SAAconst,Vopt_RANDOMIZEDconst,AllTrials_IDEALpos,AllTrials_SAApos,AllTrials_RANDOMIZEDpos,AllTrials_SAAconstpos,AllTrials_RANDOMIZEDconstpos)

