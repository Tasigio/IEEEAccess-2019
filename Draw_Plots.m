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
function []= Draw_Plots(T,delta,total_trials,N_R,Vopt_IDEAL,Vopt_SAA,Vopt_RANDOMIZED,Vopt_SAAconst,Vopt_RANDOMIZEDconst,AllTrials_IDEALpos,AllTrials_SAApos,AllTrials_RANDOMIZEDpos,AllTrials_SAAconstpos,AllTrials_RANDOMIZEDconstpos)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

close all;
figuresdir = 'Results\';

savePDF=0;
savePNG=0;

CHECK_IDEAL=mean(Vopt_IDEAL);
plot(10*log10(CHECK_IDEAL),'-','LineWidth',1.5)
hold on
CHECK_SAA=mean(Vopt_SAA);
plot(10*log10(CHECK_SAA),'-','LineWidth',1.5)
hold on
CHECK_SAAconst=mean(Vopt_SAAconst);
plot(10*log10(CHECK_SAAconst),'-','LineWidth',1.5)
hold on
CHECK_RANDOMIZED=mean(Vopt_RANDOMIZED);
plot(10*log10(CHECK_RANDOMIZED),':','LineWidth',1.5)
hold on
CHECK_RANDOMIZEDconst=mean(Vopt_RANDOMIZEDconst);
plot(10*log10(CHECK_RANDOMIZEDconst),':','LineWidth',1.5)

xlabel('Time Slot','FontWeight','bold','Interpreter','latex')
ylabel('SINR(dB)','FontWeight','bold','Interpreter','latex')

%title([num2str(total_trials), ' trials,' num2str(SAA_trials), ' SAA samples'])
legend('Ideal','SAA','SAA Constrained','Randomized','Randomized Constrained','Location','southeast','FontWeight','bold','Interpreter','latex')

COL=get(gca,'colororder'); %default color order
p1=yline(mean(10*log10(CHECK_IDEAL(1:49))),'-','Color', COL(1,:),'LineWidth',1.5);
set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p2=yline(mean(10*log10(CHECK_SAA(1:49))),'-','Color', COL(2,:),'LineWidth',1.5);
set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p3=yline(mean(10*log10(CHECK_SAAconst(1:49))),'-','Color', COL(3,:),'LineWidth',1.5);
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p4=yline(mean(10*log10(CHECK_RANDOMIZED(1:49))),':','Color', COL(4,:),'LineWidth',1.5);
set(get(get(p4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p5=yline(mean(10*log10(CHECK_RANDOMIZEDconst(1:49))),':','Color', COL(5,:),'LineWidth',1.5);
set(get(get(p5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(gca,'Ylim',[-145 -95]);
set(gca,'Xlim',[0 T]);
%ylim([min(10*log10(CHECK_RANDOMIZED)) -100])
grid on;

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

% saveas(gcf,strcat(figuresdir, 'Allpolicies'), 'fig');
if savePNG==1
    saveas(gcf,strcat(figuresdir, 'Allpolicies'), 'png');
end
if savePDF==1
    print('-painters','-dpdf','Allpolicies')
end


%%
%Histograms
for r=1:N_R
    Heatmap_IDEAL=zeros(delta,T);
    Heatmap_SAA=zeros(delta,T);
    Heatmap_RANDOMIZED=zeros(delta,T);
    Heatmap_SAAconst=zeros(delta,T);
    Heatmap_RANDOMIZEDconst=zeros(delta,T);
    
    for i=1:total_trials
        for t=1:T
            Heatmap_SAA(AllTrials_SAApos(i,t,r),t)= Heatmap_SAA(AllTrials_SAApos(i,t,r),t)+1;
            Heatmap_IDEAL(AllTrials_IDEALpos(i,t,r),t)= Heatmap_IDEAL(AllTrials_IDEALpos(i,t,r),t)+1;
            Heatmap_RANDOMIZED(AllTrials_RANDOMIZEDpos(i,t,r),t)= Heatmap_RANDOMIZED(AllTrials_RANDOMIZEDpos(i,t,r),t)+1;
            Heatmap_SAAconst(AllTrials_SAAconstpos(i,t,r),t)= Heatmap_SAAconst(AllTrials_SAAconstpos(i,t,r),t)+1;
            Heatmap_RANDOMIZEDconst(AllTrials_RANDOMIZEDconstpos(i,t,r),t)= Heatmap_RANDOMIZEDconst(AllTrials_RANDOMIZEDconstpos(i,t,r),t)+1;
        end
    end
    
    Heatmap_IDEAL=(Heatmap_IDEAL/total_trials);
    Heatmap_SAA=(Heatmap_SAA/total_trials);
    Heatmap_RANDOMIZED=(Heatmap_RANDOMIZED/total_trials);
    Heatmap_SAAconst=(Heatmap_SAAconst/total_trials);
    Heatmap_RANDOMIZEDconst=(Heatmap_RANDOMIZEDconst/total_trials);
    
    figure
    b1=bar3(Heatmap_IDEAL(:,2:end));
    %     colorbar('northoutside')
    %colorbar
    for k = 1:length(b1)
        zdata = b1(k).ZData;
        b1(k).CData = zdata;
        b1(k).FaceColor = 'interp';
    end
    caxis([0 1]);
    title(['Cluster ' num2str(r) ' (Ideal)'],'Interpreter','latex','FontSize',14)
    set(gca,'Xlim',[2 T]);
    set(gca,'Ylim',[1 delta]);
    set(gca,'Zlim',[0 1]);
    xlabel('Time Slot','Interpreter','latex','FontSize',14);
    ylabel('Relay Location','Interpreter','latex','FontSize',14);
    view([-36 14])
    saveas(gcf,strcat(figuresdir, ['Ideal_R' num2str(r)]), 'fig');
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    % set(gcf,'Units','inches');
    % screenposition = get(gcf,'Position');
    % set(gcf,...
    %     'PaperPosition',[0 0 screenposition(3:4)],...
    %     'PaperSize',[screenposition(3:4)]);
    
    if savePNG==1
        saveas(gcf,strcat(figuresdir,['Ideal_R' num2str(r)]), 'png');
    end
    if savePDF==1
        saveas(gcf,strcat(figuresdir, ['Ideal_R' num2str(r)]), 'pdf');
    end
    
    
    figure
    b2=bar3(Heatmap_SAA(:,2:end));
    for k = 1:length(b2)
        zdata = b2(k).ZData;
        b2(k).CData = zdata;
        b2(k).FaceColor = 'interp';
    end
    caxis([0 1])
    title(['Cluster ' num2str(r) ' (SAA)'],'Interpreter','latex','FontSize',14);
    set(gca,'Xlim',[2 T]);
    set(gca,'Ylim',[1 delta]);
    set(gca,'Zlim',[0 1]);
    xlabel('Time Slot','Interpreter','latex','FontSize',14);
    ylabel('Relay Location','Interpreter','latex','FontSize',14);
    view([-36 14])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    saveas(gcf,strcat(figuresdir, ['SAA_R' num2str(r)]), 'fig');
    if savePNG==1
        saveas(gcf,strcat(figuresdir, ['SAA_R' num2str(r)]), 'png');
    end
    if savePDF==1
        saveas(gcf,strcat(figuresdir, ['SAA_R' num2str(r)]), 'pdf');
    end
    
        figure
        b3=bar3(Heatmap_RANDOMIZED(:,2:end));
        %     colorbar
        for k = 1:length(b3)
            zdata = b3(k).ZData;
            b3(k).CData = zdata;
            b3(k).FaceColor = 'interp';
        end
        caxis([0 1])
        title(['Random for relay ' num2str(r)],'Interpreter','latex','FontSize',14);
        set(gca,'Xlim',[2 T]);
        set(gca,'Ylim',[1 delta]);
        set(gca,'Zlim',[0 1]);
        xlabel('Time slot','Interpreter','latex','FontSize',14);
        ylabel('Relay location','Interpreter','latex','FontSize',14);
        view([-36 14])
        saveas(gcf,strcat(figuresdir, ['Random_R' num2str(r)]), 'fig');
        if savePNG==1
            saveas(gcf,strcat(figuresdir,['Random_R' num2str(r)]), 'png');
        end
        if savePDF==1
            saveas(gcf,strcat(figuresdir,['Random_R' num2str(r)]), 'pdf');
        end
    
        figure
        b4=bar3(Heatmap_SAAconst(:,2:end));
        %colorbar
        for k = 1:length(b4)
            zdata = b4(k).ZData;
            b4(k).CData = zdata;
            b4(k).FaceColor = 'interp';
        end
        caxis([0 1])
        title(['SAA constrained for relay ' num2str(r)],'Interpreter','latex','FontSize',14);
        set(gca,'Xlim',[2 T]);
        set(gca,'Ylim',[1 delta]);
        set(gca,'Zlim',[0 1]);
        xlabel('Time Slot','Interpreter','latex','FontSize',14);
        ylabel('Relay Location','Interpreter','latex','FontSize',14);
        view([-36 14])
        saveas(gcf,strcat(figuresdir, ['SAACons_R' num2str(r)]), 'fig');
        if savePNG==1
            saveas(gcf,strcat(figuresdir, ['SAACons_R' num2str(r)]), 'png');
        end
        if savePDF==1
            saveas(gcf,strcat(figuresdir, ['SAACons_R' num2str(r)]), 'pdf');
        end
    
    
        figure
        b5=bar3(Heatmap_RANDOMIZEDconst(:,2:end));
        %     colorbar
        for k = 1:length(b5)
            zdata = b5(k).ZData;
            b5(k).CData = zdata;
            b5(k).FaceColor = 'interp';
        end
        caxis([0 1])
        title(['Random constrained for relay ' num2str(r)],'Interpreter','latex','FontSize',14);
        set(gca,'Xlim',[2 T]);
        set(gca,'Ylim',[1 delta]);
        set(gca,'Zlim',[0 1]);
        xlabel('Time slot','Interpreter','latex','FontSize',14);
        ylabel('Relay location','Interpreter','latex','FontSize',14);
        view([-36 14])
        saveas(gcf,strcat(figuresdir, ['RandomCons_R' num2str(r)]), 'fig');
        if savePNG==1
            saveas(gcf,strcat(figuresdir, ['RandomCons_R' num2str(r)]), 'png');
        end
        if savePDF==1
            saveas(gcf,strcat(figuresdir,['RandomCons_R' num2str(r)]), 'pdf');
        end
    
    
    
    
end

end

