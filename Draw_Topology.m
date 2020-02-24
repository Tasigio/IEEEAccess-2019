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
function [all_paths,all_segments,grid_cords]=Draw_Topology(x,y,S_cord,D_cord,R_cord,show_top,R_orient)



% steps_right=abs(D_cord(1,1)-S_cord(1,1))/C.d_full;
% steps_down=abs(D_cord(1,2)-S_cord(1,2))/C.d_full;

[X,Y] = meshgrid(x,y);
grid_cords=[X(:) Y(:)];

if show_top==1
    plot(grid_cords(:,1),grid_cords(:,2),'ro','markersize',15)
    xlabel('meters')
    ylabel('meters')
    xlim([0 600]);
    ylim([0 400]);
    xticks([0:100:600]);
    yticks([0:100:500]);
    
%     ax = gca;
% c = ax.Color;
% ax.GridAlpha = 0.8;
    grid on
   
    hold on
    
    for i=1:length(grid_cords)
        text(grid_cords(i,1),grid_cords(i,2),num2str(i),'HorizontalAlignment','center')
    end
    
    text(S_cord(1,1),S_cord(1,2),'S','HorizontalAlignment','center','color','#D95319','FontWeight','bold','BackgroundColor' ,'#EDB120')
    text(D_cord(1,1),D_cord(1,2),'D','HorizontalAlignment','center','color','#D95319','FontWeight','bold','BackgroundColor' ,'#EDB120')
    for r=1:length(R_cord)
        plot(R_cord(r,1),R_cord(r,2),'b-','Markersize',14)
        if R_orient(r,1)==2
            text(R_cord(r,1)+10,R_cord(r,2),num2str(r),'HorizontalAlignment','center','color','b')
        else
           text(R_cord(r,1),R_cord(r,2)+20,num2str(r),'HorizontalAlignment','center','color','b')
        end
        %     text(R_cord(i,1),R_cord(i,2),['R_' num2str(i)],'HorizontalAlignment','center','color','b')
    end
end

hold on
line([grid_cords(9,1)+13,grid_cords(12,1)-13],[grid_cords(9,2),grid_cords(12,2)],'Color','blue','LineStyle','-','LineWidth',1.5)
hold on
line([grid_cords(7,1)+13,grid_cords(10,1)-13],[grid_cords(7,2),grid_cords(10,2)],'Color','blue','LineStyle','-','LineWidth',1.5)
hold on
line([grid_cords(6,1),grid_cords(5,1)],[grid_cords(6,2)-13,grid_cords(5,2)+13],'Color','blue','LineStyle','-','LineWidth',1.5)
hold on
line([grid_cords(4,1)+13,grid_cords(7,1)-13],[grid_cords(4,2),grid_cords(7,2)],'Color','blue','LineStyle','-','LineWidth',1.5)
hold on
% line([grid_cords(9,1)+13,grid_cords(12,1)-13],[grid_cords(9,2),grid_cords(12,2)],'Color','blue','LineStyle','-','LineWidth',1.5)
% hold on
line([grid_cords(5,1)+13,grid_cords(8,1)-13],[grid_cords(5,2),grid_cords(8,2)],'Color','blue','LineStyle','-','LineWidth',1.5)
hold on
line([grid_cords(12,1),grid_cords(11,1)],[grid_cords(12,2)-13,grid_cords(11,2)+13],'Color','blue','LineStyle','-','LineWidth',1.5)
hold on



%  set(gca,'Layer','top');
% line([S_cord(1,1)+13,grid_cords(6,1)-13],[S_cord(1,2),grid_cords(6,2)],'Color','green','LineStyle',':','LineWidth',1.5)
% hold on
% line([S_cord(1,1)+13,grid_cords(6,1)-13],[S_cord(1,2)+5,grid_cords(6,2)+5],'Color','magenta','LineStyle',':','LineWidth',1.5)
% hold on
% line([grid_cords(6,1),grid_cords(5,1)],[grid_cords(6,2)-13,grid_cords(5,2)+13],'Color','green','LineStyle',':','LineWidth',1.5)
% hold on
% line([grid_cords(7,1)+13,grid_cords(10,1)-13],[grid_cords(7,2)-5,grid_cords(10,2)-5],'Color','magenta','LineStyle',':','LineWidth',1.5)
% hold on
% 
% line([grid_cords(5,1)+13,grid_cords(8,1)-13],[grid_cords(5,2)-5,grid_cords(8,2)-5],'Color','magenta','LineStyle',':','LineWidth',1.5)
% hold on
% line([grid_cords(5,1)+13,grid_cords(8,1)-13],[grid_cords(5,2)-5,grid_cords(8,2)-5],'Color','magenta','LineStyle',':','LineWidth',1.5)
% hold on
% line([grid_cords(5,1)+13,grid_cords(8,1)-13],[grid_cords(5,2)-9,grid_cords(8,2)-9],'Color','green','LineStyle',':','LineWidth',1.5)
% hold on
% line([grid_cords(8,1)+13,grid_cords(11,1)-13],[grid_cords(8,2),grid_cords(11,2)],'Color','green','LineStyle',':','LineWidth',1.5)
% hold on
% line([grid_cords(6,1)+5,grid_cords(5,1)+5],[grid_cords(6,2)-13,grid_cords(5,2)+13],'Color','magenta','LineStyle',':','LineWidth',1.5)
% hold on
% line([grid_cords(8,1),grid_cords(7,1)],[grid_cords(8,2)-13,grid_cords(7,2)+13],'Color','magenta','LineStyle',':','LineWidth',1.5)
% hold on
% line([grid_cords(11,1),grid_cords(10,1)],[grid_cords(11,2)-13,grid_cords(10,2)+13],'Color','green','LineStyle',':','LineWidth',1.5)
% hold on
% line([grid_cords(10,1)+13,D_cord(1,1)-13],[grid_cords(10,2),D_cord(1,2)],'Color','magenta','LineStyle',':','LineWidth',1.5)
% hold on
% line([grid_cords(10,1)+13,D_cord(1,1)-13],[grid_cords(10,2)-5,D_cord(1,2)-5],'Color','green','LineStyle',':','LineWidth',1.5)
% hold on


%black
% line([grid_cords(13,1)-13,D_cord(1,1)+13],[grid_cords(13,2),D_cord(1,2)],'Color','black','LineStyle','-','LineWidth',1)
% hold on
% line([grid_cords(3,1)+13,S_cord(1,1)-13],[grid_cords(3,2),S_cord(1,2)],'Color','black','LineStyle','-','LineWidth',1)
% hold on
% line([grid_cords(3,1),grid_cords(2,1)],[grid_cords(3,2)-13,grid_cords(2,2)+13],'Color','black','LineStyle','-','LineWidth',1)
% hold on
% line([grid_cords(2,1),grid_cords(1,1)],[grid_cords(2,2)-13,grid_cords(1,2)+13],'Color','black','LineStyle','-','LineWidth',1)
% hold on
% line([grid_cords(11,1),grid_cords(10,1)],[grid_cords(11,2)-13,grid_cords(10,2)+13],'Color','black','LineStyle','-','LineWidth',1)
% hold on
% line([grid_cords(14,1),grid_cords(13,1)],[grid_cords(14,2)-13,grid_cords(13,2)+13],'Color','black','LineStyle','-','LineWidth',1)
% hold on
% line([grid_cords(15,1),grid_cords(14,1)],[grid_cords(15,2)-13,grid_cords(14,2)+13],'Color','black','LineStyle','-','LineWidth',1)
% hold on
% line([grid_cords(1,1)+13,grid_cords(4,1)-13],[grid_cords(1,2),grid_cords(4,2)],'Color','black','LineStyle','-','LineWidth',1)
% hold on
% line([grid_cords(2,1)+13,grid_cords(5,1)-13],[grid_cords(2,2),grid_cords(5,2)],'Color','black','LineStyle','-','LineWidth',1)
% hold on
% line([grid_cords(12,1)+13,grid_cords(15,1)-13],[grid_cords(12,2),grid_cords(15,2)],'Color','black','LineStyle','-','LineWidth',1)
% hold on
% line([grid_cords(11,1)+13,grid_cords(14,1)-13],[grid_cords(11,2),grid_cords(14,2)],'Color','black','LineStyle','-','LineWidth',1)
% hold on
% line([grid_cords(8,1)+13,grid_cords(11,1)-13],[grid_cords(8,2),grid_cords(11,2)],'Color','black','LineStyle','-','LineWidth',1)
% hold on
%%
%Topology
% from_node=[11,11,17,17,10,10,16,16,9,9,23,23,22,22,15,15,8,8,21,21,14,14,20,20,7,13,19,29,28,27,26];
% to_node=[17,10,16,23,16,9,22,15,15,8,22,29,28,21,21,14,14,7,27,20,20,13,26,19,13,19,25,28,27,26,25];
% from_node=[11,11,17,17,10,10,16,16,9,23,22,15];
% to_node  =[17,10,16,23,16,9,22,15,15,22,21,21];
from_node=[6,5,6,5,4,9,8,9, 8,7,12,11];
to_node  =[5,4,9,8,7,8,7,12,11,10,11,10];
all_segments=[from_node',to_node'];
T=digraph(all_segments(:,1),all_segments(:,2));
[A] = getpaths(T);
%all_paths=cell2mat(A(1:35));
all_paths=cell2mat(A(1:6)); %all paths from source to destination
seg_orient=[2;2;1;1;1;2;2;1;1;1;2;2];
%seg_orient=[2;1;1;2;2;1;2;1;2;1;1;2;2;1;2;1;2;1;2;1;2;1;2;1;2;2;2;1;1;1;1];

%ADD SEGMENT NUMBERS
% hold on
% for i=1:length(all_segments)
%     if seg_orient(i)==2
%         text(  (( grid_cords(all_segments(i,1),1)+grid_cords(all_segments(i,2),1))/2), grid_cords(all_segments(i,1),2) ,['(' num2str(i) ')'],'HorizontalAlignment','center','color','b')
%     elseif seg_orient(i)==1
%         text( grid_cords(all_segments(i,1),1),(( grid_cords(all_segments(i,1),2)+grid_cords(all_segments(i,2),2))/2) ,['(' num2str(i) ')'],'HorizontalAlignment','center','color','b')
%     end
% end
end