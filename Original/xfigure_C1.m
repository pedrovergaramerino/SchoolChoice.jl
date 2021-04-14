% --------------------------------------------------------------------------------- %
% FIGURE C1 - Monte Carlo Simulations: Spatial Distribution of Students and Schools %
% --------------------------------------------------------------------------------- %

clear

if exist('set_path.m', 'file') == 2
    set_path
else 
    error('ERROR: specify path to folder containing replication files in change_directory.m and execute script')  
end

% Parameters
short_sim = 1;
propor = 5; % 500 students
M = 500; % 500 MC samples
J = 6; % 6 schools
A = 4; % Maximum size of ROL = 4
unit_cost = 0; % No cost of ranking more schools

MC_data_simulation;

% Print figure

% Map of simulated school district
h = figure('Name','Geographic Coordinates','NumberTitle','off');
set(h,'Units','Inches');
axis equal;
axis([-1 1 -1 1]);
FIG.MAP.legend_names_g=cell(1,2);
FIG.MAP.legend_names_g{1} = 'Students';
FIG.MAP.legend_names_g{2} = 'Schools';
hold on;
set(gca,'fontsize',6)
plot(MC2.check_x,MC2.check_y,'.','Color','blue','MarkerSize',8);
plot(MC.school_x,MC.school_y,'o','MarkerFaceColor','red','MarkerSize',10,'MarkerEdgeColor','black');  
plot(MC.circle_x,MC.circle_y,'Color','black','LineWidth',1);
FIG.MAP.g_legend = legend(FIG.MAP.legend_names_g,'Location','EastOutside');
set(FIG.MAP.g_legend,'FontSize',14); 
hold off;

set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[13, 8.5])
print(h, '-dpdf', '-r0',fullfile(dir_fig,'xfigure_C1.pdf'));

close(h);

