clear,clc
%Michael Gabalis
%% HW 3: Principal Component Analysis (PCA)

% NOTE: convert uint8 to double using double() before processing!
% Each frame of video should only produce a single timepoint
%% Part 1: Ideal Case
% load data
load('cam1_1.mat')
load('cam2_1.mat')
load('cam3_1.mat')

% Camera 1 
plots = [0 0 0 1 0 0]; % which plots to show
close all; clc; 
video = vidFrames1_1;
xrange = [300,400];
yrange = [200,400];
var = 1;
max_pval = 240;
[x1_1, y1_1] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);
pause(1)

% Camera 2
close all; clc; 

video = vidFrames2_1;
xrange = [250, 350];
yrange = [100, 350];
var = 1;
max_pval = 240;
[x2_1, y2_1] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);

% Camera 3 
close all; clc; 

video = vidFrames3_1;
xrange = [250, 500];
yrange = [250, 350];
var = 1;
max_pval = 240;
[x3_1, y3_1] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);

%% Principal Component Analysis part 1
close all; clc;
rank_approx = 1;
% offset accounts by aligning the frames
offset = [11, 20, 11]; 
offset = offset - (min(offset) - 1);
pc = 2; % number of principal components to plot
yrange = [100, 600];
A = my_pca(rank_approx, pc, offset, yrange, x1_1, y1_1, x2_1, y2_1, x3_1, y3_1);
%% Part 2: Nosiy Case

close all;clc
load('cam1_2.mat')
load('cam2_2.mat')
load('cam3_2.mat')

% Camera 1 
close all

video = vidFrames1_2;
xrange = [300, 400];
yrange = [250, 400];
var = 0.5;
max_pval = 230;
plots = [0 0 0 0 0 0];
[x1_2, y1_2] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);

% Camera 2
close all 

video = vidFrames2_2;
xrange = [175, 450];
yrange = [50, 450];
var = 0.5;
max_pval = 240;
plots = [0 0 0 0 0 0];
[x2_2, y2_2] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);
% Camera 3
close all

video = vidFrames3_2;
xrange = [250, 500];
yrange = [225, 300];
var = 0.8;
max_pval = 230;
plots = [0 0 0 0 0 0];
[x3_2, y3_2] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);

%% Principal Component Analysis part 2
close all

rank_approx = 3;
offset = [15, 1, 17];
offset = offset - (min(offset) - 1);
pc = 4;
yrange = [100, 600];
A = my_pca(rank_approx, pc, offset, yrange, x1_2, y1_2, x2_2, y2_2, x3_2, y3_2);

%% Part 3: Horizontal Displacement

close all
load('cam1_3.mat')
load('cam2_3.mat')
load('cam3_3.mat')

% Camera 1 
close all

video = vidFrames1_3;
xrange = [250, 400];
yrange = [200, 400];
var = 1;
max_pval = 250;
plots = [0 0 0 0 0 0];
[x1_3, y1_3] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);

% Camera 2
 close all 

video = vidFrames2_3;
xrange = [200, 400];
yrange = [175, 400];
var = 1;
max_pval = 240;
plots = [0 0 0 0 0 0];
[x2_3, y2_3] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);

% Camera 3
close all

video = vidFrames3_3;
xrange = [250, 450];
yrange = [175, 325];
var = 1;
max_pval = 245;
plots = [0 0 0 0 0 0];
[x3_3, y3_3] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);

%% Principal Component Analysis part 3
close all; clc;

rank_approx = 2;
% frame bucket moves
offset = [18 44 9];
offset = offset - (min(offset) - 1);
pc = 3;
yrange = [100,600];
my_pca(rank_approx, pc, offset, yrange, x1_3, y1_3, x2_3, y2_3, x3_3, y3_3);

%% Part 4: Horizontal Displacement and Rotation
close all; clc;

load('cam1_4.mat')
load('cam2_4.mat')
load('cam3_4.mat')

% Camera 1

close all; clc; 
video = vidFrames1_4;
xrange = [300, 450];
yrange = [225, 400];
var = 1;
max_pval = 245;
plots = [0 0 0 0 0 0];
[x1_4, y1_4] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);

% Camera 2
close all;

video = vidFrames2_4;
xrange = [210, 400];
yrange = [100, 350];
var = 1;
max_pval = 250;
plots = [0 0 0 0 0 0];
[x2_4, y2_4] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);

% Camera 3 
close all

video = vidFrames3_4;
xrange = [300, 500];
yrange = [175, 250];
var = 0.7;
max_pval = 235;
plots = [0 0 0 0 0 0];
[x3_4, y3_4] = get_xy_coords(video, xrange, yrange, var, max_pval, plots);

%% Principal Component Analysis part 4
close all

rank_approx = 2;
offset = [11, 17, 9];
offset = offset - (min(offset) - 1);
pc = 3;
yrange = [100,600];
my_pca(rank_approx, pc, offset, yrange, x1_4, y1_4, x2_4, y2_4, x3_4, y3_4);