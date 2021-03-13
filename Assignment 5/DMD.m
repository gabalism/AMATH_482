%%Michael Gabalis Code
%%Dynamic Mode Decomposition

clear,clc
%% Loading the video clips and processing

%choose which video to read

video = VideoReader('monte_carlo.mov');
%video = VideoReader('ski_drop.mov');

skip_frame = 2; % to reduce size of data
res_Red = 0.5; % reduce resolution by this factor.
imheight = video.Height*res_Red;
imwidth = video.Width*res_Red;
num_frames = ceil(video.Duration * video.FrameRate / skip_frame);
dt = 1;
t = 1:num_frames;
v = zeros(num_frames, imheight, imwidth); % stores frames as rows of 2D matrices
X = zeros(imheight*imwidth, num_frames); % stores frames as 1D columns

frame = 1;
index = 1;
while hasFrame(video)
   nframe = readFrame(video);
   if mod(index, skip_frame) == 0
        x = imresize(nframe, res_Red);
        v(frame, :, :) = imcomplement(rgb2gray(x)); % imcomplement
        X(:, frame) = reshape(v(frame,:,:), [imheight*imwidth, 1]);
        frame = frame + 1;
   end 
   index = index + 1;
end

X = X(:,1:num_frames);
%% Perform DMD

%calculate X1 and X2 matrices
X1 = X(:,1:end-1); X2 = X(:,2:end);

%find SVD and plot singular values
[U2,S2,V2] = svd(X1, 'econ');
threshold = 0.85;
r = find(cumsum(diag(S2) ./ sum(diag(S2))) > threshold  ,1);

figure(1)
plot(diag(S2)./max(diag(S2)), 'r.', 'markersize', 15);
title('Normalized Singular Values')
xlabel('Index ( j )')
ylabel('\sigma_j')
yticks(0:0.1:1)
set(gca, 'fontsize', 20);

% find constants for DMD reconstruction
U=U2(:,1:r); 
Sigma=S2(1:r,1:r);
V=V2(:,1:r);
Atilde = U'*X2*V/Sigma;
[W,D] = eig(Atilde);
Psi=X2*V/Sigma*W;

mu=diag(D);
omega=log(mu)/dt;
y0 = Psi\X(:, 1);  % pseudo-inverse initial conditions

[omega_min, mindex] = min(abs(omega));
foreground_modes_indices = find(abs(omega) > omega_min);

%%
clear all
%find low rank and sparse matrices

X_lowrank = y0(mindex).*Psi(:, mindex).*exp(omega(mindex).*t);
X_sparse = X - abs(X_lowrank);
R = X_sparse .* (X_sparse < 0); % places all negative entries in R


X_lowrank2 = abs(X_lowrank) + R; 
X_sparse2 = X_sparse - R; 

% Plot the foreground and background plots as well as R and the original

fig = figure;
clc, clear all
for frame = 1:size(X_lowrank2, 2)

   lowrank = imcomplement(reshape(X_lowrank2(:, frame), [imheight, imwidth]));
   sparse = imcomplement(reshape(X_sparse2(:, frame), [imheight, imwidth] ));
   rj = imcomplement(reshape(R(:, frame), [imheight, imwidth]));

   subplot(221)
   pcolor(flipud(lowrank)), shading interp, colormap(gray);
   title('Background (X\_lowrank)')
   set(gca, 'fontsize', 20);
  
   subplot(222)
   pcolor(flipud(sparse)), shading interp, colormap(gray);
   title('Foreground (X\_sparse)')
   set(gca, 'fontsize', 20);
    
   subplot(223)
   pcolor(flipud(sparse + lowrank)), shading interp, colormap(gray);
   title('Sum of both')
   set(gca, 'fontsize', 20);
   
   subplot(224)
   pcolor(flipud(rj)), shading interp, colormap(gray);
   title('R matrix')
   set(gca, 'fontsize', 20);
   
   drawnow;
   f(frame) = getframe(fig);
end
%%
implay(f) 