clear all; close all; clc

load('subdata.mat')
%% Initial Conditions
L=10;%Spatial Domain
n=64;%Fourier Domain
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);
t = zeros(n,n,n);
%% Setting up unflitered signal to be averaged

%we have 49 samples of the time data
%use these different realizations of the data to
%average out zero mean noise
for j = 1:49
    sub(:,:,:)=reshape(subdata(:,j),n,n,n);
    t = t +fftn(sub); 
end
%normalize data and shift in the frequency domain
ave = abs(fftshift(t))/49; %averaging unfiltered signal
%find max to normalize
shape = reshape(ave,1,n^3);
M = max(abs(shape));
ave = ave/M;
%find frequency index where max occurs
index = find(ave == 1); %Finding the max to filter off of
%central frequencies
kx0 = Kx(index);
ky0 = Ky(index);
kz0 = Kz(index);

%% Filter of Data

tau = 0.0525; % filtering bandwidth
filter = exp(-tau * ((Kx - kx0).^2 +  (Ky - ky0).^2 + (Kz - kz0).^2));

% unshift filter as it is currently centered at zero.
filtershift = ifftshift(filter);


% apply this filter to frequency domain at each timestep
%to find submarine at each step
coords = zeros(49, 3);

figure(2)

for j=1:49 % for each realization of data
   sub(:,:,:)=reshape(subdata(:,j),n,n,n);
   subt = fftn(sub);
   subft = subt.*filtershift;
   subf = ifftn(subft);
   
   [M, I] = max(reshape(abs(subf), 1,n^3));
   coords(j, :) = [X(I), Y(I), Z(I)];
   
   isosurface(X,Y,Z, abs(subf), 0.6)
   axis([-10 10 -10 10 -10 10]),grid on, drawnow
   hold on;
   title('Filtered Signal in Spatial Domain', 'fontsize' , 20)
   xlabel('X');
   ylabel('Y');
   zlabel('Z');
   if j == 49
       saveas(gcf,strcat('submarine_isoplot.jpg'))
   end
end


%% Plotting the coordinates of the submarine
final_submarine_coordinate_xyz = coords(end, :);

figure (3)
plot3(coords(:, 1), coords(:,2), coords(:, 3), 'ro'), grid on;
hold on;
plot3(coords(:, 1), coords(:,2), coords(:, 3), 'b')
hold on;
txt = strcat('\leftarrow time 49: ', mat2str(final_submarine_coordinate_xyz));
text(final_submarine_coordinate_xyz(1), final_submarine_coordinate_xyz(2), final_submarine_coordinate_xyz(3), txt, 'fontsize', 20);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Submarine Position in Spatial Domain');
set(gca, 'fontsize', 15);
set(gcf, 'position', [100, 100, 600, 500]);
saveas(gcf, strcat('submarine_position.jpg'));


