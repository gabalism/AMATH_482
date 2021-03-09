clear, clc;
%% GNR code loading data
[y, Fs] = audioread('GNR.m4a');
S = y' ;
n=length(y); L=n/Fs;  
t=(1:n)/Fs;
k=(2*pi)/(L)*[0:(n/2-1) -(n)/2:-1];
ks=fftshift(k);
%% Floyd code loading data
[y,Fs] = audioread('Floyd_trim.m4a');
%[y,Fs]=audioread('Floyd.m4a');
%y = y(1:end-1);
n=length(y); L = n/Fs;
t=(1:n)/Fs;
%rescale [-pi pi] domain fft assumes to fit our domain [-L L]
k=(2*pi)/(L)*[0:(n/2-1) -(n)/2:-1];
ks=fftshift(k);
S = y';
St = fft(S);
%% Filtering out high frequency noise
rect= @(x,a)  ones(1,numel(x)).*(abs(x)<a) ; %Creates a rectangular filter around a
%gaus = @(x,a,b) exp(-b*(x-a).^2);
filter = rect(ks,130*2*pi) ; %use 130 Hz as the highest frequency allowed
%filter = gaus(ks,200,0.01);
Sf = St.*fftshift(filter); % apply filter 
Sfil = ifft(Sf);

%% Filtering for Guitar solo
rect= @(x,a,b)  ones(1,numel(x)).*( abs(x) < b & abs(x) > a) ; % creates a band pass rectangle filter
filter = rect(ks, 230*2*pi, 600*2*pi ); %frequencies between 230 and 600
Sf = St.*fftshift(filter);
Sfil = ifft(Sf); % apply filter
plot(t,filter) %plot filter
%% Plot waveform and fft
subplot(2, 1, 1);
plot(t, S);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Waveform Signal');

subplot(2, 1, 2);
plot(ks, abs(fftshift(St))/max(abs(St)), 'r');
set(gca, 'XLim', [0 2e3]);
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Fast Fourier Transform');

sgtitle('GNR Waveform and FFT');

saveas(gcf, 'GNR_Waveform_FFT.jpg');
close(gcf);
%% Testing for different window sizes
% Plot spectrogram over varying window sizes
a = [100; 150; 250; 350; 500; 600; 750; 1000];
tslide = 0:0.1:L; %creates vector to "slide" across
%initialize matrix for spectrogram and max frequencies
Sgt_spec = zeros(length(tslide), n);
max_k = zeros(1, length(tslide));

for j=1:length(tslide)
    g = exp(-a .* (t - tslide(j)) .^ 2); %gabor filter with gaussian
    %g=exp(-a*(t - tslide(j)).^10); %super gaussian filter
    %g=(1-(t - tslide(j)).^2).*exp(-a*(t - tslide(j)).^2 ./ 2); %mexican
    %hat filter
    Sg = g .* S; %Use Sfil for Floyd Songs
    for k=1:length(a)
        Sgt = fft(Sg(k, :));
        Sgt_spec(j, :, k) = fftshift(abs(Sgt));
    end
end

for j=1:length(a)
    %Plot Spectrogram
    f = figure(j);
    pcolor(tslide, ks/(2*pi), Sgt_spec(:, :, j).');
    shading interp;
    set(gca,'Ylim',[0 500]);
    colormap(hot);
    xlabel('Time [sec]');
    ylabel('Frequency [Hz]');
    title(sprintf('GNR Spectrogram with Window Size %d', a(j)));
    set(gca , 'fontsize' ,25);
    saveas(f, sprintf('GNR_a=%d.jpg', a(j)));
    close(f);
end
%% Testing for different Time Increments
a = 150;
for dt=[0.05, 0.075, 0.1, 0.2, 0.5]

    tslide = 0:dt:L;   
    Sgt_spec = zeros(length(tslide), n);
    notefreq = zeros(1,length(tslide));
    for j=1:length(tslide)
        g = exp(-a * (t - tslide(j)).^ 2);
        %g=exp(-a*(t - tslide(j)).^10);
        %g= abs(t - tslide(j)) <= 1/(2*a);
        %g = (abs(t-tslide(j)) < 1);
        Sg = g .* S; %Use Sfil for Floyd songs
        Sgt = fft(Sg);
        [Max,Ind] = max(fftshift(abs(Sgt)));
        Sgt_spec(j, :) = fftshift(abs(Sgt));
        notefreq(j) = abs(ks(Ind)/(2*pi));
    end
    f=figure(1)
    pcolor(tslide, ks/(2*pi), Sgt_spec(:, :).');
    shading interp;
    %set(gca,'Ylim', [0 260]);
    set(gca,'Ylim',[220 630]);
    colormap(hot);
    xlabel('Time [sec]');
    ylabel('Frequency [Hz]');
    title(sprintf('Floyd Spectrogram (Only Bass)'));
    %gtext({'C5#';'D#';'F#';'G#';'#C';'F';'#F'})
    saveas(f, sprintf('Floyd_dt=%f.jpg', dt));
    close(f);
end

%% Code for GNR 
%same code as above, but tests only one dt and window
a = 150;
for dt=0.1

    tslide = 0:dt:L;   
    Sgt_spec = zeros(length(tslide), n);
    notefreq = zeros(1,length(tslide));
    for j=1:length(tslide)
        g = exp(-a * (t - tslide(j)).^ 2);
        Sg = g .*S;
        Sgt = fft(Sg);
        [Max,Ind] = max(fftshift(abs(Sgt)));
        Sgt_spec(j, :) = fftshift(abs(Sgt));
        notefreq(j) = abs(ks(Ind)/(2*pi));
    end
    f=figure(1)
    %create spectrogram for GNR
    pcolor(tslide, ks/(pi), Sgt_spec(:, :).');
    shading interp;
    %Create axis for y for better readability
    set(gca,'Ylim',[500 1500]);
    set(gca , 'fontsize' ,15);
    colormap(hot);
    xlabel('Time [sec]');
    ylabel('Frequency [Hz]');
    title(sprintf('GNR Spectrogram with Notes'));
    %Showcase where what each hotspot frequency is
    gtext({'C5#';'D#';'F#';'G#';'#C';'F';'#F'})
    saveas(f, sprintf('GNR Spectrogram with Notes.jpg'), dt);
    close(f);
end
%% Code for Floyd Song
%same code as above, but with the filter applied to Floyd 
for dt=0.075

    tslide = 0:dt:L;   
    Sgt_spec = zeros(length(tslide), n);
    notefreq = zeros(1,length(tslide));
    for j=1:length(tslide)
        g = exp(-a * (t - tslide(j)).^ 2);
        %apply filterered frequency to gabor
        Sg = g .* Sfil;
        Sgt = fft(Sg);
        [Max,Ind] = max(fftshift(abs(Sgt)));
        Sgt_spec(j, :) = fftshift(abs(Sgt));
        notefreq(j) = abs(ks(Ind)/(2*pi));
    end
    f=figure(1)
    pcolor(tslide, ks/(2*pi), Sgt_spec(:, :).');
    shading interp;
    %set(gca,'Ylim', [0 260]); %Use this range if we are doing the bass
    set(gca,'Ylim',[220 630]); %Use this range for the guitar solo
    set(gca ,'fontsize' ,15);
    colormap(hot);
    xlabel('Time [sec]');
    ylabel('Frequency [Hz]');
    %title(sprintf('Floyd Spectrogram (Only Guitar Solo)'));
    %title(sprintf('Floyd Spectrogram (Bass only)'));
    gtext({'B';'D';'E';'F#';'A';'B';'D'},'fontsize',15)
    saveas(f, sprintf('Floyd Spectrogram Guitar Solo.jpg'));
    close(f);
end
%% Make all Notes
%makes a vector of all the notes 
max_freqs = notefreq ;
notes = ["A", "A#", "B", "C", "C#","D", "D#", "E", "F", "F#", "G", "G#"];
%designates frequencies to each of these notes
fund_freqs = [27.5, 29.135, 30.863, 32.703, 34.648,36.708, 38.891, 41.203, 43.654,46.249, 48.99, 51.913];

allnotes = [];
allfreqs = [];

for j = 1:8
    for note = notes
        allnotes = [allnotes  strcat(note, num2str(j))];
    end
    for freq = fund_freqs
        %find frequencies at each octave
        % freqs double each time you move up scale
        allfreqs = [allfreqs freq*(2^(j-1))];
    end
end
%%
%create a final music score for the songs
music_score = [];

for freq = max_freqs
    %matches the frequencies we got to playable notes
    [min_val, index] = min(abs(freq - allfreqs));
    music_score = [music_score allnotes(index)];
end


