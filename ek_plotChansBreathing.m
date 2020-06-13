batchFile = 'batch.txt.keep';
fid=fopen(batchFile,'r');

% fn = 'gr98bu75_120819_183052.1805.cbin';

i = 1;
j = 0;
while (i <3)
    i = i + 1;

    
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end
    
%     if i ~= 2
%         continue
%     end

    [dat,Fs]=ReadCbinFile(fn);
    song = dat(:, 1);
    chan1 = dat(:, 2);
    if max(chan1) < 1000
        continue
    end
    % stim2 = dat(:, 3);
    
    % for looking at breathing data
    breathFilt = highpass(chan1, 1, Fs);
    breathFilt = lowpass(breathFilt, 100, Fs);
    z = zeros(1, length(breathFilt));

    figure; 
    bigFig(1) = subplot(2, 1, 1); 
%     plot(song); hold on
    
    
%     plot(chan1); hold on
    
%     N = 1024; OVERLAP = 1020; sigma = 3; F_low = 500; F_high = 10000; filter_type = 'hanningffir';
%     SPTH = 10;
%     target = 'b';
%     stimWin = [-500 500]; % in ms
%     w=hamming(N);
% 
%     filtsong=pj_bandpass(song,Fs,500,10000,'butterworth');%'hanningffir'
%     
%     [S,F,T] = spectrogram(filtsong,w,OVERLAP,N,Fs);
%     T = linspace(-.5,.5,length(T));
%     pp = find(abs(S)<=SPTH); %Find entries with very low power + scale up. Makes display nicer.
%     S(pp) = SPTH; S = log(abs(S));
%     hold on
%     imagesc(T,F,S)
%     ylim([0 10000])
%     colormap(jet);
%     caxis([8 15])
   
    
    bigFig(2) = subplot(2, 1, 2); 
    plot(breathFilt); hold on; plot(z, 'r')
    
%     bigFig(2) = subplot(3, 1, 3); 
%     plot(cumsum(breathFilt))
%     plot(chan1);
    % plot(stim2)
    
    
    linkaxes([bigFig(1), bigFig(2)], 'x');
end

fclose(fid)

% need breath segmentation algorithm


