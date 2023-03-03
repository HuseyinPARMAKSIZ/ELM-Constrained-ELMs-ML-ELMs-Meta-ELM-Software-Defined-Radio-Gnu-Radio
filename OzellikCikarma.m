clear all; close all; clc; warning off
tic
% Specify the folder where the files live.
%myFolder = 'F:\Doktora_My_RF_Dataset\RFson\'; % pwd
myFolder = 'D:\RFson\'; % pwd

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '**/*.*'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for a = 1 : length(theFiles)
    baseFileName = theFiles(a).name;
   if (string(baseFileName) ~= '.') && (string(baseFileName) ~= '..' )
           dosyaAd(a)=string(baseFileName);
    end 
end
%dosya(1:3397)=dosyaAd(3:3399);
dosya(1:7)=dosyaAd(3:9);
for x=1:7
pat='.'; dosya(x); N = split(dosya(x),pat);
Label(x)=N(1);   
 
seconds_to_read = 1;
Fs = 1E6*2;          %1 megabit/s data rate for basic DPSK ?
%Fs = 20E6*2;          %1 megabit/s data rate for basic DPSK ?
samples_to_read = floor(seconds_to_read * Fs);
%filename1 = 'F:\Doktora_My_RF_Dataset\RFson\'+dosya(x);   %as appropriate
filename1 = 'D:\RFson\'+dosya(x);   %as appropriate
[fid, msg] = fopen(filename1, 'r');
if fid < 0
    error('Failed to open file "%s" because "%s"', filename1, msg);
end
%data is interleaved real then complex. When we fread into two rows
%then the top row becomes the reals and the bottom row becomes the imag
data1 = fread(fid, [2 samples_to_read], '*float32');
%data1 = fread(fid,[2 samples_to_read],'uint8=>double');
fclose(fid);
%reformulate as complex
data1 = complex(data1(1, :), data1(2,:));

%converting abs 2D signal to 1D
mgA=abs(double(data1));

%removing zero values
mg = mgA(mgA~=0);

mg1= mg.';
mg2=mg1(:);
%converting complex 2D signal to 1D
cx1= data1.';
cx2=cx1(:);
% abs signal max point
k=max(mg2);
%determine the peaks that bigger than threshold we define
pks=findpeaks(mg2,Fs,'MinPeakHeight',k-50);
transient=[];
for i=1:length(mg2)
    %finding first peak 
    if mg2(i)==pks(1)
        b=i;
        firstpeak=b;
    end
end
for j=1:200
    %take 200 point before first peak
    transient(j)=cx2(firstpeak-200+j);
end

 f=figure(1);
 f.Position = [100 200 1700 750];
 subplot(2,1,1)
 plot (abs(data1(:)))
 plot (abs(data1(1:firstpeak)))
 title("Captured RF Signal Class " + Label(x))
 subplot(2,1,2)
 plot (abs(transient))
 title("Captured RF Signal Transient")

%obtained signal normalized
gol = ([transient-min(transient)]/[max(transient)-min(transient)]);
% %upsampling process
gol_up=interp(gol,10);
%features matrices create
a_mn=[]; a_var=[];  a_skw=[]; a_krts=[];
%divide the transient response 16 pieces
for i=1:16
    %every 125 sample take the values
    gol_up_d=gol_up(((i*125)-124):i*125); %take the value every 125 sample
    hill=hilbert(abs(gol_up_d));
    inamp1=abs(hill);
    inp1=unwrap(angle(hill));
    %instf1=diff(unwrap(angle(hill)))/((1/Fs)*2*pi);
    instf1 = Fs/(2*pi)*diff(unwrap(angle(hill)));
       
   %instantaneous amplitude
    a_mn(i)=mean(inamp1);          a_var(i)=var(inamp1);
    a_skw(i)=skewness(inamp1);     a_krts(i)=kurtosis(inamp1);
    %instantaneous phase
    p_mn(i)=mean(inp1);            p_var(i)=var(inp1);
    p_skw(i)=skewness(inp1);       p_krts(i)=kurtosis(inp1);
    %instantaneous frequency
    f_mn(i)=mean(instf1);          f_var(i)=var(instf1);
    f_skw(i)=skewness(instf1);     f_krts(i)=kurtosis(instf1);
end
x
 
RF=[Label(x),a_mn,a_var,a_skw,a_krts,p_mn,p_var,p_skw,p_krts,f_mn,f_var,f_skw,f_krts];
writematrix(RF,'RF_kayit.csv','Delimiter','tab','WriteMode','append')
end

-----------------------------------------------------------------------------------------
Version 2
clear all; close all; clc; warning off
tic
% Specify the folder where the files live.
%myFolder = 'F:\Doktora_My_RF_Dataset\RFson\'; % pwd
myFolder = 'D:\RF\'; % pwd

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '**/*.*'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for a = 1 : length(theFiles)
    baseFileName = theFiles(a).name;
   if (string(baseFileName) ~= '.') && (string(baseFileName) ~= '..' )
           dosyaAd(a)=string(baseFileName);
    end 
end
%dosya(1:3397)=dosyaAd(3:3399);
aaa=3;
dosya(1:aaa)=dosyaAd(3:aaa+2);
for x=1:aaa
pat='.'; dosya(x); N = split(dosya(x),pat);
Label(x)=N(1);   
 
seconds_to_read = 1;
Fs = 10E6*2;          %1 megabit/s data rate for basic DPSK ?
samples_to_read = floor(seconds_to_read * Fs);
%filename1 = 'F:\Doktora_My_RF_Dataset\RFson\'+dosya(x);   %as appropriate
filename1 = 'D:\RF\'+dosya(x);   %as appropriate
[fid, msg] = fopen(filename1, 'r');
if fid < 0
    error('Failed to open file "%s" because "%s"', filename1, msg);
end
%data is interleaved real then complex. When we fread into two rows
%then the top row becomes the reals and the bottom row becomes the imag
data1 = fread(fid, [2 samples_to_read], '*float32');
%data1 = fread(fid,[2 samples_to_read],'uint8=>double');
fclose(fid);
%reformulate as complex
data1 = complex(data1(1, :), data1(2,:));

%converting abs 2D signal to 1D
mgA=abs(double(data1));

%removing zero values
mg = mgA(mgA~=0);
%mg=medfilt1(mgilk,10,'omitnan','truncate');
mg1= mg.';
mg2=mg1(:);
%converting complex 2D signal to 1D
cx1= data1.';
cx2=cx1(:);

% abs signal max point
k=max(mg2);
%determine the peaks that bigger than threshold we define
pks=findpeaks(mg2,Fs,'MinPeakHeight',k-50);
transient=[];
for i=1:length(mg2)
    %finding first peak 
    if mg2(i)==pks(1)
        b=i;
        firstpeak=b;
    end
end
for j=1:200
    %take 200 point before first peak
    transient(j)=cx2(firstpeak-200+j);
end

 f=figure(1);
 f.Position = [200 200 1700 750];
 subplot(2,1,1)
 plot (abs(data1(:)))
 plot (abs(data1(1:firstpeak)))
 title("Captured RF Signal Class:  " + Label(x))
 subplot(2,1,2)
 plot (abs(transient))
 title("Transient Response of Captured RF Signal:  " + x)

%obtained transient signal normalized
gol = ([transient-min(transient)]/[max(transient)-min(transient)]);
% %upsampling process
gol_up=interp(gol,5);
%features matrices create
a_mn=[]; a_var=[];  a_skw=[]; a_krts=[];
%divide the transient response 16 pieces/subregion
for i=1:8
    %every 125 sample take the values
    gol_up_d=gol_up(((i*125)-124):i*125); %take the value every 125 sample
    % performing the HILBERT TRANSFORM
    hill=hilbert(abs(gol_up_d));
    % calculating the INSTANTANEOUS AMPLITUDE (ENVELOPE)
    inamp1=abs(hill);
    % calculating the INSTANTANEOUS PHASE
    inp1=unwrap(angle(hill));
    % calculating the INSTANTANEOUS FREQUENCY
    instf1=diff(unwrap(angle(hill)))/((1/Fs)*2*pi);
    %instf1 = Fs/(2*pi)*diff(unwrap(angle(hill)));
       
   %instantaneous amplitude
   n = length(inamp1);
   sorted_data = sort(inamp1);
   a_mdnH(i) = sorted_data(n/2 + 1.5);
   a_mdnL(i) = sorted_data(((n+1)/2)-1);
      
   
   % median = L + interval * (N / 2 - CF) / F 
   %L = lower limit of the median interval
   %N = total number of data points
   %CF = number of data points below the median interval
   %F = number of data points in the median interval
   interval=0.01;
   edges = [min(sorted_data):+interval:max(sorted_data)];
   frequency = [1:length(edges)/1.25];
   N = sum(frequency);
   cumfreq = cumsum(frequency);
   median_index = (N+1)/2; 
   group_index = find(cumfreq >= median_index, 1);
   class_width = edges(2) - edges(1);
   a_mdnG(i)= edges(group_index) + ((median_index - cumfreq(group_index-1)) / frequency(group_index)) * class_width;
   
   
   
   
   
    a_mdn(i)=median(inamp1,'omitnan') ;         
    a_hmn(i)=harmmean(inamp1);                  a_gmn(i)=geomean(inamp1);
    a_mn(i)=mean(inamp1);                       a_var(i)=var(inamp1);
    a_skw(i)=skewness(inamp1);                  a_krts(i)=kurtosis(inamp1);
    a_std(i)=std(inamp1);
    %instantaneous phase
    p_mn(i)=mean(inp1);                         p_var(i)=var(inp1);
    p_skw(i)=skewness(inp1);                    p_krts(i)=kurtosis(inp1);
    %instantaneous frequency
    f_mn(i)=mean(instf1);                       f_var(i)=var(instf1);
    f_skw(i)=skewness(instf1);                  f_krts(i)=kurtosis(instf1);
end
x
%16 subregion * 3 amplitude/phase/freq * 4 statistical_method  = 192Features 
%RF=[Label(x),a_mn,a_var,a_skw,a_krts,p_mn,p_var,p_skw,p_krts,f_skw,f_krts];
% stats.mean(_params), stats.geometric_mean(_params), stats.harmonic_mean(_params), stats.median(_params), stats.median_low(_params), stats.median_high(_params), stats.median_grouped(_params), stats.variance(_params), stats.stdev(_params))

RF=[Label(x) ,a_mdn,a_mdnL,a_mdnH,a_mdnG, a_std,a_hmn,a_gmn,a_mn,a_var,a_skw,a_krts,p_mn,p_var,p_skw,p_krts,f_mn,f_var,f_skw,f_krts];
writematrix(RF,'RF_010323_all.csv','Delimiter','tab','WriteMode','append')
%pause(3);
end

