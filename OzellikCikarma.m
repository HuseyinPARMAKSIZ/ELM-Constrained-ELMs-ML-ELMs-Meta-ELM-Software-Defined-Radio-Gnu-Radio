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

