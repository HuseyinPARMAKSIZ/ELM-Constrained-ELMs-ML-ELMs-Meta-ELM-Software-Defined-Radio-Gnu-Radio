
All codes are in the testing phase yet!!!
It can be used after the test phase is completed with the permission of the author.

clear all; close all; clc; warning off
tic
% Specify the folder where the files live.
myFolder = 'F:\Doktora_My_RF_Dataset\RFson\'; % pwd
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '**/*.*'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for a = 1 : length(theFiles)
    baseFileName = theFiles(a).name;
   if (string(baseFileName) ~= '.') && (string(baseFileName) ~= '..' )
           dosyaAd(a)=string(baseFileName);
    end 
end
dosya(1:3397)=dosyaAd(3:3399);
for x=1:3397
pat='.'; dosya(x); N = split(dosya(x),pat);
Label(x)=N(1);   
 
seconds_to_read = 1;
Fs = 20E6*2;          %1 megabit/s data rate for basic DPSK ?
samples_to_read = floor(seconds_to_read * Fs);
filename1 = 'F:\Doktora_My_RF_Dataset\RFson\'+dosya(x);   %as appropriate
%filename1 = 'D:\RFson\'+dosya(x);   %as appropriate
[fid, msg] = fopen(filename1, 'r');
if fid < 0
    error('Failed to open file "%s" because "%s"', filename1, msg);
end
%data is interleaved real then complex. When we fread into two rows
%then the top row becomes the reals and the bottom row becomes the imag
data1 = fread(fid, [2 samples_to_read], '*float32');
fclose(fid);
%reformulate as complex
data1 = complex(data1(1, :), data1(2,:));
%converting abs 2D signal to 1D
mg=abs(double(data1));
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

% f=figure(1);
% f.Position = [100 200 1700 750];
% % subplot(4,1,1)
% %plot (abs(data1(:)))
% plot (abs(data1(1:10240)))
% title("Class " + Label(x))

%obtained signal normalized
gol = ([transient-min(transient)]/[max(transient)-min(transient)]);
% %upsampling process
gol_up=interp(gol,10);
%----------------- IMFs'lerden Faz, Frekans ve Genlik---------------------%
% IMFs=emd(real(gol));
  IMFs_i=emd(real(gol_up));
  IMFs_i2=(IMFs_i(:,2));
  a_p2p = max(abs(IMFs_i2))*2;
  a_rms = sqrt(mean(IMFs_i2.^2));
% VMFs=vmd(real(gol));
% VMFs_i=vmd(real(gol_up));

% %----------------- IMFs'lerden Faz, Frekans ve Genlik---------------------%
% for ii=2:4
%     IMFs_p=unwrap(angle(IMFs(:,ii)))';
%     IMFs_f=diff(unwrap(angle(IMFs(:,ii))))/((1/Fs)*2*pi)'; 
%     IMFs_amp=abs(IMFs(:,ii))'; 
%     
%     VMFs_p=unwrap(angle(VMFs(:,ii)))';
%     VMFs_f=diff(unwrap(angle(VMFs(:,ii))))/((1/Fs)*2*pi)';
%     VMFs_amp=abs(VMFs(:,ii))'; 
% 
%     IMFs_a_mn(ii-1)=mean(IMFs_amp); 
%     IMFs_a_var(ii-1)=var(IMFs_amp); 
%     IMFs_a_skw(ii-1)=skewness(IMFs_amp); 
%     IMFs_a_krts(ii-1)=kurtosis(IMFs_amp);
% 
%     VMFs_a_mn(ii-1)=mean(VMFs_amp); 
%     VMFs_a_var(ii-1)=var(VMFs_amp); 
%     VMFs_a_skw(ii-1)=skewness(VMFs_amp); 
%     VMFs_a_krts(ii-1)=kurtosis(VMFs_amp);
%  
% %----------------- IMFs'lerden Faz, Frekans ve Genlik---------------------%
% end
%features matrices create
a_mn=[];
a_var=[];
a_skw=[];
a_krts=[];

%divide the transient response 16 pieces

for i=1:16
    %every 125 sample take the values
    gol_up_d=gol_up(((i*125)-124):i*125); %take the value every 125 sample
    hill=hilbert(abs(gol_up_d));
    inamp1=abs(hill);
    inp1=unwrap(angle(hill));
    instf1=diff(unwrap(angle(hill)))/((1/Fs)*2*pi);

   %instantaneous amplitude
    a_std(i)=std(inamp1);
    a_mn(i)=mean(inamp1);
    a_gmn(i)=geomean(inamp1);
    a_hmn(i)=harmmean(inamp1);
    a_var(i)=var(inamp1);
    a_skw(i)=skewness(inamp1);
    a_krts(i)=kurtosis(inamp1);
    %instantaneous phase
    p_std(i)=std(inp1);
    p_mn(i)=mean(inp1);
    p_hmn(i)=harmmean(inp1);
    p_var(i)=var(inp1);
    p_skw(i)=skewness(inp1);
    p_krts(i)=kurtosis(inp1);
    %instantaneous frequency
    f_std(i)=std(instf1);
    f_mn(i)=mean(instf1);
    f_hmn(i)=harmmean(instf1);
    f_var(i)=var(instf1);
    f_skw(i)=skewness(instf1);
    f_krts(i)=kurtosis(instf1);
end

x



% a_std =mean(a_std,'all'); a_mn = mean(a_mn,'all'); a_gmn = mean(a_gmn,'all');a_hmn = mean(a_hmn,'all');a_var = mean(a_var,'all'); a_skw = mean(a_skw,'all'); a_krts = mean(a_krts,'all');
% p_std =mean(p_std,'all');p_mn= mean(p_mn,'all');p_hmn = mean(p_hmn,'all'); p_var= mean(p_var,'all');p_skw= mean(p_skw,'all'); p_krts= mean(p_krts,'all');
% f_std =mean(f_std,'all');f_mn= mean(f_mn,'all');f_hmn = mean(f_hmn,'all'); f_var= mean(f_var,'all'); f_skw= mean(f_skw,'all'); f_krts= mean(f_krts,'all');

 
% IMFs_a_mn = mean(IMFs_a_mn,'all'); IMFs_a_var = mean(IMFs_a_var,'all'); IMFs_a_skw = mean(IMFs_a_skw,'all'); IMFs_a_krts = mean(IMFs_a_krts,'all');
% VMFs_a_mn = mean(VMFs_a_mn,'all'); VMFs_a_var = mean(VMFs_a_var,'all'); VMFs_a_skw = mean(VMFs_a_skw,'all'); VMFs_a_krts = mean(VMFs_a_krts,'all');

%F=[Label(x),PSD_mn,PSD_var,PSD_skw,PSD_krts,a_mn,a_var,a_skw,a_krts,p_mn,p_var,p_skw,p_krts,f_mn,f_var,f_skw,f_krts,IMFs_a_mn,IMFs_a_var,IMFs_a_skw,IMFs_a_krts,VMFs_a_mn,VMFs_a_var,VMFs_a_skw,VMFs_a_krts];
F=[Label(x),a_p2p,a_rms,a_std,a_mn,a_gmn,a_hmn,a_var,a_skw,a_krts,p_std,p_mn,p_hmn,p_var,p_skw,p_krts,f_std,f_hmn,f_skw,f_krts];

writematrix(F,'RF_kayit.csv','Delimiter','tab','WriteMode','append')
s = seconds(toc) ;
end
s.Format = 'hh:mm'

