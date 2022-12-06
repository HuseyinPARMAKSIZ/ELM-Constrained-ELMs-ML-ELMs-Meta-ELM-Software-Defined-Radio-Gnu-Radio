clear all
RFdata9F = readtable('C:\Users\matlab\Desktop\h4ck3r-elm\F_38.csv');
label9 = RFdata9F {:,1};
X=RFdata9F {:,2:38};
Y=label9;
NormalizeX=X;

%% Normalizasyon
nmin=0; nmax=1; % Veriler nmin-nmax değerleri arasına uysun
NEdilecekSutunlar = [1 2 3 4 5 6 7 8 9 10 11 12 13] ;
for i = 1:length(NEdilecekSutunlar)
sutun = NEdilecekSutunlar(i);
data=X(:,sutun);
NormalizeX(:,sutun)=nmin + ((data - min(data))*(nmax-nmin) / (max(data)-min(data)));
end
clear i nmax nmin sutun NormalizeEdilecekSutunlar;
%% Verinin ne kadarını alırsak, orjinal veriyi % kaç temsil etmiş oluruz?
[PCAkatsayilar,PCANormalizeX, eigvalues] = pca(NormalizeX);
% Kaç bileşen alırsak veriyi % kaç temsil edebiliyoruz?
BTemsilAgirliklari=cumsum(eigvalues)./sum(eigvalues);
YuzdeBilesenTA=BTemsilAgirliklari.*100;
BilesenNo=[1:1:size(BTemsilAgirliklari,1)]';
fprintf('PCA sonucu oluşanlardan kaç bileşen alırsak,orijinal veri %% kaç temsil edilebilir?\n')
fprintf('%0.0f %0.3f\n',[BilesenNo,YuzdeBilesenTA]')
% Grafikle görelim
time=[1:1:size(BTemsilAgirliklari,1)];
plot(YuzdeBilesenTA,time,'--gs','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor',[0.5,0.5,0.5]);
for i = 1:length(BilesenNo)
text(YuzdeBilesenTA(i,:), i, sprintf('..:%s',num2str(YuzdeBilesenTA(i,:),'%.2f')))
end
set(gca,'YTick',time(1):1:time(end));
title('Percentage Representation of Data with PCA'); ylabel('Number of Principal Components'); xlabel('Representing the Original Data %'); grid on; grid minor;
