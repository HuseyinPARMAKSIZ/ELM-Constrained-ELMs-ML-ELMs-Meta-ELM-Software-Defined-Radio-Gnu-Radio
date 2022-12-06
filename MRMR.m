%% Rank features for classification using minimum redundancy maximum relevance (MRMR) algorithm %%
clear all; close all; clc; warning off
%------------------Loading datasets ----------------% 
load('Our_RF_tr-te_X_37F_70e30.mat')
Yd_train=train_data(:,1);
X=(train_data(:,2:size(train_data,2)));
Yd_test=test_data(:,1);
Xt=(test_data(:,2:size(test_data,2)));
SecIMFSay=8;
IMFOzellik=X; IMFClass= Yd_train;
%% feature selection - MRMR (Minimum Redundancy Maximum Relevance) algorithm
[index, etki] = fscmrmr(IMFOzellik,IMFClass); % Feature Selection Classification MRMR
ONEtkileri = [index(1:size(IMFOzellik,2));etki]'; % Etki miktari ve oznitelik numarasi bir tabloda
SiraliONEtkileri = sortrows(ONEtkileri,2,'descend');
x=SiraliONEtkileri (:,1); vals = SiraliONEtkileri (:,2); 
EtkiliSUTUNsiralamasi=x';
SecIMF = EtkiliSUTUNsiralamasi(1:SecIMFSay)
