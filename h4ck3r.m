
All codes are in the testing phase yet!!!
It can be used after the test phase is completed with the permission of the author.
  
%% Written by me last quarter of 2022 (for RF classification problems.)
clear all; close all; clc; warning off
%------------------Loading datasets ----------------% 
%-----Load Ready-RF dataset -----------------------------%
 load('Ready9Ftrvete.mat');
% load('Ilk_RF_set.mat')
%-----Load Our RF dataset -------------------------------%
% load('Our_RF_tr-te_43F.mat');
%-----Load Our RF dataset after PCA ---------------------%
%load('Our_RF_43_PCA_tr-te_26F.mat');
%-----Load Our RF dataset after PCA and normalized-------%
%load('Our_RF_43_PCANormalize1-13_tr-te_19F.mat');
%load('Our_RF_data_21F_70e30.mat');
% load('Our_RF_tr-te_X_37F_70e30.mat');
Yd_train=train_data(:,1)'; X=train_data(:,2:size(train_data,2))';
Yd_test=test_data(:,1)';   Xt=test_data(:,2:size(test_data,2))';
%clear train_data;          clear test_data; %Release raw training & testing data array
%------------------Loading datasets ----------------% 
%load('RF_3397R-274F-20C_70e30.mat');
% Yd_train=train_data(:,1)'; X=train_data(:,2:9)'; X= zscore(X);
% Yd_test=test_data(:,1)';   Xt=test_data(:,2:9)'; Xt= zscore(Xt); %size(test_data,2)
%-----Definitions--------%
ite=1; Elm_Type = 1; REGRESSION=0; CLASSIFIER=1; 
Fn=size(X,1); % attribute_number 9 (Label Mean GMean HMean Median MedianL MedianH MedianG Variance Stdev) 
N=size(X,2); % the number of samples for train data N=2790 = 3985*0,7
Nt=size(Xt,2); % the number of samples for test data Nt=1195 = 3985*0,3
C=max(eig(X'*X)); % regularization factor for meta-elm 3.312241605718222e+06
time=[]; RMSEtrain=[]; RMSEtest=[]; % cost array in terms of RMSE
Dogruluk_Meta_ELMtrain=[];Dogruluk_Meta_ELMtest=[]; Avg_Meta_ELMtrain=[];Avg_Meta_ELMtest=[];
TrFcn = 'logsig';%transfer func of the node
M=[500,600]; % Meta groups 
Nnode=[2,8]; % sub-ELMs of each Meta
Msize=size(M,2); % size of deneme Meta groups
Nnodesize=size(Nnode,2); % size of deneme Nnode
BosNnode=[]; BosNnodeT=[];
DenSay=3; % iteration size
%----- Preprocessing the data of classification
if Elm_Type~=REGRESSION
    sorted_target=sort(cat(2,Yd_train,Yd_test),2); % 1*3985 bütün sınıflar sıralı yan yana 
    label=zeros(1,1); % Find and save in 'label' class label from training and testing data sets
    label(1,1)=sorted_target(1,1);
    j=1;
    for i = 2:(N+Nt)
        if sorted_target(1,i) ~= label(1,j)
            j=j+1;
            label(1,j) = sorted_target(1,i);
        end
    end
    number_class=j; % class size
    NofOutNeurons=number_class;
    %----- Preprocessing the targets of training
    temp_T=zeros(NofOutNeurons, N);
    for i = 1:N
        for j = 1:number_class
            if label(1,j) == Yd_train(1,i)
                break;
            end
        end
    temp_T(j,i)=1;
    end
    Yd_train=temp_T*2-1;
    Yd_train=Yd_train'; % 2790*4 hangi sınıfa aitse 1 diğerleri -1
    %----- Preprocessing the targets of testing
    temp_TV_T=zeros(NofOutNeurons, Nt);
    for i = 1:Nt
        for j = 1:number_class
            if label(1,j) == Yd_test(1,i)
                break;
            end
        end
    temp_TV_T(j,i)=1;
    end
    Yd_test=temp_TV_T*2-1;
    Yd_test=Yd_test';   % 1195*4 hangi sınıfa aitse 1 diğerleri -1  
end % Classification_if
for aa=1:Msize % ELM size
    MM=M(aa);
    DS=round(N/MM);% Data size
    for bb=1:Nnodesize % Each ELM has node_size nodes/neuron
        NN=Nnode(bb);
        while ite<=DenSay %---LOOP for numerical experiments
        tstart = tic;
        Xm=[]; Nm=[]; W=[]; b=[]; Beta=[]; HM=[]; HMtest=[]; H=[]; Hm=[]; %Accuracy_Meta_ELMtrain=[]; Accuracy_Meta_ELMtest=[];
        %---Stage:1 Loop for training of sub ELMs
            for m=1:MM
            %---SLFN structural parameters
                if m==1 
                        Xm=X(:,1:m*DS); 
                    else 
                        Xm=X(:,(m-1)*DS+1:min(m*DS, max(size(X)))); 
                end %Xm=X(:,(m-1)*M+1:m*M); end
                Fnm=size(Xm,1); %the number of external inputs
                Nm=size(Xm,2);%the number of subset samples=N/M
                out_size=size(Yd_train,2);%the number of output
            %--- Construction of mth SLFN and initial assigmants
                if NN<Fnm
                    W(:,:,m)=orth(rand(Fnm,NN)); %input weights of hidden nodes
                else
                    W(:,:,m)= (get_orthonormal(NN,Fnm))'; %randn(n,Nnode)*2-1;
                end
                b(:,:,m)=(get_orthonormal(NN,1))';
            %--- Construction of H matrice
                H=feval(TrFcn,(Xm'*W(:,:,m)+repmat(b(:,:,m),Nm,1))); % DS*Nnode
            %--- ELM training,
                if m==1
                    Beta(:,:,m)=(inv(H'*H)*H')*(Yd_train(1:DS,:));% 
                else
                    Beta(:,:,m)=(inv(H'*H)*H')*Yd_train(((m-1)*DS+1):min(m*DS, max(size(X))),:); %Nnode * Class * M
                end %if
        end % for m 1 .. MM
        %---Stage:2 Training of Meta-ELMs
            % H computation of MetaELM
            for j=1:MM
                Hm(:,:,j)=feval(TrFcn,(X'*W(:,:,j)+repmat(b(:,:,j),N,1)))*Beta(:,:,j); % train * Class * M
                HM=[HM Hm(:,:,j)]; %Her M çıkışını yan yana ekleme
            end
            %--- Training of MetaELM
            BetaM=pinv(HM)*Yd_train; %implementation without regularization factor //refer to 2006 Neurocomputing paper
            %BetaM=HM'*inv(eye(size(HM * HM'))/C + HM * HM')*Yd_train; % faster method 1 //refer to 2012 IEEE TSMC-B paper
            %BetaM=((eye(size(HM * HM'))/C+HM * HM') \ HM)' * Yd_train; % faster method 2 //refer to 2012 IEEE TSMC-B paper
            
            %OutputWeight=pinv(H') * T';% implementation without regularization factor //refer to 2006 Neurocomputing paper
            %OutputWeight=inv(eye(size(H,1))/C+H * H') * H * T'; % faster method 1 //refer to 2012 IEEE TSMC-B paper
            %OutputWeight=(eye(size(H,1))/C+H * H') \ H * T'; % faster method 2 //refer to 2012 IEEE TSMC-B paper

            
            
            time(ite)= toc(tstart);
        %--- PERFORMANCE FOR TRAINING DATA with founded output weight (Beta)parameters
            Ytrain=HM*BetaM;
            Y=Ytrain'; % 4*2790 sınıfa dahil olanlar pozitif diğerleri negatif
            %--- Cost computatiton for training
            Cost_train=sqrt(mse(Ytrain,Yd_train)); %sqrt(sum((Yd_train-Ytrain).^2)/N);
            RMSEtrain(ite)=Cost_train;  
        %--- PERFORMANCE FOR TESTING DATA with founded rule parameter
            %--- H computation of MetaELM
            for k=1:MM
                Htest(:,:,k)=feval(TrFcn,(Xt'*W(:,:,k)+repmat(b(:,:,k),Nt,1)))*Beta(:,:,k); % test * Class * M
                HMtest=[HMtest Htest(:,:,k)];
            end
            Ytest=HMtest*BetaM;
            TY=Ytest'; % 4*1195 sınıfa dahil olanlar pozitif diğerleri negatif
            %--- Cost computatiton for training
            Cost_test=sqrt(mse(Ytest,Yd_test));%sqrt(sum((Yd_test-Ytest).^2)/Nt);
            RMSEtest(ite)=Cost_test; %Cost saving in an array
            %---  Calculate training & testing classification accuracy
            if Elm_Type == CLASSIFIER
                WrongGuessClass=0;
                WrongGuessClassT=0;
                Yd_tr=Yd_train'; Yd_te=Yd_test';
                %--- Calculate training classification accuracy
                    for i = 1 : size(Yd_tr, 2)
                        [x, i_Labeld]=max(Yd_tr(:,i));
                        [x, i_labela]=max(Y(:,i));
                        if i_labela~=i_Labeld
                            WrongGuessClass=WrongGuessClass+1;
                        end
                    end
                TrainingAccuracy=1-WrongGuessClass/size(Yd_tr,2);
                Accuracy_Meta_ELMtrain(ite)=TrainingAccuracy;
                %--- Calculate training testing classification accuracy
                    for i = 1 : size(Yd_te, 2)
                        [x, i_Labeld]=max(Yd_te(:,i));
                        [x, i_labela]=max(TY(:,i)); % TY: the actual output of the testing data
                        if i_labela~=i_Labeld
                            WrongGuessClassT=WrongGuessClassT+1;
                        end
                    end
                TestingAccuracy=1-WrongGuessClassT/size(Yd_te,2);
                Accuracy_Meta_ELMtest(ite)=TestingAccuracy;
            end % Classifier_if
            ite=ite+1; %iteration increment
        end % while ite
        %Accuracy_Meta_ELMtrain
        %Accuracy_Meta_ELMtest
        ite=1;
        % disp('the elapsed time metrics for training');
        average_time=mean(time);    best_time=min(time);    worst_time=max(time);   std_dev_time=std(time);
        % disp('The metrics for training data');
        average_cost=mean(RMSEtrain);       best_cost=min(RMSEtrain);   worst_cost=max(RMSEtrain);  std_dev_cost=std(RMSEtrain);
        % disp('The metrics for testing data');
        Taverage_cost=mean(RMSEtest);   Tbest_cost=min(RMSEtest);   Tworst_cost=max(RMSEtest);  Tstd_dev_cost=std(RMSEtest);
        %Avg_Meta_ELMtrain(aa,bb)=mean(Accuracy_Meta_ELMtrain);  Avg_Meta_ELMtest(aa,bb)=mean(Accuracy_Meta_ELMtest);
            if M~=0
                BosNnode(aa,bb)=mean(Accuracy_Meta_ELMtrain); BosNnodeT(aa,bb)=mean(Accuracy_Meta_ELMtest);
                %EE=nonzeros(BosNnode); TT=nonzeros(BosNnodeT);
            end
         disp('---------------------------------------------------------------------------------------------------------');
         disp('[  RegFac.    M       Nnode     TrAcc.       TeAcc.   ActivationFunc.  ]');
         XX=[      C     MM       NN    BosNnode(aa,bb)  BosNnodeT(aa,bb)  string(TrFcn)];
         disp(XX); 
    end
end




