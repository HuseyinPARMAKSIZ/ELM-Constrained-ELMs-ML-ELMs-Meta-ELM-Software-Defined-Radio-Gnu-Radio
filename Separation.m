
All codes are in the testing phase yet!!!
It can be used after the test phase is completed with the permission of the author.
%-- training and test data separation
clear all
load('Ready6Ftrvete.mat');

% randomly select indexes to split data into 70% 
    % training set, 0% validation set and 30% test set.
    [train_idx, ~, test_idx] = dividerand(3397, 0.7, 0, 0.3);
    % slice training data with train indexes 
    %(take training indexes in all 10 features)
    train_data = RFkayit(train_idx, :);
    % select test data
    test_data = RFkayit(test_idx, :);

-------------------------------------------------------------
clear all
load('RF_020323.mat');

% randomly select indexes to split data into 70% 
    % training set, 0% validation set and 30% test set.
    [train_idx, ~, test_idx] = dividerand(2443, 0.7, 0, 0.3);
    % slice training data with train indexes 
    %(take training indexes in all 10 features)
    train_data = RF010323all(train_idx, :);
    % select test data
    test_data = RF010323all(test_idx, :);
