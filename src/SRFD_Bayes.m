function [train_reference, SRFD_Bayes_results] = SRFD_Bayes(train_data, test_data)

global param;

if strcmp(param.dataset_dir,'simulation_dataset')
    train_labels = train_data(end,:);
    [train_labels,sorted_train_index] = sort(train_labels);
    train_data = train_data(1:(end-1),sorted_train_index);

    test_labels = test_data(end,:);
    [test_labels,sorted_test_index] = sort(test_labels);
    test_data = test_data(1:(end-1),sorted_test_index);

    train_class_index = unique(train_labels);
    param.class_num = size(train_class_index,2);
    param.train_sample_num = [];
    for i = 1:size(train_class_index,2)
        param.train_sample_num = [param.train_sample_num, sum(train_labels==train_class_index(i))];
    end

    test_class_index = unique(test_labels);
    param.test_sample_num = [];
    for i = 1:size(test_class_index,2)
        param.test_sample_num = [param.test_sample_num, sum(test_labels==test_class_index(i))];
    end

elseif strcmp(param.dataset_name,'Xu_data')
    param.class_num = 2;
    param.train_sample_num = [418,704];
    param.test_sample_num = [417,346];
    
elseif strcmp(param.dataset_name,'Chen_data')
    param.class_num = 6;
    param.train_sample_num = [207,5,20,39,45,35];
    param.test_sample_num = [207,2,3,29,11,34];

end

param.train_gt_label = [];
for i = 1:size(param.train_sample_num,2)
    param.train_gt_label = [param.train_gt_label, i*ones(1,param.train_sample_num(i))];
end
param.test_gt_label = [];
for i = 1:size(param.test_sample_num,2)
    param.test_gt_label = [param.test_gt_label, i*ones(1,param.test_sample_num(i))];
end

%% Semi-reference-free deconvolution (SRFD)
fprintf('1. Performing semi-reference-free deconvolution (SRFD)\n');
fprintf('This may take a while...\n');
[train_W,train_R,~] = SRFD(train_data,param);
train_reference = train_W;

%% Oversampling for Chen_data using borderline-SMOTE
% load('D:\bioinformatics\code\data\TCGA_new\NC_data\diagnosis\late_training\result/train_U_0.mat');
% load('D:\bioinformatics\code\data\TCGA_new\NC_data\diagnosis\late_training\result/train_V_0.mat');
% train_W = train_U;
% train_R = train_V;

if param.oversampling_flag
    [os_R,os_labels] = borderline_smote(train_R, param.train_gt_label, 40,'fraction', param);
    reconstruction_err = train_data-train_W*train_R;
    os_data = mix_WH(train_W, os_R, os_labels,reconstruction_err,param);

    train_data = [train_data, os_data];
    param.train_gt_label = [param.train_gt_label, os_labels];
    [param.train_gt_label,sorted_index] = sort(param.train_gt_label);
    train_data = train_data(:,sorted_index);
end

%% Compute source fraction (deconcolution without structural constraints)
fprintf('2. Estimating source fractions\n');
fprintf('This may take a while...\n');
train_sf =  deconvolution_with_reference(train_data, train_W, zeros(1,size(train_data,2)), param);
test_sf =  deconvolution_with_reference(test_data, train_W, zeros(1,size(test_data,2)), param);

%% Bayesian diagnosis
fprintf('3. Performing Bayesian inference\n');
bayes_pred_prob_all = zeros(param.class_num, size(test_data,2),10);
warning off;
for j = 1:10
    [classifier_prob, classifier_prediction, Bayes_prob, Bayes_prediction] = Bayes_diagnosis(train_data, param.train_gt_label,...
        test_data, param.prior, train_sf, test_sf, param);
    bayes_pred_prob_all(:,:,j) = Bayes_prob;
end
mean_bayes_pre_prob = mean(bayes_pred_prob_all,3);
[SRFD_Bayes_prob, SRFD_Bayes_prediction] = max(mean_bayes_pre_prob,[],1);

% for cancer detection
SRFD_Bayes_prob(SRFD_Bayes_prediction==1) = 1- SRFD_Bayes_prob(SRFD_Bayes_prediction==1);

%% final results
% Three rows: the first row is cancerous probability; the second row is predicted labels; the third row is predicted tumor fraction.
tumor_fraction = sum(test_sf((param.healthy_pattern_num+1):end,:),1);
SRFD_Bayes_results = [SRFD_Bayes_prob;SRFD_Bayes_prediction;tumor_fraction];
fprintf('Finish!\n');

end
