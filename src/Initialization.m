global param;
seed = 666;
rand('seed',seed);

%% parameter configuration
param.p = 0.5;  % 2-p norm
param.eta = 1000;  % coefficient of structural penalty
param.cancer_pattern_num = 2;
param.healthy_pattern_num = 7;
param.convergence_threshold = 1e-2;
param.maximum_iterations = 1000;
param.prior = 'SVM';  %RF,MLP, SVM

if strcmp(param.dataset_dir,'simulation_dataset')
    save_path = strcat('../results/',param.dataset_dir);
else
    save_path = strcat('../results/',param.dataset_dir, '/',param.dataset_name);
end
mkdir(save_path);

%% load training and  test data
file_dir = strcat('../data/',param.dataset_dir,'/');

if strcmp(param.dataset_dir,'simulation_dataset')
    
    disp(['Dataset: simulation dataset']);
    disp(['CNV event probability = ', param.dataset_name, '%']);
    train_data_dir = strcat(file_dir,'/train_data.mat');
    test_data_dir = strcat(file_dir,param.dataset_name,'/test_data.mat');
    load(train_data_dir);
    load(test_data_dir);
    
    load(strcat(file_dir,'/train_theta.mat'));
    load(strcat(file_dir,param.dataset_name,'/test_theta.mat'));
    evaluate_deconvolution_flag = true;
    
elseif strcmp(param.dataset_dir,'real_dataset')
    disp(['Dataset: real dataset']);
%     disp(['To be updated...']);
%     pause;
    if strcmp(param.dataset_name,'validation_real_data')
        data_dir = strcat(file_dir,'validation_real_data.mat');
        load(data_dir);
        test_data = validation_real_data;
        reference_dir = strcat(file_dir,'train_reference.mat');
        load(reference_dir);
    else
        data_dir = strcat(file_dir, param.dataset_name, '/');
        load(strcat(data_dir, 'train_data.mat'));
        load(strcat(data_dir, 'test_data.mat'));
        evaluate_deconvolution_flag = false;
        if strcmp(param.dataset_name,'Chen_data')
            param.oversampling_flag = true;
        else
            param.oversampling_flag = false;
        end
    end
end
