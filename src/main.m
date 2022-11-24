clc;
clear;

%% initialization
global param;
param.dataset_dir = 'real_dataset';  % choose the dataset type: 'simulation_dataset' and 'real_dataset'
param.dataset_name = 'validation_real_data';  % choose a dataset: 
% simulation_dataset: 'CNV10', 'CNV30', 'CNV50'
% real_dataset: 'validation_real_data (GSE122126, GSE108462, GSE129373)', 'Xu_data', and Chen_data''


Initialization;  % data loading, parameter setting, etc.

if strcmp(param.dataset_dir, 'real_dataset') && strcmp(param.dataset_name,'validation_real_data')
    validation_tf =  deconvolution_with_reference(validation_real_data, train_reference, zeros(1,size(validation_real_data,2)), param);
    save(strcat(save_path,'/validation_tf.mat'), 'validation_tf');
else

    %% SRFD-Bayes algorithm
    %-------------------------------------------------------------------%
    %----Input:
    %------training data matrix, test data matrix, parameters
    %----Output:
    %------diagnostic results:
    %------The first row represents the predicted class clabels
    %------The second row represents the predicted tumor fraction
    [train_reference, SRFD_Bayes_results] = SRFD_Bayes(train_data, test_data);
    %-------------------------------------------------------------------%

    %% save results
    save(strcat(save_path,'/SRFD_Bayes_results.mat'), 'SRFD_Bayes_results');

    if strcmp(param.dataset_dir, 'simulation_dataset')
        save(strcat('../data/real_dataset/train_reference.mat'), 'train_reference');

        %% evaluate deconvolution (simulation dataset only)
        if evaluate_deconvolution_flag
            fprintf('Evaluating the deconvolution performance\n');
            evaluate_deconvolution(SRFD_Bayes_results(3,:), test_theta, param);
        end
    end

    %% evaluate diagnosis
    if param.class_num ==2
        [sorted_prob, sorted_index] = sort(SRFD_Bayes_results(1,:));
        sorted_label = param.test_gt_label(sorted_index)-1;

        [FPR,TPR] = cal_roc(sorted_prob, sorted_label);
        auc = abs(trapz(FPR,TPR));
        disp(['AUC=' num2str(auc)]);
    else
        disp('>>>>>>>>Evaluating the diagnostic performance...');
        average_ACC = evaluate_diagnosis(SRFD_Bayes_results(2,:), param);
        disp(['average ACC=' num2str(average_ACC)]);
    end
end
