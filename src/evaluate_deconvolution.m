function evaluate_deconvolution(test_tf,test_theta, param)

%% healthy evaluation
disp('>>>>>>>>Healthy evaluation<<<<<<<<<<');
hRMSE = healthy_eval(test_tf,test_theta, param);
disp(['RMSE=' num2str(hRMSE)]);
%% tumor fraction evaluation
disp('>>>>>>>>Tumor fraction evaluation<<<<<<<<<<');
[tRMSE,tPCC] = tumor_eval(test_tf,test_theta, param);
disp(['RMSE=' num2str(tRMSE)]);
disp(['PCC=' num2str(tPCC)]);
%
end

function RMSE = healthy_eval(test_tf,test_theta, param)
%% evaluate healthy samples
healthy_theta_pre = test_tf(1:param.test_sample_num(1));
healthy_theta_gt =1 -  test_theta(1:param.test_sample_num(1));
RMSE = sqrt(sum((healthy_theta_pre - healthy_theta_gt).^2)/size(healthy_theta_pre,2));
end

function [RMSE,PCC] = tumor_eval(test_tf,test_theta,param)
%% evaluate tumor samples
predict_theta = test_tf((param.test_sample_num(1)+1):end);
theta_gt = test_theta((param.test_sample_num(1)+1):end);
RMSE = sqrt(sum((predict_theta - theta_gt).^2)/size(predict_theta,2));
PCC = corr(predict_theta', theta_gt');
end