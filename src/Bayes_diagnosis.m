function [posterior, pre_labels, bayes_prob, bayes_labels] = Bayes_diagnosis(train_data, gt_label, test_data, method, train_V_without_label, test_V_without_label, param)


if strcmp(method, 'RF') % Random forest
    model = TreeBagger(100, train_data', gt_label', 'Method','classification');
    [pre_results,posterior] = predict(model,test_data');
    pre_labels = zeros(size(pre_results,1),1);
    for i =1:size(pre_results,1)
        pre_labels(i) = str2num(pre_results{i,1});
    end
    
elseif strcmp(method, 'SVM')  % Support vector machine
    t = templateSVM('Standardize',1,'KernelFunction','linear');
    options = statset('UseParallel',true);
    SVMModel = fitcecoc(train_data',gt_label,'Options',options,'Learners',t, 'FitPosterior',true,'Holdout',0.15);
    Mdl = SVMModel.Trained{1};

    [pre_labels,neg_loss, pre_score, posterior] = predict(Mdl,test_data','Options',options);
    
elseif strcmp(method, 'MLP')  % Neural networks.
    
    gt_onehot = zeros(param.class_num,size(train_data,2));
    for i = 1:size(train_data,2)
        gt_onehot(gt_label(i),i) = 1;
    end
    trainFcn = 'trainscg';  % Scaled conjugate gradient backpropagation.
    
    % Create a Pattern Recognition Network
    hiddenLayerSize = 10;
    net = patternnet(hiddenLayerSize, trainFcn);
    
    % Setup Division of Data for Training, Validation, Testing
    net.divideParam.trainRatio = 70/100;
    net.divideParam.valRatio = 15/100;
    net.divideParam.testRatio = 15/100;
    
    % Train the Network
    [net,tr] = train(net,train_data,gt_onehot);
    
    posterior = sim(net,test_data);
    [~,pre_labels] = max(posterior,[],1);
    posterior = posterior';
    pre_labels = pre_labels';
    
    save_Model = net;
    
elseif strcmp(method, 'Naive')
    
    posterior = repmat(1/param.class_num,size(test_data,2),param.class_num);%ones(size(test_data,2),param.class_num)/param.class_num; 
    pre_labels = randi(param.class_num,  size(test_data,2), 1);
    
end

%% Bayes 
class_type = unique(gt_label);
theta_dist = {};
for i = 1:size(class_type,2)
    train_burden = train_V_without_label((param.healthy_pattern_num+1):end, gt_label == class_type(i));
%     tumor_fraction = sum(train_burden,1);
%     train_burden = [train_burden; tumor_fraction];
    for j = 1:size(train_burden,1)
        train_theta_one = train_burden(j,:);
        beta_ab = betafit(train_theta_one);
        theta_dist{i}(j,:) = beta_ab;
    end
end

bayes_prob = posterior';
test_burden = test_V_without_label((param.healthy_pattern_num+1):end,:);
% test_burden = [test_burden; sum(test_burden,1)];
for i = 1:size(test_V_without_label,2)
    for j = 1:size(class_type,2)
        for k = 1:size(theta_dist{j},1)
            beta_param = theta_dist{j}(k,:);
            likelyhood = betapdf(test_burden(k,i),beta_param(1),beta_param(2));
            bayes_prob(j,i) = bayes_prob(j,i)*likelyhood;
        end
    end
    bayes_prob(:,i) = bayes_prob(:,i)/sum(bayes_prob(:,i));
end

[bayes_p, bayes_labels] = max(bayes_prob,[],1);

end
% 
% function [classifier_prob, classifier_prediction, SRFD_Bayes_prob, SRFD_Bayes_prediction] = Bayes_diagnosis(train_data, test_data, method, train_sf, test_sf, param)
% 
% gt_label = param.train_gt_label;
% 
% if strcmp(method, 'RF')  % Random forest
%     model = TreeBagger(100, train_data', gt_label', 'Method','classification');
%     [pre_results,classifier_prob] = predict(model,test_data');
%     classifier_prediction = zeros(size(pre_results,1),1);
%     for i =1:size(pre_results,1)
%         classifier_prediction(i) = str2num(pre_results{i,1});
%     end
%     
% elseif strcmp(method, 'SVM')  % Support vector machine
%     t = templateSVM('Standardize',1,'KernelFunction','linear');
%     options = statset('UseParallel',true);
%     SVMModel = fitcecoc(train_data',gt_label,'Options',options,'Learners',t, 'FitPosterior',true,'Holdout',0.15);
%     Mdl = SVMModel.Trained{1};
%     [classifier_prediction,neg_loss, pre_score, classifier_prob] = predict(Mdl,test_data','Options',options);
%     
% elseif strcmp(method, 'MLP')  % Neural networks.
%     
%     gt_onehot = zeros(param.class_num,size(train_data,2));
%     for i = 1:size(train_data,2)
%         gt_onehot(gt_label(i),i) = 1;
%     end
%     trainFcn = 'trainscg';  % Scaled conjugate gradient backpropagation.
%     
%     % Create a Pattern Recognition Network
%     hiddenLayerSize = 10;
%     net = patternnet(hiddenLayerSize, trainFcn);
%     
%     % Setup Division of Data for Training, Validation, Testing
%     net.divideParam.trainRatio = 70/100;
%     net.divideParam.valRatio = 15/100;
%     net.divideParam.testRatio = 15/100;
%     
%     % Train the Network
%     [net,tr] = train(net,train_data,gt_onehot);
%     
%     classifier_prob = sim(net,test_data);
%     [~,classifier_prediction] = max(classifier_prob,[],1);
%     classifier_prob = classifier_prob';
%     classifier_prediction = classifier_prediction';
% end
% classifier_prediction = classifier_prediction';
% 
% %% Bayesian inference
% sample_num = [0,param.train_sample_num];
% theta_dist = {};
% for i = 1:param.class_num
%     train_burden = train_sf((param.healthy_pattern_num+1):end, (sum(sample_num(1:i)) + 1):sum(sample_num(1:(i+1))));
%     for j = 1:size(train_burden,1)
%         train_theta_one = train_burden(j,:);
%         beta_ab = betafit(train_theta_one);
%         theta_dist{i}(j,:) = beta_ab;
%     end
% end
% 
% SRFD_Bayes_prob = classifier_prob';
% for i = 1:size(test_sf,2)
%     for j = 1:param.class_num
%         for k = 1:size(theta_dist{j},1)
%             beta_param = theta_dist{j}(k,:);
%             likelyhood = betapdf(test_sf(param.healthy_pattern_num+k,i),beta_param(1),beta_param(2));
%             SRFD_Bayes_prob(j,i) = SRFD_Bayes_prob(j,i)*likelyhood;
%         end
%     end
%     SRFD_Bayes_prob(:,i) = SRFD_Bayes_prob(:,i)/sum(SRFD_Bayes_prob(:,i));
% end
% 
% [bayes_p, SRFD_Bayes_prediction] = max(SRFD_Bayes_prob,[],1);
% 
% end