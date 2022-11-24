function [synthetic_data,synthetic_labels] = borderline_smote(train_data, train_labels, post_num,data_type, param)
% seed = 666;
% rand('seed',seed);
% global param;
param.ratio = 1;
param.knn = 10;
%
% train_sample_num = param.train_sample_num;
class_index = 1:param.class_num;

train_data_without_noise = [];
train_labels_without_noise = [];
danger_data = [];
danger_labels = [];
for i = 1:size(class_index,2)
    class_data = train_data(:,train_labels == class_index(i));
    randindex = randperm(size(class_data,2));
    class_data = class_data(:,randindex);

    nearest_in_num = zeros(1,size(class_data,2));
    for j = 1:size(class_data,2)
        single_h = class_data(:,j);
        dis_out = compute_tumor_dis(single_h, train_data, param, data_type);%sum((repmat(single_h,1,size(train_data,2))-train_data).^2);
        dis_out(dis_out==0) = inf;
        [sorted_dis,sorted_index] = sort(dis_out);
        
        nearest_labels = train_labels(sorted_index(1:param.knn));
        nearest_in_num(j) = sum(nearest_labels==class_index(i));

    end
    nearest_in_num(nearest_in_num==0) = inf;  %% noise data
    [sorted_num,sorted_index] = sort(nearest_in_num);
    sorted_index(sorted_num==inf) = [];
    
    if isempty(sorted_index)
        continue;   %sorted_index = randi(size(class_data,2),1,2);
    end
    
    danger_data = [danger_data,class_data(:,sorted_index(1:min(param.knn,size(sorted_index,2))))];
    danger_labels = [danger_labels,repmat(class_index(i), 1, min(param.knn,size(sorted_index,2)))];
    
    train_data_without_noise = [train_data_without_noise, class_data(:,sorted_index)];
    train_labels_without_noise = [train_labels_without_noise,i*ones(1,size(sorted_index,2))];

end

extend_index = unique(danger_labels);

if size(extend_index,2)<param.class_num
    test = 1;
end

extend_num = [];
for i = 1:size(extend_index,2)
    danger_num = sum(danger_labels==extend_index(i));
    train_num = sum(train_labels_without_noise==extend_index(i));
    temp = (max([post_num, danger_num])-train_num)/danger_num;
    extend_num = [extend_num, temp];
end
extend_num(extend_num<0) = 0;

synthetic_data = [];
synthetic_labels = [];

for i = 1:size(extend_index,2)
    class_data = danger_data(:,danger_labels == extend_index(i));
    class_in_data = train_data_without_noise(:,train_labels_without_noise == extend_index(i));
    for j = 1:size(class_data,2)
        single_h = class_data(:,j);
        dis_in = compute_tumor_dis(single_h, class_in_data, param, data_type);%sum((repmat(single_h,1,size(class_in_data,2))-class_in_data).^2);
        
        if strcmp(data_type, 'methylation')  %% avoid same data
            dis_in(dis_in==0) = inf;
            [sorted_dis,sorted_index] = sort(dis_in);
            sorted_index(sorted_dis==inf) = [];
        elseif strcmp(data_type, 'fraction')
            [sorted_dis,sorted_index] = sort(dis_in);
        end
            
        max_sample = min(param.knn,size(sorted_index,2));
        nearest_H = class_in_data(:,sorted_index(1:max_sample));  %compute the nearest neighbors
        
        if extend_num(i) == 0
            continue;
        elseif extend_num(i) < 1
            if rand(1) > extend_num(i)
                continue;
            else
                extend_num_i = 1;
            end
        else
            extend_num_i = round(extend_num(i));
        end
        
        for k = 1:extend_num_i
            rand_neighbor = nearest_H(:,randi(max_sample));  %random neighbor
            mixed_coff = rand(1); %random mixed coefficient
            synthetic_data = [synthetic_data, single_h + mixed_coff*(rand_neighbor -single_h)];
            synthetic_labels = [synthetic_labels, extend_index(i)];
        end
        
    end
end

end

function dis_out = compute_tumor_dis(single_h, train_data, param, data_type)
if strcmp(data_type, 'fraction')
    tumor_h = kron(eye(param.class_num-1),ones(1,param.cancer_pattern_num))*single_h(param.healthy_pattern_num+1:end,:);
    tumor_train = kron(eye(param.class_num-1),ones(1,param.cancer_pattern_num))*train_data(param.healthy_pattern_num+1:end,:);
    dis_out = sum((repmat(tumor_h,1,size(tumor_train,2)) - tumor_train).^2);
elseif strcmp(data_type, 'methylation')
    dis_out = sum((repmat(single_h,1,size(train_data,2)) - train_data).^2);
end
end