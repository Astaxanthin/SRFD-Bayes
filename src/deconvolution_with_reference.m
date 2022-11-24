function test_H = deconvolution_with_reference(test_data, train_W, predicted_label, param)
seed = 666;
rand('seed',seed);

class_num=param.class_num;
healthy_pattern_num = param.healthy_pattern_num;
cancer_pattern_num = param.cancer_pattern_num;
R=healthy_pattern_num + (class_num-1)*cancer_pattern_num;
N2=size(test_data,2);

test_H = zeros(R,N2);

parfor i = 1:size(test_data,2)
%     i
    k = predicted_label(i);
    if k ==0  % no structural constraints
        mask =zeros(1,R);
    elseif k ==1  % predicted healthy controls
        mask = [zeros(1,healthy_pattern_num),ones(1,cancer_pattern_num*(class_num-1))];
    else  % predicted tumor types
        mask = [zeros(1,healthy_pattern_num),ones(1,cancer_pattern_num*(class_num-1))];
        mask((healthy_pattern_num+cancer_pattern_num*(k-2)+1):(healthy_pattern_num+cancer_pattern_num*(k-1))) = 0;
    end
    
    fun = @(x)norm(test_data(:,i)-train_W*x,'fro');
    x0 = rand([R,1]);
    A = [];
    b = [];
    Aeq = [ones(1,R);mask];
    beq = [1,0];
    lb = zeros(1,R);
    ub = ones(1,R);
    options = optimoptions('fmincon','Display','none');
    test_H(:,i) = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options);
end
end
