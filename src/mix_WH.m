function extend_data = mix_WH(train_W, extend_H, extend_labels,reconstruction_err,param)
extend_WH = train_W*extend_H;
extend_id = unique(extend_labels);

extend_sample_num = [];
for i = 1:size(extend_id,2)
    extend_sample_num = [extend_sample_num,sum(extend_labels==extend_id(i))];
end
extend_sample_num = [0,extend_sample_num];
train_sample_num = [0,param.train_sample_num];

extend_data = [];
for i = 1:size(extend_id,2)
    extend_class_WH = extend_WH(:,sum(extend_sample_num(1:i))+1:sum(extend_sample_num(1:i+1)));
    rec_class_err = reconstruction_err(:,sum(train_sample_num(1:i))+1:sum(train_sample_num(1:i+1)));
    rec_class_err = rec_class_err(:);
    [mu,sigma] = normfit(rec_class_err);
    
    extend_err = zeros(size(extend_class_WH,1),size(extend_class_WH,2));
    for j = 1:size(extend_class_WH,2)
        extend_err(:,j) = normrnd(mu,sigma,[size(extend_class_WH,1),1]);
    end
    extend_WH_err = extend_class_WH + extend_err;
    extend_WH_err(extend_WH_err>1) = 1;
    extend_WH_err(extend_WH_err<0) = 0;
    extend_data = [extend_data, extend_WH_err];
    
end

end