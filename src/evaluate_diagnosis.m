function average_ACC = evaluate_diagnosis(diagosis_results, param)
disp(['confusion_matrix=']);
confusion_matrix = compute_confusion_matrix(diagosis_results,param.test_sample_num)
ACC = diag(confusion_matrix)'./param.test_sample_num;
average_ACC = mean(ACC);
end

%% compute confusion matrix
 function confusion_matrix=compute_confusion_matrix(predict_label,num_in_class)
num_class=length(num_in_class);
num_in_class=[0 num_in_class];

confusion_matrix=size(num_class,num_class);
 
for ci=1:num_class
    for cj=1:num_class
        summer=0;
        c_start=sum(num_in_class(1:ci))+1;
        c_end=sum(num_in_class(1:ci+1));
        summer=size(find(predict_label(c_start:c_end)==cj),2);  
        %confusion_matrix(ci,cj)=summer/num_in_class(ci+1);
        confusion_matrix(ci,cj)=summer;
    end
end
 
end