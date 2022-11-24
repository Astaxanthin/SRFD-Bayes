function [W,R,err]=SRFD(X, param)

K=size(X,1);
N=size(X,2);

class_num=param.class_num;
healthy_pattern_num = param.healthy_pattern_num;
cancer_pattern_num = param.cancer_pattern_num;
C=param.healthy_pattern_num + (class_num-1)*cancer_pattern_num; %
train_sample_num = param.train_sample_num;

%% Initialization

W=rand([K,C]); %
R=rand([C,N]);
norm_coff=sum(R,1);
R=R./repmat(norm_coff,size(R,1),1);%normalization

%% structural mask
healthy_top_mask = zeros(healthy_pattern_num, N);
healthy_left_mask = ones(C-healthy_pattern_num, train_sample_num(1));
cancer_mask = [];
for i = 1:size(train_sample_num,2)-1
    cancer_mask =  blkdiag(cancer_mask, ones(cancer_pattern_num, train_sample_num(i+1)));
end
cancer_mask = 1-cancer_mask;
train_mask = [healthy_left_mask,cancer_mask];
mask = [healthy_top_mask;train_mask];

%% parameters
p=param.p; %L2-p
eta=param.eta;

%% objective function
Omega = 1/2*eta*norm(mask.*R,'fro')^2;
% calculate_2p_norm(X-W*H,p) 
err(1)=1/2*calculate_2p_norm(X-W*R,p)  + Omega;

iter=1;
while iter<param.maximum_iterations 
    iter;
    Z=X-W*R;
    L2_p= p./(2*sum(Z.^2,1).^(1-p/2));
    D=diag(L2_p);
    
    % update W
    factorW_numerator=X*D*R';
    factorW_denominator=W*(R*D*R');
    W=W.*(factorW_numerator./factorW_denominator);
    W(W>1) = 1.0-1e-4;
    
    % update H
    factorH_numerator=W'*X*D ;
    factorH_denominator=W'*W*R*D+eta*mask.*R;
    R=R.*(factorH_numerator./factorH_denominator);
    norm_coff=sum(R,1);
    R=R./repmat(norm_coff,size(R,1),1);%normalization
    
    iter=iter+1;
    Omega = 1/2*eta*norm(mask.*R,'fro')^2;
    err(iter) = 1/2*calculate_2p_norm(X-W*R,p)  + Omega;
%     err(iter)
    
    if abs(err(iter-1)-err(iter))<param.convergence_threshold
        break;
    end
end

end

function x_norm = calculate_2p_norm(X,p)
 x_norm = sum(sum(X.^2,1).^(p/2));
end