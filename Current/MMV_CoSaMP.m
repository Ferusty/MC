function [S]=MMV_CoSaMP(A,B,K)

%Subspace-Pursuit 
%Parameters : We are given the matrix A, the vector B, and the sparsity level K.

N=size(A,2);
L=size(B,2);

done = 0;
corr=zeros(1,N);

%Initialization
for j=1:N
    corr(j) = norm(A(:,j)'*B,2)/norm(A(:,j),2);
end
[sortedValues,sortIndex] = sort(abs(corr(:)),'descend');
T_cap = sortIndex(1:K);

x2=zeros(N,L);
x2(T_cap,:)=A(:,T_cap)\B;
res=B-A*x2;
   
if norm(res,'fro')==0
    done=1;
end

%Iteration
while ~done
   
    for j=1:N
        corr(j) = norm(A(:,j)'*res,2)/norm(A(:,j),2);
    end
    [sortedValues,sortIndex] = sort(abs(corr(:)),'descend');
    
    T_dash=union(T_cap,sortIndex(1:2*K));
    
    xp_dash=zeros(N,L);
    xp_dash(T_dash,:)=A(:,T_dash)\B;
    
    [sortedValues1,sortIndex1] = sort(sqrt(sum(abs(xp_dash).^2,2)),'descend');
    T_tilda=sortIndex1(1:K);
    
    x2 = zeros(N,L);
    x2(T_tilda,:)=A(:,T_tilda)\B;
    
    res_tilda=B-A*x2;
    
    if (norm(res_tilda,'fro') < norm(res,'fro')) 
        T_cap=T_tilda;
        res=res_tilda;
    else
        done=1;
    end

end

% x=zeros(N,L);
% x(T_cap,:)=pinv(A(:,T_cap))*B;
%T_cap=T_cap';
S=pinv(A(:,T_cap))*B;