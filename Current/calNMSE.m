function nmse=calNMSE(orgSig,recSig,varargin)

if isempty(varargin)
    boun = 0;
else boun = varargin{1};
end

if size(orgSig,2)==1       % if signal is 1-D
    orgSig = orgSig(boun+1:end-boun,:);
    recSig = recSig(boun+1:end-boun,:);
else                       % if signal is 2-D or 3-D
    orgSig = orgSig(boun+1:end-boun,boun+1:end-boun,:);
    recSig = recSig(boun+1:end-boun,boun+1:end-boun,:);
end

mse=norm(orgSig(:)-recSig(:),2)^2/length(orgSig(:));
sigEner=norm(orgSig(:))^2;
nmse=(mse/sigEner);