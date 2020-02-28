function xhat_sm = H_Fuse(A, y, Hsm)
% reconstruction using smoothness regularization (usin pinv for final reconstruction)
        
        [~, N] = size(Hsm);
        Asm = [A;Hsm];
        ysm = [y; zeros(N-2+1,1)];
        xhat_sm = (pinv(Asm)*ysm).';
        
end

