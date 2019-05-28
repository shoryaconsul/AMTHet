% MixClone generates a different order of subclones. This function finds
% the order of the columns
function [ord_ind,L_shuffle] = L_perm(L,L_in)
    [m,n] = size(L);
    ind = perms(1:n);
    err = zeros(size(ind,1),1);
    for i=1:size(ind,1)
        err(i) = sum(sum(abs(L-L_in(:,ind(i,:)))));
    end
    
    [~,pos] = min(err);
    ord_ind = ind(pos,:);
    L_shuffle = L_in(:,ord_ind);
    
end