%% Function to find mu from the read counts and computed copy number aberrations (CNAs(
% r: vector of read counts (normalized by coverage)
% L: Matrix with each row corresponding to CNAs for associated region
% res: Determined mu

function res = proj_simplex(r,L)
    y = pinv(L)*r;
    u = sort(y,'descend');
    ulen = length(u);
    rho_ser = (1-cumsum(u))./((1:ulen)');
    rho = sum((u + rho_ser)>0);
    lambda = rho_ser(rho);
    res = max(y+lambda,0);
end