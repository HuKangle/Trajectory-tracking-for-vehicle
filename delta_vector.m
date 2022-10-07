function deltas = delta_vector(k)
    values = [-1, 1];
    
    n = numel(values);
    combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
    combs = reshape(values(combs),[],k);
    
    deltas = [];
    for i = 1 : size(combs, 1)
        deltas = [deltas; unique(perms(combs(i, :)), 'rows')];
    end
end