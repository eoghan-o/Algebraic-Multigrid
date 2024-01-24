starting_size = 4;

L = eye(starting_size);

max_level = 5;
coarse_size = 3;

L_i = L;

L_list = cell(1, max_level);
P_list = cell(1, max_level);

n = starting_size;

for i = 1:(max_level - 1)
    if(size(L_i) <= coarse_size)
        [aggregate1, aggregate2] = form_aggregate(L_i, n);
        P_l = build_prolongation(L_i, aggregate1, aggregate2);
        P_list{i} = P_l;
        L_i = P_l'*L_i*P_l;
        L_list{i} = L_i;
    end
end


function P = build_prolongation(L_i, aggregate1, aggregate2)
    n = size(L_i);
    P = zeros(n,n-1);
    for i = 1:n
        for j = 1:n
            if i==aggregate2 && j== aggregate1
                P(i,j) = 1;
            elseif i < aggregate2 && i == j  
                P(i,j) = 1;
            elseif i > aggregate2 && i == j + 1
                P(i,j) = 1;
            end
        end
    end
end

function [aggregate1, aggregate2] = form_aggregate(L_i, n)
    max_value = max(L_i, [], 'all');
    for i = 1:n
        for j = 1:n
            if(L_i(i,j) == max_value)
                aggregate1 = j;
                aggregate2 = i;
                return
            end
        end
    end
end

function x_l = v_cycle(L_l, x_l, b_l, l)
    if l == max_level
        x_l = pinv(L_l)*b_l;
        return
    else
        x_l = pre_smoothing(L_l, x_l, b_l, l);
        r_l = b_l - L_l * x_l;
        r_l = P_list{l}'*r_l;
        e_l = v_cycle(L_list{l+1}, 0, r_l, l + 1);
        x_l = x_l + P_list{l}*e_l;
        x_l = post_smoothing(L_l, x_l, b_l, l);
    end
    return
end