function y = isnilpotent(A)
    % by Lukas Brunke, lukas.brunke@mail.utoronto.ca
    y = abs(det(eye(size(A)) + A) - 1) <= 1e-8;
end

