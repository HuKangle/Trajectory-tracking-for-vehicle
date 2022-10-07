function F_s = rpi_nilpotent(A, W)
    % by Lukas Brunke, lukas.brunke@mail.utoronto.ca
    % Compute the mRPI set for the case that A is nilpotent
    eps = 1e-8;
    state_dim = size(A, 1);
    
    % Initialiize F_s = {0} and A^0 = I, at s = 0
    F_s = Polyhedron(zeros(2, state_dim));
    A_k = eye(state_dim);
    
    % state_dim is the maximal nilpotent degree
    for i = 1 : state_dim
        % Stopping when the nilpotent degree is reached
        if all(abs(A_k(:)) < eps)
            break;
        end
        
        % Compute the Minkowski sum for F_s from eq. (2) from Rakovic et al. (2005).
        F_s = F_s + A_k * W;
        F_s.minHRep(); 
        
        % Compute next power of A
        A_k = A_k * A;
    end
end
