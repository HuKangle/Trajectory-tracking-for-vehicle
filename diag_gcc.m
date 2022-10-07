function [K, P, e] = diag_gcc(A, B, E, Cy, Dy, Q, R, N)
    % gcc Generates the infinite horizon Guaranteed Cost Control for a
    %     parametric uncertain linear system
    %    Inputs: A  - System matrix
    %            B  - Input matrix
    %            E  - Disturbance matrix (LHS)
    %            Cy - System disturbance matrix (RHS)
    %            Dy - Input disturbance matrix (RHS)
    %            Q  - State cost matrix
    %            R  - Input cost matrix
    %            N  - State/Input cross cost matrix
    %
    %    Outputs: K - GCC gain matrix
    %             P - GCC cost matrix
    %             e - optimal epsilon
    %
    %    Author: Carlos M. Massera
    %    Instituition: University of São Paulo

    % Get system size
    Nx = size(A,2);
    Nu = size(B,2);
    Nw = size(E,2);
    E  = 1000 * E;
    Cy = Cy / 1000;
    Dy = Dy / 1000;
    if nargin <= 7
        N = zeros(Nx, Nu);
    end

    % Calculate cost decomposition
    W = [Q N; N' R];
    [~, s, vt] = svd(W);
    mask = (diag(s) >= 1e-7);

    Whalf = sqrt(s(mask, mask)) * vt(:, mask)';
    Cz = Whalf(:,1:Nx);
    Dz = Whalf(:,Nx+1:end);
    
    % Get cost size
    Nz = size(Whalf,1);
    
    % Generate Quadratic Cost LMI
    Pinv = sdpvar(Nx, Nx);
    KPinv = sdpvar(Nu, Nx, 'full');
    S = sdpvar(Nx, Nx);

    M = blkvar;
    M22 = - Pinv;
    M(1,1) = - eye(Nz);
    M(1,3) = Cz * Pinv - Dz * KPinv;
%     M(2,2) = - Pinv;
    M(2,3) = A * Pinv - B * KPinv;
    M(3,3) = - Pinv;

    N = blkvar;
    N(1,1) = - Pinv;
    N(1,2) = eye(Nx);
    N(2,2) = - S;
    
    % Add robustness to M
    e = sdpvar(Nw, 1, 'full');
    
    for i = 1:Nw
        Ei = E(:,i);
        Cyi = Cy(i,:);
        Dyi = Dy(i,:);
        
        idx = 3 + i;
        M22 = M22 + e(i) * (Ei * Ei');
        M(3, idx) = (Cyi * Pinv - Dyi * KPinv)';
        M(idx, idx) = - e(i);
    end

    M(2,2) = M22;

    constraints = [Pinv >= 0;
                   M <= 0;
                   N <= 0];

    opt = sdpsettings('solver', '+sedumi', 'verbose', 0);
    sol = optimize(constraints, trace(S), opt);

    if sol.problem ~= 0
        warning('gcc:solver_failed','Solver did not converge');
    end

    K = value(KPinv) / value(Pinv);
    P = inv(value(Pinv));
    e = value(e);

end