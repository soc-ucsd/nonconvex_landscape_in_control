clearvars;

A = -1;
B2 = 1;
C2 = 1;

[nx, nu] = size(B2);
[ny, ~] = size(C2);

nw = nx + ny;
nz = nx + nu;

B1 = [eye(nx) zeros(nx, ny)];
C1 = [eye(nx); zeros(nu, nx)];
D11 = zeros(nz, nw);
D12 = [zeros(nx, nu); eye(nu)];
D21 = [zeros(ny, nx) eye(ny)];

sys_pl = ss(A, [B1 B2], [C1; C2], [D11 D12; D21 zeros(ny, nu)]);
% [K_opt, sys_cl_opt, gamma_opt] = hinfsyn(sys_pl, ny, nu);

n_ax = 101;
n_ay = 101;
n_az = 61;

Ak_vals = linspace(-2,2,n_ax);
Bk_vals = linspace(-4,4,n_ay);
Dk_vals = linspace(-1.5,1.5,n_az);

[Ak_grid,Bk_grid,Dk_grid] = meshgrid(Ak_vals, Bk_vals, Dk_vals);

P12vals = zeros(n_ax, n_ay, n_az);
Jvals = zeros(n_ax, n_ay, n_az);

eps = 1e-9;

for kk = 1:n_az
    for ii = 1:n_ax
        for jj = 1:n_ay
            AK = Ak_grid(ii, jj, kk);
            BK = Bk_grid(ii, jj, kk);
            CK = 1;
            DK = Dk_grid(ii, jj, kk);

            Acl = [A + B2*DK*C2, B2*CK;
                BK*C2, AK];
            Bcl = [B1 + B2*DK*D21; BK*D21];
            Ccl = [C1 + D12*DK*C2, D12*CK];
            Dcl = D11 + D12*DK*D21;
            
            sys_cl = ss(Acl, Bcl, Ccl, Dcl);
            hinf_norm = hinfnorm(sys_cl, eps);
            
            Jvals(ii, jj, kk) = hinf_norm;
            
            if hinf_norm == Inf || hinf_norm > 1e10
                P12vals(ii, jj, kk) = Inf;
                Jvals(ii, jj, kk) = Inf;
                continue;
            end
            
            % Find P by Riccati-based approach
            P = findP_Riccati(Acl, Bcl, Ccl ,Dcl, hinf_norm/(1-eps));
            minP = min(eig(P));
            
            if isempty(P) || minP <= 0
                % Find P by LMI-based approach
                [P, tmin] = findP_LMI(Acl, Bcl, Ccl, Dcl, hinf_norm);
                minP = min(eig(P));
            end
            
            test_M = [Acl'*P + P*Acl, P*Bcl, Ccl';
                Bcl'*P, -hinf_norm*eye(nw), Dcl';
                Ccl, Dcl, -hinf_norm*eye(nz)];
            
            maxM = max(eig(test_M));
            
            if maxM >= 1e-6
                % Flag there exists no P>=0 satisfying the LMI
                P12vals(ii, jj, kk) = NaN;
            end
            
            if minP >= 1e-4
                P12vals(ii, jj, kk) = det(P(1:nx,nx+1:end));
            else
                % Flag P may not be positive definite
                P12vals(ii, jj, kk) = NaN;
            end
        end
    end
end

% save data.mat;