clear;clc
load("H_yNoise_signa2_iter_Nt");
% Convert complex H and y_noise to real-valued representations
Hr = [real(H), -imag(H); imag(H), real(H)];
yr = [real(y_noise); imag(y_noise)];

%Hr = [Hr,sum(Hr,2)]; %Add ABFT checksums
% Initialize variables
A = Hr' * Hr;
A = A + sigma2 * eye(2 * Nt ); %Regularize A for numerical stability by adding a small value proportional to 𝜎2
b = Hr' * yr;
Ainv = diag(1./diag(A));  %%% Initialize A0 inverse
E = 2 *eye(2 * Nt);


A = [A;sum(A,1)];
Ainv = [Ainv,sum(Ainv,2)];
E = [E,sum(E,2);sum(E,1),sum(E(:))];

% Newton Iteration
for i = 1:iter
    t1 = (E - A * Ainv);
    Ainv(:,end)=0;
    Ainv = Ainv * t1;
end

% Compute the estimated x_hat (real-valued)
x_hat_real = Ainv(:,1:end-1) * b;

% Convert x_hat back to complex form
x_hat = x_hat_real(1:Nt) + 1i * x_hat_real(Nt+1:end)

%x_hat =
%   5.6885 - 4.2631i
%   5.4215 + 1.7480i
%   0.4144 - 5.3477i
%   3.7611 + 5.2414i
%   5.5210 + 1.1007i
%   6.3614 - 2.5327i
%  -4.0187 + 3.4213i
%  -2.3859 + 5.3711i

