% Load the data
load('data_tsb.mat');

% Reshape data into meshgrid standard
X_raw = axeX_Z0';
Y_raw = flipud(axeY_Z0-axeY_Z0(1,1))';
U_raw = permute(U_Z0_VER_Run1, [2,1,3]);
V_raw = permute(V_Z0_VER_Run1, [2,1,3]);

% Undersample every 3 points
X = X_raw(1:3:end,1:3:end);
Y = Y_raw(1:3:end,1:3:end);
U = U_raw(1:3:end,1:3:end,:);
V = V_raw(1:3:end,1:3:end,:);

% Reshape the data into 2D matrices
[nx, ny, nt] = size(U);
U_2d = reshape(U, nx*ny, nt)';
V_2d = reshape(V, nx*ny, nt)';

% Subtract the mean from the data
meanU = mean(U_2d);
U_2d = U_2d - meanU;
meanV = mean(V_2d);
V_2d = V_2d - meanV;

% Compute the covariance matrix and its eigenvalues/eigenvectors
C = U_2d * V_2d';
[U_s, S, V_s] = svd(C);
lambda = diag(S).^2;

% Compute the POD modes
Phi = V_2d' * U_s * diag(1./sqrt(lambda));

% % Plot the first POD mode
% mode1 = reshape(Phi(:,1), nx, ny);
% contourf(X, Y, mode1');
% colorbar;
% title('POD mode 1');
% xlabel('x');
% ylabel('y');

% Reshape the first POD mode into a 2D matrix
mode1_2d = reshape(Phi(:,1), nx, ny);

% Plot the first POD mode
figure;
contourf(X, Y, mode1_2d);
colorbar;
title('POD mode 1');
xlabel('x');
ylabel('y');
