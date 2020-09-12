%% L=0.5, W = 0.2

% XZ
centers = [ 0.7     1;  
            1.3     1; 
            ];
radia = [0.4 0.4];
figure;
subplot(131)
viscircles( centers, radia);axis equal; xlim([0.5 1.5]); ylim([0.5 1.5]); title('XZ'); set(gca, 'FontSize', 16)

% XY
centers = [ 0.7     1;  
            1.3     1; 
            1       0.7; 
            1       1.3;
%             0.788   0.788; 
%             1.212   0.788;
%             0.788   1.212; 
%             1.212   1.212;
            ];
radia = 0.4*ones(1, size(centers,1));
subplot(132)
viscircles( centers, radia); axis equal; xlim([0.7 1.3]); ylim([0.7 1.3]); title('XY'); set(gca, 'FontSize', 16)

subplot(133)
scatter( centers(:,1), centers(:,2), 'ro', 'LineWidth', 3);axis equal;xlim([0 2]); ylim([0 2]); title('Centers XY'); set(gca, 'FontSize', 16)
suptitle('L = 0.5, W = 0.2')
%% L=1, W = 0.4

% XZ
centers = [ 0.5     1;  
            1.5     1; 
            ];
radia = [0.7 0.7];
figure; subplot(131)

viscircles( centers, radia); axis equal;xlim([0.5 1.5]); ylim([0.5 1.5]); title('XZ'); set(gca, 'FontSize', 16)

% XY
centers = [ 0.5     1;  
            1.5     1; 
            1       0.5; 
            1       1.5;
            0.65    0.65; 
            0.65    1.35;
            1.35    0.65; 
            1.35    1.35;
            ];
radia = 0.7*ones(1, size(centers,1));
subplot(132)
viscircles( centers, radia); axis equal;xlim([0.7 1.3]); ylim([0.7 1.3]); title('XY'); set(gca, 'FontSize', 16)


subplot(133)
scatter( centers(:,1), centers(:,2), 'ro', 'LineWidth', 3);axis equal;xlim([0 2]); ylim([0 2]); title('Centers XY'); set(gca, 'FontSize', 16)
suptitle('L = 1.0, W = 0.4')

%% L=2, W = 0.8

% XZ
centers = [ 0     1;  
            2     1; 
            ];
radia = [1.4 1.4];
figure;subplot(131)

viscircles( centers, radia); axis equal;xlim([0 2]); ylim([0 2]); title('XZ'); set(gca, 'FontSize', 16)

% XY
centers = [ 0     1;  
            2     1; 
            1     0; 
            1     2;
            0.29    0.29; 
            0.29    1.71;
            1.71    0.29; 
            1.71    1.71;
            ];
radia = 1.4*ones(1, size(centers,1));
subplot(132)
viscircles( centers, radia); axis equal;xlim([0.5 1.5]); ylim([0.5 1.5]); title('XY'); set(gca, 'FontSize', 16)


subplot(133)
scatter( centers(:,1), centers(:,2), 'ro', 'LineWidth', 3);axis equal;xlim([0 2]); ylim([0 2]); title('Centers XY'); set(gca, 'FontSize', 16)
suptitle('L = 2.0, W = 0.8')
