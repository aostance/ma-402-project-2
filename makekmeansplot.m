function f = makekmeansplot(X,k,visibility,C0)
    % INPUTS: X, the data
    %         k, the number of clusters
    % OPTIONAL: visibility, to display or hide figures
    %           C0, a starting guess for the means, size (k x d)
    % OUTPUTS: A figure with the k-means clustering
    % In order to produce a single plot, 
    % just input X and k and leave out 'visibility'
    
    if nargin < 3, visibility = 'on'; end
    
    f = figure('Visible',visibility);
    
    % Run the k-means algorithm
    if nargin == 4
        [idx,C] = kmeans(X,k,'start',C0);
    else
        [idx,C] = kmeans(X,k); 
    end
    
    % Define a grid on the plot
    x1 = min(X(:,1)):0.01:max(X(:,1));
    x2 = min(X(:,2)):0.01:max(X(:,2));
    [x1G,x2G] = meshgrid(x1,x2);
    XGrid = [x1G(:),x2G(:)]; 

    idx2Region = kmeans(XGrid,k,'MaxIter',1,'Start',C);
    
    % Plot the Voronoi cells
    COL = [.75,0.85,0.85;0.85,.75,0.85;0.85,0.85,.75]; 
    if k > 3
        COL = [COL; [.5,.5,.5] + .5*rand(k-3,3)]; 
    end      
   
    gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
        COL,'..');
    hold on;
    axis equal

    gscatter(X(:,1),X(:,2),idx);
    p = plot(C(:,1),C(:,2),'k*','MarkerSize',15,'LineWidth',2);

    lgd = legend(p,'Means'); 
    lgd.Location = 'northwest'; 
    lgd.FontSize = 14; 
    xlim([min(x1),max(x1)]); 
    ylim([min(x2),max(x2)]); 

end
