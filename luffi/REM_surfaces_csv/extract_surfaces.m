mo = matfile('BEST_fitting_models.mat');

v = who(mo);

for i=1:length(v)
    %eval(v{i})
    ax = plot(eval(v{i}));
    
    %X = ax.XData(1,:)
    %Y = ax.XData(:,1)
    %Z = ax.ZData
    csvwrite([v{i},'.csv'], [ax.XData(:), ax.YData(:), ax.ZData(:)]);
end

