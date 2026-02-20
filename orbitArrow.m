function h = orbitArrow(x, y, method, pos, lArr, col)

if method == "ind"
    ind1 = pos;
elseif method == "rel_ind"
    ind1 = round(pos*length(x));
elseif method == "horiz"
    ind1 = find((x-x(1))./(x(end)-x(1)) > pos, 1, 'first')-1;
elseif method == "vert"
    ind1 = find((y-y(1))./(y(end)-y(1)) > pos, 1, 'first')-1;
else
    ind1 = [];
end

if ~isempty(ind1)
    % Find last point in the series that is less than length from start point
    dist2 = (x-x(ind1)).^2 + (y-y(ind1)).^2;
    ind2 = find(dist2 < lArr^2, 1, 'last');

    % If this returns the same index as the starting point, just use the next
    % point in the series
    if ind2 == ind1
        ind2 = ind1+1;
    end

    % Set arrow's start and end points
    x_start = x(ind1);
    x_end = x(ind2);
    y_start = y(ind1);
    y_end = y(ind2);
    
    h = annotation('arrow');
    h.Parent = gca;
    
    h.X = [x_start, x_end]; 
    h.Y = [y_start, y_end];
    h.Color = col;
else
    fprintf('   Warning: Invalid method used in call to orbitArrow, exiting without doing anything\n')
end