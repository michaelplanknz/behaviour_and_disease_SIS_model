function h = orbitArrow(x, y, method, pos, col)

if method == 'ind'
    ind = pos;
elseif method == 'rel_ind'
    ind = round(pos*length(x));
elseif method == 'horiz'
    ind = find((x-x(1))./(x(end)-x(1)) > pos, 1, 'first')-1;
elseif method == 'vert'
    ind = find((y-y(1))./(y(end)-y(1)) > pos, 1, 'first')-1;
else
    ind = [];
end

if ~isempty(ind)
    % Set arrow's start and end points
    x_start = x(ind);
    x_end = x(ind+1);
    y_start = y(ind);
    y_end = y(ind+1);
    
    h = annotation('arrow');
    h.Parent = gca;
    
    h.X = [x_start, x_end]; 
    h.Y = [y_start, y_end];
    h.Color = col;
else
    fprintf('   Warning: Invalid method used in call to orbitArrow, exiting without doing anything\n')
end