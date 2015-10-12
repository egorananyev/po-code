function [ imageWindow ] = renderWindow ( m, rMulti )
% Renders a nuttallwin / flattopwin window based on matrix dimensions of m.
% rMulti = radius multiplier

[r, c] = size(m);
wc = window(@nuttallwin,c).^.4; % TODO: nuttallwin or flattopwin?
wr = window(@nuttallwin,r).^.4;
[maskr,maskc] = meshgrid(wr,wc);
imageWindow = maskr.*maskc;
imageWindow(imageWindow<0)=0; % not really needed with nuttallwin

% Applying additional, "hard" circular window on top of the "smooth" one above.
n2 = floor(r/2); % assuming that r (rows) = c, i.e., image is square
[x,y] = meshgrid(-n2:n2); % a series of x and y coordinates
M = (sqrt(x.^2 + y.^2)) < (rMulti*n2); % 1 if true that the coordinates are inside
    % ... the specified radius.
imageWindow = imageWindow .* M;

end