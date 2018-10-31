function yret = downsample(y,numpoints)
    if nargin<2
        numpoints = 100;
    end
	x = 1:length(y);
    xq = linspace(1,length(y),numpoints);
    yret = interp1(x,y,xq);        