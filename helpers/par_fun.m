function results = par_fun(fun,range)

tic;

n = length(range);

prob = 1/ceil(n/75);

reps = round(n*prob);
fprintf('%s\n\n',repmat('-',1,reps)) %make the "background"

parfor i = 1:n   
    results(i,:) = fun(range(i));
    if (rand < prob)
        fprintf(1,'\b*\n'); % \b is backspace
    end
end

toc;
