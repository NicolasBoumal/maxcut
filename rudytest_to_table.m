function rudytest_to_table()
% Transform data from rudytest.m in a Latex table for paper

data = load('rudytest.mat');
if isfield(data, 'record')
    data = data.record;
else
    data = data.data;
end

fid = fopen('rudytable.tex', 'w+');

for graphid = [1:67, 70, 72, 77, 81]
    
    graph = load(sprintf('Gset/g%d.mat', graphid));
    
    dat = squeeze(data(graphid, :, :));
    
    fprintf(fid, [...
		'Graph %d & Cut bound & %s & %s & %s & %s & %s \\\\\\nopagebreak\r\n'...
		'\\quad %d nodes & $\\lambdamin(S)$ & %s & %s & %s & %s & %s \\\\\\nopagebreak\r\n'...
		'\\quad %d edges & Time [s] & %s & %s & %s & %s & %s \\\\\r\n'...
		'\\hline\r\n\r\n'], ...
        graphid, metric1(dat(:, 1)), metric1(dat(:, 2)), metric1(dat(:, 3)), metric1(dat(:, 4)), metric1(dat(:, 5)), ...
        graph.n, metric2(dat(:, 1)), metric2(dat(:, 2)), metric2(dat(:, 3)), metric2(dat(:, 4)), metric2(dat(:, 5)), ...
        graph.m, metric3(dat(:, 1)), metric3(dat(:, 2)), metric3(dat(:, 3)), metric3(dat(:, 4)), metric3(dat(:, 5)) ...
        );
    
end

fclose(fid);


    % Cut bound
    function s = metric1(dat)
        if ~isnan(dat(1)) && ~isnan(dat(2))
            s = sprintf('%.1f', dat(1)-dat(2)*graph.n/4);
        else
            s = '-';
        end
    end

    % lambdamin(S)
    function s = metric2(dat)
        if ~isnan(dat(2))
            s = ['$' double2latex(dat(2), 0) '$'];
        else
            s = '-';
        end
    end

    % Computation time
    function s = metric3(dat)
        if ~isnan(dat(4))
            s = sprintf('%.1f', dat(4));
        else
            s = '-';
        end
    end

end
