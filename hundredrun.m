Nruns  = 100;
outDir = fullfile(pwd, 'logs_tv_varT');
if ~exist(outDir,'dir'), mkdir(outDir); end

summary = table('Size',[Nruns 7], ...
    'VariableTypes', {'double','double','double','double','double','double','double'}, ...
    'VariableNames', {'run','seed','k_end','t_final','min_h','goal_dist','runtime_s'});

for r = 1:Nruns
    seed = r;  % or any list you want
    fprintf('Run %d/%d (seed=%d)...\n', r, Nruns, seed);

    tRun = tic;

    out = tv_cbf_varT_run(seed, struct('doViz',false,'doPlots',false));

    runtime_s = toc(tRun);

    % Save ONLY logs (MAT is best for big arrays)
    out.runtime_s = runtime_s;
    save(fullfile(outDir, sprintf('run_%03d.mat', r)), '-struct', 'out', '-v7.3');

    % Summary row for quick paper tables/plots
    t_final   = out.thist(out.k_end+1);
    min_h     = out.h_min;
    goal_dist = norm(out.X(1:2, out.k_end+1) - out.goal(:));

    summary{r,:} = [r, seed, out.k_end, t_final, min_h, goal_dist, runtime_s];
end

writetable(summary, fullfile(outDir, 'summary.csv'));
save(fullfile(outDir, 'summary.mat'), 'summary');
