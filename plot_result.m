function [] = plot_result(res, theory, d, n, k, gamma, xi, num_rep, outpath, ylab, name, fontsize)
    savefigs = 1; 
    figure, hold on
    me = mean(res, 1);
    re = std(res, 1);
    h1 = errorbar(d, me, re,'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
    h2 = plot(d, theory, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
    if k > 1
        xlabel('$$d_i$$', 'Interpreter', 'LaTex');
    else
        xlabel('$$d_1$$', 'Interpreter', 'LaTex');
    end
    if ylab == 'cos' & k > 1
        ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
    elseif  ylab == 'cos' & k == 1
        ylabel('$$|\langle u_1,\tilde\xi_1\rangle|^2$$', 'Interpreter', 'LaTeX');
    elseif ylab == 'lbd' & k > 1
        ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
    else
        ylabel('$$\tilde\lambda_1$$', 'Interpreter', 'LaTeX');
    end
    if ~exist('fontsize', 'var')
        fontsize = 30;
    end
    set(gca,'FontSize', fontsize)
    if exist('name', 'var')
        title(name, 'FontSize', fontsize, 'Interpreter', 'LaTex');
    end
    legend('location', 'best');
    grid on;
    if savefigs==1
        filename = append(outpath, sprintf('gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.pdf', gamma, xi, num_rep, n, k));
        saveTightFigure(gcf, filename);
        fprintf(['Saved Results to ' filename '\n']);
        close(gcf)
    end
end