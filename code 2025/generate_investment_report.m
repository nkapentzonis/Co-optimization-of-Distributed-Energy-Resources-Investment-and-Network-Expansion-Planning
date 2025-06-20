% Report new transmission lines
fprintf(fid, '\nNew Transmission Lines:\n');
selected_lines = find(uline > 0.9); % Only consider lines with values very close to 1
num_new_lines = length(selected_lines);
fprintf(fid, 'Number of new lines: %d\n', num_new_lines);
if num_new_lines > 0
    fprintf(fid, 'Lines: %s\n', num2str(selected_lines'));
end

% Investment Summary
fprintf(fid, '\nInvestment Summary:\n');
fprintf(fid, 'Total DER Investment: %.2f M€\n', total_der_investment/1e6);
fprintf(fid, 'Total Line Investment: %.2f M€\n', total_line_investment/1e6);
fprintf(fid, 'Total Investment: %.2f M€\n', (total_der_investment + total_line_investment)/1e6);
fprintf(fid, 'Annualized DER Investment: %.2f M€/year\n', annualized_der_investment/1e6);
fprintf(fid, 'Annualized Line Investment: %.2f M€/year\n', annualized_line_investment/1e6);
fprintf(fid, 'Total Annualized Investment: %.2f M€/year\n', (annualized_der_investment + annualized_line_investment)/1e6);

% Add new metrics
fprintf(fid, '\nInvestment Ratios:\n');
fprintf(fid, 'Lines to DER Ratio (LtDR): %.2f%%\n', (total_line_investment/total_der_investment)*100);
fprintf(fid, 'Lines to Total Ratio (LtTR): %.2f%%\n', (total_line_investment/(total_der_investment + total_line_investment))*100);

fclose(fid);
fprintf('Investment report generated: %s\n', filename); 