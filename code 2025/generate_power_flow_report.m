function generate_power_flow_report(ppf, x_opt)
    % Create timestamp for filename
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('power_flow_results_%s.csv', timestamp);
    
    % Get power flow results
    power_flow = evaluate(ppf, x_opt);
    
    % Round to 4 decimal places
    power_flow = round(power_flow, 4);
    
    % Replace values outside [-0.1500, 0.1500] with 'X'
    power_flow_modified = power_flow;
    power_flow_modified(power_flow > 0.1500 | power_flow < -0.1500) = 'X';
    
    % Save to CSV file
    csvwrite(filename, power_flow_modified);
    
    fprintf('Power flow results saved to: %s\n', filename);
end 