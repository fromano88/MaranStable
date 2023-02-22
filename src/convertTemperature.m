function convertedTemperature = convertTemperature(h, temperature, direction)
% this function turns the displayed function of the edit text box in a
% handle function (string format) that can be used in the code

if isempty(temperature)
    convertedTemperature = temperature;
else
    % (x,y) --> (r,z)
        temperature = strrep(temperature,'x','r');
        temperature = strrep(temperature,'y','z');

    % coordinate transformation
        if length(h.blocks) == 2
            temperature = strrep(temperature,'z','(d_1+d/2-z)');
        else
            temperature = strrep(temperature,'z','(d/2-z)');
        end
        
    % insert numbers instead of geometry parameters
        temperature = strrep(temperature,'r_c',num2str(h.r_c));
        temperature = strrep(temperature,'r_i',num2str(h.r_i));
        if isfield(h,'r_o')
            temperature = strrep(temperature,'r_o',num2str(h.r_o));
            temperature = strrep(temperature,'d_1',num2str(h.l_d1));
            temperature = strrep(temperature,'d_2',num2str(h.l_d2));
        end
        temperature = strrep(temperature,'d',num2str(h.l_lb));

    % to avoid errors: in front of every '^' put a '.'
        temperature = strrep(temperature,'^','.^');
        
    % update temperature profile
        convertedTemperature = temperature;
        if strcmp(direction,'r')
            convertedTemperature = ['@(r)' convertedTemperature];
        else
            convertedTemperature = ['@(z)' convertedTemperature];
        end
end
        
end