function convertedVelocity = convertVelocity(h,velocity)
% this function turns the displayed function of the edit text box in a
% handle function (string format) that can be used in the code

if isempty(velocity)
    convertedVelocity = velocity;
else
    % (x,y) --> (r,z)
        velocity = strrep(velocity,'x','r');

    % change 'ln(' into 'log('
        velocity = strrep(velocity,'ln(','log(');

    % insert handles.r_i(r_o) instead of r_i(r_o)
        velocity = strrep(velocity,'r_i',num2str(h.r_i));
        velocity = strrep(velocity,'r_o',num2str(h.r_o));

    % to avoid errors: in front of every '^' put a '.'
        velocity = strrep(velocity,'^','.^');

    % coordinate change: positive z-axis of GUI points in opposite direction of positive z-axis of the code
        velocity = ['@(r)-1*(' velocity ')'];
        convertedVelocity = velocity;
end

end