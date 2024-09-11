% Creates string given an input number and a desired length, padding the
% front of the string with zeroes
%
% Nathan Laxague, 2015
%
function outstr = leading_zeroes(in_num,str_length)

% Generates string from input number
in_str = num2str(in_num);

% Check to see if the input number string is longer than the input length
if str_length > length(in_str)
    
    % Create an empty 'holder' string
    zero_str = [];
    
    % Compute the number of loop iterations
    str_length = str_length - length(in_str);

    % Loop and append the 'holder' string with zero strings
    for i = 1:str_length
    
        zero_str = [zero_str '0'];
    
    end
    
    % Concatenate the zero string block and the input string
    outstr = [zero_str in_str];

else
    
    outstr = in_str;
    
end

