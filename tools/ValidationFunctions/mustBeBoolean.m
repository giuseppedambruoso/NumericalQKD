function mustBeBoolean(x)
% mustBeBoolean Validation function to ensure input is logical (boolean)
%
% Throws an error if the input is not a logical scalar or array.

if ~islogical(x)
    error("mustBeBoolean:NotLogical", ...
        "Value must be of type logical (boolean).");
end
end