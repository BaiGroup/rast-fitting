% process_variable_arguments Utility function for processing variable length input argument
% list
%
% Syntax
%   arguments = process_variable_arguments(names, default_values, varArg)
%   [arguments, varargout] = process_variable_arguments(names, default_values, varArg)
%
% Description
%   arguments = process_variable_arguments(names, default_values, var_args)
%   generates a commaseparated list and call struct, arguments, with k attribute
%   names defined in names (k by 1 matrix) and k default values defined in
%   default_values (k by 1 matrix). The property/value pairs that are explictly
%   given in var_args overwrite those in default_values. var_args is usually the
%   variable length input argument list, varargin, to the function calling
%   process_variable_arguments.
%
%   [arguments, varargout] = process_variable_arguments(names, default_values, var_args)
%   stores parameters from var_args that are not found in names in varargout.
%   Currently, vargargout support only one variable that will store all the
%   parameters not in names.
%
% See also varargin varargout

function [arguments, varargout] = process_variable_arguments(names, default_values, var_args)

structInput = cell(2, length(names));
structInput(1, :) = names';
structInput(2, :) = default_values';
arguments = struct(structInput{:});

nVArgOut = max(nargout, 1) - 1;
if nVArgOut > 1
    error('process_variable_arguments:TooManyOutputVar','Too many output variables. Currently support at most two output variables.');
end
iVarArgOut = 0;

nVarArg = size(var_args, 2);
iVarArg = 1;
flag = 2;

while iVarArg <= nVarArg
    arg = var_args{iVarArg};
    if flag == 1
        if ischar(arg)
            arg=lower(deblank(arg));
        end
        arguments.(names{iName,:})=arg;
        flag=2;
    elseif flag==2
        if ~ischar(arg)
           error('process_variable_arguments:ParamNotString','Parameter name must be a string.',iVarArg);
        end
        iName=find(strncmpi(arg,names,length(arg)));
        if isempty(iName)
            flag=0;
            continue
        elseif length(iName)>1
            jName=find(strcmpi(arg,names));
            if length(jName)==1
                iName=jName;
            else
                allMatchedNames=['(',names{iName(1),:}];
                for jName=iName(2:length(iName))'
                    allMatchedNames=[allMatchedNames,', ',names{jName,:}];
                end
                allMatchedNames=[allMatchedNames,')'];
                error('process_variable_arguments:AmbiguousParamName','Ambiguous parameter name: ',arg,allMatchedNames);
            end
        end
        flag=1;
    else
        iVarArgOut=iVarArgOut+1;
        varargout{iVarArgOut}=var_args{iVarArg};
    end
    iVarArg=iVarArg+1;
end