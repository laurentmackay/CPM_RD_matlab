function ref = nameref(nm)

code=' \t\f\*\+\-\/\,\=\(\)\[\]\>\<\&\~\;\:\|\{\}\^\.\n';
ref = strcat(['(?<=(?:^|[' code  ']))('],nm,[')(?=(?:$|[' code ']))'])';
end

