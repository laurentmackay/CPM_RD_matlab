function ref = nameref(nm)

code=' \t\f\*\+\-\/\,\=\(\)\[\]\>\<\&\~\;\:\|\{\}\^\.';
ref = strcat(['(?<=(?:^|[' code  ']))('],nm,[')(?=(?:$|[' code ']))'])';
end

