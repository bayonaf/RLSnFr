function [row, col] = calculateGridFromStacksNumber(stacksNumber, caseType)
%CALCULATEGRIDFROMSTACKSNUMBER Summary of this function goes here

    switch caseType
        case "SynapsesPlots"
            if mod(stacksNumber,2) == 0
                fact = factor(stacksNumber + 4);
            else
                fact = factor(stacksNumber + 3);
            end
        case "DFMontage"
            fact = factor(stacksNumber);
    end
    if length(fact) == 6
        m = fact(2) * fact(3) * fact(4);
        M = fact(1) * fact(5) * fact(6);
        row = max([m M]);
        col = min([m M]);
    elseif length(fact) == 5
        m = fact(3) * fact(4);
        M = fact(1) * fact(2) * fact(5);
        row = max([m M]);
        col = min([m M]);
    elseif length(fact) == 4
        m = fact(1) * fact(4);
        M = fact(2) * fact(3);
        row = max([m M]);
        col = min([m M]);
    elseif length(fact) == 3
        M = (fact(1) * fact(2));
        m = fact(3);
        row = max([m M]);
        col = min([m M]);
    elseif length(fact) == 2
        row = ceil(sqrt(prod(fact)));
        col = row;
    end
end


