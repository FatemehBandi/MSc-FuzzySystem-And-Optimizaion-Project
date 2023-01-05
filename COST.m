function C=COST(y,yh)

C=(sum(yh-y)).^2/numel(y);

end