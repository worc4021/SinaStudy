function [Aout,bout] = bigReduce(Ain,bin)

if ~isempty(Ain)
    idx = true(size(bin));

    for i = 1:length(bin)
        idx(i) = false;
        curNum = numel(find(idx));
        model = struct('obj',ones(curNum,1),...
            'A',sparse([Ain(idx,:)';bin(idx)']),...
            'rhs', [Ain(i,:)';bin(i)],...
            'sense', [char(ones(size(Ain,2),1)*'=');char(ones(1)*'<')],...
            'lb',zeros(curNum,1));
        param = struct('OutputFlag', 0);
        res = gurobi(model,param);
        idx(i) = ~or(strcmp(res.status,'OPTIMAL'),strcmp(res.status,'SUBOPTIMAL'));
    end

    Aout = Ain(idx,:);
    bout = bin(idx);
else
    Aout = [];
    bout = [];
end