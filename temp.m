s = cell(2,8,1);
each_cell={[inf]};
% Time for each delayed reaction is scheduled in these vectors.
s(:,:)=each_cell;
Td5 = [1;2];
s{1,5}
s{1,5} = s{1,5} (1:Td5(1)-1);
s{1,5} = [s{1,5}, 72, inf];
Td5(1) = Td5(1) + 1;
s{1,5}
Td5
s{1,5}=s{1,5}(2:Td5(1));
s{1,5}
Td5

% %% s5(ck,:) = [s5(ck,1:Td5(ck) - 1), nph6, inf]; 
% %s(1,5) = {[Inf]}
% %first get {[]}
% %how to index the cell array?
% if s{1,5} == inf
%     s{1,5} = {nph6, inf};
% else
%     s{1,5} = {s{1,5}(2:end), nph6, inf};
% end
% % should become s(1,5) = {[nph6value, inf]}
% % without s(2,5) changing
% 
% %% s5(ck,:)=s5(ck,2:Td5(ck));
% % s(1,5) = {[nph6value, inf]}
% % should become s(1,5) = {[inf]} 
% % without s(2,5) changing
% 
% s{1,5} = s{1,5}(2:end);