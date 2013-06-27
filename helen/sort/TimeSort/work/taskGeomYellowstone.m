function [taskList] = taskGeomYellowstone(ntasks, tasksPerNode, ensSize)
%TASKGEOM
%
% if ensSize > ntasks this fails.  Of course, you have no need to run it
% if ensSize > ntasks

%ntasks = 128;
%tasksPerNode = 16;
%ensSize = 56;
% EnsIndex = index of which tasks have ensemble members

if rem(ntasks,tasksPerNode) ~= 0
   
    error(' ntasks/tasksPerNode must be a whole number')
    
end

nodes = ntasks/tasksPerNode;


% make task list for each node
tasks = 0:1:ntasks-1;

% reshape tasks
tasksOnNodes = reshape(tasks, tasksPerNode, nodes);

% blank copy of this
taskList = zeros(size(tasksOnNodes));
taskList(:) = NaN;

% split ensemble members across nodes
perNode = zeros(nodes,1);

% need to account for ensSize/nodes having a remainder
%perNode(:) = ensSize/nodes;
perNode(:) = idivide(int32(ensSize),int32(nodes));

% distibute the leftovers to end nodes since task zero
% on node 1, already has the most memory
leftOvers = ensSize - sum(perNode);

for m = 0:leftOvers - 1
   
    perNode(end - m) = perNode(end - m) + 1;
    
end

count = 0;

for col = 1:nodes
    
    for row = 1:perNode(col)
        
        taskList(row, col) = count;
        
        count = count + 1;
        
    end
end

tasksFromOne=taskList(:) + 1;

EnsIndex = find(~isnan(tasksFromOne));

% split rest of tasks

count = ensSize;

for col = 1:nodes
    
    for row = perNode(col)+1:tasksPerNode
        
        taskList(row, col) = count;
        
        count = count + 1;
        
    end
    
end

% print task geometry
fid = fopen('task_list', 'w');

fprintf(fid, '%d \n', ntasks);

for i=1:ntasks    
fprintf(fid, '%d \n', taskList(i));
end
    

fclose(fid);

end
