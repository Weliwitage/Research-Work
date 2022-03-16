function [result] = repcolumn(A, n)
    %n - how many times each column from A should be repeated

    [rows columns] = size(A);
    result = repmat(A(:,1),1,n);

    for i = 2:columns
        result = [result,repmat(A(:,i),1,n)];
    end
end

