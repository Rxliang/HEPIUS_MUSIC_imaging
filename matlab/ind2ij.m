function [i,j] = ind2ij(ind, rows)

j = ceil(ind/rows);
i = ind-(j-1)*rows;
