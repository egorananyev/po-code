function outnum = c2n(cellTable, varName, curRow)

outnum = cell2mat( cellTable( curRow+1, strcmp(cellTable(1,:),varName) ) );

end