function  simtable=DLEM_AG_read(outnamenew)
    % read output
    lines2=readlines(outnamenew);
    tf=strncmpi(lines2,"year:",5);
    allday=lines2(tf);

      values=[];
for i=1:length(allday)
dd=split(allday(i,1));
dd2=split(dd([1:62,64,66:140,146:151,159:181],1),':');%pheonology
values(i,:)=str2double(dd2(:,2))';
end

simtable=array2table(values,'VariableNames',dd2(:,1));