for i=1:512;
echo(i,1)=exp(-time(i,1)*10.0/5.0);
end;
plot(echo);