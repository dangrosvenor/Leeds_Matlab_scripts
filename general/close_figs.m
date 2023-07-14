function close_figs(range)

for i=range
    try
        close(i)
    catch
    end
end