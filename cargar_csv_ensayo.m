function T = cargar_csv_ensayo(fname)
T = readtable(fname, "Delimiter",";", "TextType","string", "VariableNamingRule","preserve");

% quitar fila de unidades si la primera fila contiene texto típico
row1 = lower(string(T{1,:}));
if any(contains(row1, ["mm","kn","mpa","(","s)","(mm)","(kn)"]))
    T(1,:) = [];
end

% convertir columnas a double (coma->punto) si hace falta
for k = 1:width(T)
    vn = T.Properties.VariableNames{k};
    col = T.(vn);
    if isstring(col) || iscellstr(col) || ischar(col)
        T.(vn) = str2double(strrep(string(col), ",", "."));
    end
end
end
