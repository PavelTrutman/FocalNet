function gfield(hObject,eventdata,handles)
d = str2double(get(hObject,'string'))
if isnan(d)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
