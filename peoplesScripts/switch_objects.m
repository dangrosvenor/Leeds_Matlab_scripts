% This function is used to bring the original axes
% on the foreground and the transformed axis to the background
function switch_objects_depth(parenthandle,obj1,obj2);
children = get(parenthandle,'Children');
obj1pos = find(children==obj1);
obj2pos = find(children==obj2);
children(obj1pos) = obj2;
children(obj2pos) = obj1;
set(parenthandle,'Children',children);