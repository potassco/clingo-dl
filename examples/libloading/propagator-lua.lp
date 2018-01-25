#script (lua)

package.cpath  = package.cpath .. ';' .. '.' .. '/?.so'

function main(prg)
    require("libluaclingoDL")
    loadPropagator(prg:get_c_control())
    prg:ground({{"base", {}}})
    prg:solve{}
end

#end.

{a;b}.

% both can't be true at the same time
&diff{x-y} <= 5 :- a.
&diff{y-x} <= -6 :- b.

