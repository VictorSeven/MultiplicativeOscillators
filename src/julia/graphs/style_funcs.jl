function cm_to_inches(cm)
    return cm * 0.393701
end

function one_col_size(ratio=1.618)
    width_pt = 246
    heigth_pt = width_pt / ratio
    return (width_pt, heigth_pt)
end

function two_col_size(ratio=1.618)
    width_pt = 510 
    heigth_pt = width_pt / ratio
    return (width_pt, heigth_pt)
end


function create_legend(axis, pos)
    axislegend(axis, position=pos, framevisible=false, rowgap=0, patchlabelgap=1, patchsize=(10,10)) 
end

function label_axes(axes; pos=[0.05, 0.9])
    if length(size(pos)) > 1
        alphabet = "abcdefghijklmnopqrstuvwxyz"
        for (i,ax) in enumerate(axes)
            text!(ax, pos[i,1], textsize=10, pos[i,2]; text="($(alphabet[i:i]))", space=:relative)
        end
    else
        alphabet = "abcdefghijklmnopqrstuvwxyz"
        for (i,ax) in enumerate(axes)
            text!(ax, pos[1], pos[2]; textsize=10, text="($(alphabet[i:i]))", space=:relative)
        end
    end
end

