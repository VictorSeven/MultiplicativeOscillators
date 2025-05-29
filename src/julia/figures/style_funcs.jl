module StyleFuncs
export  one_col_figure, two_col_figure, label_axes, cm_to_inches, mm_to_inches

using CairoMakie

function cm_to_inches(cm)
    return cm * 0.393701
end

function mm_to_inches(mm)
    return mm * 3.93701
end

function label_axes(axes; pos = [0.05, 0.9])
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    if length(size(pos)) == 1
        for (i,ax) in enumerate(axes)
            text!(ax, pos[i], 0.9; fontsize=10, text="($(alphabet[i:i]))", space=:relative)
        end
    else
        for (i,ax) in enumerate(axes)
            text!(ax, pos[i,1], pos[i,2]; fontsize=10, text="($(alphabet[i:i]))", space=:relative)
        end
    end
end



# --- A theme adequate for a paper figure

paper_theme = Theme(
    fontsize=9,
    figure_padding=5,
    
    #Default: one col figure
    size=(246., 170.), 

    Axis=(

        #Disable grid
        xgridvisible=false,
        ygridvisible=false,

        #Set the label size
        xlabelsize=10,
        ylabelsize=10,

        #The way I like the spines
        spinewidth=1.2,
        rightspinevisible=false,
        topspinevisible=false,

        #Reduced tick size
        xminorticksize=2.0,
        yminorticksize=2.0,
        xticksize=3.5,
        yticksize=3.5,

        #Reduce label pads
        xticklabelpad=0.0,
        yticklabelpad=1.0,

        xlabelpadding=  1.0,
        ylabelpadding=  2.0,
    ), 

    #Set markersize to real size
    Scatter=(markersize=4.0,),

    #Customize legend
    Legend=(
        framevisible=false,
        rowgap = 0.0,
        patchlabelgap = 2.0,
        patchsize=(7.5, 10.0),
        labelsize = 8,
        backgroundcolor =:transparent
    )
)

#Return a theme for a figure with single-column size
function one_col_figure(ratio=1.618, width_pt=246)
    heigth_pt = width_pt / ratio
    paper_theme.size = (width_pt, heigth_pt)
    return paper_theme
end

#Return a theme for a figure with double-column size
function two_col_figure(ratio=1.618, width_pt=510)
    heigth_pt = width_pt / ratio
    paper_theme.size = (width_pt, heigth_pt)
    return paper_theme
end




# --- Slide theme, more adequate for presentations --- 
slide_theme = Theme(
    fontsize=24,
    figure_padding=5,
    
    #Default: one col figure
    size=(800,600), 

    Axis=(

        #Disable grid
        xgridvisible=false,
        ygridvisible=false,

        #Set the label size
        xlabelsize=24,
        ylabelsize=24,

        #The way I like the spines
        spinewidth=1.5,
        rightspinevisible=false,
        topspinevisible=false,

        #Reduced tick size
        xminorticksize=2.0,
        yminorticksize=2.0,
        xticksize=10,
        yticksize=10,

        #Reduce label pads
        xticklabelpad=0.0,
        yticklabelpad=1.0,

        xlabelpadding=  2.0,
        ylabelpadding=  5.0,
    ), 

    #Set markersize to real size
    Scatter=(markersize=12.0,),

    Lines=(linewidth=3.0,),

    StepHist=(linewidth=3.0,),

    #Customize legend
    Legend=(
        framevisible=false,
        rowgap = 0.0,
        patchlabelgap = 3.0,
        patchsize=(17.0, 30.0),
        labelsize = 24,
    )
)

function slide_figure(ratio=1.618, width_pt=800.)
    heigth_pt = width_pt / ratio
    slide_theme.size = (width_pt, heigth_pt)
    return slide_theme 
end


end


