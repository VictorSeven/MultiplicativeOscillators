module ArtsyPalettes
using CairoMakie 

export list_available, is_colorblindfriendly, met_brew


met_brewer = Dict(
    "Morgenstern" => Dict("order"=>[7, 5, 4, 6, 3, 2, 1], "colors"=>["#7c668c", "#b08ba5", "#dfbbc8", "#ffc680", "#ffb178", "#db8872", "#a56457"], "colorblind"=>true),
    "Moreau" => Dict("order"=>[2, 5, 3, 4, 7, 1, 6], "colors"=>["#421600", "#792504", "#bc7524", "#8dadca", "#527baa", "#104839", "#082844"], "colorblind"=>false),
    "Tara" => Dict("order"=>[1, 3, 2, 5, 4], "colors"=>["#eab1c6", "#d35e17", "#e18a1f", "#e9b109", "#829d44"], "colorblind"=>false),
    "Java" => Dict("order"=>[1, 4, 2, 5, 3], "colors"=>["#663171", "#cf3a36", "#ea7428", "#e2998a", "#0c7156"], "colorblind"=>true),
    "Kandinsky" => Dict("order"=>[1, 2, 3, 4], "colors"=>["#3b7c70", "#ce9642", "#898e9f", "#3b3a3e"], "colorblind"=>true),
    "Klimt" => Dict("order"=>[5, 2, 3, 4, 6, 1], "colors"=>["#df9ed4", "#c93f55", "#eacc62", "#469d76", "#3c4b99", "#924099"], "colorblind"=>false),
    "Degas" => Dict("order"=>[5, 2, 1, 3, 4, 7, 6], "colors"=>["#591d06", "#96410e", "#e5a335", "#556219", "#418979", "#2b614e", "#053c29"], "colorblind"=>false),
    "Greek" => Dict("order"=>[2, 3, 5, 1, 4], "colors"=>["#3c0d03", "#8d1c06", "#e67424", "#ed9b49", "#f5c34d"], "colorblind"=>true),
    "Monet" => Dict("order"=>[2, 5, 8, 3, 4, 9, 1, 6, 7], "colors"=>["#4e6d58", "#749e89", "#abccbe", "#e3cacf", "#c399a2", "#9f6e71", "#41507b", "#7d87b2", "#c2cae3"], "colorblind"=>false),
    "Homer2" => Dict("order"=>[3, 7, 1, 4, 6, 2, 5], "colors"=>["#bf3626", "#e9724c", "#e9851d", "#f9c53b", "#aeac4c", "#788f33", "#165d43"], "colorblind"=>false),
    "Redon" => Dict("order"=>[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], "colors"=>["#5b859e", "#1e395f", "#75884b", "#1e5a46", "#df8d71", "#af4f2f", "#d48f90", "#732f30", "#ab84a5", "#59385c", "#d8b847", "#b38711"], "colorblind"=>false),
    "Ingres" => Dict("order"=>[4, 5, 3, 6, 2, 7, 1, 8], "colors"=>["#041d2c", "#06314e", "#18527e", "#2e77ab", "#d1b252", "#a97f2f", "#7e5522", "#472c0b"], "colorblind"=>true),
    "Hokusai2" => Dict("order"=>[5, 2, 4, 1, 6, 3], "colors"=>["#abc9c8", "#72aeb6", "#4692b0", "#2f70a1", "#134b73", "#0a3351"], "colorblind"=>true),
    "Manet" => Dict("order"=>[8, 3, 10, 4, 7, 9, 11, 2, 6, 1, 5], "colors"=>["#3b2319", "#80521c", "#d29c44", "#ebc174", "#ede2cc", "#7ec5f4", "#4585b7", "#225e92", "#183571", "#43429b", "#5e65be"], "colorblind"=>false),
    "Peru2" => Dict("order"=>[4, 1, 3, 5, 2, 7, 6], "colors"=>["#65150b", "#961f1f", "#c0431f", "#b36c06", "#f19425", "#c59349", "#533d14"], "colorblind"=>false),
    "VanGogh3" => Dict("order"=>[7, 5, 1, 4, 8, 2, 3, 6], "colors"=>["#e7e5cc", "#c2d6a4", "#9cc184", "#669d62", "#447243", "#1f5b25", "#1e3d14", "#192813"], "colorblind"=>true),
    "Derain" => Dict("order"=>[4, 2, 5, 7, 1, 3, 6], "colors"=>["#efc86e", "#97c684", "#6f9969", "#aab5d5", "#808fe1", "#5c66a8", "#454a74"], "colorblind"=>true),
    "Isfahan1" => Dict("order"=>[5, 2, 4, 6, 1, 7, 3], "colors"=>["#4e3910", "#845d29", "#d8c29d", "#4fb6ca", "#178f92", "#175f5d", "#1d1f54"], "colorblind"=>true),
    "Signac" => Dict("order"=>[13, 3, 2, 1, 11, 5, 8, 14, 12, 10, 7, 4, 6, 9], "colors"=>["#fbe183", "#f4c40f", "#fe9b00", "#d8443c", "#9b3441", "#de597c", "#e87b89", "#e6a2a6", "#aa7aa1", "#9f5691", "#633372", "#1f6e9c", "#2b9b81", "#92c051"], "colorblind"=>false),
    "Nattier" => Dict("order"=>[1, 6, 3, 4, 7, 2, 5], "colors"=>["#52271c", "#944839", "#c08e39", "#7f793c", "#565c33", "#184948", "#022a2a"], "colorblind"=>false),
    "Veronese" => Dict("order"=>[5, 1, 7, 2, 3, 6, 4], "colors"=>["#67322e", "#99610a", "#c38f16", "#6e948c", "#2c6b67", "#175449", "#122c43"], "colorblind"=>true),
    "Thomas" => Dict("order"=>[3, 2, 8, 6, 1, 4, 7, 5], "colors"=>["#b24422", "#c44d76", "#4457a5", "#13315f", "#b1a1cc", "#59386c", "#447861", "#7caf5c"], "colorblind"=>false),
    "Cassatt2" => Dict("order"=>[7, 3, 9, 1, 5, 6, 2, 10, 4, 8], "colors"=>["#2d223c", "#574571", "#90719f", "#b695bc", "#dec5da", "#c1d1aa", "#7fa074", "#466c4b", "#2c4b27", "#0e2810"], "colorblind"=>true),
    "Pillement" => Dict("order"=>[4, 3, 2, 5, 1, 6], "colors"=>["#a9845b", "#697852", "#738e8e", "#44636f", "#2b4655", "#0f252f"], "colorblind"=>true),
    "Gauguin" => Dict("order"=>[2, 5, 4, 3, 1, 6], "colors"=>["#b04948", "#811e18", "#9e4013", "#c88a2c", "#4c6216", "#1a472a"], "colorblind"=>false),
    "Paquin" => Dict("order"=>[10, 6, 1, 8, 4, 3, 5, 9, 2, 7, 11], "colors"=>["#831818", "#c62320", "#f05b43", "#f78462", "#feac81", "#f7dea3", "#ced1af", "#98ab76", "#748f46", "#47632a", "#275024"], "colorblind"=>false),
    "Cross" => Dict("order"=>[4, 7, 1, 8, 2, 6, 3, 5, 9], "colors"=>["#c969a1", "#ce4441", "#ee8577", "#eb7926", "#ffbb44", "#859b6c", "#62929a", "#004f63", "#122451"], "colorblind"=>false),
    "Peru1" => Dict("order"=>[3, 1, 5, 2, 4, 6], "colors"=>["#b5361c", "#e35e28", "#1c9d7c", "#31c7ba", "#369cc9", "#3a507f"], "colorblind"=>false),
    "Archambault" => Dict("order"=>[2, 7, 5, 1, 6, 4, 3], "colors"=>["#88a0dc", "#381a61", "#7c4b73", "#ed968c", "#ab3329", "#e78429", "#f9d14a"], "colorblind"=>true),
    "Cassatt1" => Dict("order"=>[3, 6, 1, 8, 4, 5, 2, 7], "colors"=>["#b1615c", "#d88782", "#e3aba7", "#edd7d9", "#c9c9dd", "#9d9dc7", "#8282aa", "#5a5a83"], "colorblind"=>true),
    "Benedictus" => Dict("order"=>[9, 5, 11, 1, 7, 3, 13, 4, 8, 2, 12, 6, 10], "colors"=>["#9a133d", "#b93961", "#d8527c", "#f28aaa", "#f9b4c9", "#f9e0e8", "#ffffff", "#eaf3ff", "#c5daf6", "#a1c2ed", "#6996e3", "#4060c8", "#1a318b"], "colorblind"=>false),
    "Tsimshian" => Dict("order"=>[6, 1, 7, 4, 1, 5, 3], "colors"=>["#582310", "#aa361d", "#82c45f", "#318f49", "#0cb4bb", "#2673a3", "#473d7d"], "colorblind"=>false),
    "Nizami" => Dict("order"=>[5, 2, 6, 8, 3, 7, 4, 1], "colors"=>["#dd7867", "#b83326", "#c8570d", "#edb144", "#8cc8bc", "#7da7ea", "#5773c0", "#1d4497"], "colorblind"=>false),
    "Robert" => Dict("order"=>[2, 5, 3, 1, 6, 4], "colors"=>["#11341a", "#375624", "#6ca4a0", "#487a7c", "#18505f", "#062e3d"], "colorblind"=>false),
    "Navajo" => Dict("order"=>[1, 2, 3, 4, 5], "colors"=>["#660d20", "#e59a52", "#edce79", "#094568", "#e1c59a"], "colorblind"=>false),
    "Hiroshige" => Dict("order"=>[6, 2, 9, 3, 7, 5, 1, 10, 4, 8], "colors"=>["#e76254", "#ef8a47", "#f7aa58", "#ffd06f", "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", "#376795", "#1e466e"], "colorblind"=>true),
    "Wissing" => Dict("order"=>[2, 3, 5, 4, 1], "colors"=>["#4b1d0d", "#7c291e", "#ba7233", "#3a4421", "#2d5380"], "colorblind"=>false),
    "OKeeffe2" => Dict("order"=>[7, 1, 6, 4, 2, 5, 3], "colors"=>["#fbe3c2", "#f2c88f", "#ecb27d", "#e69c6b", "#d37750", "#b9563f", "#92351e"], "colorblind"=>true),
    "Isfahan2" => Dict("order"=>[4, 2, 3, 5, 1], "colors"=>["#d7aca1", "#ddc000", "#79ad41", "#34b6c6", "#4063a3"], "colorblind"=>true),
    "Tiepolo" => Dict("order"=>[1, 2, 8, 4, 3, 5, 7, 6], "colors"=>["#802417", "#c06636", "#ce9344", "#e8b960", "#646e3b", "#2b5851", "#508ea2", "#17486f"], "colorblind"=>false),
    "Troy" => Dict("order"=>[2, 7, 4, 5, 1, 8, 3, 6], "colors"=>["#421401", "#6c1d0e", "#8b3a2b", "#c27668", "#7ba0b4", "#44728c", "#235070", "#0a2d46"], "colorblind"=>true),
    "Egypt" => Dict("order"=>[1, 2, 3, 4], "colors"=>["#dd5129", "#0f7ba2", "#43b284", "#fab255"], "colorblind"=>true),
    "Stevens" => Dict("order"=>[4, 2, 3, 5, 1, 6], "colors"=>["#042e4e", "#307d7f", "#598c4c", "#ba5c3f", "#a13213", "#470c00"], "colorblind"=>false),
    "Homer1" => Dict("order"=>[6, 3, 2, 7, 4, 8, 5, 1], "colors"=>["#551f00", "#a62f00", "#df7700", "#f5b642", "#fff179", "#c3f4f6", "#6ad5e8", "#32b2da"], "colorblind"=>false),
    "Pissaro" => Dict("order"=>[6, 2, 4, 1, 7, 5, 3], "colors"=>["#134130", "#4c825d", "#8cae9e", "#8dc7dc", "#508ca7", "#1a5270", "#0e2a4d"], "colorblind"=>false),
    "OKeeffe1" => Dict("order"=>[8, 6, 1, 4, 10, 3, 11, 5, 2, 7, 9], "colors"=>["#6b200c", "#973d21", "#da6c42", "#ee956a", "#fbc2a9", "#f6f2ee", "#bad6f9", "#7db0ea", "#447fdd", "#225bb2", "#133e7e"], "colorblind"=>true),
    "Juarez" => Dict("order"=>[1, 2, 3, 4, 5, 6], "colors"=>["#a82203", "#208cc0", "#f1af3a", "#cf5e4e", "#637b31", "#003967"], "colorblind"=>false),
    "Renoir" => Dict("order"=>[2, 5, 9, 12, 3, 8, 7, 10, 4, 1, 6, 11], "colors"=>["#17154f", "#2f357c", "#6c5d9e", "#9d9cd5", "#b0799a", "#f6b3b0", "#e48171", "#bf3729", "#e69b00", "#f5bb50", "#ada43b", "#355828"], "colorblind"=>false),
    "VanGogh2" => Dict("order"=>[1, 5, 8, 2, 7, 4, 6, 3], "colors"=>["#bd3106", "#d9700e", "#e9a00e", "#eebe04", "#5b7314", "#c3d6ce", "#89a6bb", "#454b87"], "colorblind"=>false),
    "Lakota" => Dict("order"=>[1, 2, 3, 4, 5, 6], "colors"=>["#04a3bd", "#f0be3d", "#931e18", "#da7901", "#247d3f", "#20235b"], "colorblind"=>false),
    "Hokusai1" => Dict("order"=>[2, 7, 4, 6, 5, 1, 3], "colors"=>["#6d2f20", "#b75347", "#df7e66", "#e09351", "#edc775", "#94b594", "#224b5e"], "colorblind"=>false),
    "Demuth" => Dict("order"=>[9, 5, 1, 7, 3, 4, 8, 2, 6, 10], "colors"=>["#591c19", "#9b332b", "#b64f32", "#d39a2d", "#f7c267", "#b9b9b8", "#8b8b99", "#5d6174", "#41485f", "#262d42"], "colorblind"=>true),
    "Tam" => Dict("order"=>[3, 8, 1, 6, 2, 7, 4, 5], "colors"=>["#ffd353", "#ffb242", "#ef8737", "#de4f33", "#bb292c", "#9f2d55", "#62205f", "#341648"], "colorblind"=>true),
    "VanGogh1" => Dict("order"=>[3, 5, 7, 4, 6, 2, 1], "colors"=>["#2c2d54", "#434475", "#6b6ca3", "#969bc7", "#87bcbd", "#89ab7c", "#6f9954"], "colorblind"=>false),
    "NewKingdom" => Dict("order"=>[2, 1, 3, 4, 5], "colors"=>["#e1846c", "#9eb4e0", "#e6bb9e", "#9c6849", "#735852"], "colorblind"=>false),
    "Johnson" => Dict("order"=>[3, 1, 4, 2, 5], "colors"=>["#a00e00", "#d04e00", "#f6c200", "#0086a8", "#132b69"], "colorblind"=>true),
    "Austria" => Dict("order"=>[1, 2, 3, 4, 6, 5, 7], "colors"=>["#a40000", "#16317d", "#007e2f", "#ffcd12", "#b86092", "#721b3e", "#00b7a7"], "colorblind"=>false),
    "Hokusai3" => Dict("order"=>[4, 2, 5, 3, 1, 6], "colors"=>["#d8d97a", "#95c36e", "#74c8c3", "#5a97c1", "#295384", "#0a2e57"], "colorblind"=>true),
)

    """
        met_brew(name::String; n::Int=128, categorical::Bool=true)

    Generate a color palette from the MetBrewer collection, invoking it by 
    its `name`. If `categorical` is false, generates a continuous colormap
    with `n` colors.
    """
    function met_brew(name::String; n=128, categorical=true)
        palette = met_brewer[name]

        #Generate the color palette
        if categorical
            return palette["colors"]
        else
            return cgrad(palette["colors"], n) 
        end
    end

    """
        list_available(only_colorblindfriendly=true)

    Shows a list with the names of the available palettes and the 
    collection they come from. Set `only_colorblindfriendly=false` to 
    see the full list.
    """
    function list_available(;only_colorblindfriendly=true)
        #Filter the colorblind friendly ones
        if only_colorblindfriendly
            #Function performing a filter on key,value 
            f((k,v)) = is_colorblindfriendly(k) 
            #Use the builtin filter with our f
            filtered_dict = filter(f, met_brewer)

            println(keys(filtered_dict))
        else
            #Show all
            println(keys(met_brewer))
        end
    end

    """
        is_colorblindfriendly(name::String)
    
    Checks if the color palette is adequate for colorblind people. 
    """
    function is_colorblindfriendly(name::String)
        return met_brewer[name]["colorblind"]
    end
    
    
end