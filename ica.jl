
module WZS
    export ICA
    module ICA
        export
            Empire

        mutable struct Empire
            imp_index::Int64
            colonies_index::Array{Int64}
        end
    end
end
