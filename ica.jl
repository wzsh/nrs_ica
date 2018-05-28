
module WZS
    export ICA
    module ICA
        export
            Imperialist,
            Colony,
            Empire

        mutable struct Imperialist
            index::Int64
            nf::Int64
            #power::Float64
        end

        mutable struct Colony
            index::Int64
            nf::Int64
        end

        mutable struct Empire
            imp_index::Int64
            colonies_index::Array{Int64}
            total_fitness::Float32
            power::Float32
        end
    end
end
