abstract type AbstractAnalysis end

abstract type AbstractPlotSpecification <: AbstractAnalysis end
abstract type AbstractSpaceTimePlotSpecification <: AbstractPlotSpecification end

# Default that should work for most
function output_name(aps::AbstractPlotSpecification)
	aps.output_name
end
