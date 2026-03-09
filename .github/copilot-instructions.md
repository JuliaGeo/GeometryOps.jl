# GeometryOps.jl Developer Instructions

GeometryOps.jl is a Julia package for geometric calculations on (primarily 2D) geometries, built on top of GeoInterface.jl. The package uses a literate programming approach with documentation embedded directly in source code.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively

### Environment Setup and Build Process
**CRITICAL**: GeometryOps.jl consists of two packages that must be developed together:
- Main package: `GeometryOps` (in repository root)
- Core package: `GeometryOpsCore` (in `./GeometryOpsCore` subdirectory)

Bootstrap, build, and test the repository:
1. `cd /path/to/GeometryOps.jl`
2. `julia --project=. -e "using Pkg; Pkg.develop(path=joinpath(\".\", \"GeometryOpsCore\")); Pkg.instantiate()"` -- takes 3-5 minutes. NEVER CANCEL. Set timeout to 10+ minutes.
3. `julia --project=. -e "using Pkg; Pkg.build()"` -- takes 1-2 minutes. NEVER CANCEL. Set timeout to 5+ minutes.
4. Verify installation: `julia --project=. -e "using GeometryOps; println(\"Package loaded successfully\")"` -- takes < 1 minute

### Testing
**WARNING**: Full test suite currently has dependency conflicts in some environments. Core functionality works correctly.

Test core functionality manually:
- `julia --project=. -e "using GeometryOps, GeoInterface; poly = GeoInterface.Polygon([[(0,0), (1,0), (1,1), (0,1), (0,0)]]); println(\"Area: \", GeometryOps.area(poly)); println(\"Centroid: \", GeometryOps.centroid(poly))"`
- Expected output: `Area: 1.0` and `Centroid: (0.5, 0.5)`

For full test suite (if dependencies allow):
- `julia --project=. -e "using Pkg; Pkg.test()"` -- takes 15-30 minutes if successful. NEVER CANCEL. Set timeout to 45+ minutes.

### Documentation Build
**WARNING**: Documentation build requires extensive dependencies and may fail in limited network environments.

Build documentation (if network allows):
1. `julia --project=docs -e "using Pkg; Pkg.develop([Pkg.PackageSpec(name=\"GeometryOpsCore\", path=\"./GeometryOpsCore\"), Pkg.PackageSpec(name=\"GeometryOps\", path=\".\")]); Pkg.instantiate()"` -- takes 10-20 minutes. NEVER CANCEL. Set timeout to 30+ minutes.
2. `julia --project=docs docs/make.jl` -- takes 5-15 minutes. NEVER CANCEL. Set timeout to 25+ minutes.

## Validation Scenarios

**ALWAYS test these scenarios after making changes:**

### Core Functionality Validation
Run these exact commands to verify basic operations work:
```julia
using GeometryOps, GeoInterface

# Test area calculation
poly = GeoInterface.Polygon([[(0,0), (1,0), (1,1), (0,1), (0,0)]])
area_result = GeometryOps.area(poly)  # Should be 1.0

# Test centroid calculation  
centroid_result = GeometryOps.centroid(poly)  # Should be (0.5, 0.5)

# Test distance calculation
point = (0.5, 0.5)
dist_result = GeometryOps.distance(point, poly)  # Should be 0.0

# Test signed area
signed_area_result = GeometryOps.signed_area(poly)  # Should be 1.0

println("Area: $area_result, Centroid: $centroid_result, Distance: $dist_result, Signed Area: $signed_area_result")
```

### Package Loading Test
Always verify the package loads without errors:
```bash
julia --project=. -e "using GeometryOps; println(\"SUCCESS: GeometryOps loaded\")"
```

## Repository Structure and Navigation

### Key Directories
- `src/` - Main package source code (124 Julia files)
  - `src/GeometryOps.jl` - Main module file
  - `src/methods/` - Core geometric operations (area, distance, centroid, etc.)
  - `src/transformations/` - Geometry transformation functions
  - `src/utils/` - Utility functions and spatial tree interfaces
- `GeometryOpsCore/` - Core package with fundamental types and interfaces
- `test/` - Test suite organized to match `src/` structure
- `docs/` - Documentation using Literate.jl
- `ext/` - Package extensions for Proj, LibGEOS, etc.

### Important Files
- `Project.toml` - Main package dependencies and metadata
- `GeometryOpsCore/Project.toml` - Core package dependencies
- `docs/make.jl` - Documentation build script
- `.github/workflows/CI.yml` - CI configuration with build and test jobs

### Coding Standards
- **Literate Programming**: Each source file includes extensive documentation at the top with examples
- **Example-Driven**: Every file should have visual examples using plots (typically Makie)
- **GeoInterface Compatible**: All functions work with any GeoInterface.jl compatible geometry
- **Algorithm Types**: Define Algorithm types to control function behavior
- **Internal Functions**: Begin with `_` and generally lack docstrings

## Common Development Tasks

### Adding New Functionality
1. Follow the pattern in `src/methods/area.jl` or `src/methods/distance.jl`
2. Include comprehensive documentation with examples at the top of the file
3. Add visual examples using Makie when possible
4. Define appropriate Algorithm types for different behaviors
5. Export only functions intended for end users
6. Add corresponding tests in `test/` directory

### Key Methods to Understand
Review these files as templates for new functionality:
- `src/methods/area.jl` - Calculate area and signed area
- `src/methods/distance.jl` - Distance calculations
- `src/methods/centroid.jl` - Centroid computation
- `src/methods/contains.jl` - Spatial relationship tests

### Working with Extensions
The package includes extensions for:
- `Proj` - Coordinate reference system transformations
- `LibGEOS` - GEOS library integration
- `FlexiJoins` - Table joining operations
- `TGGeometry` - Alternative geometry representations

## Julia Version and Compatibility

- **Minimum Julia Version**: 1.10
- **CI Tests**: Julia 1.10, 1.11, and nightly
- **Platform Support**: Ubuntu, Windows, macOS (ARM64 for macOS)
- **Architecture**: Primarily x64, ARM64 for macOS

## Validation Commands Reference

### Quick Health Check (< 1 minute)
```bash
julia --project=. -e "using GeometryOps; println(\"OK\")"
```

### Functionality Test (< 1 minute)  
```bash
julia --project=. -e "using GeometryOps, GeoInterface; poly = GeoInterface.Polygon([[(0,0), (1,0), (1,1), (0,1), (0,0)]]); println(\"Area: \", GeometryOps.area(poly))"
```

### Package Status Check
```bash
julia --project=. -e "using Pkg; Pkg.status()"
```

### Check GeometryOpsCore Development
```bash
julia --project=. -e "using Pkg; println(Pkg.project().dependencies)"
```

## Performance Considerations

- Use statically sized, immutable types where possible
- Avoid unnecessary allocations by using GeoInterface constructs directly
- Propagate compile-time information for better performance
- Consider using TimerOutputs.jl for performance analysis
- Profile with ProfileView.jl and diagnose with Cthulhu.jl for type stability

## Documentation and Examples

- All source files serve as documentation via literate programming
- Examples should be visual when dealing with geometric operations
- Documentation is generated from source code using Literate.jl
- Browse generated docs at: https://juliageo.org/GeometryOps.jl/stable

## Troubleshooting

### Common Issues
- **"GeometryOpsCore not found"**: Run the development setup commands above
- **Test failures**: Often due to missing optional dependencies; core functionality should still work
- **Documentation build failures**: Usually due to network issues or missing dependencies
- **Package loading errors**: Ensure both main package and GeometryOpsCore are properly developed

### Emergency Recovery
If the environment gets corrupted:
1. `rm -rf ~/.julia/compiled` (clears compiled cache)
2. Re-run the bootstrap commands from the "Working Effectively" section

Always commit your changes before attempting recovery procedures.