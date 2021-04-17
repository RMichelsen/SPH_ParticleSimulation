using System.Linq;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using UnityEngine.Rendering;
using UnityEngine.Profiling;
using UnityEditor;
using System;
using System.IO;

struct IndirectShaderData {
    float4x4 LocalToWorld;
    float4x4 WorldToLocal;
}

public struct Particle {
    public int bin;
    public int index;
}

public struct SimulateParticlesConstantBuffer {
    public int3 grid_dimensions;
    public int bin_count;
    public float cell_size;
    public float rest_density;
    public float particle_mass;
    public float radius;
    public float stiffness;
    public float viscosity;
    public float timestep;
}

public enum ScenarioName {
    Cube,
    DamBreak,
    DoubleDamBreak,
    DamBreakWithBunny
}

public class ParticleSpawner : MonoBehaviour {
    public ScenarioName scenarioName;
    public bool record;
    public Mesh particleMesh;
    public Material particleMaterial;
    public Material surfaceMaterial;
    public Material mcGridMaterial;
    public bool pauseSimulation;
    public bool drawSurface;
    public bool drawParticles;

    private bool recordingDone = false;
    private const float visualScale = 100.0f;
    private const float particleRadius = 0.01f;
    private const float restDensity = 1000.0f;
    private const float particleMass = smoothingLength * smoothingLength * smoothingLength * restDensity;
    private const float stiffness = 7.0f;
    private const float viscosity = 0.01f;
    private const float timeStep = 0.001f;
    private const float smoothingLength = particleRadius * 2.0f;
    private static readonly float3 bounds = new float3(3.0f, 3.0f, 3.0f);
    private const float cellSize = smoothingLength * 2.0f;
    private static readonly int3 gridDimensions = new int3(
        (int)(bounds.x / cellSize),
        (int)(bounds.y / cellSize),
        (int)(bounds.z / cellSize)
    );
    private const float marchingCubesCellSize = particleRadius;
    private static readonly int3 marchingCubesGridDimensions = new int3(
        (int)math.ceil(bounds.x / marchingCubesCellSize) + 1,
        (int)math.ceil(bounds.y / marchingCubesCellSize) + 1,
        (int)math.ceil(bounds.z / marchingCubesCellSize) + 1
    );
    private static readonly int marchingCubeGridPointCount = 
        marchingCubesGridDimensions.x * marchingCubesGridDimensions.y * marchingCubesGridDimensions.z;
    private static readonly int binCount = gridDimensions.x * gridDimensions.y * gridDimensions.z;
    private int particleCount;

    List<float3> positionsToDraw;

    private ComputeBuffer boundaryParticles;
    private ComputeBuffer boundaryPositions;
    private ComputeBuffer boundaryVolume;
    private ComputeBuffer boundaryBins;
    private ComputeBuffer boundaryPrefixSums;
    private ComputeBuffer positions;
    private ComputeBuffer sortedPositions;
    private ComputeBuffer particles;
    private ComputeBuffer sortedParticles;
    private ComputeBuffer velocities;
    private ComputeBuffer sortedVelocities;
    private ComputeBuffer bins;
    private ComputeBuffer prefixSums;
    private ComputeBuffer pressures;
    private ComputeBuffer densities;
    private ComputeBuffer colorGradients;
    private ComputeBuffer forces;
    private ComputeBuffer surfaceIdentifiers;
    private ComputeBuffer cellTags;
    private ComputeBuffer cellsToTriangulate;
    private ComputeBuffer cellsToTriangulateArgs;
    private ComputeBuffer scalarField;
    private ComputeBuffer triangulatedCells;
    private ComputeBuffer triangleConnectionTableBuffer;
    private ComputeBuffer edgeTableBuffer;
    private ComputeBuffer drawProceduralIndirectArgs;
    private ComputeBuffer particleShaderArgs;
    private ComputeBuffer particleShaderData;
    private ComputeBuffer colors;

    private ComputeShader sortParticles;
    private ComputeShader simulateParticles;
    private ComputeShader surfaceReconstruction;
    private ComputeShader positionsToShaderData;

    void InitializeComputeBuffers(List<float3> initialParticlePositions) {
        positions = new ComputeBuffer(particleCount, sizeof(float) * 3);                                                positions.name = "Positions";
        positions.SetData(initialParticlePositions);
        sortedPositions = new ComputeBuffer(particleCount, sizeof(float) * 3);                                          sortedPositions.name = "Sorted Positions";
        particles = new ComputeBuffer(particleCount, sizeof(int) * 2);                                                  particles.name = "Particles";
        sortedParticles = new ComputeBuffer(particleCount, sizeof(int) * 2);                                            sortedParticles.name = "Sorted Particles";
        velocities = new ComputeBuffer(particleCount, sizeof(float) * 3);                                               velocities.name = "Velocities";
        velocities.SetData(new float3[particleCount]);
        sortedVelocities = new ComputeBuffer(particleCount, sizeof(float) * 3);                                         sortedVelocities.name = "Sorted Velocities";
        bins = new ComputeBuffer(binCount, sizeof(int));                                                                bins.name = "Bins";
        prefixSums = new ComputeBuffer(binCount, sizeof(int));                                                          prefixSums.name = "Prefix Sums";
        boundaryBins = new ComputeBuffer(binCount, sizeof(int));                                                        boundaryBins.name = "Boundary Bins";
        boundaryPrefixSums = new ComputeBuffer(binCount, sizeof(int));                                                  boundaryPrefixSums.name = "Boundary Prefix Sums";
        pressures = new ComputeBuffer(particleCount, sizeof(float));                                                    pressures.name = "Pressures";
        densities = new ComputeBuffer(particleCount, sizeof(float));                                                    densities.name = "Densities";
        colorGradients = new ComputeBuffer(particleCount, sizeof(float) * 3);                                           colorGradients.name = "Color Gradients";
        forces = new ComputeBuffer(particleCount, sizeof(float) * 3);                                                   forces.name = "Forces";
        surfaceIdentifiers = new ComputeBuffer(particleCount, sizeof(int));                                             surfaceIdentifiers.name = "Surface Identifiers";
        cellTags = new ComputeBuffer(marchingCubeGridPointCount, sizeof(uint) * 4);                                     cellTags.name = "Cell Tags";
        cellsToTriangulate = new ComputeBuffer(marchingCubeGridPointCount, sizeof(uint) * 3, ComputeBufferType.Append); cellsToTriangulate.name = "Cells to Triangulate";
        cellsToTriangulateArgs = new ComputeBuffer(1, 3 * sizeof(uint), ComputeBufferType.IndirectArguments);           cellsToTriangulateArgs.name = "Cells to Triangulate Args";
        scalarField = new ComputeBuffer(marchingCubeGridPointCount, sizeof(float));                                     scalarField.name = "Scalar Field";
        triangulatedCells = new ComputeBuffer(marchingCubeGridPointCount, sizeof(float) * 6, ComputeBufferType.Append); triangulatedCells.name = "Triangulated Cells";
        triangleConnectionTableBuffer = new ComputeBuffer(MarchingCubesData.triangleConnectionTable.Length, sizeof(int)); triangleConnectionTableBuffer.name = "Triangle Connection Table Buffer";
        triangleConnectionTableBuffer.SetData(MarchingCubesData.triangleConnectionTable);
        edgeTableBuffer = new ComputeBuffer(MarchingCubesData.edgeTable.Length, sizeof(int));                           edgeTableBuffer.name = "Edge Table Buffer";
        edgeTableBuffer.SetData(MarchingCubesData.edgeTable);
        drawProceduralIndirectArgs = new ComputeBuffer(4, sizeof(uint), ComputeBufferType.IndirectArguments);           drawProceduralIndirectArgs.name = "Draw Procedural Indirect Args";
        particleShaderArgs = new ComputeBuffer(1, 5 * sizeof(uint), ComputeBufferType.IndirectArguments);               particleShaderArgs.name = "Particle Shader Args";
        uint[] args = { particleMesh.GetIndexCount(0), (uint)particleCount, particleMesh.GetIndexStart(0), particleMesh.GetBaseVertex(0), 0 };
        particleShaderArgs.SetData(args);
        particleShaderData = new ComputeBuffer(particleCount, sizeof(float) * 32);                                      particleShaderData.name = "Particle Shader Data";
        colors = new ComputeBuffer(particleCount, sizeof(float) * 3);                                                   colors.name = "Colors";
    }

    void SetSortParameters() {
        sortParticles = Resources.Load<ComputeShader>("Shaders/SortParticles");
        sortParticles.SetInt("bin_count", binCount);
        sortParticles.SetFloat("cell_size", cellSize);
        sortParticles.SetInts("grid_dimensions", new int[] { gridDimensions.x, gridDimensions.y, gridDimensions.z });
    }

    void SetSimulationParameters(List<float3> initialParticlePositions) {
        simulateParticles = Resources.Load<ComputeShader>("Shaders/SimulateParticles");
        simulateParticles.SetInts("grid_dimensions", new int[] { gridDimensions.x, gridDimensions.y, gridDimensions.z });
        simulateParticles.SetFloats("bounds", new float[] { bounds.x, bounds.y, bounds.z });
        simulateParticles.SetInt("bin_count", binCount);
        simulateParticles.SetInt("particle_count", initialParticlePositions.Count);
        simulateParticles.SetFloat("cell_size", cellSize);
        simulateParticles.SetFloat("mc_cell_size", marchingCubesCellSize);
        simulateParticles.SetInts("mc_grid_dimensions", new int[] { marchingCubesGridDimensions.x, marchingCubesGridDimensions.y, marchingCubesGridDimensions.z });
        simulateParticles.SetFloat("rest_density", restDensity);
        simulateParticles.SetFloat("particle_mass", particleMass);
        simulateParticles.SetFloat("stiffness", stiffness);
        simulateParticles.SetFloat("viscosity", viscosity);
        simulateParticles.SetFloat("dt", timeStep);
        simulateParticles.SetFloat("h", smoothingLength);

        particleCount = initialParticlePositions.Count;
    }

    void SetSurfaceReconstructionParameters() {
        surfaceReconstruction = Resources.Load<ComputeShader>("Shaders/SurfaceReconstruction");
        surfaceReconstruction.SetFloat("timestep", timeStep);
        surfaceReconstruction.SetFloat("particle_radius", particleRadius);
        surfaceReconstruction.SetFloat("cell_size", cellSize);
        surfaceReconstruction.SetInts("grid_dimensions", new int[] { gridDimensions.x, gridDimensions.y, gridDimensions.z });
        surfaceReconstruction.SetFloat("mc_cell_size", marchingCubesCellSize);
        surfaceReconstruction.SetInts("mc_grid_dimensions", new int[] { marchingCubesGridDimensions.x, marchingCubesGridDimensions.y, marchingCubesGridDimensions.z });
        surfaceReconstruction.SetInt("mc_num_grid_points", marchingCubeGridPointCount);
    }

    void InitializeBoundaries(List<float3> initBoundaryPositions) {
        int3 dims = new int3(
            (int)math.ceil(bounds.x / smoothingLength),
            (int)math.ceil(bounds.y / smoothingLength),
            (int)math.ceil(bounds.z / smoothingLength)
        );
        for (int z = 0; z < dims.z; ++z) {
            for (int y = 0; y < dims.y; ++y) {
                for (int x = 0; x < dims.x; ++x) {
                    if (x != 0 && y != 0 && z != 0 &&
                       x != dims.x - 1 && y != dims.y - 1 && z != dims.z - 1) continue;
                    initBoundaryPositions.Add((new float3(x, y, z) + new float3(0.5f)) * smoothingLength);
                }
            }
        }
        ComputeBuffer initialBoundaryParticles = new ComputeBuffer(initBoundaryPositions.Count, sizeof(int) * 2);
        ComputeBuffer initialBoundaryPositions = new ComputeBuffer(initBoundaryPositions.Count, sizeof(float) * 3);
        initialBoundaryPositions.SetData(initBoundaryPositions);
        boundaryParticles = new ComputeBuffer(initBoundaryPositions.Count, sizeof(int) * 2);    boundaryParticles.name = "Boundary Particles";
        boundaryPositions = new ComputeBuffer(initBoundaryPositions.Count, sizeof(float) * 3);  boundaryPositions.name = "Boundary Positions";
        boundaryVolume = new ComputeBuffer(initBoundaryPositions.Count, sizeof(float));         boundaryVolume.name = "Boundary Volume";
        SortBoundaryParticles(initialBoundaryParticles, initialBoundaryPositions);
        initialBoundaryParticles.Dispose();
        initialBoundaryPositions.Dispose();

        int kernel = simulateParticles.FindKernel("ComputeBoundaryVolume");
        simulateParticles.SetBuffer(kernel, "sorted_boundary_positions", boundaryPositions);
        simulateParticles.SetBuffer(kernel, "boundary_bins", boundaryBins);
        simulateParticles.SetBuffer(kernel, "boundary_prefix_sums", boundaryPrefixSums);
        simulateParticles.SetBuffer(kernel, "boundary_volume", boundaryVolume);
        simulateParticles.Dispatch(kernel, (int)math.ceil(initBoundaryPositions.Count / 32.0f), 1, 1);
    }

    void SetParticleRenderParameters() {
        positionsToShaderData = Resources.Load<ComputeShader>("Shaders/PositionsToShaderData");
        positionsToShaderData.SetMatrix("_LocalToWorld", transform.localToWorldMatrix);
        positionsToShaderData.SetMatrix("_WorldToLocal", transform.worldToLocalMatrix);
        positionsToShaderData.SetFloat("particle_radius", particleRadius);
        positionsToShaderData.SetFloat("visual_scale", visualScale);
    }

    void SortBoundaryParticles(ComputeBuffer initialParticles, ComputeBuffer initialPositions) {
        Profiler.BeginSample("Reset Prefix Sums and Bins");
        int kernel = sortParticles.FindKernel("ResetPrefixSumsAndBins");
        sortParticles.SetBuffer(kernel, "bins", boundaryBins);
        sortParticles.SetBuffer(kernel, "prefix_sums", boundaryPrefixSums);
        sortParticles.Dispatch(kernel, (int)math.ceil(binCount / 32.0), 1, 1);
        Profiler.EndSample();

        Profiler.BeginSample("Distribute Particles");
        kernel = sortParticles.FindKernel("DistributeParticles");
        sortParticles.SetBuffer(kernel, "particles", initialParticles);
        sortParticles.SetBuffer(kernel, "positions", initialPositions);
        sortParticles.SetBuffer(kernel, "bins", boundaryBins);
        sortParticles.Dispatch(kernel, (int)math.ceil(initialParticles.count / 32.0f), 1, 1);
        Profiler.EndSample();

        Profiler.BeginSample("Calculate Prefix Sums");
        kernel = sortParticles.FindKernel("CalculatePrefixSums");
        sortParticles.SetBuffer(kernel, "bins", boundaryBins);
        sortParticles.SetBuffer(kernel, "prefix_sums", boundaryPrefixSums);
        sortParticles.Dispatch(kernel, 1, 1, 1);
        Profiler.EndSample();

        Profiler.BeginSample("Copy Boundary Particles");
        kernel = sortParticles.FindKernel("CopyBoundaryParticles");
        sortParticles.SetBuffer(kernel, "prefix_sums", boundaryPrefixSums);
        sortParticles.SetBuffer(kernel, "particles", initialParticles);
        sortParticles.SetBuffer(kernel, "positions", initialPositions);
        sortParticles.SetBuffer(kernel, "sorted_particles", boundaryParticles);
        sortParticles.SetBuffer(kernel, "sorted_positions", boundaryPositions);
        sortParticles.Dispatch(kernel, (int)math.ceil(initialParticles.count / 32.0f), 1, 1);
        Profiler.EndSample();
    }

    void SortParticles() {
        Profiler.BeginSample("Reset Prefix Sums and Bins");
        int kernel = sortParticles.FindKernel("ResetPrefixSumsAndBins");
        sortParticles.SetBuffer(kernel, "bins", bins);
        sortParticles.SetBuffer(kernel, "prefix_sums", prefixSums);
        sortParticles.Dispatch(kernel, (int)math.ceil(binCount / 32.0), 1, 1);
        Profiler.EndSample();

        Profiler.BeginSample("Distribute Particles");
        kernel = sortParticles.FindKernel("DistributeParticles");
        sortParticles.SetBuffer(kernel, "particles", particles);
        sortParticles.SetBuffer(kernel, "positions", positions);
        sortParticles.SetBuffer(kernel, "bins", bins);
        sortParticles.Dispatch(kernel, (int)math.ceil(particleCount / 32.0f), 1, 1);
        Profiler.EndSample();

        Profiler.BeginSample("Calculate Prefix Sums");
        kernel = sortParticles.FindKernel("CalculatePrefixSums");
        sortParticles.SetBuffer(kernel, "bins", bins);
        sortParticles.SetBuffer(kernel, "prefix_sums", prefixSums);
        sortParticles.Dispatch(kernel, 1, 1, 1);
        Profiler.EndSample();

        Profiler.BeginSample("Copy Particles");
        kernel = sortParticles.FindKernel("CopyParticles");
        sortParticles.SetBuffer(kernel, "prefix_sums", prefixSums);
        sortParticles.SetBuffer(kernel, "particles", particles);
        sortParticles.SetBuffer(kernel, "positions", positions);
        sortParticles.SetBuffer(kernel, "velocities", velocities);
        sortParticles.SetBuffer(kernel, "sorted_particles", sortedParticles);
        sortParticles.SetBuffer(kernel, "sorted_positions", sortedPositions);
        sortParticles.SetBuffer(kernel, "sorted_velocities", sortedVelocities);
        sortParticles.Dispatch(kernel, (int)(particleCount / 32.0f), 1, 1);
        Profiler.EndSample();
    }

    void StepSimulation() {
        int kernel = simulateParticles.FindKernel("CalculateDensityPressure");
        simulateParticles.SetBuffer(kernel, "sorted_positions", sortedPositions);
        simulateParticles.SetBuffer(kernel, "sorted_boundary_positions", boundaryPositions);
        simulateParticles.SetBuffer(kernel, "boundary_volume", boundaryVolume);
        simulateParticles.SetBuffer(kernel, "bins", bins);
        simulateParticles.SetBuffer(kernel, "boundary_bins", boundaryBins);
        simulateParticles.SetBuffer(kernel, "prefix_sums", prefixSums);
        simulateParticles.SetBuffer(kernel, "boundary_prefix_sums", boundaryPrefixSums);
        simulateParticles.SetBuffer(kernel, "densities", densities);
        simulateParticles.SetBuffer(kernel, "pressures", pressures);
        Profiler.BeginSample("Calculate Density and Pressure");
        simulateParticles.Dispatch(kernel, (int)math.ceil(particleCount / 32.0f), 1, 1);
        Profiler.EndSample();

        kernel = simulateParticles.FindKernel("CalculateForces");
        simulateParticles.SetBuffer(kernel, "sorted_particles", sortedParticles);
        simulateParticles.SetBuffer(kernel, "sorted_positions", sortedPositions);
        simulateParticles.SetBuffer(kernel, "sorted_boundary_positions", boundaryPositions);
        simulateParticles.SetBuffer(kernel, "sorted_velocities", sortedVelocities);
        simulateParticles.SetBuffer(kernel, "boundary_volume", boundaryVolume);
        simulateParticles.SetBuffer(kernel, "bins", bins);
        simulateParticles.SetBuffer(kernel, "boundary_bins", boundaryBins);
        simulateParticles.SetBuffer(kernel, "prefix_sums", prefixSums);
        simulateParticles.SetBuffer(kernel, "boundary_prefix_sums", boundaryPrefixSums);
        simulateParticles.SetBuffer(kernel, "densities", densities);
        simulateParticles.SetBuffer(kernel, "pressures", pressures);
        simulateParticles.SetBuffer(kernel, "forces", forces);
        simulateParticles.SetBuffer(kernel, "color_gradients", colorGradients);
        simulateParticles.SetBuffer(kernel, "surface_identifiers", surfaceIdentifiers);
        Profiler.BeginSample("Calculate Forces");
        simulateParticles.Dispatch(kernel, (int)math.ceil(particleCount / 32.0f), 1, 1);
        Profiler.EndSample();

        kernel = simulateParticles.FindKernel("CalculateSurfaceTensionForce");
        simulateParticles.SetBuffer(kernel, "sorted_particles", sortedParticles);
        simulateParticles.SetBuffer(kernel, "sorted_positions", sortedPositions);
        simulateParticles.SetBuffer(kernel, "bins", bins);
        simulateParticles.SetBuffer(kernel, "prefix_sums", prefixSums);
        simulateParticles.SetBuffer(kernel, "densities", densities);
        simulateParticles.SetBuffer(kernel, "forces", forces);
        simulateParticles.SetBuffer(kernel, "color_gradients", colorGradients);
        Profiler.BeginSample("Calculate Surface Tension Force");
        simulateParticles.Dispatch(kernel, (int)math.ceil(particleCount / 32.0f), 1, 1);
        Profiler.EndSample();

        kernel = simulateParticles.FindKernel("PerformTimeIntegration");
        simulateParticles.SetBuffer(kernel, "sorted_particles", sortedParticles);
        simulateParticles.SetBuffer(kernel, "sorted_positions", sortedPositions);
        simulateParticles.SetBuffer(kernel, "sorted_velocities", sortedVelocities);
        simulateParticles.SetBuffer(kernel, "positions", positions);
        simulateParticles.SetBuffer(kernel, "forces", forces);
        simulateParticles.SetBuffer(kernel, "velocities", velocities);
        simulateParticles.SetBuffer(kernel, "surface_identifiers", surfaceIdentifiers);
        simulateParticles.SetBuffer(kernel, "cell_tags", cellTags);
        Profiler.BeginSample("Perform Time Integration");
        simulateParticles.Dispatch(kernel, (int)math.ceil(particleCount / 32.0f), 1, 1);
        Profiler.EndSample();
    }

    void DrawSurface() {
        triangulatedCells.SetCounterValue(0);
        cellsToTriangulate.SetCounterValue(0);

        int kernel = surfaceReconstruction.FindKernel("ConvertCellTags");
        surfaceReconstruction.SetBuffer(kernel, "cell_tags", cellTags);
        surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate", cellsToTriangulate);
        Profiler.BeginSample("Convert Cell Tags");
        surfaceReconstruction.Dispatch(kernel, (int)math.ceil(marchingCubeGridPointCount / 1024.0f), 1, 1);
        Profiler.EndSample();

        ComputeBuffer.CopyCount(cellsToTriangulate, cellsToTriangulateArgs, 0);
        kernel = surfaceReconstruction.FindKernel("SetTriangulateArgs");
        surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate_args", cellsToTriangulateArgs);
        Profiler.BeginSample("Set Triangulate Args");
        surfaceReconstruction.Dispatch(kernel, 1, 1, 1);
        Profiler.EndSample();

        kernel = surfaceReconstruction.FindKernel("CalculateScalarField");
        surfaceReconstruction.SetBuffer(kernel, "cell_tags", cellTags);
        surfaceReconstruction.SetBuffer(kernel, "bins", bins);
        surfaceReconstruction.SetBuffer(kernel, "prefix_sums", prefixSums);
        surfaceReconstruction.SetBuffer(kernel, "positions", positions);
        surfaceReconstruction.SetBuffer(kernel, "scalar_field", scalarField);
        surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate_read", cellsToTriangulate);
        surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate_args", cellsToTriangulateArgs);
        Profiler.BeginSample("Calculate Scalar");
        surfaceReconstruction.DispatchIndirect(kernel, cellsToTriangulateArgs, 0);
        Profiler.EndSample();

        kernel = surfaceReconstruction.FindKernel("TriangulateCells");
        surfaceReconstruction.SetBuffer(kernel, "cell_tags", cellTags);
        surfaceReconstruction.SetBuffer(kernel, "scalar_field", scalarField);
        surfaceReconstruction.SetBuffer(kernel, "edge_table", edgeTableBuffer);
        surfaceReconstruction.SetBuffer(kernel, "triangle_connection_table", triangleConnectionTableBuffer);
        surfaceReconstruction.SetBuffer(kernel, "triangulated_cells", triangulatedCells);
        surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate_read", cellsToTriangulate);
        surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate_args", cellsToTriangulateArgs);
        Profiler.BeginSample("Triangulate Cells");
        surfaceReconstruction.DispatchIndirect(kernel, cellsToTriangulateArgs, 0);
        Profiler.EndSample();

        ComputeBuffer.CopyCount(triangulatedCells, drawProceduralIndirectArgs, 0);

        kernel = surfaceReconstruction.FindKernel("FillIndirectArgs");
        surfaceReconstruction.SetBuffer(kernel, "draw_procedural_indirect_args", drawProceduralIndirectArgs);
        Profiler.BeginSample("Fill Indirect Args");
        surfaceReconstruction.Dispatch(kernel, 1, 1, 1);
        Profiler.EndSample();

        MaterialPropertyBlock properties = new MaterialPropertyBlock();
        properties.SetFloat("visualScale", visualScale);
        properties.SetBuffer("vertices", triangulatedCells);

        Graphics.DrawProceduralIndirect(
            surfaceMaterial,
            new Bounds(transform.position, transform.lossyScale * 10000),
            MeshTopology.Triangles,
            drawProceduralIndirectArgs,
            0,
            null,
            properties,
            ShadowCastingMode.Off
        );
    }

    void DrawParticles() {
        int kernel = positionsToShaderData.FindKernel("PositionsToShaderData");
        positionsToShaderData.SetBuffer(kernel, "positions", positions);
        positionsToShaderData.SetBuffer(kernel, "velocities", velocities);
        positionsToShaderData.SetBuffer(kernel, "colors", colors);
        positionsToShaderData.SetBuffer(kernel, "_InstancedIndirectShaderData", particleShaderData);
        Profiler.BeginSample("Positions To Shader Data");
        positionsToShaderData.Dispatch(kernel, (int)math.ceil(particleCount / 32.0f), 1, 1);
        Profiler.EndSample();

        MaterialPropertyBlock properties = new MaterialPropertyBlock();
        properties.SetBuffer("_InstancedIndirectShaderData", particleShaderData);
        properties.SetBuffer("colors", colors);
        properties.SetMatrix("_LocalToWorld", transform.localToWorldMatrix);

        Graphics.DrawMeshInstancedIndirect(
            particleMesh,
            0,
            particleMaterial,
            new Bounds(transform.position, transform.lossyScale * 10000),
            particleShaderArgs,
            0,
            properties,
            ShadowCastingMode.Off
        );

    }

    void Start() {
        SetSortParameters();
        SetSurfaceReconstructionParameters();
        SetParticleRenderParameters();

        Scenario scenario;
        switch(scenarioName) {
            case ScenarioName.Cube:
                scenario = Scenarios.Cube(smoothingLength);
                break;
            case ScenarioName.DamBreak:
                scenario = Scenarios.DamBreak(smoothingLength);
                break;
            case ScenarioName.DoubleDamBreak:
                scenario = Scenarios.DoubleDamBreak(smoothingLength);
                break;
            case ScenarioName.DamBreakWithBunny:
                scenario = Scenarios.DoubleDamBreak(smoothingLength);
                break;
            default:
                scenario = Scenarios.DamBreak(smoothingLength);
                break;
        }

        SetSimulationParameters(scenario.particlePositions);
        InitializeComputeBuffers(scenario.particlePositions);
        InitializeBoundaries(scenario.boundaryPositions);

    }

    void RecordFrames() {
        ComputeBuffer countBuffer = new ComputeBuffer(1, sizeof(int), ComputeBufferType.Raw);

        string destinationFolder = Environment.GetFolderPath(Environment.SpecialFolder.Desktop) + @"\" + scenarioName.ToString();
        Directory.CreateDirectory(destinationFolder);

        int frame = 0;
        // Every second is 1000 steps
        for (int i = 0; i < 10000; ++i) {
            SortParticles();
            StepSimulation();

            triangulatedCells.SetCounterValue(0);
            cellsToTriangulate.SetCounterValue(0);

            int kernel = surfaceReconstruction.FindKernel("ConvertCellTags");
            surfaceReconstruction.SetBuffer(kernel, "cell_tags", cellTags);
            surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate", cellsToTriangulate);
            Profiler.BeginSample("Convert Cell Tags");
            surfaceReconstruction.Dispatch(kernel, (int)math.ceil(marchingCubeGridPointCount / 1024.0f), 1, 1);
            Profiler.EndSample();

            ComputeBuffer.CopyCount(cellsToTriangulate, cellsToTriangulateArgs, 0);
            kernel = surfaceReconstruction.FindKernel("SetTriangulateArgs");
            surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate_args", cellsToTriangulateArgs);
            Profiler.BeginSample("Set Triangulate Args");
            surfaceReconstruction.Dispatch(kernel, 1, 1, 1);
            Profiler.EndSample();

            kernel = surfaceReconstruction.FindKernel("CalculateScalarField");
            surfaceReconstruction.SetBuffer(kernel, "cell_tags", cellTags);
            surfaceReconstruction.SetBuffer(kernel, "bins", bins);
            surfaceReconstruction.SetBuffer(kernel, "prefix_sums", prefixSums);
            surfaceReconstruction.SetBuffer(kernel, "positions", positions);
            surfaceReconstruction.SetBuffer(kernel, "scalar_field", scalarField);
            surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate_read", cellsToTriangulate);
            surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate_args", cellsToTriangulateArgs);
            Profiler.BeginSample("Calculate Scalar");
            surfaceReconstruction.DispatchIndirect(kernel, cellsToTriangulateArgs, 0);
            Profiler.EndSample();

            kernel = surfaceReconstruction.FindKernel("TriangulateCells");
            surfaceReconstruction.SetBuffer(kernel, "cell_tags", cellTags);
            surfaceReconstruction.SetBuffer(kernel, "scalar_field", scalarField);
            surfaceReconstruction.SetBuffer(kernel, "edge_table", edgeTableBuffer);
            surfaceReconstruction.SetBuffer(kernel, "triangle_connection_table", triangleConnectionTableBuffer);
            surfaceReconstruction.SetBuffer(kernel, "triangulated_cells", triangulatedCells);
            surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate_read", cellsToTriangulate);
            surfaceReconstruction.SetBuffer(kernel, "cells_to_triangulate_args", cellsToTriangulateArgs);
            Profiler.BeginSample("Triangulate Cells");
            surfaceReconstruction.DispatchIndirect(kernel, cellsToTriangulateArgs, 0);
            Profiler.EndSample();

            if (i == 700) {
                ComputeBuffer.CopyCount(triangulatedCells, countBuffer, 0);
                int[] count = new int[1];
                countBuffer.GetData(count);
                float3x2[] vertices = new float3x2[count[0] * 3];
                triangulatedCells.GetData(vertices);
                ParticleMesher.WriteParticlesToDisk(vertices, destinationFolder + @"\" + string.Format("Particles_Frame{0:D3}", frame));
                ++frame;
                print("Recorded Frame " + frame);
            }
        }

        countBuffer.Dispose();
    }

    void Update() {
        if(recordingDone) {
            return;
        }

        if(record) {
            RecordFrames();
            recordingDone = true;
        }

        simulateParticles.SetFloat("dt", pauseSimulation ? 0.0f : timeStep);

        SortParticles();
        StepSimulation();

        if (drawSurface) {
            DrawSurface();
        }

        if (drawParticles) {
            DrawParticles();
        }
    }

    private void OnDrawGizmos() {
        Gizmos.color = Color.white;
        if (positionsToDraw != null) {
            foreach (float3 p in positionsToDraw) {
                Gizmos.DrawSphere(p * visualScale, particleRadius * visualScale);
            }
        }

        //int steps = 5;        
        //float step = bounds.x / 5.0f;

        //for(int z = 0; z < steps; ++z) {
        //    for(int y = 0; y < steps; ++y) {
        //        for(int x = 0; x < steps; ++x) {
        //            Vector3 cube = new Vector3(x, y, z) * step + new Vector3(step / 2.0f, step / 2.0f, step / 2.0f);
        //            Gizmos.DrawWireCube(cube * visualScale, new Vector3(step * visualScale, step * visualScale, step * visualScale));
        //        }
        //    }
        //}

        Vector3 center = new Vector3(bounds.x / 2.0f, bounds.y / 2.0f, bounds.z / 2.0f);
        Gizmos.DrawWireCube(center * visualScale, (bounds - new float3(smoothingLength) * 2.0f) * visualScale);
    }

    private void OnDestroy() {
        positions.Dispose();
        sortedPositions.Dispose();
        particles.Dispose();
        sortedParticles.Dispose();
        velocities.Dispose();
        sortedVelocities.Dispose();
        bins.Dispose();
        prefixSums.Dispose();
        boundaryBins.Dispose();
        boundaryPrefixSums.Dispose();
        pressures.Dispose();
        densities.Dispose();
        colorGradients.Dispose();
        forces.Dispose();
        surfaceIdentifiers.Dispose();
        cellTags.Dispose();
        cellsToTriangulate.Dispose();
        cellsToTriangulateArgs.Dispose();
        scalarField.Dispose();
        triangulatedCells.Dispose();
        triangleConnectionTableBuffer.Dispose();
        edgeTableBuffer.Dispose();
        drawProceduralIndirectArgs.Dispose();
        particleShaderArgs.Dispose();
        particleShaderData.Dispose();
        colors.Dispose();
        boundaryParticles.Dispose();
        boundaryPositions.Dispose();
        boundaryVolume.Dispose();
    }
}
