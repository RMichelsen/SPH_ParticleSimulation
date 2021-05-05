using System.Collections.Generic;
using System.Runtime.Versioning;
using Unity.Mathematics;
using UnityEngine;
using Random = Unity.Mathematics.Random;

public struct Scenario {
    public List<float3> particlePositions;
    public List<float3> boundaryPositions;
}

public class Scenarios {
    static Random rnd = new Random(0xdeadbeef);

    public static Scenario Cube(float smoothingLength) {
        List<float3> initialPositions = new List<float3>();
        float step = smoothingLength;
        int xSteps = (int)(1.38f / step);
        int ySteps = (int)(1.40f / step);
        int zSteps = (int)(1.38f / step);
        for (int z = 0; z < zSteps; ++z) {
            for (int y = 0; y < ySteps; ++y) {
                for (int x = 0; x < xSteps; ++x) {
                    initialPositions.Add((new float3(x, y, z) + new float3(1.5f)) * step + new float3(0.81f, 0.30f, 0.81f));
                }
            }
        }

        int remainingParticles = (32 - (initialPositions.Count % 32)) % 32;
        for (int i = 0; i < remainingParticles; ++i) {
            initialPositions.Add(
                new float3(0.81f, 0.30f, 0.81f) + rnd.NextFloat3(new float3(smoothingLength * 1.5f), new float3(0.81f, 0.30f, 0.81f) - new float3(smoothingLength * 1.5f))
            );
        }

        return new Scenario {
            particlePositions = initialPositions,
            boundaryPositions = new List<float3>()
        };
    }

    public static Scenario CubeWithBunny(float smoothingLength)
    {
        List<float3> initialPositions = new List<float3>();
        float step = smoothingLength;
        int xSteps = (int)(1.38f / step);
        int ySteps = (int)(0.7f / step);
        int zSteps = (int)(1.38f / step);
        for (int z = 0; z < zSteps; ++z)
        {
            for (int y = 0; y < ySteps; ++y)
            {
                for (int x = 0; x < xSteps; ++x)
                {
                    initialPositions.Add((new float3(x, y, z) + new float3(1.5f)) * step + new float3(0.81f, 1.0f, 0.81f));
                }
            }
        }

        int remainingParticles = (32 - (initialPositions.Count % 32)) % 32;
        for (int i = 0; i < remainingParticles; ++i)
        {
            initialPositions.Add(
                new float3(0.81f, 1.0f, 0.81f) + rnd.NextFloat3(new float3(smoothingLength * 1.5f), new float3(1.38f, 0.7f, 1.38f) - new float3(smoothingLength * 1.5f))
            );
        }

        GameObject bunny = Resources.Load<GameObject>("Meshes/bunny");
        MeshFilter meshFilter = bunny.GetComponentInChildren<MeshFilter>();

        return new Scenario
        {
            particlePositions = initialPositions,
            boundaryPositions = Obstacles.GetBoundaryPositionsFromMesh(0.5f, new float3(1.4f, 0.0f, 1.4f), meshFilter.sharedMesh, smoothingLength)
        };
    }

    public static Scenario DamBreak(float smoothingLength) {
        List<float3> initialPositions = new List<float3>();
        float step = smoothingLength;
        int xSteps = (int)(1.0f / step);
        int ySteps = (int)(2.0f / step);
        int zSteps = (int)(1.0f / step);
        for (int z = 0; z < zSteps; ++z) {
            for (int y = 0; y < ySteps; ++y) {
                for (int x = 0; x < xSteps; ++x) {
                    initialPositions.Add((new float3(x, y, z) + new float3(1.5f)) * step);
                }
            }
        }

        int remainingParticles = (32 - (initialPositions.Count % 32)) % 32;
        for (int i = 0; i < remainingParticles; ++i) {
            initialPositions.Add(
                rnd.NextFloat3(new float3(smoothingLength * 1.5f), new float3(1.0f, 2.0f, 1.0f) - new float3(smoothingLength * 1.5f))
            );
        }

        return new Scenario {
            particlePositions = initialPositions,
            boundaryPositions = new List<float3>()
        };
    }

    public static Scenario DoubleDamBreak(float smoothingLength) {
        List<float3> initialPositions = new List<float3>();
        float step = smoothingLength;
        int xSteps = (int)(1.0f / step);
        int ySteps = (int)(2.0f / step);
        int zSteps = (int)(1.0f / step);
        for (int z = 0; z < zSteps; ++z) {
            for (int y = 0; y < ySteps; ++y) {
                for (int x = 0; x < xSteps; ++x) {
                    initialPositions.Add((new float3(x, y, z) + new float3(1.5f)) * step);
                }
            }
        }
        xSteps = (int)(0.75f / step);
        ySteps = (int)(1.5f / step);
        zSteps = (int)(0.75f / step);
        for (int z = 0; z < zSteps; ++z) {
            for (int y = 0; y < ySteps; ++y) {
                for (int x = 0; x < xSteps; ++x) {
                    initialPositions.Add((new float3(x, y, z) + new float3(1.5f)) * step + new float3(2.25f - (smoothingLength * 1.5f), 0.0f, 2.25f - (smoothingLength * 1.5f)));
                }
            }
        }

        int remainingParticles = (32 - (initialPositions.Count % 32)) % 32;
        for (int i = 0; i < remainingParticles; ++i) {
            initialPositions.Add(
                rnd.NextFloat3(new float3(smoothingLength * 1.5f), new float3(1.0f, 2.0f, 1.0f) - new float3(smoothingLength * 1.5f))
            );
        }

        return new Scenario {
            particlePositions = initialPositions, 
            boundaryPositions = new List<float3>()
        };
    }

    public static Scenario FloodWithCubes(float smoothingLength) {
        List<float3> initialPositions = new List<float3>();
        float step = smoothingLength;
        int xSteps = (int)(2.9f / step);
        int ySteps = (int)(2.0f / step);
        int zSteps = (int)(0.3f / step);
        for (int z = 0; z < zSteps; ++z) {
            for (int y = 0; y < ySteps; ++y) {
                for (int x = 0; x < xSteps; ++x) {
                    initialPositions.Add((new float3(x, y, z) + new float3(1.5f)) * step);
                }
            }
        }
        int remainingParticles = (32 - (initialPositions.Count % 32)) % 32;
        for (int i = 0; i < remainingParticles; ++i)
        {
            initialPositions.Add(
                rnd.NextFloat3(new float3(smoothingLength * 1.5f), new float3(2.9f, 2.0f, 0.3f) - new float3(smoothingLength * 1.5f))
            );
        }

        List<float3> boundaryPositions = new List<float3>();
        Obstacles.AddCubeToBoundaries(new float3(1.0f, 0.0f, 1.0f), new float3(0.2f, 1.2f, 0.2f), boundaryPositions, smoothingLength);
        Obstacles.AddCubeToBoundaries(new float3(2.0f, 0.0f, 2.0f), new float3(0.3f, 0.7f, 0.4f), boundaryPositions, smoothingLength);
        Obstacles.AddCubeToBoundaries(new float3(2.0f, 0.0f, 0.5f), new float3(0.15f, 0.2f, 0.4f), boundaryPositions, smoothingLength);
        Obstacles.AddCubeToBoundaries(new float3(2.5f, 0.0f, 1.3f), new float3(0.1f, 2.2f, 0.1f), boundaryPositions, smoothingLength);
        Obstacles.AddCubeToBoundaries(new float3(1.2f, 0.0f, 1.4f), new float3(0.15f, 0.4f, 0.1f), boundaryPositions, smoothingLength);
        Obstacles.AddCubeToBoundaries(new float3(0.1f, 0.0f, 2.3f), new float3(1.0f, 0.2f, 0.2f), boundaryPositions, smoothingLength);

        return new Scenario {
            particlePositions = initialPositions,
            boundaryPositions = boundaryPositions
        };
    }

    //    MeshFilter obstacleMesh = obstacle.GetComponentInChildren<MeshFilter>();

    //    List<float3> initBoundaryPositions = new List<float3>();
    //    positionsToDraw = new List<float3>();
    //        for (int i = 0; i<obstacleMesh.mesh.vertexCount; i += 180) {
    //            if (i > obstacleMesh.mesh.vertexCount) break;
    //            if (i + 1 > obstacleMesh.mesh.vertexCount) break;
    //            if (i + 2 > obstacleMesh.mesh.vertexCount) break;
    //            Vector3 v0 = obstacleMesh.mesh.vertices[i + 0];
    //    Vector3 v1 = obstacleMesh.mesh.vertices[i + 1];
    //    Vector3 v2 = obstacleMesh.mesh.vertices[i + 2];

    //    // Formula to calculate centroid
    //    float x = (v0.x + v1.x + v2.x) / 3.0f;
    //    float y = (v0.y + v1.y + v2.y) / 3.0f;
    //    float z = (v0.z + v1.z + v2.z) / 3.0f;

    //    //initBoundaryPositions.Add((new float3(x, y, z) + new float3(3.5f, 0.0f, 3.5f)) / visualScale * 40.0f);
    //    //positionsToDraw.Add((new float3(x, y, z) + new float3(3.5f, 0.0f, 3.5f)) / visualScale * 40.0f);
    //}


}

