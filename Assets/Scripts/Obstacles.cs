using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;
using Random = Unity.Mathematics.Random;

public class Obstacles {
    static Random rnd = new Random(0xdeadbeef);

    static bool SampleRandomPoint(List<float3> currentPoints, List<Vector3> vertices, float allowedDistance) {
        int randomFace = (int)(rnd.NextFloat() * (vertices.Count / 3));
        float3 v0 = vertices[randomFace + 0];
        float3 v1 = vertices[randomFace + 1];
        float3 v2 = vertices[randomFace + 2];
        float u0 = rnd.NextFloat();
        float u1 = rnd.NextFloat();
        float3 randomPoint = (1 - math.sqrt(u0)) * v0 + ((1 - u1) * math.sqrt(u0)) * v1 + (u1 * math.sqrt(u0)) * v2;

        float allowedDistanceSquared = allowedDistance * allowedDistance;
        foreach(float3 p in currentPoints) {
            if(math.dot(p, p) < allowedDistanceSquared) {
                return false;
            }
        }

        currentPoints.Add(randomPoint);
        return true;
    }

    public static List<float3> GetBoundaryPositionsFromMesh(Mesh mesh, float smoothingLength) {
        List<Vector3> vertices = new List<Vector3>();
        mesh.GetVertices(vertices);
        List<float3> currentPoints = new List<float3>();

        int consecutiveFails = 0;
        while(consecutiveFails < 50) {
            bool sampled = SampleRandomPoint(currentPoints, vertices, smoothingLength);
            consecutiveFails = sampled ? consecutiveFails + 1 : 0;
        }

        return currentPoints;
    }
}