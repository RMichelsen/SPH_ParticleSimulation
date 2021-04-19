using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using Unity.Mathematics;
using Unity.VisualScripting.FullSerializer;
using UnityEngine;
using Random = Unity.Mathematics.Random;

public class Obstacles {
    static Random rnd = new Random(0xdeadbeef);

    static List<float> GetTriangleCDFTable(List<Vector3> vertices) {
        List<float> triangleCDFTable = new List<float>();
        float totalSurfaceArea = 0.0f;

        for(int i = 0; i < (vertices.Count / 3); i += 3) {
            Vector3 p0 = vertices[i + 0];
            Vector3 a = vertices[i + 1] - p0;
            Vector3 b = vertices[i + 2] - p0;
            float area = 0.5f * math.length(math.cross(a, b));
            totalSurfaceArea += area;
            triangleCDFTable.Add(totalSurfaceArea);
        }

        for (int i = 0; i < triangleCDFTable.Count; ++i) {
            triangleCDFTable[i] /= totalSurfaceArea;
        }

        return triangleCDFTable;
    }


    public static List<float3> GetBoundaryPositionsFromMesh(float scale, float3 offset, Mesh mesh, float smoothingLength) {
        float allowedDistanceSquared = smoothingLength * smoothingLength;

        List<Vector3> vertices = new List<Vector3>();
        mesh.GetVertices(vertices);

        List<float> triangleCDFTable = GetTriangleCDFTable(vertices);
        List<float3> points = new List<float3>();
        for (int i = 0; i < 50000; ++i) {
            int index = triangleCDFTable.BinarySearch(rnd.NextFloat());
            if(index < 0) {
                index = ~index;
            }
            index *= 3;

            float3 v0 = vertices[index + 0];
            float3 v1 = vertices[index + 1];
            float3 v2 = vertices[index + 2];
            float u0 = rnd.NextFloat();
            float u1 = rnd.NextFloat();
            float3 randomPoint = (1 - math.sqrt(u0)) * v0 + ((1 - u1) * math.sqrt(u0)) * v1 + (u1 * math.sqrt(u0)) * v2;
            float3 candidate = scale * randomPoint + offset;

            foreach (float3 p in points) {
                if (math.dot(candidate, p) < allowedDistanceSquared) {
                    --i;
                    continue;
                }
            }

            points.Add(scale * randomPoint + offset);
        }

        return points;
    }
}